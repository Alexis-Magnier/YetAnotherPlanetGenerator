use std::collections::VecDeque;
use glam::Vec3;
use super::{GeneratorStage, PlanetData};

#[derive(Debug, Default)]
pub struct Stress{
    pub transmission_threshold: f32,
    pub base_strength: f32,
    pub strength_mult: f32,
    pub base_stress_mult: f32,
    pub reference_distance: f32,
    pub diffuse_weight: f32,
    pub diffuse_steps: u32,
    pub noise_seed: u32,
    pub noise_octaves: u32,
    pub noise_frequency: f32,
    pub noise_amplitude: f32,
    pub noise_lacunarity: f32,
    pub noise_persistence: f32,
    pub noise_strength_weight: f32,
    pub noise_max: f32,
}

impl Stress{
    fn transfer_stress(&self, data: &mut PlanetData) -> anyhow::Result<()> {
        let stresses = &mut data.cells.stress;

        // Each queued packet now carries the physical distance it has covered
        // since it last went through a full directional split. This lets the
        // fan-out logic below be driven by real distance instead of by the
        // number of mesh edges crossed — which is what previously made the
        // effective transmission radius shrink in denser regions of the mesh
        // (more, shorter hops there meant more splitting events per unit of
        // physical distance travelled).
        let mut queue: VecDeque<(u32, Vec3, f32)> = data
            .convergent_cells
            .iter()
            .chain(&data.divergent_cells)
            .chain(&data.transform_cells)
            .copied()
            .map(|a| (a, stresses[a as usize] * self.base_stress_mult, 0.0))
            .collect();

        stresses.fill(Vec3::ZERO);

        let ref_distance = self.reference_distance;

        // Reused across iterations instead of reallocated every time.
        let mut neighbors_info: Vec<(u32, f32, f32)> = Vec::with_capacity(10);

        while let Some((c, stress, dist_since_split)) = queue.pop_front() {
            let c_idx = c as usize;
            let pos_c = data.cells.position[c_idx];

            let stress_len = stress.length();
            if stress_len < self.transmission_threshold {
                continue;
            }

            let force_dir = stress / stress_len;

            let c_strength = (data.cells.yeild_strength[c_idx] * self.strength_mult + self.base_strength)
                .clamp(0.0, 1.0);

            neighbors_info.clear();
            let mut total_weight = 0.0;
            let mut best_neighbor: Option<u32> = None;
            let mut best_alignment = f32::NEG_INFINITY;
            let mut best_distance = 0.0f32;

            for n in data.mesh.r_circulate_r(c) {
                let n_idx = n as usize;
                let pos_n = data.cells.position[n_idx];

                let vec_to_neighbor = pos_n - pos_c;
                let distance = vec_to_neighbor.length();

                if distance > 0.0001 {
                    let dir_to_neighbor = vec_to_neighbor / distance;
                    let alignment: f32 = force_dir.dot(dir_to_neighbor);
                    let forward_alignment = alignment.max(0.0);

                    if alignment > best_alignment {
                        best_alignment = alignment;
                        best_neighbor = Some(n);
                        best_distance = distance;
                    }

                    // Directional weight for splitting force among forward neighbors
                    let weight = (forward_alignment + 0.1) / distance;

                    total_weight += weight;
                    neighbors_info.push((n, weight, distance));
                }
            }

            if total_weight > 0.0 {
                let best = best_neighbor.expect("total_weight > 0 implies a best neighbor exists");
                let mut total_absorbed = Vec3::ZERO;

                // Ramp from 0 (send everything undiluted to the single
                // best-aligned neighbor) to 1 (the original fully-normalized
                // split across every forward neighbor) as the distance covered
                // since the last split approaches `ref_distance`. A chain of
                // many short hops through dense cells now behaves like one
                // longer hop of the same total physical length.
                let cumulative_distance = dist_since_split + best_distance;
                let split_factor = (cumulative_distance / ref_distance).clamp(0.0, 1.0);
                let next_dist_since_split = if split_factor >= 1.0 { 0.0 } else { cumulative_distance };

                for (n, weight, distance) in neighbors_info.iter().copied() {
                    let split_weight_fraction = weight / total_weight;
                    let pass_through_fraction = if n == best { 1.0 } else { 0.0 };

                    let weight_fraction = split_factor * split_weight_fraction
                        + (1.0 - split_factor) * pass_through_fraction;

                    if weight_fraction <= 0.0 {
                        continue;
                    }

                    let directed_stress = stress * weight_fraction;

                    // Scale transmission fraction exponentially with step distance: T = strength ^ (d / ref_d)
                    let step_transmission = c_strength.powf(distance / ref_distance);

                    let transmitted_to_n = directed_stress * step_transmission;
                    let absorbed_from_n = directed_stress * (1.0 - step_transmission);

                    total_absorbed += absorbed_from_n;

                    if transmitted_to_n.length() >= self.transmission_threshold {
                        queue.push_back((n, transmitted_to_n, next_dist_since_split));
                    } else {
                        // Previously this remainder was silently dropped instead of
                        // absorbed, quietly leaking stress out of the system (worse
                        // wherever splitting produced many small fragments). Absorb
                        // it here instead so total stress is conserved.
                        total_absorbed += transmitted_to_n;
                    }
                }

                stresses[c_idx] += total_absorbed;
            } else {
                // If there are no valid neighbors, absorb all remaining stress
                stresses[c_idx] += stress;
            }
        }

        Ok(())
    }

    fn diffuse(&self, data: &mut PlanetData) -> anyhow::Result<()>{

        let stress = &mut data.cells.stress;
        let weight = self.diffuse_weight;

        for c in 0..data.cells.count{
            let c_idx = c as usize;

            let c_s = stress[c_idx];
            let c_center = data.cells.position[c_idx];

            let mut delta_t = Vec3::ZERO;

            for n in data.mesh.r_circulate_r(c){
                let n_idx = n as usize;

                let n_s = stress[n_idx];

                let n_center = data.cells.position[n_idx];
                let dist = n_center.distance(c_center);

                delta_t += (n_s - c_s) * weight / (1.0 + dist);
            }

            stress[c_idx] = c_s + delta_t;
        }

        Ok(())
    }

    fn diffuse_stress(&self, data: &mut PlanetData) -> anyhow::Result<()>{
        
        for _ in 0..self.diffuse_steps{
            self.diffuse(data)?;
        }

        Ok(())
    }

    fn add_noise(&self, data: &mut PlanetData) -> anyhow::Result<()>{
        use noise::{Fbm, MultiFractal, NoiseFn, Simplex};

        let noise = Fbm::<Simplex>::new(self.noise_seed)
            .set_octaves(self.noise_octaves as usize)           // Number of noise layers stacked together
            .set_frequency(self.noise_frequency as f64)       // The initial scale of the noise
            .set_lacunarity(self.noise_lacunarity as f64)      // How much the frequency increases per octave
            .set_persistence(self.noise_persistence as f64);

        data.cells.stress.iter_mut()
            .enumerate()
            .for_each(|(idx, temp)|{

                let position = &data.cells.position[idx];
                let strength = data.cells.yeild_strength[idx];
                
                let n = noise.get(position.as_dvec3().to_array()) as f32 * self.noise_amplitude;

                *temp += (n / (1.0 + strength * self.noise_strength_weight)).min(self.noise_max);
        });

        Ok(())
    }
}

impl GeneratorStage for Stress{
    fn name(&self) -> &'static str {
        "Stress"
    }

    fn run(&self, planet: &mut PlanetData) -> anyhow::Result<()> {

        self.transfer_stress(planet)?;
        self.diffuse_stress(planet)?;
        self.add_noise(planet)?;

        Ok(())
    }
}