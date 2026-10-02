use glam::Vec3;
use super::{GeneratorStage, PlanetData};

#[derive(Debug, Default)]
pub struct Isostasy{
    pub mantle_density: f32,
    pub phase_iterations: u32,
    pub fold_wavelength_km: f32,
    pub fold_gain: f32,
    pub fault_spacing_km: f32,
    pub fault_gain: f32,
    pub yield_stress: f32,
    pub reference_thickness_km: f32,
    pub smoothing_iterations: u32,
    pub seed: u32,
}

impl Isostasy{

    fn stress_divergence(data: &PlanetData) -> Vec<f32>{
        (0..data.cells.count).map(|c|{
            let idx = c as usize;

            let p = data.cells.position[idx];
            let s = data.cells.stress[idx];

            let mut acc = 0.0;
            let mut count = 0.0;

            for neighbor in data.mesh.r_circulate_r(c as u32) {
                let j = neighbor as usize;
                let q = data.cells.position[j];

                let d = q - p;
                let len = d.length();

                if len < 1e-5 {
                    continue;
                }

                // directional derivative of S along the edge, projected on the edge
                acc += (data.cells.stress[j] - s).dot(d / len) / len;
                count += 1.0;
            }

            2.0 * acc / count

        }).collect()
    }

    fn isostatic_elevation(&self, data: &PlanetData) -> Vec<f32>{
        (0..data.cells.count).map(|c|{
            let idx = c as usize;

            let rho = data.cells.density[idx].min(self.mantle_density - 1.0);
            let buoyancy = (self.mantle_density - rho) / self.mantle_density;

            data.cells.thickness[idx] * buoyancy
        }).collect()
    }

    fn shortening_phase(&self, data: &PlanetData) -> Vec<f32> {
        let n = data.cells.position.len();

        // tangential unit stress direction per cell
        let dir: Vec<Vec3> = (0..n)
            .map(|c| {
                let up = data.cells.position[c]; // unit sphere -> surface normal
                let s = data.cells.stress[c];
                let t = s - up * s.dot(up);
                if t.length_squared() > 1e-12 {
                    t.normalize()
                } else {
                    Vec3::ZERO
                }
            })
            .collect();

        let mut phi = vec![0.0f32; n];
        let mut next = vec![0.0f32; n];

        for _ in 0..self.phase_iterations {
            for c in 0..n {
                let p = data.cells.position[c];

                let mut sum = 0.0f32;
                let mut count = 0.0f32;

                for neighbor in data.mesh.r_circulate_r(c as u32) {
                    let j = neighbor as usize;
                    let q = data.cells.position[j];

                    // phi(q) ~ phi(p) + g . (q - p)   =>   phi(p) ~ phi(q) - g . d
                    let g = (dir[c] + dir[j]) * 0.5;
                    sum += phi[j] - g.dot(q - p);
                    count += 1.0;
                }

                next[c] = if count > 0.0 { sum / count } else { phi[c] };
            }

            std::mem::swap(&mut phi, &mut next);
        }

        phi
    }


    fn add_deformation(
        &self,
        data: &PlanetData,
        phase: &[f32],
        div: &[f32],
        h: &mut [f32],
    ) {
        let n = data.cells.position.len();
        let tau = std::f32::consts::TAU;

        for c in 0..n {
            let stress = data.cells.stress[c].length();
            if stress <= 1e-6 {
                continue;
            }

            let strength = data.cells.thickness[c] / self.reference_thickness_km;
            let strain = stress / strength;

            let phi = phase[c];
            let compressive = div[c] < 0.0;

            let fold = (tau * phi / self.fold_wavelength_km).sin();
            h[c] += fold * strain * self.fold_gain;

            let excess = stress - self.yield_stress * strength;
            if excess <= 0.0 {
                continue;
            }

            let u = phi / self.fault_spacing_km;
            let block = u.floor();
            let frac = u - block;

            // random throw per block
            let r = hash_f32(block as i32 as u32 ^ self.seed) * 2.0 - 1.0;

            // In-block ramp. Compression stacks sheets upward with a steep
            // leading edge; extension tilts blocks and drops them.
            let ramp = if compressive {
                frac // sawtooth: gentle back-limb, sharp front
            } else {
                0.5 - frac // tilted block, down-dropped
            };

            let throw = excess * self.fault_gain / strength;
            h[c] += throw * (r * 0.6 + ramp * 0.8);
        }
    }

    fn relax(data: &PlanetData, h: &mut [f32], iterations: u32, rate: f32) {
        if iterations == 0 {
            return;
        }
        let n = h.len();
        let mut tmp = vec![0.0f32; n];

        for _ in 0..iterations {
            for c in 0..n {
                let mut sum = 0.0f32;
                let mut count = 0.0f32;
                for neighbor in data.mesh.r_circulate_r(c as u32) {
                    sum += h[neighbor as usize];
                    count += 1.0;
                }
                tmp[c] = if count > 0.0 {
                    h[c] + (sum / count - h[c]) * rate
                } else {
                    h[c]
                };
            }
            h.copy_from_slice(&tmp);
        }
    }


    fn compute(&self, data: &PlanetData) -> anyhow::Result<Vec<f32>>{

        let div = Self::stress_divergence(data);
        let mut h = self.isostatic_elevation(data);

        let phase = self.shortening_phase(data);

        self.add_deformation(
            data,
            &phase,
            &div,
            &mut h
        );

        Self::relax(data, &mut h, self.smoothing_iterations, 0.25);

        Ok(h)
    }
}

#[inline]
fn hash_f32(mut x: u32) -> f32 {
    x ^= x >> 16;
    x = x.wrapping_mul(0x7feb_352d);
    x ^= x >> 15;
    x = x.wrapping_mul(0x846c_a68b);
    x ^= x >> 16;
    (x >> 8) as f32 / 16_777_216.0
}

impl GeneratorStage for Isostasy{
    fn name(&self) -> &'static str {
        "Isostasy"
    }

    fn run(&self, planet: &mut PlanetData) -> anyhow::Result<()> {
        planet.cells.height = self.compute(planet)?;

        Ok(())
    }
}