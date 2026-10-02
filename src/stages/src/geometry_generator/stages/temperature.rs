use std::collections::VecDeque;
use bitvec::order::Msb0;
use super::{GeneratorStage, PlanetData};

#[derive(Debug, Default)]
pub struct Temperatures{
    pub divergent_temp: f32,
    pub convergent_temp: f32,
    pub transform_temp: f32,
    pub cooldown_factor: f32,
    pub diffuse_step: u32,
    pub diffuse_weight: f32,
    pub diffusion_coef: f32,
    pub noise_seed: u32,
    pub noise_octaves: u32,
    pub noise_frequency: f32,
    pub noise_amplitude: f32,
    pub noise_lacunarity: f32,
    pub noise_persistence: f32,
}

impl Temperatures{
    fn set_temperature(cells: &Vec<u32>, target: f32, temperatures: &mut Vec<f32>){
        for cell in cells{
            temperatures[*cell as usize] = target; 
        }
    }

    fn assign_base_temperature(&self, data: &mut PlanetData) -> anyhow::Result<Vec<f32>>{
        let mut temperatures = vec![0.0; data.cells.count as usize];

        Self::set_temperature(
            &data.divergent_cells,
            self.divergent_temp,
            &mut temperatures
        );

        Self::set_temperature(
            &data.convergent_cells,
            self.convergent_temp,
            &mut temperatures
        );

        Self::set_temperature(
            &data.transform_cells,
            self.transform_temp,
            &mut temperatures
        );

        Ok(temperatures)
    }

    fn transfer_tempetature(&self, data: &PlanetData, temperatures: &mut Vec<f32>) -> anyhow::Result<()>{

        use bitvec::bitvec;

        let mut queue: VecDeque<u32> = data
            .convergent_cells
            .iter()
            .chain(&data.divergent_cells)
            .chain(&data.transform_cells)
            .copied()
            .collect();
        
        let mut visited = bitvec![u64, Msb0; 0; data.cells.count as usize];
        
        while let Some(cell) = queue.pop_back(){

            let c_idx = cell as usize;
            
            let c_t = temperatures[c_idx];
            let c_center = data.cells.position[c_idx];

            if visited[c_idx]{
                continue;
            }

            visited.set(c_idx, true);

            for n in data.mesh.r_circulate_r(cell){
                let n_idx = n as usize;
                
                let n_temperature = temperatures[n_idx];

                let plate_idx = data.cells.plate_id[n_idx] as usize;
                let cooldown_factor = data.plates.thermal_cooldown[plate_idx];

                let n_center = data.cells.position[n_idx];
                let dist = n_center.distance(c_center);

                let n_t = c_t / (1.0 + cooldown_factor * dist * self.cooldown_factor);
                
                temperatures[n_idx] = (n_temperature + n_t) * 0.5;

                queue.push_front(n);
            }
        }

        temperatures
            .iter_mut()
            .for_each(|temp|{
                if temp.is_nan(){
                    *temp = 0.0;
                }
            });

        Ok(())
    }

    fn diffuse(&self, data: &PlanetData, temperatures: &mut Vec<f32>, forward: bool){
        let range = 
            if forward {
                0..data.cells.count
            } else {
                data.cells.count..0
            };

        let mut next_temperatures = temperatures.clone();

        for c in range {
            let c_idx = c as usize;
            let c_t = temperatures[c_idx];
            let c_center = data.cells.position[c_idx];

            let mut total_weight = 0.0;
            let mut weighted_heat_flux = 0.0;

            for n in data.mesh.r_circulate_r(c as u32) {
                let n_idx = n as usize;
                let n_t = temperatures[n_idx];
                let n_center = data.cells.position[n_idx];

                let dist = n_center.distance_squared(c_center);

                if dist > 0.0001 {
                    let weight = 1.0 / dist;

                    weighted_heat_flux += (n_t - c_t) * weight * self.diffuse_weight;
                    total_weight += weight;
                }
            }

            let plate_id = data.cells.plate_id[c as usize] as usize;
            let diffusion_rate = data.plates.thermal_diffusion_rate[plate_id];

            if total_weight > 0.0 {
                let average_gradient = weighted_heat_flux / total_weight;
                next_temperatures[c_idx] = c_t + self.diffusion_coef * diffusion_rate * average_gradient;
            }
        }

        *temperatures = next_temperatures;
    }

    fn diffuse_tempetature(&self, data: &PlanetData, temperatures: &mut Vec<f32>) -> anyhow::Result<()>{

        for i in 0..self.diffuse_step{
            self.diffuse(data, temperatures, i % 2 == 0);
        }

        Ok(())
    }

    fn add_noise(&self, data: &PlanetData, temperatures: &mut Vec<f32>) -> anyhow::Result<()>{
        use noise::{Fbm, MultiFractal, NoiseFn, Simplex};

        let noise = Fbm::<Simplex>::new(self.noise_seed)
            .set_octaves(self.noise_octaves as usize)           // Number of noise layers stacked together
            .set_frequency(self.noise_frequency as f64)       // The initial scale of the noise
            .set_lacunarity(self.noise_lacunarity as f64)      // How much the frequency increases per octave
            .set_persistence(self.noise_persistence as f64);

        temperatures.iter_mut()
            .enumerate()
            .for_each(|(idx, temp)|{
                let position = &data.cells.position[idx];

                *temp += noise.get(position.as_dvec3().to_array()) as f32 * self.noise_amplitude;
        });

        Ok(())
    }

}

impl GeneratorStage for Temperatures{
    fn name(&self) -> &'static str {
        "Temperatures"
    }

    fn run(&self, planet: &mut PlanetData) -> anyhow::Result<()> {
        let mut temperatures = self.assign_base_temperature(planet)?;
        self.transfer_tempetature(&planet, &mut temperatures)?;
        self.diffuse_tempetature(&planet, &mut temperatures)?;
        self.add_noise(&planet, &mut temperatures)?;

        planet.cells.temperature =  temperatures;

        Ok(())
    }
}