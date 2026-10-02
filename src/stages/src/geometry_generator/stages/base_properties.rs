use noise::{Fbm, MultiFractal, NoiseFn, Simplex};

use super::{GeneratorStage, PlanetData};

#[derive(Debug, Default)]
pub struct AssignBaseProperties{
    pub density_noise_seed: u32,
    pub height_noise_seed: u32,
    pub octaves: u32,
    pub frequency: f32,
    pub lacunarity: f32,
    pub persistence: f32,
}

impl AssignBaseProperties{
    fn assign_density(&self, data: &mut PlanetData) -> anyhow::Result<()>{
        let noise = Fbm::<Simplex>::new(self.density_noise_seed)
            .set_octaves(self.octaves as usize)
            .set_frequency(self.frequency as f64)
            .set_lacunarity(self.lacunarity as f64)
            .set_persistence(self.persistence as f64);
    
        data.cells.density = (0..data.cells.count).map(|cell|{
            let cell_idx = cell as usize;

            let plate_idx = data.cells.plate_id[cell_idx] as usize;

            let base_noise = {
                let base = noise.get(data.cells.position[cell_idx].as_dvec3().to_array()) as f32;
                
                let ampl = data.plates.density_ampl[plate_idx];
                let mean = data.plates.density_mean[plate_idx];

                base * ampl + mean
            };

            base_noise
        }).collect();

        Ok(())
    }

    fn assign_width(&self, data: &mut PlanetData) -> anyhow::Result<()>{
        let noise = Fbm::<Simplex>::new(self.height_noise_seed)
            .set_octaves(self.octaves as usize)
            .set_frequency(self.frequency as f64)
            .set_lacunarity(self.lacunarity as f64)
            .set_persistence(self.persistence as f64);

        data.cells.thickness = (0..data.cells.count).map(|cell|{
            let cell_idx = cell as usize;

            let plate_idx = data.cells.plate_id[cell_idx] as usize;

            let base_noise = {
                let base = noise.get(data.cells.position[cell_idx].as_dvec3().to_array()) as f32;
                
                let ampl = data.plates.thickness_ampl[plate_idx];
                let mean = data.plates.thickness_mean[plate_idx];

                base * ampl + mean
            };

            base_noise
        }).collect();

        Ok(())
    }
}

impl GeneratorStage for AssignBaseProperties{
    fn name(&self) -> &'static str {
        "Assign properties"
    }

    fn run(&self, planet: &mut PlanetData) -> anyhow::Result<()> {
        self.assign_density(planet)?;
        self.assign_width(planet)?;

        Ok(())
    }
}