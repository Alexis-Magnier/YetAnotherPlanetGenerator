use super::{GeneratorStage, PlanetData};

#[derive(Debug, Default)]
pub struct Strength{
    pub noise_seed: u32,
    pub noise_octaves: u32,
    pub noise_frequency: f32,
    pub noise_lacunarity: f32,
    pub noise_persistence: f32,
    pub noise_amplitude: f32,
    pub noise_mean: f32,
    pub noise_weight: f32,
    pub thermal_weight: f32,
}

impl Strength{

    fn add_noise(&self, data: &mut PlanetData) -> anyhow::Result<()>{

        use noise::{Fbm, MultiFractal, NoiseFn, Simplex};

        let noise = Fbm::<Simplex>::new(self.noise_seed)
            .set_octaves(self.noise_octaves as usize)
            .set_frequency(self.noise_frequency as f64)
            .set_lacunarity(self.noise_lacunarity as f64)
            .set_persistence(self.noise_persistence as f64);

        data.cells.yeild_strength = (0..data.cells.count).map(|c|{
            let c_idx = c as usize;

            let p_idx = data.cells.plate_id[c_idx] as usize;

            let position = data.cells.position[c_idx];

            let base = noise.get(position.as_dvec3().to_array()) as f32;

            let amplitude = data.plates.strenght_ampl[p_idx] * self.noise_amplitude;
            let mean = data.plates.strenght_mean[p_idx] + self.noise_mean;

            (base * amplitude + mean).clamp(0., 1.0) * self.noise_weight
        }).collect();

        Ok(())
    }

    fn compute_strength(&self, data: &mut PlanetData) -> anyhow::Result<()>{
        
        data.cells.yeild_strength
            .iter_mut()
            .enumerate()
            .for_each(|(idx, strength)|{
                let t = data.cells.temperature[idx];

                *strength /= 1. + t * self.thermal_weight;
            });
            

        Ok(())
    }
}

impl GeneratorStage for Strength{
    fn name(&self) -> &'static str {
        "Strength"
    }

    fn run(&self, planet: &mut PlanetData) -> anyhow::Result<()> {
        self.add_noise(planet)?;
        self.compute_strength(planet)?;
        Ok(())
    }
}