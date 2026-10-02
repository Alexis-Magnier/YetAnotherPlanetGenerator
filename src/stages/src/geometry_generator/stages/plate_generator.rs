use core::f32;

use anyhow::Context;

use rand::rngs::StdRng;
use rand::{Rng, RngExt, SeedableRng};
use glam::Vec3;

use priority_queue::PriorityQueue;

use crate::geometry_generator::PlanetData;
use crate::geometry_generator::stage::GeneratorStage;

pub struct WeightedSampler {
    cumulative_weights: Vec<f32>,
    total_weight: f32,
}

impl WeightedSampler {
    pub fn new(weights: &[f32]) -> Option<Self> {
        if weights.is_empty() {
            return None;
        }

        let mut cumulative_weights = Vec::with_capacity(weights.len());
        let mut total = 0.0;

        for &weight in weights {
            // Reject invalid weights to prevent silent logic bugs
            if weight < 0.0 || weight.is_nan() {
                return None;
            }
            total += weight;
            cumulative_weights.push(total);
        }

        if total <= 0.0 {
            return None; 
        }

        Some(Self {
            cumulative_weights,
            total_weight: total,
        })
    }

    pub fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        let target = rng.random_range(0.0..self.total_weight);
        self.cumulative_weights.partition_point(|&w| w <= target)
    }
}

#[derive(Debug, Default)]
pub struct PlateBaseMaterial{
    pub probabilistic_weight: f32,
    pub density_mean: f32,
    pub density_dev: f32,
    pub height_mean: f32,
    pub height_dev: f32,
    pub thermal_cooldown: f32,
    pub thermal_diffusion_rate: f32,
    pub strength_mean: f32,
    pub strength_dev: f32,
}

#[derive(Debug, Default)]
pub struct PlateGenerator{
    pub seed: u64,
    pub plate_count: u32,
    pub priority_mean: f32,
    pub priority_dev: f32,
    pub directionnal_preference_mean: f32,
    pub directionnal_preference_dev: f32,
    pub velocity_mean: f32,
    pub velocity_dev: f32,
    pub materials: Vec<PlateBaseMaterial>
}

impl PlateGenerator{
    fn make_rng(&self) -> impl Rng{
        Box::new(StdRng::seed_from_u64(self.seed as u64))
    }
}

impl PlateGenerator{
    /// Generates a list of `plate_count` seed cells.
    /// They are uniformally sampled and may form awkardly close seeds
    /// TODO : Maybe implement blue noise
    fn generate_seeds(&self, rng: &mut impl Rng, data: &PlanetData) -> Vec<u32>{
        let cell_count = data.cells.position.len() as u32;

        (0..self.plate_count).map(|_| {
            rng.next_u32() % cell_count
        }).collect()
    }

    fn assign_priorities(&self, rng: &mut impl Rng) -> anyhow::Result<Vec<f32>>{
        use rand_distr::{Normal, Distribution};

        let distribution = Normal::new(self.priority_mean, self.priority_dev)?;
        
        Ok((0..self.plate_count).map(|_| {
            distribution.sample(rng)
        }).collect())
    }

    fn assign_directions(&self, rng: &mut impl Rng) -> anyhow::Result<Vec<(f32, Vec3)>>{
        use rand_distr::{Normal, Distribution};

        let distribution = Normal::new(self.directionnal_preference_mean, self.directionnal_preference_dev)?;
        
        Ok((0..self.plate_count).map(|_| {
            (distribution.sample(rng), Vec3::from_array(rng.random()).normalize())
        }).collect())
    }

    fn generate_plates(&self, rng: &mut impl Rng, data: &mut PlanetData) -> anyhow::Result<()>{
        let cell_count = data.cells.position.len() as u32;

        let mut cells_plate = vec![u32::MAX; cell_count as usize];
        let mut seed_dist = vec![0.0; cell_count as usize];

        let mut queue: PriorityQueue<u32, u64> = PriorityQueue::new();

        let seeds = self.generate_seeds(rng, data);
        
        // Push the seeds into the priority queue
        seeds.iter().enumerate().for_each(|(id, seed)|{

            // Assigned with highest priority to ensure that they are treated first
            queue.push(*seed, u64::MAX);
            cells_plate[*seed as usize] = id as u32;
        });

        let priorities = self.assign_priorities(rng)?;
        let directions = self.assign_directions(rng)?;

        while let Some((cell, current_cost)) = queue.pop(){
            let current_plate = cells_plate[cell as usize];
            let current_pos = data.cells.position[cell as usize];

            let (direction_weight, direction) = directions[current_plate as usize];
            let priority = priorities[current_plate as usize];

            for neighbor in data.mesh.r_circulate_r(cell) {
                let n_plate = cells_plate[neighbor as usize];
                
                // Skip already assigned plates
                if n_plate != u32::MAX{
                    
                    // We could also find neighboring plates here, but not for now
                    continue;
                }

                let neighbor_position = data.cells.position[neighbor as usize];
                let dif = neighbor_position - current_pos;

                let dist = dif.length();

                let distance_cost = 1. / dif.length();

                let direction_cost = dif.dot(direction) * direction_weight;

                let cost = current_cost - (distance_cost * direction_cost * priority * 100.0).round() as u64;

                cells_plate[neighbor as usize] = current_plate;
                seed_dist[neighbor as usize] = seed_dist[cell as usize] + dist;

                queue.push(neighbor, cost);
            }
        }

        data.cells.plate_id = cells_plate;
        data.cells.seed_dist = seed_dist;
        data.plates.seeds = seeds;
        data.plates.count = self.plate_count;

        Ok(())
    }

    fn assign_plate_velocity(&self, rng: &mut impl Rng, data: &mut PlanetData) -> anyhow::Result<()>{
        use rand_distr::{Normal, Distribution};

        let distribution = Normal::new(self.velocity_mean, self.velocity_dev)?;

        data.plates.velocity = (0..self.plate_count).map(|p| {
            let seed = data.plates.seeds[p as usize];
            let center = data.cells.position[seed as usize];

            let neighbor = data.mesh.r_circulate_r(seed).next().expect("Seed region has no neighbor");
            
            (center - data.cells.position[neighbor as usize]).normalize() * distribution.sample(rng)
        }).collect();

        Ok(())
    }

    fn assign_plate_properties(&self, rng: &mut impl Rng, data: &mut PlanetData) -> anyhow::Result<()>{

        let index_weights:Vec<f32> = self.materials.iter().map(|m|{
                m.probabilistic_weight
            }).collect();

        let sampler = WeightedSampler::new(&index_weights)
            .context("Failed to create sampler")?;

        let materials:Vec<u32> = (0..self.plate_count).map(|_| {
            sampler.sample(rng) as u32
        }).collect();
        
        let (density_mean, density_ampl) = self.make_density(rng, &materials)?;
        let (height_mean, height_ampl) = self.make_height(rng, &materials)?;
        let (strength_mean, strength_ampl) = self.make_strength(rng, &materials)?;

        data.plates.density_mean = density_mean;
        data.plates.density_ampl = density_ampl;

        data.plates.thickness_mean = height_mean;
        data.plates.thickness_ampl = height_ampl;

        data.plates.strenght_mean = strength_mean;
        data.plates.strenght_ampl = strength_ampl;

        data.plates.thermal_diffusion_rate = materials
            .iter()
            .map(|material_id|{
                let material = &self.materials[*material_id as usize];
                material.thermal_diffusion_rate
            }).collect();

        data.plates.thermal_cooldown = materials
            .iter()
            .map(|material_id| {
                let material = &self.materials[*material_id as usize];
                material.thermal_cooldown
            }).collect();

        Ok(())
    }

    fn make_density(&self, mut rng: &mut impl Rng, materials: &Vec<u32>) -> anyhow::Result<(Vec<f32>, Vec<f32>)> {
        use rand_distr::{Normal, Distribution};

        Ok((0..self.plate_count)
            .map(|p| {
    
                let material_id = materials[p as usize] as usize;
    
                let material = &self.materials[material_id];
    
                let normal = Normal::new(
                    material.density_mean,
                    material.density_dev,
                )?;
    
                Ok((normal.sample(&mut rng), material.density_dev))
            })
            .collect::<anyhow::Result<Vec<_>>>()?
            .into_iter()
            .unzip()
        )
    }

    fn make_height(&self, mut rng: &mut impl Rng, materials: &Vec<u32>) -> anyhow::Result<(Vec<f32>, Vec<f32>)> {
        use rand_distr::{Normal, Distribution};

        Ok((0..self.plate_count)
            .map(|p| {
    
                let material_id = materials[p as usize] as usize;
    
                let material = &self.materials[material_id];
    
                let normal = Normal::new(
                    material.height_mean,
                    material.height_dev,
                )?;
    
                Ok((normal.sample(&mut rng), material.density_dev))
            })
            .collect::<anyhow::Result<Vec<_>>>()?
            .into_iter()
            .unzip()
        )
    }

    fn make_strength(&self, mut rng: &mut impl Rng, materials: &Vec<u32>) -> anyhow::Result<(Vec<f32>, Vec<f32>)> {
        use rand_distr::{Normal, Distribution};

        Ok((0..self.plate_count)
            .map(|p| {
    
                let material_id = materials[p as usize] as usize;
    
                let material = &self.materials[material_id];
    
                let normal = Normal::new(
                    material.strength_mean,
                    material.strength_dev,
                )?;
    
                Ok((normal.sample(&mut rng).clamp(0.0, 1.0), material.strength_dev))
            })
            .collect::<anyhow::Result<Vec<_>>>()?
            .into_iter()
            .unzip()
        )
    }
}

impl GeneratorStage for PlateGenerator{
    fn name(&self) -> &'static str {
        "Plate Generator"
    }

    fn run(&self, data: &mut PlanetData) -> anyhow::Result<()> {
        let mut rng = self.make_rng();
        
        self.generate_plates(&mut rng, data)?;
        self.assign_plate_velocity(&mut rng, data)?;
        self.assign_plate_properties(&mut rng, data)?;

        Ok(())
    }
}