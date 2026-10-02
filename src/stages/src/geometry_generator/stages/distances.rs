use std::collections::VecDeque;

use super::{GeneratorStage, PlanetData};

#[derive(Debug, Default)]
pub struct AssignDistances{}

impl AssignDistances{
    fn assign_distances(&self, data: &mut PlanetData) -> anyhow::Result<()>{
        data.cells.convergent_dist = self.assign_distance_field(&data.convergent_cells, data);
        data.cells.divergent_dist = self.assign_distance_field(&data.divergent_cells, data);
        Ok(())
    }

    fn assign_distance_field(
        &self,
        seeds: &Vec<u32>,
        data: &PlanetData
    ) -> Vec<f32> {
        let cell_count = data.cells.position.len() as u32;

        let mut distances = vec![f32::INFINITY; cell_count as usize];
        let mut queue: VecDeque<u32> = VecDeque::new();


        for c in seeds{
            queue.push_back(*c);
            distances[*c as usize] = 0.0;
        }

        while let Some(cell) = queue.pop_front(){
            let cell_idx = cell as usize;

            let current_distance = distances[cell_idx];
            let current_center = data.cells.position[cell_idx];
            
            for neighbor in data.mesh.r_circulate_r(cell) {

                let d = &mut distances[neighbor as usize];
                
                if *d != f32::INFINITY {
                    continue;
                }

                let n_pos = data.cells.position[neighbor as usize];

                *d = current_distance + (current_center - n_pos).length();
                queue.push_back(neighbor as u32);
            }
        }

        distances
    }
}

impl GeneratorStage for AssignDistances{
    fn name(&self) -> &'static str {
        "Assign Distances"
    }

    fn run(&self, planet: &mut PlanetData) -> anyhow::Result<()> {
        self.assign_distances(planet)
    }
}