use tracing::{instrument, warn};
use yapg_core::{PlanetComponent, PropertyMap};
use crate::geometry_generator::stage::GeneratorStage;

use super::Mesh;
use super::data::PlanetData;
use glam::Vec3;

#[derive(Debug, Default)]
pub struct PlanetGeometry{
    pub mesh: Mesh,

    pub cell_count: u32,
    pub cell_position: Vec<Vec3>,
    pub cell_height: Vec<f32>,
}

#[derive(Debug)]
struct PlanetGeometryGeneratorComponent{
    pub stages: Vec<Box<dyn GeneratorStage>>,
}

impl PlanetComponent for PlanetGeometryGeneratorComponent{
    fn name(&self) -> &'static str {
        "Planet Goemetry Generator"
    }

    #[instrument]
    fn on_generate(&mut self, props: &mut PropertyMap) -> anyhow::Result<()> {
        if props.has::<PlanetGeometry>(){
            warn!("The planet already contains a PlanetGeometry component and will be overwritten");
        }

        let mut planet_data = PlanetData::default();

        for stage in &self.stages{
            stage.run(&mut planet_data)?;
        }

        Ok(())
    }
}