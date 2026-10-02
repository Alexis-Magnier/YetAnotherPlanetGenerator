
use super::Mesh;
use glam::Vec3;


#[derive(Debug, Clone, Default)]
pub struct Layer{
    pub thickness: u16,
    pub density: f32,
}

#[derive(Debug, Clone, Default)]
pub struct CellLayers{
    pub sedimentary: Layer,
    pub igneous: Layer,
}

#[derive(Debug, Clone, Default)]
pub struct CellData{
    pub plate_id: Vec<u32>,

    pub position: Vec<Vec3>,

    pub stress: Vec<Vec3>,
    pub temperature: Vec<f32>,
    
    /// Distance from plate seed
    pub seed_dist: Vec<f32>,
    
    pub convergent_dist: Vec<f32>,
    pub divergent_dist: Vec<f32>,
    
    pub layers: Vec<CellLayers>,

    pub yeild_strength: Vec<f32>,
    pub height: Vec<f32>,

    pub density: Vec<f32>,
    pub thickness: Vec<f32>,

    pub count: u32,
}

#[derive(Debug, Clone, Default)]
pub struct PlateData{
    pub seeds: Vec<u32>,
    pub velocity: Vec<Vec3>,

    /// Density base noise
    pub density_mean: Vec<f32>,
    pub density_ampl: Vec<f32>,

    /// Strenght base noise
    pub strenght_mean: Vec<f32>,
    pub strenght_ampl: Vec<f32>,

    /// thickness base noise
    pub thickness_mean: Vec<f32>,
    pub thickness_ampl: Vec<f32>,

    /// How much distance from plate seed affects strenght 
    pub strenght_craton_factor: Vec<f32>,
    /// How much distance from divergent collision affect strenght
    pub strenght_thermal_factor: Vec<f32>,

    /// How quickly the temperature will decrease over time
    pub thermal_cooldown: Vec<f32>,
    pub thermal_diffusion_rate: Vec<f32>,

    pub count: u32,
}

#[derive(Debug, Clone, Default)]
pub struct PlanetData {
    pub mesh: Mesh,
    pub cells: CellData,
    pub plates: PlateData,

    pub convergent_cells: Vec<u32>,
    pub divergent_cells: Vec<u32>,
    pub transform_cells: Vec<u32>,
}