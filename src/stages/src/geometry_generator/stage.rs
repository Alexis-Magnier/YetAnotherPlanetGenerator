use super::PlanetData;

pub trait GeneratorStage{
    fn name(&self) -> &'static str;
    fn run(&self, data: &mut PlanetData) -> anyhow::Result<()>;
}

impl std::fmt::Debug for dyn GeneratorStage{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("Geometry generator stage \"{}\"", self.name()))
    }
}