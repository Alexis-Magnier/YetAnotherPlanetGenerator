use std::any::{self, Any, TypeId};
use std::collections::HashMap;

use tracing::instrument;

/// Type-safe dynamic property storage for the planet.
#[derive(Default, Debug)]
pub struct PropertyMap {
    data: HashMap<TypeId, Box<dyn Any>>,
}

impl PropertyMap {
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert or overwrite a property struct of type T.
    pub fn insert<T: 'static>(&mut self, value: T) {
        self.data.insert(TypeId::of::<T>(), Box::new(value));
    }

    /// Read a property struct of type T.
    pub fn get<T: 'static>(&self) -> Option<&T> {
        self.data
            .get(&TypeId::of::<T>())
            .and_then(|boxed| boxed.downcast_ref::<T>())
    }

    /// Mutably access a property struct of type T.
    pub fn get_mut<T: 'static>(&mut self) -> Option<&mut T> {
        self.data
            .get_mut(&TypeId::of::<T>())
            .and_then(|boxed| boxed.downcast_mut::<T>())
    }

    /// Remove a property of type T from the planet.
    pub fn remove<T: 'static>(&mut self) -> Option<T> {
        self.data
            .remove(&TypeId::of::<T>())
            .and_then(|boxed| boxed.downcast::<T>().ok().map(|b| *b))
    }

    pub fn has<T: 'static>(&self) -> bool {
        self.data.contains_key(&TypeId::of::<T>())
    }
}

/// Lifecycle trait implemented by all planet subsystems.
pub trait PlanetComponent {
    fn name(&self) -> &'static str;

    /// Called once when the planet generates.
    fn on_generate(&mut self, _props: &mut PropertyMap) -> anyhow::Result<()> {Ok(())}

    /// Called every simulation tick.
    fn on_tick(&mut self, _props: &mut PropertyMap, _delta_time: f32) -> anyhow::Result<()> {Ok(())}
}

impl std::fmt::Debug for dyn PlanetComponent{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("PlanetComponent \"{}\"", self.name()))
    }
}

/// The core Planet object holding state properties and component behaviors.
#[derive(Debug)]
pub struct Planet {
    pub properties: PropertyMap,
    components: Vec<Box<dyn PlanetComponent>>,
}

impl Planet {
    pub fn new() -> Self {
        Self {
            properties: PropertyMap::new(),
            components: Vec::new(),
        }
    }

    pub fn add_component<C: PlanetComponent + 'static>(&mut self, component: C) {
        self.components.push(Box::new(component));
    }

    #[instrument]
    pub fn generate(&mut self) {
        for component in &mut self.components {
            component.on_generate(&mut self.properties);
        }
    }

    /// Advance simulation by delta_time across all registered components.
    #[instrument]
    pub fn tick(&mut self, delta_time: f32) {
        for component in &mut self.components {
            component.on_tick(&mut self.properties, delta_time);
        }
    }
}