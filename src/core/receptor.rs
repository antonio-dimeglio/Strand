use nalgebra::Vector3;
use crate::core::atom::Atom;
use crate::core::spatial_grid::SpatialGrid;

pub struct Receptor {
    pub atoms: Vec<Atom>,
    pub grid_size: f32,
    pub grid: SpatialGrid
}

impl Receptor {
    pub fn new(atoms: Vec<Atom>, grid_size: f32) -> Receptor {
        let grid = SpatialGrid::new(&atoms, grid_size);
        Receptor {
            atoms,
            grid_size,
            grid
        }
    }

    pub fn find_nearby_atoms(&self, position: &Vector3<f32>, cutoff: f32) -> Vec<usize> {
        self.grid.query(position, &self.atoms, cutoff)
    }
}

