//! Internal representation of an atom, and of an atom type.
//! For docking simulation the most important information relating to an atom are its coordinates
//! and the "type" of atom is. The type of atom is defined by its Van der Waals
//! radius (for repulsion in case of clashes), the wall depth it has (for attractive strength),
//! partial charges and its role in a hydrogen bond.

use nalgebra::Vector3;

#[derive(PartialEq)]
pub enum HydrogenBond {
    Donor,
    Acceptor,
    Both,
    Neither
}

// Force field parameter data type, obtained from AMBER, OPLS, ...
pub struct AtomType {
    pub vdw_radius: f32, // Van der Waals radius
    pub well_depth: f32,  // Attraction strength to other atoms
    pub partial_charge: f32,
    pub hydrogen_bond: HydrogenBond
}

pub struct Atom {
    pub element: String,
    pub coords: Vector3<f32>,
    pub atom_type: AtomType,
}

impl Atom {
    pub fn new(element: impl Into<String>, coords: Vector3<f32>, atom_type: AtomType) -> Self {
        Self {
            element: element.into(),
            coords,
            atom_type,
        }
    }

    /// Returns a copy of the atom's coordinates.
    pub fn copy_coords(&self) -> Vector3<f32> {
        self.coords
    }
}