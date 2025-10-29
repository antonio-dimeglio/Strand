//! Internal representation of a ligand.
//! A ligand is a small molecule that binds to a protein receptor.
//! Here, a ligand is represented as a list of atoms, and a list of bonds.
//! Such a list is implemented as an adjacency list where each outer index represents an atom
//! and the list associated with it represents a bond, which is just an index referring to the atom it binds
//! to and a boolean flag stating if the bond is rotational or not.

use nalgebra::Vector3;
use crate::core::atom::Atom;

pub struct Bond {
    pub(crate) target: usize,
    pub(crate) is_rotatable: bool,
}

pub struct Ligand {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Vec<Bond>>,
    pub rotatable_bonds: Vec<(usize, usize)>
}

impl Ligand {
    pub fn new(atoms: Vec<Atom>, bonds: Vec<Vec<Bond>>, rotatable_bonds: Vec<(usize, usize)>) -> Ligand {
        Ligand {
            atoms,
            bonds,
            rotatable_bonds
        }
    }

    pub fn is_bond_rotational(&self, src: usize, dest: usize) -> bool {
        self.bonds[src]
            .iter()
            .find(|b| b.target == dest)
            .map(|b| b.is_rotatable)
            .unwrap_or(false)
    }

    pub fn copy_coordinates(&self) -> Vec<Vector3<f32>> {
        self.atoms
        .iter()
        .map(|v| v.copy_coords()).
        collect()
    }
}
