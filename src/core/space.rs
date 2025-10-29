//! Internal representation of the search space for a docking simulation problem.
//! A search space in docking simulation is defined by three main components:
//! 1.  Spatial boundary: where is the binding site of a protein, and how far can
//!     a ligand move within said space in each direction?
//! 2.  Rotational space: can the ligand rotate along all axes by 360 degrees, or are there
//!     some constraints on how it can move?
//! 3.  Conformational space: some bonds are rotatable in a ligand, if so, which ones, and what 
//!     is the rotational range?
use nalgebra::{Point3, Vector3};

pub struct Boundaries {
    center: Point3<f32>,
    extents: Vector3<f32> // +- dx, +- dy, +- dz
}

pub struct RotationConstraints {
    roll_range: Option<(f32, f32)>,   // None = no constraint on this axis
    pitch_range: Option<(f32, f32)>,
    yaw_range: Option<(f32, f32)>,
}

pub enum RotationalSpace {
    Constrained(RotationConstraints),
    Unconstrained
}

pub struct RotatableBond {
    atom1_idx: usize,
    atom2_idx: usize,
    angle_range: (f32, f32)
}
pub struct ConformationalSpace {
    // Which bonds are rotatable in a ligand?
    rotatable_bonds: Vec<RotatableBond>
}

pub struct Space {
    spatial_boundaries: Boundaries,
    rotational_space: RotationalSpace,
    conformational_space: ConformationalSpace
}