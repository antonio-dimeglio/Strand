//! Internal representation of a pose in docking simulation.
//! A pose is a specific configuration of the ligand relative to the
//! receptor at one moment in the search space.
//! Said configuration contains:
//!     1. Translation: the 3d position of the ligand's center offset to the protein
//!     2. Rotation: Quaternion representing the angular orientation of the ligand w.r.t the protein orientation
//!     3. Conformation: The torsion angles for rotatable bonds within the ligand (internal flexibility)
//!     4. Binding energy: the energy of the binding, lower is better.
//!     5. Rmsd: An optional value that can be used to compare to literature results.
//! Such a struct makes it easier to separate the actual space searching task - which consists of only
//! generating/iterating over coordinates/angles -  vs the actual realization of said
//! transformations. I.e. instead of actually storing all the rotations of a ligand that the search
//! space allows, we only store the mathematical parameters that define how said transformation is performed
//! which is much more efficient.

use std::collections::HashSet;
use nalgebra::{DVector, Rotation3, Translation3, UnitQuaternion, UnitVector3, Vector3};
use rand::{thread_rng, Rng};
use rand_distr::{Distribution, Normal};
use crate::core::space::RotationalSpace;
use crate::core::ligand::{Ligand};
use crate::core::space::{Boundaries, ConformationalSpace, Space};

type Coordinates = Vec<Vector3<f32>>;

/// Step sizes for pose perturbation in Monte Carlo search
pub struct StepSizes {
    pub translation: f32,        // Step size in Angstroms
    pub rotation_degrees: f32,   // Step size in degrees
    pub torsion_degrees: f32,    // Step size in degrees for torsion angles
}

impl StepSizes {
    /// Default step sizes optimized for typical docking scenarios
    pub fn default() -> Self {
        StepSizes {
            translation: 1.0,      // 1.0 Å
            rotation_degrees: 20.0, // 20°
            torsion_degrees: 20.0,  // 20°
        }
    }

    /// Create custom step sizes
    pub fn new(translation: f32, rotation_degrees: f32, torsion_degrees: f32) -> Self {
        StepSizes {
            translation,
            rotation_degrees,
            torsion_degrees,
        }
    }
}

#[derive(Clone)]
pub struct Pose {
    translation: Translation3<f32>,
    rotation: UnitQuaternion<f32>,
    conformation: DVector<f32>,
    energy: f32,  // binding energy
    rmsd: Option<f32>,  // RMSD from reference
}

impl Pose {
    fn sample_angle(rng: &mut impl Rng, range: Option<(f32, f32)>) -> f32 {
        match range {
            Some((min, max)) => rng.gen_range(min..max),
            None => rng.gen_range(0.0..2.0 * std::f32::consts::PI),
        }
    }

    pub fn generate_random_translation(boundaries: &Boundaries) -> Translation3<f32> {
        let center = boundaries.center;
        let extents = boundaries.extents;
        let mut rng = thread_rng();

        let random_x = center.x + rng.gen_range(-extents.x..extents.x);
        let random_y = center.y + rng.gen_range(-extents.y..extents.y);
        let random_z = center.z + rng.gen_range(-extents.z..extents.z);

        Translation3::new(random_x, random_y, random_z)
    }

    pub fn generate_random_rotation(rotation: &RotationalSpace) -> UnitQuaternion<f32>{
        let mut rng = thread_rng();

        match rotation {
            RotationalSpace::Unconstrained => {
                // Sample 4D unit sphere via Gaussian method
                let q = nalgebra::Quaternion::new(
                    rng.gen_range(-1.0..1.0),
                    rng.gen_range(-1.0..1.0),
                    rng.gen_range(-1.0..1.0),
                    rng.gen_range(-1.0..1.0),

                );
                UnitQuaternion::from_quaternion(q)
            }
            RotationalSpace::Constrained(constraints) => {
                let roll = Self::sample_angle(&mut rng, constraints.roll_range);
                let pitch = Self::sample_angle(&mut rng, constraints.pitch_range);
                let yaw = Self::sample_angle(&mut rng, constraints.yaw_range);

                UnitQuaternion::from_euler_angles(roll, pitch, yaw)
            }
        }
    }

    pub fn generate_random_conformation(conformational_space: &ConformationalSpace) -> DVector<f32> {
        let mut conformation = DVector::zeros(conformational_space.rotatable_bonds.len());
        let mut rng = thread_rng();

        for (i, bond) in conformational_space.rotatable_bonds.iter().enumerate() {
            let min_angle = bond.angle_range.0;
            let max_angle = bond.angle_range.1;
            conformation[i] = rng.gen_range(min_angle..max_angle);
        }
        return conformation
    }
    pub fn random_pose(space: &Space, _ligand: &Ligand) -> Pose {
        let translation = Pose::generate_random_translation(&space.spatial_boundaries);
        let rotation = Pose::generate_random_rotation(&space.rotational_space);
        let conformation = Pose::generate_random_conformation(&space.conformational_space);

        Pose {
            translation,
            rotation,
            conformation,
            energy: 0.0,
            rmsd: None
        }
    }

    /// Generate a random unit vector uniformly distributed on the unit sphere
    /// Uses normal distribution to ensure uniform sampling
    fn random_unit_vector() -> UnitVector3<f32> {
        let mut rng = thread_rng();
        let normal = Normal::new(0.0, 1.0).unwrap();

        let x = normal.sample(&mut rng);
        let y = normal.sample(&mut rng);
        let z = normal.sample(&mut rng);

        UnitVector3::new_normalize(Vector3::new(x, y, z))
    }

    /// Perturb translation by a small random displacement
    fn perturb_translation(
        current: &Translation3<f32>,
        boundaries: &Boundaries,
        step_size: f32
    ) -> Translation3<f32> {
        let mut rng = thread_rng();

        // Random displacement in each dimension
        let dx = rng.gen_range(-step_size..step_size);
        let dy = rng.gen_range(-step_size..step_size);
        let dz = rng.gen_range(-step_size..step_size);

        let new_x = current.vector.x + dx;
        let new_y = current.vector.y + dy;
        let new_z = current.vector.z + dz;

        // Clamp to boundaries
        let center = boundaries.center;
        let extents = boundaries.extents;

        let clamped_x = new_x.clamp(center.x - extents.x, center.x + extents.x);
        let clamped_y = new_y.clamp(center.y - extents.y, center.y + extents.y);
        let clamped_z = new_z.clamp(center.z - extents.z, center.z + extents.z);

        Translation3::new(clamped_x, clamped_y, clamped_z)
    }

    /// Perturb rotation by applying a small random rotation
    fn perturb_rotation(
        current: &UnitQuaternion<f32>,
        max_angle_degrees: f32
    ) -> UnitQuaternion<f32> {
        let mut rng = thread_rng();

        // Random axis
        let random_axis = Self::random_unit_vector();

        // Random angle in range [-max_angle, +max_angle]
        let angle_degrees = rng.gen_range(-max_angle_degrees..max_angle_degrees);
        let angle_radians = angle_degrees.to_radians();

        // Create small rotation
        let small_rotation = UnitQuaternion::from_axis_angle(&random_axis, angle_radians);

        // Compose with current rotation
        small_rotation * current
    }

    /// Perturb conformation by adjusting one or more torsion angles
    fn perturb_conformation(
        current: &DVector<f32>,
        conformational_space: &ConformationalSpace,
        max_angle_degrees: f32
    ) -> DVector<f32> {
        let mut rng = thread_rng();
        let mut new_conformation = current.clone();

        if conformational_space.rotatable_bonds.is_empty() {
            return new_conformation;
        }

        // Randomly select one torsion angle to perturb
        let bond_idx = rng.gen_range(0..conformational_space.rotatable_bonds.len());
        let bond = &conformational_space.rotatable_bonds[bond_idx];

        // Random angle change
        let angle_delta = rng.gen_range(-max_angle_degrees..max_angle_degrees);
        let angle_delta_radians = angle_delta.to_radians();

        // Apply perturbation
        new_conformation[bond_idx] += angle_delta_radians;

        // Clamp to bond's allowed range
        let min_angle = bond.angle_range.0;
        let max_angle = bond.angle_range.1;
        new_conformation[bond_idx] = new_conformation[bond_idx].clamp(min_angle, max_angle);

        new_conformation
    }

    /// Perturb the current pose by making small random changes to translation, rotation, and conformation
    /// Returns a new perturbed pose
    pub fn perturb(&self, space: &Space, step_sizes: &StepSizes) -> Pose {
        let new_translation = Self::perturb_translation(
            &self.translation,
            &space.spatial_boundaries,
            step_sizes.translation
        );

        let new_rotation = Self::perturb_rotation(
            &self.rotation,
            step_sizes.rotation_degrees
        );

        let new_conformation = Self::perturb_conformation(
            &self.conformation,
            &space.conformational_space,
            step_sizes.torsion_degrees
        );

        Pose {
            translation: new_translation,
            rotation: new_rotation,
            conformation: new_conformation,
            energy: 0.0,
            rmsd: None
        }
    }

    /// DFS algorithm to explore the ligand molecule and find atoms that are effected downstream
    /// by rotation of a bond
    fn find_atoms_downstream(&self, ligand: &Ligand, bond: (usize, usize)) -> Vec<usize> {
        let mut visited: HashSet<usize> = HashSet::new();
        let mut queue = Vec::new();
        let mut moving: Vec<usize> = Vec::new();

        visited.insert(bond.0);
        visited.insert(bond.1);

        for neighbor in &ligand.bonds[bond.1] {
            queue.push(neighbor.target);
        }

        while !queue.is_empty() {
            let curr = queue.pop().unwrap();

            if curr == bond.0 {
                continue;
            }

            if !visited.contains(&curr) {
                visited.insert(curr);
                moving.push(curr);

                for neighbor in &ligand.bonds[curr] {
                    queue.push(neighbor.target);
                }
            }
        }

        return moving
    }


    /// Applying a conformational translation is the first step to apply a pose onto a ligand.
    /// To do so, first the rotatable bonds of the ligand are fetched, alongside with the atoms that,
    /// downstream, are rotated if the source atom is rotated.
    /// Then, all are rotated on the correct rotational axis
    fn apply_conformational_translation(&self, ligand: &Ligand, coords: &mut Coordinates) {
        for i in 0..self.conformation.len() {
            let (atom1, atom2) = ligand.rotatable_bonds[i];
            let angle = self.conformation[i];
            let moving_downstream = self.find_atoms_downstream(ligand, (atom1, atom2));
            let axis = UnitVector3::new_normalize(coords[atom2] - coords[atom1]);
            let pivot = coords[atom1];

            for atom_idx in moving_downstream {
                let pos = coords[atom_idx] - pivot;
                let rotation = Rotation3::from_axis_angle(&axis, angle);
                let rotated_pos = rotation * pos;
                coords[atom_idx] = rotated_pos + pivot;
            }
        }
    }

    fn apply_global_rotation(&self, ligand: &Ligand, coords: &mut Coordinates) {
        for i in 0..ligand.atoms.len() {
            coords[i] = self.rotation * coords[i];
        }
    }

    fn apply_translation(&self, ligand: &Ligand, coords: &mut Coordinates) {
        for i in 0..ligand.atoms.len() {
            coords[i] += self.translation.vector;
        }
    }

    // Applies a pose transformation to all ligand's atoms
    pub fn generate_coordinates(&self, ligand: &Ligand) -> Vec<Vector3<f32>> {
        let mut coords = ligand.copy_coordinates();

        self.apply_conformational_translation(ligand, &mut coords);
        self.apply_global_rotation(ligand, &mut coords);
        self.apply_translation(ligand, &mut coords);

        return coords;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Translation3, UnitQuaternion, Vector3, DVector};
    use crate::core::ligand::{Ligand, Bond};
    use crate::core::atom::{Atom, AtomType, HydrogenBond};
    use std::f32::consts::PI;

    /// Helper to create a dummy AtomType for testing (values don't matter for coordinate generation)
    fn dummy_atom_type() -> AtomType {
        AtomType {
            vdw_radius: 1.5,
            well_depth: 0.1,
            partial_charge: 0.0,
            hydrogen_bond: HydrogenBond::Neither,
        }
    }

    /// Helper to create a test atom at specific coordinates
    fn make_test_atom(x: f32, y: f32, z: f32) -> Atom {
        Atom::new("C", Vector3::new(x, y, z), dummy_atom_type())
    }

    /// Helper to check if two coordinates are approximately equal (within tolerance)
    fn coords_approx_equal(a: &Vector3<f32>, b: &Vector3<f32>, tolerance: f32) -> bool {
        (a - b).norm() < tolerance
    }

    #[test]
    fn test_pure_translation() {
        // Create a simple 3-atom linear ligand: A(0,0,0), B(1,0,0), C(2,0,0)
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),
            make_test_atom(1.0, 0.0, 0.0),
            make_test_atom(2.0, 0.0, 0.0),
        ];

        // No bonds needed for this test
        let bonds = vec![vec![], vec![], vec![]];
        let rotatable_bonds = vec![];

        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        // Create pose with only translation (5, 3, 2)
        let pose = Pose {
            translation: Translation3::new(5.0, 3.0, 2.0),
            rotation: UnitQuaternion::identity(), // No rotation
            conformation: DVector::from_vec(vec![]), // No torsions
            energy: 0.0,
            rmsd: None,
        };

        // Generate coordinates
        let result = pose.generate_coordinates(&ligand);

        // Expected results: all atoms shifted by (5, 3, 2)
        let expected = vec![
            Vector3::new(5.0, 3.0, 2.0),
            Vector3::new(6.0, 3.0, 2.0),
            Vector3::new(7.0, 3.0, 2.0),
        ];

        // Verify
        assert_eq!(result.len(), expected.len());
        for i in 0..result.len() {
            assert!(
                coords_approx_equal(&result[i], &expected[i], 0.001),
                "Atom {} mismatch: got {:?}, expected {:?}",
                i, result[i], expected[i]
            );
        }
    }

    #[test]
    fn test_pure_rotation() {
        // Create a simple 3-atom linear ligand along X-axis
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),
            make_test_atom(1.0, 0.0, 0.0),
            make_test_atom(2.0, 0.0, 0.0),
        ];

        let bonds = vec![vec![], vec![], vec![]];
        let rotatable_bonds = vec![];

        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        // Create pose with 90° rotation around Z-axis
        let pose = Pose {
            translation: Translation3::new(0.0, 0.0, 0.0), // No translation
            rotation: UnitQuaternion::from_axis_angle(&Vector3::z_axis(), PI / 2.0), // 90° around Z
            conformation: DVector::from_vec(vec![]),
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // After 90° rotation around Z, X becomes Y
        // A(0,0,0) stays at origin
        // B(1,0,0) -> (0,1,0)
        // C(2,0,0) -> (0,2,0)
        let expected = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 2.0, 0.0),
        ];

        assert_eq!(result.len(), expected.len());
        for i in 0..result.len() {
            assert!(
                coords_approx_equal(&result[i], &expected[i], 0.001),
                "Atom {} mismatch: got {:?}, expected {:?}",
                i, result[i], expected[i]
            );
        }
    }

    #[test]
    fn test_single_torsion() {
        // Create a 3-atom bent chain where C is OFF the A-B axis
        // A at origin, B along X-axis, C offset in Y so we can see rotation
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),  // A
            make_test_atom(1.0, 0.0, 0.0),  // B
            make_test_atom(2.0, 1.0, 0.0),  // C - offset in Y direction
        ];

        // Adjacency list: A-B, B-A and B-C, C-B
        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: true }],  // A connects to B
            vec![
                Bond { target: 0, is_rotatable: true },   // B connects to A
                Bond { target: 2, is_rotatable: false },  // B connects to C
            ],
            vec![Bond { target: 1, is_rotatable: false }], // C connects to B
        ];

        // A-B is rotatable (atoms 0-1)
        let rotatable_bonds = vec![(0, 1)];

        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        // Create pose with 90° torsion around A-B bond (X-axis)
        let pose = Pose {
            translation: Translation3::new(0.0, 0.0, 0.0),
            rotation: UnitQuaternion::identity(),
            conformation: DVector::from_vec(vec![PI / 2.0]), // 90° torsion
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // After rotation around A-B (X-axis) by 90°:
        // A(0,0,0) - fixed (upstream of bond)
        // B(1,0,0) - fixed (pivot point)
        // C(2,1,0) - relative to B: (1,1,0)
        //          - rotate (1,1,0) by 90° around X-axis
        //          - Y and Z components rotate: (1,0,1)
        //          - absolute position: B + (1,0,1) = (2,0,1)
        let expected = vec![
            Vector3::new(0.0, 0.0, 0.0), // A unchanged
            Vector3::new(1.0, 0.0, 0.0), // B unchanged (pivot)
            Vector3::new(2.0, 0.0, 1.0), // C rotated: was (2,1,0), now (2,0,1)
        ];

        assert_eq!(result.len(), expected.len());
        for i in 0..result.len() {
            assert!(
                coords_approx_equal(&result[i], &expected[i], 0.001),
                "Atom {} mismatch: got {:?}, expected {:?}",
                i, result[i], expected[i]
            );
        }
    }

    #[test]
    fn test_combined_transformation() {
        // Test all three transformations together
        // Use bent geometry so torsion is visible
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),  // A
            make_test_atom(1.0, 0.0, 0.0),  // B
            make_test_atom(2.0, 1.0, 0.0),  // C - offset in Y
        ];

        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: true }],
            vec![
                Bond { target: 0, is_rotatable: true },
                Bond { target: 2, is_rotatable: false },
            ],
            vec![Bond { target: 1, is_rotatable: false }],
        ];

        let rotatable_bonds = vec![(0, 1)];
        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        // Pose with torsion, rotation, and translation
        let pose = Pose {
            translation: Translation3::new(5.0, 0.0, 0.0),
            rotation: UnitQuaternion::from_axis_angle(&Vector3::z_axis(), PI / 4.0), // 45° around Z
            conformation: DVector::from_vec(vec![PI / 2.0]), // 90° torsion
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // Step-by-step calculation:
        // After torsion (90° around X-axis from A to B):
        //   A: (0,0,0) - fixed
        //   B: (1,0,0) - pivot
        //   C: (2,1,0) -> rotate (1,1,0) rel to B by 90° -> (1,0,1) -> abs (2,0,1)
        //
        // After 45° rotation around Z-axis:
        //   A: (0,0,0) -> (0,0,0)
        //   B: (1,0,0) -> (cos45°, sin45°, 0) = (0.707, 0.707, 0)
        //   C: (2,0,1) -> (2*cos45°, 2*sin45°, 1) = (1.414, 1.414, 1)
        //
        // After translation by (5,0,0):
        //   A: (5, 0, 0)
        //   B: (5.707, 0.707, 0)
        //   C: (6.414, 1.414, 1)

        let sqrt2_2 = (2.0_f32).sqrt() / 2.0;
        let sqrt2 = (2.0_f32).sqrt();
        let expected = vec![
            Vector3::new(5.0, 0.0, 0.0),
            Vector3::new(5.0 + sqrt2_2, sqrt2_2, 0.0),
            Vector3::new(5.0 + sqrt2, sqrt2, 1.0),
        ];

        assert_eq!(result.len(), expected.len());
        for i in 0..result.len() {
            assert!(
                coords_approx_equal(&result[i], &expected[i], 0.001),
                "Atom {} mismatch: got {:?}, expected {:?}",
                i, result[i], expected[i]
            );
        }
    }

    #[test]
    fn test_water_molecule_rigid() {
        // Test with a realistic water molecule (no rotatable bonds)
        // O at origin, H atoms forming ~104.5° angle
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),           // O
            make_test_atom(0.757, 0.586, 0.0),       // H1
            make_test_atom(0.757, -0.586, 0.0),      // H2
        ];

        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: false }, Bond { target: 2, is_rotatable: false }],
            vec![Bond { target: 0, is_rotatable: false }],
            vec![Bond { target: 0, is_rotatable: false }],
        ];

        let rotatable_bonds = vec![]; // Water is rigid
        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        // Move water to (2.8, 0, 0)
        let pose = Pose {
            translation: Translation3::new(2.8, 0.0, 0.0),
            rotation: UnitQuaternion::identity(),
            conformation: DVector::from_vec(vec![]),
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        let expected = vec![
            Vector3::new(2.8, 0.0, 0.0),
            Vector3::new(3.557, 0.586, 0.0),
            Vector3::new(3.557, -0.586, 0.0),
        ];

        assert_eq!(result.len(), expected.len());

        // Verify all atoms shifted correctly
        for i in 0..result.len() {
            assert!(
                coords_approx_equal(&result[i], &expected[i], 0.001),
                "Atom {} mismatch: got {:?}, expected {:?}",
                i, result[i], expected[i]
            );
        }

        // Also verify O-H bond lengths are preserved (should be ~0.96 Å)
        let oh1_original = ((0.757_f32 * 0.757) + (0.586 * 0.586)).sqrt();
        let oh1_after = (result[1] - result[0]).norm();
        assert!((oh1_original - oh1_after).abs() < 0.001, "O-H1 bond length changed");

        let oh2_original = ((0.757_f32 * 0.757) + (0.586 * 0.586)).sqrt();
        let oh2_after = (result[2] - result[0]).norm();
        assert!((oh2_original - oh2_after).abs() < 0.001, "O-H2 bond length changed");
    }

    #[test]
    fn test_multiple_rotatable_bonds() {
        // Test a chain with TWO rotatable bonds: A-B-C-D
        // Both B-C rotations should be applied sequentially
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),  // A
            make_test_atom(1.0, 0.0, 0.0),  // B
            make_test_atom(2.0, 0.0, 0.0),  // C
            make_test_atom(3.0, 1.0, 0.0),  // D - offset in Y
        ];

        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: true }],   // A-B
            vec![
                Bond { target: 0, is_rotatable: true },    // B-A
                Bond { target: 2, is_rotatable: true },    // B-C rotatable
            ],
            vec![
                Bond { target: 1, is_rotatable: true },    // C-B
                Bond { target: 3, is_rotatable: false },   // C-D
            ],
            vec![Bond { target: 2, is_rotatable: false }],  // D-C
        ];

        // Two rotatable bonds: (0,1) and (1,2)
        let rotatable_bonds = vec![(0, 1), (1, 2)];
        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        // Rotate first bond by 90°, second bond by 90°
        let pose = Pose {
            translation: Translation3::new(0.0, 0.0, 0.0),
            rotation: UnitQuaternion::identity(),
            conformation: DVector::from_vec(vec![PI / 2.0, PI / 2.0]),
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // First rotation (A-B): affects C and D
        // Second rotation (B-C): affects only D
        // Expected after first: C and D rotate around X-axis
        // Expected after second: D rotates around new B-C axis

        // A and B don't move
        assert!(coords_approx_equal(&result[0], &Vector3::new(0.0, 0.0, 0.0), 0.001));
        assert!(coords_approx_equal(&result[1], &Vector3::new(1.0, 0.0, 0.0), 0.001));

        // C should stay on X-axis (it's on both rotation axes)
        assert!(coords_approx_equal(&result[2], &Vector3::new(2.0, 0.0, 0.0), 0.001));

        // D should have moved due to both rotations
        // Complex to calculate exactly, but it should NOT be at original position
        assert!(!coords_approx_equal(&result[3], &Vector3::new(3.0, 1.0, 0.0), 0.1));
    }

    #[test]
    fn test_negative_rotation_angle() {
        // Test rotation in opposite direction (negative angle)
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),
            make_test_atom(1.0, 0.0, 0.0),
            make_test_atom(2.0, 1.0, 0.0),
        ];

        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: true }],
            vec![
                Bond { target: 0, is_rotatable: true },
                Bond { target: 2, is_rotatable: false },
            ],
            vec![Bond { target: 1, is_rotatable: false }],
        ];

        let rotatable_bonds = vec![(0, 1)];
        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        // -90° rotation (opposite of previous test)
        let pose = Pose {
            translation: Translation3::new(0.0, 0.0, 0.0),
            rotation: UnitQuaternion::identity(),
            conformation: DVector::from_vec(vec![-PI / 2.0]),
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // After -90° rotation around X-axis:
        // C should rotate from (2,1,0) to (2,0,-1) instead of (2,0,1)
        let expected = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 0.0, -1.0),  // Negative Z
        ];

        for i in 0..result.len() {
            assert!(
                coords_approx_equal(&result[i], &expected[i], 0.001),
                "Atom {} mismatch: got {:?}, expected {:?}",
                i, result[i], expected[i]
            );
        }
    }

    #[test]
    fn test_180_degree_rotation() {
        // Test 180° rotation (flips atom to opposite side)
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),
            make_test_atom(1.0, 0.0, 0.0),
            make_test_atom(2.0, 1.0, 0.0),
        ];

        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: true }],
            vec![
                Bond { target: 0, is_rotatable: true },
                Bond { target: 2, is_rotatable: false },
            ],
            vec![Bond { target: 1, is_rotatable: false }],
        ];

        let rotatable_bonds = vec![(0, 1)];
        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        let pose = Pose {
            translation: Translation3::new(0.0, 0.0, 0.0),
            rotation: UnitQuaternion::identity(),
            conformation: DVector::from_vec(vec![PI]),  // 180°
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // After 180° rotation: (2,1,0) -> (2,-1,0)
        let expected = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, -1.0, 0.0),  // Flipped Y
        ];

        for i in 0..result.len() {
            assert!(
                coords_approx_equal(&result[i], &expected[i], 0.001),
                "Atom {} mismatch: got {:?}, expected {:?}",
                i, result[i], expected[i]
            );
        }
    }

    #[test]
    fn test_360_degree_rotation() {
        // Test full 360° rotation (should return to original position)
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),
            make_test_atom(1.0, 0.0, 0.0),
            make_test_atom(2.0, 1.0, 0.0),
        ];

        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: true }],
            vec![
                Bond { target: 0, is_rotatable: true },
                Bond { target: 2, is_rotatable: false },
            ],
            vec![Bond { target: 1, is_rotatable: false }],
        ];

        let rotatable_bonds = vec![(0, 1)];
        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        let pose = Pose {
            translation: Translation3::new(0.0, 0.0, 0.0),
            rotation: UnitQuaternion::identity(),
            conformation: DVector::from_vec(vec![2.0 * PI]),  // 360°
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // Should be back to original positions
        let expected = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 1.0, 0.0),  // Back to original
        ];

        for i in 0..result.len() {
            assert!(
                coords_approx_equal(&result[i], &expected[i], 0.01),  // Slightly larger tolerance for accumulated error
                "Atom {} mismatch: got {:?}, expected {:?}",
                i, result[i], expected[i]
            );
        }
    }

    #[test]
    fn test_branched_molecule() {
        // Test Y-shaped molecule where rotation affects one branch but not another
        //     D
        //     |
        // A-B-C
        //     |
        //     E
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),   // A
            make_test_atom(1.0, 0.0, 0.0),   // B
            make_test_atom(2.0, 0.0, 0.0),   // C (branch point)
            make_test_atom(2.0, 1.0, 0.0),   // D (branch 1)
            make_test_atom(2.0, -1.0, 0.0),  // E (branch 2)
        ];

        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: true }],   // A-B
            vec![
                Bond { target: 0, is_rotatable: true },    // B-A
                Bond { target: 2, is_rotatable: false },   // B-C
            ],
            vec![
                Bond { target: 1, is_rotatable: false },   // C-B
                Bond { target: 3, is_rotatable: false },   // C-D
                Bond { target: 4, is_rotatable: false },   // C-E
            ],
            vec![Bond { target: 2, is_rotatable: false }],  // D-C
            vec![Bond { target: 2, is_rotatable: false }],  // E-C
        ];

        // Only A-B is rotatable
        let rotatable_bonds = vec![(0, 1)];
        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        // 90° rotation around A-B (X-axis)
        let pose = Pose {
            translation: Translation3::new(0.0, 0.0, 0.0),
            rotation: UnitQuaternion::identity(),
            conformation: DVector::from_vec(vec![PI / 2.0]),
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // A and B unchanged
        assert!(coords_approx_equal(&result[0], &Vector3::new(0.0, 0.0, 0.0), 0.001));
        assert!(coords_approx_equal(&result[1], &Vector3::new(1.0, 0.0, 0.0), 0.001));

        // C stays on X-axis
        assert!(coords_approx_equal(&result[2], &Vector3::new(2.0, 0.0, 0.0), 0.001));

        // D and E both rotate 90° around X-axis
        // D: (2,1,0) -> (2,0,1)
        // E: (2,-1,0) -> (2,0,-1)
        assert!(coords_approx_equal(&result[3], &Vector3::new(2.0, 0.0, 1.0), 0.001));
        assert!(coords_approx_equal(&result[4], &Vector3::new(2.0, 0.0, -1.0), 0.001));
    }

    #[test]
    fn test_zero_rotation() {
        // Edge case: rotation angle is 0 (should not change anything)
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),
            make_test_atom(1.0, 0.0, 0.0),
            make_test_atom(2.0, 1.0, 0.0),
        ];

        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: true }],
            vec![
                Bond { target: 0, is_rotatable: true },
                Bond { target: 2, is_rotatable: false },
            ],
            vec![Bond { target: 1, is_rotatable: false }],
        ];

        let rotatable_bonds = vec![(0, 1)];
        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        let pose = Pose {
            translation: Translation3::new(0.0, 0.0, 0.0),
            rotation: UnitQuaternion::identity(),
            conformation: DVector::from_vec(vec![0.0]),  // 0° rotation
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // Everything should stay in original position
        let expected = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(2.0, 1.0, 0.0),
        ];

        for i in 0..result.len() {
            assert!(
                coords_approx_equal(&result[i], &expected[i], 0.001),
                "Atom {} mismatch: got {:?}, expected {:?}",
                i, result[i], expected[i]
            );
        }
    }

    #[test]
    fn test_longer_chain() {
        // Test a longer chain with 6 atoms to catch indexing bugs
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),
            make_test_atom(1.0, 0.0, 0.0),
            make_test_atom(2.0, 0.0, 0.0),
            make_test_atom(3.0, 1.0, 0.0),
            make_test_atom(4.0, 1.0, 0.0),
            make_test_atom(5.0, 0.0, 0.0),
        ];

        let bonds = vec![
            vec![Bond { target: 1, is_rotatable: true }],
            vec![
                Bond { target: 0, is_rotatable: true },
                Bond { target: 2, is_rotatable: true },
            ],
            vec![
                Bond { target: 1, is_rotatable: true },
                Bond { target: 3, is_rotatable: false },
            ],
            vec![
                Bond { target: 2, is_rotatable: false },
                Bond { target: 4, is_rotatable: false },
            ],
            vec![
                Bond { target: 3, is_rotatable: false },
                Bond { target: 5, is_rotatable: false },
            ],
            vec![Bond { target: 4, is_rotatable: false }],
        ];

        let rotatable_bonds = vec![(0, 1), (1, 2)];
        let ligand = Ligand::new(atoms, bonds, rotatable_bonds);

        let pose = Pose {
            translation: Translation3::new(0.0, 0.0, 0.0),
            rotation: UnitQuaternion::identity(),
            conformation: DVector::from_vec(vec![PI / 4.0, -PI / 4.0]),  // 45° and -45°
            energy: 0.0,
            rmsd: None,
        };

        let result = pose.generate_coordinates(&ligand);

        // Just verify we got 6 atoms back and none are NaN
        assert_eq!(result.len(), 6);
        for i in 0..result.len() {
            assert!(!result[i].x.is_nan(), "Atom {} X is NaN", i);
            assert!(!result[i].y.is_nan(), "Atom {} Y is NaN", i);
            assert!(!result[i].z.is_nan(), "Atom {} Z is NaN", i);
        }

        // A and B should be unchanged
        assert!(coords_approx_equal(&result[0], &Vector3::new(0.0, 0.0, 0.0), 0.001));
        assert!(coords_approx_equal(&result[1], &Vector3::new(1.0, 0.0, 0.0), 0.001));
    }
}

