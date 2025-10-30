use crate::server::protein::{Atom, MoleculeData, DockingRequest};
use crate::core::atom::{Atom as CoreAtom, AtomType, HydrogenBond};
use crate::core::ligand::{Ligand, Bond};
use crate::core::receptor::Receptor;
use crate::core::space::{Space, Boundaries, RotationalSpace, ConformationalSpace};
use crate::core::pose::StepSizes;
use crate::docking::search::{random_search, mc_simulated_annealing};
use nalgebra::{Point3, Vector3};

/// Creates a water molecule centered at the origin
/// Water molecule (H2O):
/// - 1 Oxygen atom at (0, 0, 0)
/// - 2 Hydrogen atoms forming ~104.5° angle
pub fn create_water_molecule_1() -> MoleculeData {
    MoleculeData {
        name: "Water_1".to_string(),
        molecule_type: "small_molecule".to_string(),
        file_content: None,
        file_format: None,
        atoms: vec![
            // Oxygen atom
            Atom {
                element: "O".to_string(),
                x: 0.0,
                y: 0.0,
                z: 0.0,
                charge: Some(-0.834),  // Typical partial charge for TIP3P water
                atom_type: Some("O.3".to_string()),
                atom_name: Some("O".to_string()),
                residue_name: Some("HOH".to_string()),
                residue_number: Some(1),
                chain_id: Some("A".to_string()),
                serial_number: Some(1),
            },
            // Hydrogen atom 1 (at ~0.957 Å distance, typical O-H bond length)
            Atom {
                element: "H".to_string(),
                x: 0.757,
                y: 0.586,
                z: 0.0,
                charge: Some(0.417),
                atom_type: Some("H".to_string()),
                atom_name: Some("H1".to_string()),
                residue_name: Some("HOH".to_string()),
                residue_number: Some(1),
                chain_id: Some("A".to_string()),
                serial_number: Some(2),
            },
            // Hydrogen atom 2
            Atom {
                element: "H".to_string(),
                x: 0.757,
                y: -0.586,
                z: 0.0,
                charge: Some(0.417),
                atom_type: Some("H".to_string()),
                atom_name: Some("H2".to_string()),
                residue_name: Some("HOH".to_string()),
                residue_number: Some(1),
                chain_id: Some("A".to_string()),
                serial_number: Some(3),
            },
        ],
        bonds: vec![
            [0, 1],  // O-H1 bond
            [0, 2],  // O-H2 bond
        ],
    }
}

/// Creates a second water molecule displaced from the first
/// Positioned to form a hydrogen bond with the first water molecule
/// Expected optimal distance: ~2.8 Å between oxygen atoms
pub fn create_water_molecule_2() -> MoleculeData {
    MoleculeData {
        name: "Water_2".to_string(),
        molecule_type: "small_molecule".to_string(),
        file_content: None,
        file_format: None,
        atoms: vec![
            // Oxygen atom (positioned ~3.5 Å away initially, docking should find ~2.8 Å)
            Atom {
                element: "O".to_string(),
                x: 3.5,
                y: 0.0,
                z: 0.0,
                charge: Some(-0.834),
                atom_type: Some("O.3".to_string()),
                atom_name: Some("O".to_string()),
                residue_name: Some("HOH".to_string()),
                residue_number: Some(2),
                chain_id: Some("B".to_string()),
                serial_number: Some(4),
            },
            // Hydrogen atom 1 (oriented towards first water)
            Atom {
                element: "H".to_string(),
                x: 2.743,
                y: 0.586,
                z: 0.0,
                charge: Some(0.417),
                atom_type: Some("H".to_string()),
                atom_name: Some("H1".to_string()),
                residue_name: Some("HOH".to_string()),
                residue_number: Some(2),
                chain_id: Some("B".to_string()),
                serial_number: Some(5),
            },
            // Hydrogen atom 2
            Atom {
                element: "H".to_string(),
                x: 2.743,
                y: -0.586,
                z: 0.0,
                charge: Some(0.417),
                atom_type: Some("H".to_string()),
                atom_name: Some("H2".to_string()),
                residue_name: Some("HOH".to_string()),
                residue_number: Some(2),
                chain_id: Some("B".to_string()),
                serial_number: Some(6),
            },
        ],
        bonds: vec![
            [0, 1],  // O-H1 bond
            [0, 2],  // O-H2 bond
        ],
    }
}

/// Creates a complete water dimer docking request
pub fn create_water_dimer_request() -> DockingRequest {
    DockingRequest {
        receptor: create_water_molecule_1(),
        ligand: create_water_molecule_2(),
    }
}

/// Convert server Atom to core Atom
fn convert_atom(server_atom: &Atom) -> CoreAtom {
    let atom_type = AtomType {
        vdw_radius: if server_atom.element == "O" { 1.52 } else { 1.2 },
        well_depth: if server_atom.element == "O" { 0.21 } else { 0.02 },
        partial_charge: server_atom.charge.unwrap_or(0.0) as f32,
        hydrogen_bond: if server_atom.element == "O" {
            HydrogenBond::Acceptor
        } else {
            HydrogenBond::Donor
        },
    };

    CoreAtom::new(
        &server_atom.element,
        Vector3::new(server_atom.x as f32, server_atom.y as f32, server_atom.z as f32),
        atom_type,
    )
}

/// Create receptor from water molecule
fn create_receptor() -> Receptor {
    let water = create_water_molecule_1();
    let atoms: Vec<CoreAtom> = water.atoms.iter().map(|a| convert_atom(a)).collect();
    Receptor::new(atoms, 5.0)
}

/// Create ligand from water molecule
fn create_ligand() -> Ligand {
    let water = create_water_molecule_2();
    let atoms: Vec<CoreAtom> = water.atoms.iter().map(|a| convert_atom(a)).collect();

    // Water is rigid - no rotatable bonds
    let bonds = vec![
        vec![Bond { target: 1, is_rotatable: false }, Bond { target: 2, is_rotatable: false }],
        vec![Bond { target: 0, is_rotatable: false }],
        vec![Bond { target: 0, is_rotatable: false }],
    ];

    Ligand::new(atoms, bonds, vec![])
}

/// Create search space for water dimer
fn create_search_space() -> Space {
    // Allow translation within 5 Å cube around origin
    let boundaries = Boundaries {
        center: Point3::new(0.0, 0.0, 0.0),
        extents: Vector3::new(5.0, 5.0, 5.0),
    };

    // Allow full rotation
    let rotational_space = RotationalSpace::Unconstrained;

    // No conformational flexibility (rigid water)
    let conformational_space = ConformationalSpace {
        rotatable_bonds: vec![],
    };

    Space {
        spatial_boundaries: boundaries,
        rotational_space,
        conformational_space,
    }
}

/// CLI runner function (deprecated - use cargo test instead)
pub fn run_water_dimer_test() {
    println!("=== Water Dimer Docking Test ===");
    println!("Note: This CLI runner is deprecated. Use 'cargo test water_dimer' instead.");
    println!("\nRunning tests programmatically...\n");

    // test_water_dimer_random_search();
    println!();
    // test_water_dimer_monte_carlo();
}

#[test]
fn test_water_dimer_random_search() {
    println!("\n=== Water Dimer Random Search Test ===");

    let receptor = create_receptor();
    let ligand = create_ligand();
    let space = create_search_space();

    println!("Running random search with 100 samples...");
    let (_best_pose, best_score) = random_search(100, &space, &ligand, &receptor);

    println!("\nResults:");
    println!("  Total Energy: {:.3} kcal/mol", best_score.total);
    println!("  VDW:          {:.3} kcal/mol", best_score.vdw);
    println!("  Electrostatic: {:.3} kcal/mol", best_score.esi);
    println!("  H-bond:       {:.3} kcal/mol", best_score.h_bond);
    println!("  Torsional:    {:.3} kcal/mol", best_score.torsional);

    // Validate that we found a reasonable energy
    assert!(best_score.total < 0.0, "Best score should be negative (favorable binding)");
    println!("\n✓ Random search completed successfully");
}

#[test]
fn test_water_dimer_monte_carlo() {
    println!("\n=== Water Dimer Monte Carlo Test ===");

    let receptor = create_receptor();
    let ligand = create_ligand();
    let space = create_search_space();
    let step_sizes = StepSizes::default();

    println!("Running Monte Carlo with 1000 iterations...");
    println!("  Temperature: 1.0 kcal/mol");
    println!("  Step sizes: trans={:.1}Å, rot={:.0}°, tors={:.0}°",
             step_sizes.translation, step_sizes.rotation_degrees, step_sizes.torsion_degrees);

    let (_best_pose, best_score) = mc_simulated_annealing(
        1000,
        &space,
        &ligand,
        &receptor,
        &step_sizes,
        1.0,
    );

    println!("\nResults:");
    println!("  Total Energy: {:.3} kcal/mol", best_score.total);
    println!("  VDW:          {:.3} kcal/mol", best_score.vdw);
    println!("  Electrostatic: {:.3} kcal/mol", best_score.esi);
    println!("  H-bond:       {:.3} kcal/mol", best_score.h_bond);
    println!("  Torsional:    {:.3} kcal/mol", best_score.torsional);

    // Validate that we found a favorable binding energy
    assert!(best_score.total < 0.0, "Best score should be negative (favorable binding)");

    // For water dimer, expect significant hydrogen bonding contribution
    assert!(best_score.h_bond < -2.0, "Should have hydrogen bond contribution");

    println!("\n✓ Monte Carlo search completed successfully");
}
