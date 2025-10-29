use crate::server::protein::{Atom, MoleculeData, DockingRequest};

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

/// Test runner for water dimer docking
/// Expected result: Optimal O-O distance of ~2.8 Å with hydrogen bond formed
pub fn run_water_dimer_test() {
    let request = create_water_dimer_request();

    println!("=== Water Dimer Docking Test ===");
    println!("Receptor: {}", request.receptor.name);
    println!("  Atoms: {}", request.receptor.atoms.len());
    println!("Ligand: {}", request.ligand.name);
    println!("  Atoms: {}", request.ligand.atoms.len());
    println!("\nExpected outcome:");
    println!("  - Optimal O-O distance: ~2.8 Å");
    println!("  - Hydrogen bond: O...H distance ~1.8 Å");
    println!("  - Linear O-H...O arrangement");
    println!("================================\n");

    // TODO: Implement docking algorithm
    // Commented out - implement these functions:

    // let result = dock(&request.receptor, &request.ligand, &DockingConfig::default());
    // println!("Docking complete!");
    // println!("  Best score: {}", result.score);
    // println!("  Best pose: {:?}", result.pose);
    // println!("  O-O distance: {:.2} Å", result.metrics.distance);

    // Validate result
    // assert!((result.metrics.distance - 2.8).abs() < 0.5,
    //         "Water dimer O-O distance should be ~2.8 Å");
}
