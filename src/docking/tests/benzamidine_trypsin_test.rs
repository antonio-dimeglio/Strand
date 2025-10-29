use crate::server::protein::{Atom, MoleculeData, DockingRequest};

/// Creates benzamidine ligand (C7H9N2+)
/// Benzamidine is a small competitive inhibitor of trypsin
/// Structure: Benzene ring with C(NH2)2+ group
/// Molecular weight: ~120 Da
/// Known to bind in the S1 specificity pocket of trypsin
pub fn create_benzamidine() -> MoleculeData {
    MoleculeData {
        name: "Benzamidine".to_string(),
        molecule_type: "small_molecule".to_string(),
        file_content: None,
        file_format: Some("pdb".to_string()),
        atoms: vec![
            // Benzene ring carbons (6 carbons in ring)
            Atom {
                element: "C".to_string(),
                x: 0.0,
                y: 0.0,
                z: 0.0,
                charge: Some(-0.06),
                atom_type: Some("C.ar".to_string()),
                atom_name: Some("C1".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(1),
            },
            Atom {
                element: "C".to_string(),
                x: 1.21,
                y: 0.70,
                z: 0.0,
                charge: Some(-0.06),
                atom_type: Some("C.ar".to_string()),
                atom_name: Some("C2".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(2),
            },
            Atom {
                element: "C".to_string(),
                x: 1.21,
                y: 2.10,
                z: 0.0,
                charge: Some(-0.06),
                atom_type: Some("C.ar".to_string()),
                atom_name: Some("C3".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(3),
            },
            Atom {
                element: "C".to_string(),
                x: 0.0,
                y: 2.80,
                z: 0.0,
                charge: Some(0.35),
                atom_type: Some("C.cat".to_string()),  // Cationic carbon
                atom_name: Some("C4".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(4),
            },
            Atom {
                element: "C".to_string(),
                x: -1.21,
                y: 2.10,
                z: 0.0,
                charge: Some(-0.06),
                atom_type: Some("C.ar".to_string()),
                atom_name: Some("C5".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(5),
            },
            Atom {
                element: "C".to_string(),
                x: -1.21,
                y: 0.70,
                z: 0.0,
                charge: Some(-0.06),
                atom_type: Some("C.ar".to_string()),
                atom_name: Some("C6".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(6),
            },
            // Amidinium group (C(NH2)2+)
            Atom {
                element: "C".to_string(),
                x: 0.0,
                y: 4.28,
                z: 0.0,
                charge: Some(0.64),
                atom_type: Some("C.cat".to_string()),
                atom_name: Some("C7".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(7),
            },
            // First amino group
            Atom {
                element: "N".to_string(),
                x: 1.15,
                y: 4.95,
                z: 0.0,
                charge: Some(-0.80),
                atom_type: Some("N.pl3".to_string()),
                atom_name: Some("N1".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(8),
            },
            Atom {
                element: "H".to_string(),
                x: 1.15,
                y: 5.95,
                z: 0.0,
                charge: Some(0.40),
                atom_type: Some("H".to_string()),
                atom_name: Some("H1A".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(9),
            },
            Atom {
                element: "H".to_string(),
                x: 2.05,
                y: 4.50,
                z: 0.0,
                charge: Some(0.40),
                atom_type: Some("H".to_string()),
                atom_name: Some("H1B".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(10),
            },
            // Second amino group
            Atom {
                element: "N".to_string(),
                x: -1.15,
                y: 4.95,
                z: 0.0,
                charge: Some(-0.80),
                atom_type: Some("N.pl3".to_string()),
                atom_name: Some("N2".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(11),
            },
            Atom {
                element: "H".to_string(),
                x: -1.15,
                y: 5.95,
                z: 0.0,
                charge: Some(0.40),
                atom_type: Some("H".to_string()),
                atom_name: Some("H2A".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(12),
            },
            Atom {
                element: "H".to_string(),
                x: -2.05,
                y: 4.50,
                z: 0.0,
                charge: Some(0.40),
                atom_type: Some("H".to_string()),
                atom_name: Some("H2B".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(13),
            },
            // Aromatic hydrogens on benzene ring (5 hydrogens)
            Atom {
                element: "H".to_string(),
                x: 0.0,
                y: -1.08,
                z: 0.0,
                charge: Some(0.12),
                atom_type: Some("H".to_string()),
                atom_name: Some("H3".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(14),
            },
            Atom {
                element: "H".to_string(),
                x: 2.15,
                y: 0.15,
                z: 0.0,
                charge: Some(0.12),
                atom_type: Some("H".to_string()),
                atom_name: Some("H4".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(15),
            },
            Atom {
                element: "H".to_string(),
                x: 2.15,
                y: 2.65,
                z: 0.0,
                charge: Some(0.12),
                atom_type: Some("H".to_string()),
                atom_name: Some("H5".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(16),
            },
            Atom {
                element: "H".to_string(),
                x: -2.15,
                y: 2.65,
                z: 0.0,
                charge: Some(0.12),
                atom_type: Some("H".to_string()),
                atom_name: Some("H6".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(17),
            },
            Atom {
                element: "H".to_string(),
                x: -2.15,
                y: 0.15,
                z: 0.0,
                charge: Some(0.12),
                atom_type: Some("H".to_string()),
                atom_name: Some("H7".to_string()),
                residue_name: Some("BEN".to_string()),
                residue_number: Some(1),
                chain_id: Some("L".to_string()),
                serial_number: Some(18),
            },
        ],
        bonds: vec![
            // Benzene ring
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0],
            // Amidinium connection
            [3, 6],
            // Amino groups
            [6, 7], [7, 8], [7, 9],
            [6, 10], [10, 11], [10, 12],
            // Aromatic hydrogens
            [0, 13], [1, 14], [2, 15], [4, 16], [5, 17],
        ],
    }
}

/// Creates a simplified trypsin binding pocket
/// Trypsin is a serine protease with a deep S1 specificity pocket
/// Key residues in binding site:
/// - Asp189: Provides negative charge for binding positively charged ligands
/// - Ser195, His57, Asp102: Catalytic triad
/// - Gly216, Gly219: Form oxyanion hole
///
/// This is a SIMPLIFIED version with just the key binding pocket residues
/// In reality, trypsin has ~1600 atoms in the full binding site
pub fn create_trypsin_pocket() -> MoleculeData {
    MoleculeData {
        name: "Trypsin_Binding_Pocket".to_string(),
        molecule_type: "protein".to_string(),
        file_content: None,
        file_format: Some("pdb".to_string()),
        atoms: vec![
            // ASP189 - Critical for binding benzamidine (provides negative charge)
            // Backbone
            Atom {
                element: "N".to_string(),
                x: 10.5,
                y: 12.3,
                z: 8.7,
                charge: Some(-0.47),
                atom_type: Some("N.am".to_string()),
                atom_name: Some("N".to_string()),
                residue_name: Some("ASP".to_string()),
                residue_number: Some(189),
                chain_id: Some("A".to_string()),
                serial_number: Some(1),
            },
            Atom {
                element: "C".to_string(),
                x: 11.2,
                y: 11.4,
                z: 9.5,
                charge: Some(0.07),
                atom_type: Some("C.3".to_string()),
                atom_name: Some("CA".to_string()),
                residue_name: Some("ASP".to_string()),
                residue_number: Some(189),
                chain_id: Some("A".to_string()),
                serial_number: Some(2),
            },
            Atom {
                element: "C".to_string(),
                x: 10.8,
                y: 9.9,
                z: 9.2,
                charge: Some(0.51),
                atom_type: Some("C.2".to_string()),
                atom_name: Some("C".to_string()),
                residue_name: Some("ASP".to_string()),
                residue_number: Some(189),
                chain_id: Some("A".to_string()),
                serial_number: Some(3),
            },
            Atom {
                element: "O".to_string(),
                x: 11.5,
                y: 9.0,
                z: 9.6,
                charge: Some(-0.51),
                atom_type: Some("O.2".to_string()),
                atom_name: Some("O".to_string()),
                residue_name: Some("ASP".to_string()),
                residue_number: Some(189),
                chain_id: Some("A".to_string()),
                serial_number: Some(4),
            },
            // Side chain
            Atom {
                element: "C".to_string(),
                x: 12.7,
                y: 11.6,
                z: 9.3,
                charge: Some(-0.16),
                atom_type: Some("C.3".to_string()),
                atom_name: Some("CB".to_string()),
                residue_name: Some("ASP".to_string()),
                residue_number: Some(189),
                chain_id: Some("A".to_string()),
                serial_number: Some(5),
            },
            Atom {
                element: "C".to_string(),
                x: 13.4,
                y: 12.8,
                z: 9.9,
                charge: Some(-0.16),
                atom_type: Some("C.3".to_string()),
                atom_name: Some("CG".to_string()),
                residue_name: Some("ASP".to_string()),
                residue_number: Some(189),
                chain_id: Some("A".to_string()),
                serial_number: Some(6),
            },
            // Carboxylate group (negatively charged - key for benzamidine binding)
            Atom {
                element: "C".to_string(),
                x: 14.9,
                y: 12.7,
                z: 9.8,
                charge: Some(0.75),
                atom_type: Some("C.2".to_string()),
                atom_name: Some("CD".to_string()),
                residue_name: Some("ASP".to_string()),
                residue_number: Some(189),
                chain_id: Some("A".to_string()),
                serial_number: Some(7),
            },
            Atom {
                element: "O".to_string(),
                x: 15.5,
                y: 13.7,
                z: 10.3,
                charge: Some(-0.80),
                atom_type: Some("O.co2".to_string()),
                atom_name: Some("OD1".to_string()),
                residue_name: Some("ASP".to_string()),
                residue_number: Some(189),
                chain_id: Some("A".to_string()),
                serial_number: Some(8),
            },
            Atom {
                element: "O".to_string(),
                x: 15.5,
                y: 11.7,
                z: 9.3,
                charge: Some(-0.80),
                atom_type: Some("O.co2".to_string()),
                atom_name: Some("OD2".to_string()),
                residue_name: Some("ASP".to_string()),
                residue_number: Some(189),
                chain_id: Some("A".to_string()),
                serial_number: Some(9),
            },
            // SER195 - Part of catalytic triad
            Atom {
                element: "N".to_string(),
                x: 9.6,
                y: 9.7,
                z: 8.5,
                charge: Some(-0.47),
                atom_type: Some("N.am".to_string()),
                atom_name: Some("N".to_string()),
                residue_name: Some("SER".to_string()),
                residue_number: Some(195),
                chain_id: Some("A".to_string()),
                serial_number: Some(10),
            },
            Atom {
                element: "C".to_string(),
                x: 9.1,
                y: 8.4,
                z: 8.0,
                charge: Some(0.07),
                atom_type: Some("C.3".to_string()),
                atom_name: Some("CA".to_string()),
                residue_name: Some("SER".to_string()),
                residue_number: Some(195),
                chain_id: Some("A".to_string()),
                serial_number: Some(11),
            },
            Atom {
                element: "C".to_string(),
                x: 10.2,
                y: 7.8,
                z: 7.1,
                charge: Some(-0.16),
                atom_type: Some("C.3".to_string()),
                atom_name: Some("CB".to_string()),
                residue_name: Some("SER".to_string()),
                residue_number: Some(195),
                chain_id: Some("A".to_string()),
                serial_number: Some(12),
            },
            Atom {
                element: "O".to_string(),
                x: 9.7,
                y: 6.7,
                z: 6.3,
                charge: Some(-0.65),
                atom_type: Some("O.3".to_string()),
                atom_name: Some("OG".to_string()),
                residue_name: Some("SER".to_string()),
                residue_number: Some(195),
                chain_id: Some("A".to_string()),
                serial_number: Some(13),
            },
            Atom {
                element: "H".to_string(),
                x: 10.4,
                y: 6.2,
                z: 5.9,
                charge: Some(0.43),
                atom_type: Some("H".to_string()),
                atom_name: Some("HG".to_string()),
                residue_name: Some("SER".to_string()),
                residue_number: Some(195),
                chain_id: Some("A".to_string()),
                serial_number: Some(14),
            },
        ],
        bonds: vec![
            // ASP189 backbone
            [0, 1], [1, 2], [2, 3],
            // ASP189 side chain
            [1, 4], [4, 5], [5, 6], [6, 7], [6, 8],
            // SER195 backbone and side chain
            [2, 9], [9, 10], [10, 11], [11, 12], [12, 13],
        ],
    }
}

/// Creates a complete benzamidine-trypsin docking request
pub fn create_benzamidine_trypsin_request() -> DockingRequest {
    DockingRequest {
        receptor: create_trypsin_pocket(),
        ligand: create_benzamidine(),
    }
}

/// Test runner for benzamidine-trypsin docking
/// Expected result: Benzamidine binds in S1 pocket with amidinium group
/// forming salt bridge with Asp189 carboxylate
pub fn run_benzamidine_trypsin_test() {
    let request = create_benzamidine_trypsin_request();

    println!("=== Benzamidine-Trypsin Docking Test ===");
    println!("Receptor: {} (simplified binding pocket)", request.receptor.name);
    println!("  Atoms: {} (key residues: ASP189, SER195)", request.receptor.atoms.len());
    println!("Ligand: {}", request.ligand.name);
    println!("  Atoms: {}", request.ligand.atoms.len());
    println!("  Molecular formula: C7H9N2+");
    println!("\nExpected outcome:");
    println!("  - Benzamidine binds in S1 specificity pocket");
    println!("  - Salt bridge: Amidinium group <-> Asp189 carboxylate");
    println!("  - Distance: C(NH2)2+ to COO- ~3-4 Å");
    println!("  - Known Ki: ~20 µM (strong binding)");
    println!("  - Reference: PDB 3PTB");
    println!("========================================\n");

    // TODO: Implement docking algorithm
    // Commented out - implement these functions:

    // let config = DockingConfig {
    //     search_algorithm: SearchAlgorithm::GeneticAlgorithm,
    //     scoring_function: ScoringFunction::AutoDockVina,
    //     exhaustiveness: 8,
    //     num_modes: 10,
    //     ..Default::default()
    // };

    // let result = dock(&request.receptor, &request.ligand, &config);

    // println!("Docking complete!");
    // println!("  Best score: {:.2} kcal/mol", result.score);
    // println!("  Best pose: {:?}", result.pose);
    // println!("  Key interactions:");
    // for interaction in &result.interactions {
    //     println!("    - {}: {:.2} Å", interaction.type_name, interaction.distance);
    // }

    // Validate result
    // assert!(result.has_salt_bridge("BEN", "ASP189"),
    //         "Expected salt bridge between benzamidine and Asp189");
    // assert!(result.score < -5.0,
    //         "Expected favorable binding energy (< -5 kcal/mol)");
}
