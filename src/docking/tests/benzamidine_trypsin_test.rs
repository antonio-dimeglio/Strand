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

use crate::core::atom::{Atom as CoreAtom, AtomType, HydrogenBond};
use crate::core::ligand::{Ligand, Bond};
use crate::core::receptor::Receptor;
use crate::core::space::{Space, Boundaries, RotationalSpace, ConformationalSpace};
use crate::core::pose::StepSizes;
use crate::docking::search::{random_search, mc_simulated_annealing};
use nalgebra::{Point3, Vector3};

/// Convert server Atom to core Atom
fn convert_atom(server_atom: &Atom) -> CoreAtom {
    let atom_type = AtomType {
        vdw_radius: match server_atom.element.as_str() {
            "C" => 1.70,
            "N" => 1.55,
            "O" => 1.52,
            "H" => 1.20,
            _ => 1.50,
        },
        well_depth: match server_atom.element.as_str() {
            "C" => 0.15,
            "N" => 0.16,
            "O" => 0.21,
            "H" => 0.02,
            _ => 0.10,
        },
        partial_charge: server_atom.charge.unwrap_or(0.0) as f32,
        hydrogen_bond: match server_atom.element.as_str() {
            "O" => HydrogenBond::Acceptor,
            "N" => HydrogenBond::Both,
            "H" if server_atom.charge.unwrap_or(0.0) > 0.2 => HydrogenBond::Donor,
            _ => HydrogenBond::Neither,
        },
    };

    CoreAtom::new(
        &server_atom.element,
        Vector3::new(server_atom.x as f32, server_atom.y as f32, server_atom.z as f32),
        atom_type,
    )
}

/// Create receptor from trypsin pocket
fn create_receptor() -> Receptor {
    let trypsin = create_trypsin_pocket();
    let atoms: Vec<CoreAtom> = trypsin.atoms.iter().map(|a| convert_atom(a)).collect();
    Receptor::new(atoms, 8.0)
}

/// Create ligand from benzamidine
fn create_ligand() -> Ligand {
    let benzamidine = create_benzamidine();
    let atoms: Vec<CoreAtom> = benzamidine.atoms.iter().map(|a| convert_atom(a)).collect();

    // Create bonds from the molecule data
    let num_atoms = atoms.len();
    let mut bonds: Vec<Vec<Bond>> = vec![vec![]; num_atoms];

    for bond in &benzamidine.bonds {
        let atom1 = bond[0] as usize;
        let atom2 = bond[1] as usize;

        // Benzamidine is mostly rigid, only the amidinium group has some flexibility
        // For simplicity, treat all bonds as rigid for now
        bonds[atom1].push(Bond { target: atom2, is_rotatable: false });
        bonds[atom2].push(Bond { target: atom1, is_rotatable: false });
    }

    Ligand::new(atoms, bonds, vec![])
}

/// Create search space for benzamidine-trypsin
fn create_search_space() -> Space {
    // Center search around Asp189 binding pocket
    let boundaries = Boundaries {
        center: Point3::new(14.0, 12.0, 9.0),  // Near Asp189 carboxylate
        extents: Vector3::new(5.0, 5.0, 5.0),
    };

    // Allow full rotation
    let rotational_space = RotationalSpace::Unconstrained;

    // No conformational flexibility (rigid ligand for now)
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
pub fn run_benzamidine_trypsin_test() {
    println!("=== Benzamidine-Trypsin Docking Test ===");
    println!("Note: This CLI runner is deprecated. Use 'cargo test benzamidine' instead.");
    println!("\nRunning tests programmatically...\n");

    // test_benzamidine_trypsin_random_search();
    println!();
    // test_benzamidine_trypsin_monte_carlo();
}

#[test]
fn test_benzamidine_trypsin_random_search() {
    println!("\n=== Benzamidine-Trypsin Random Search Test ===");

    let receptor = create_receptor();
    let ligand = create_ligand();
    let space = create_search_space();

    println!("Receptor: Trypsin binding pocket");
    println!("  Atoms: {} (ASP189, SER195)", receptor.atoms.len());
    println!("Ligand: Benzamidine (C7H9N2+)");
    println!("  Atoms: {}", ligand.atoms.len());

    println!("\nRunning random search with 100 samples...");
    let (_best_pose, best_score) = random_search(100, &space, &ligand, &receptor);

    println!("\nResults:");
    println!("  Total Energy: {:.3} kcal/mol", best_score.total);
    println!("  VDW:          {:.3} kcal/mol", best_score.vdw);
    println!("  Electrostatic: {:.3} kcal/mol", best_score.esi);
    println!("  H-bond:       {:.3} kcal/mol", best_score.h_bond);
    println!("  Torsional:    {:.3} kcal/mol", best_score.torsional);

    // Validate that we found a reasonable energy
    assert!(best_score.total < 0.0, "Best score should be negative (favorable binding)");

    // Expect strong electrostatic attraction (positive amidinium + negative carboxylate)
    assert!(best_score.esi < -5.0, "Should have strong electrostatic attraction for salt bridge");

    println!("\n✓ Random search completed successfully");
}

#[test]
fn test_benzamidine_trypsin_monte_carlo() {
    println!("\n=== Benzamidine-Trypsin Monte Carlo Test ===");

    let receptor = create_receptor();
    let ligand = create_ligand();
    let space = create_search_space();
    let step_sizes = StepSizes::default();

    println!("Receptor: Trypsin binding pocket");
    println!("  Atoms: {} (ASP189, SER195)", receptor.atoms.len());
    println!("Ligand: Benzamidine (C7H9N2+)");
    println!("  Atoms: {}", ligand.atoms.len());

    println!("\nRunning Monte Carlo with 2000 iterations...");
    println!("  Temperature: 2.0 kcal/mol");
    println!("  Step sizes: trans={:.1}Å, rot={:.0}°, tors={:.0}°",
             step_sizes.translation, step_sizes.rotation_degrees, step_sizes.torsion_degrees);

    let (_best_pose, best_score) = mc_simulated_annealing(
        2000,
        &space,
        &ligand,
        &receptor,
        &step_sizes,
        2.0,
    );

    println!("\nResults:");
    println!("  Total Energy: {:.3} kcal/mol", best_score.total);
    println!("  VDW:          {:.3} kcal/mol", best_score.vdw);
    println!("  Electrostatic: {:.3} kcal/mol", best_score.esi);
    println!("  H-bond:       {:.3} kcal/mol", best_score.h_bond);
    println!("  Torsional:    {:.3} kcal/mol", best_score.torsional);

    // Validate favorable binding
    assert!(best_score.total < -5.0, "Expected strong binding (< -5 kcal/mol)");

    // Salt bridge should dominate (electrostatics)
    assert!(best_score.esi < -10.0, "Expected strong salt bridge electrostatics");

    println!("\n✓ Monte Carlo search completed successfully");
    println!("\nExpected: Salt bridge between amidinium (BEN) and Asp189 carboxylate");
    println!("Reference: PDB 3PTB, Ki ~20 µM");
}
