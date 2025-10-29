
use serde::{Deserialize, Serialize};

/// Docking request containing receptor and ligand molecules
#[derive(Serialize, Deserialize, Debug)]
pub struct DockingRequest {
    pub receptor: MoleculeData,
    pub ligand: MoleculeData,
}

/// Minimal molecule representation for docking simulations
/// Can represent proteins, small molecules, or complexes
#[derive(Serialize, Deserialize, Debug)]
pub struct MoleculeData {
    /// Molecule name/identifier
    pub name: String,

    /// Type of molecule: "protein", "small_molecule", "complex", "peptide", etc.
    #[serde(rename = "moleculeType", default = "default_molecule_type")]
    pub molecule_type: String,

    /// Original file content (PDB, SDF, MOL2, etc.) - PREFERRED for docking
    /// Preserves all format-specific information (charges, bond orders, etc.)
    #[serde(skip_serializing_if = "Option::is_none", rename = "fileContent")]
    pub file_content: Option<String>,

    /// File format: "pdb", "cif", "sdf", "mol", "mol2", "pdbqt", etc.
    #[serde(skip_serializing_if = "Option::is_none", rename = "fileFormat")]
    pub file_format: Option<String>,

    /// Parsed atom data (fallback if fileContent not available)
    #[serde(default)]
    pub atoms: Vec<Atom>,

    /// Bond connectivity (atom index pairs)
    #[serde(default)]
    pub bonds: Vec<[u32; 2]>,
}

/// Default molecule type if not specified
fn default_molecule_type() -> String {
    "unknown".to_string()
}

/// Minimal atom representation for docking
/// Only includes data relevant to molecular docking algorithms
#[derive(Debug, Serialize, Deserialize)]
pub struct Atom {
    /// Element symbol (C, N, O, H, etc.) - ESSENTIAL for force fields
    pub element: String,

    /// 3D coordinates - ESSENTIAL
    pub x: f64,
    pub y: f64,
    pub z: f64,

    /// Optional: Partial charge (for electrostatics in docking)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub charge: Option<f64>,

    /// Optional: Atom type (for force field parameters: sp2, sp3, aromatic, etc.)
    #[serde(skip_serializing_if = "Option::is_none", rename = "atomType")]
    pub atom_type: Option<String>,

    // Fields below are useful for context but not strictly required for docking

    /// Atom name (CA, CB, N1, etc.) - useful for protein residue identification
    #[serde(skip_serializing_if = "Option::is_none", rename = "atomName")]
    pub atom_name: Option<String>,

    /// Residue name (ALA, GLY, etc. for proteins; ligand name for small molecules)
    #[serde(skip_serializing_if = "Option::is_none", rename = "residueName")]
    pub residue_name: Option<String>,

    /// Residue number (for proteins)
    #[serde(skip_serializing_if = "Option::is_none", rename = "residueNumber")]
    pub residue_number: Option<i32>,

    /// Chain ID (for multi-chain proteins)
    #[serde(skip_serializing_if = "Option::is_none", rename = "chainID")]
    pub chain_id: Option<String>,

    /// Serial number (atom index)
    #[serde(skip_serializing_if = "Option::is_none", rename = "serialNumber")]
    pub serial_number: Option<u32>,
}

// All visualization-specific structs removed - not needed for docking!

