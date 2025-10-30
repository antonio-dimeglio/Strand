use std::cmp::Ordering;
use crate::core::atom::Atom;
use crate::core::atom::HydrogenBond::{Acceptor, Both, Donor};
use crate::core::ligand::Ligand;
use crate::core::pose::Pose;
use crate::core::receptor::Receptor;

const COULOMB_CONSTANT: f32 = 332.0;
const DIELECTRIC_CONSTANT: f32 = 4.0;
const E_0: f32 = 5.0; // Maximum H-Bond strength
const R_0: f32 = 2.8; // Optimal H-bond distance
const W_TORSION: f32 = 0.5; // KCal/Mol, energy penalization for each rotatable bond, so that simpler conformations are preferred.
#[derive(PartialEq, Clone)]
pub struct ScoringResult {
    pub vdw: f32,
    pub esi: f32,
    pub h_bond: f32,
    pub torsional: f32,
    pub total: f32
}

impl ScoringResult {
    pub fn weighted_scoring(&self,
        weights: [f32; 4]) -> ScoringResult {
        ScoringResult {
            vdw: self.vdw * weights[0],
            esi: self.esi * weights[1],
            h_bond: self.h_bond * weights[2],
            torsional: self.torsional * weights[3],
            total: self.vdw * weights[0] +
                self.esi * weights[1] +
                self.h_bond * weights[2] +
                self.torsional * weights[3]
        }
    }
}

/// Applies Van der Waals interaction scoring on ligand-protein binding.
/// This interaction only considers energy between ligand and receptor pairs, i.e., for each
/// atom of the ligand it checks which protein atoms are in the neighborhood and if
/// within significantly low enough distance (usually, cutoff is ~10 A) it computes the Lennard-Jones potential.
/// This works under the assumption that the ligand is already in an energy stable conformation and
/// that the receptor is rigid, i.e., that its internal energy is constant.
/// The same assumptions are held for all other scoring functions.
fn apply_vdw_scoring(ligand_atom: &Atom, receptor_atom: &Atom, distance: f32) -> f32 {
    let r_min = ligand_atom.atom_type.vdw_radius + receptor_atom.atom_type.vdw_radius;
    let epsilon = (ligand_atom.atom_type.well_depth * receptor_atom.atom_type.well_depth).sqrt();
    let ratio = r_min / distance;
    let e_vdw = epsilon * (ratio.powi(12) - 2.0 * ratio.powi(6));
    e_vdw
}

/// Applies Electrostatic interaction scoring on ligand-protein binding.
/// Unlike VdW scoring, this affects pairs of atoms at great distances and is much stronger in terms
/// of effects to the energy of the system. Currently, the function assume a fixed constant dielectric,
/// however, this should be computed depending on the environment (i.e. protein interior, bulk water)
/// and is often time distance dependent.
fn apply_esi_scoring(ligand_atom: &Atom, receptor_atom: &Atom, distance: f32) -> f32 {
    let q1 = ligand_atom.atom_type.partial_charge;
    let q2 = receptor_atom.atom_type.partial_charge;
    let esi = (COULOMB_CONSTANT * q1 * q2) / (DIELECTRIC_CONSTANT * distance);
    esi
}

fn is_hbond_compatible(atom1: &Atom, atom2: &Atom) -> bool {
    let type1 = &atom1.atom_type.hydrogen_bond;
    let type2 = &atom2.atom_type.hydrogen_bond;

    let donors = [Donor, Both];
    let acceptors = [Acceptor, Both];

    (donors.contains(&type1) && acceptors.contains(&type2)) ||
    (donors.contains(&type2) && acceptors.contains(&type1))
}

/// Applies Hydrogen bond scoring on ligand-protein binding.
/// This can be only considered under very short distances.
/// Intuitively, it can be thought of as a directional lennard-jones with angular terms.
/// Fully calculating this requires knowing H atom position for donor-H vector. This is simplified by
/// assuming that the hydrogen bonds is impactful if the distance between atoms is within 2.5 and 3.5 A.
fn apply_hydrogen_scoring(ligand_atom: &Atom, receptor_atom: &Atom, distance: f32) -> f32 {
    if !is_hbond_compatible(ligand_atom, receptor_atom) && (distance < 2.5 || distance > 3.5){
        return 0.0;
    }
    // TODO: Add angle check, for now the assumption is that of favourable geometry.
    -5.0
}

/// Computes the conformational entropy associated with each rotatable bond, as upon binding the
/// rotational freedom is lost. This measurement is computed solely on the ligand, and is independent
/// of pose geometry.
fn apply_torsional_penalty(ligand: &Ligand) -> f32 {
    let n_rot = ligand.rotatable_bonds.len();
    n_rot as f32 * W_TORSION
}

pub fn compute_scoring(ligand: &Ligand, receptor: &Receptor, pose: &Pose, cutoff: f32) -> ScoringResult {
    // Generate transformed coordinates from pose
    let transformed_coords = pose.generate_coordinates(ligand);

    let mut total_vdw = 0.0;
    let mut total_esi = 0.0;
    let mut total_h_bond = 0.0;
    let total_tors = apply_torsional_penalty(ligand);

    // Score using transformed coordinates
    for (i, ligand_atom) in ligand.atoms.iter().enumerate() {
        let atom_pos = transformed_coords[i];  // Use transformed position from pose
        let nearby_receptor_atoms = receptor.find_nearby_atoms(&atom_pos, cutoff);

        for idx in nearby_receptor_atoms {
            let receptor_atom = &receptor.atoms[idx];
            let distance = (atom_pos - receptor_atom.coords).norm();
            total_vdw += apply_vdw_scoring(ligand_atom, receptor_atom, distance);
            total_esi += apply_esi_scoring(ligand_atom, receptor_atom, distance);
            total_h_bond += apply_hydrogen_scoring(ligand_atom, receptor_atom, distance);
        }
    }

    ScoringResult {
        vdw: total_vdw,
        esi: total_esi,
        h_bond: total_h_bond,
        torsional: total_tors,
        total: total_vdw + total_esi + total_h_bond  + total_tors
    }
}