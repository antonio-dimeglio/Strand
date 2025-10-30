use std::cmp::{Ordering, PartialOrd};
use rand::{random, thread_rng, Rng};
use crate::core::ligand::Ligand;
use crate::core::pose::{Pose, StepSizes};
use crate::core::receptor::Receptor;
use crate::core::space::Space;
use crate::docking::scoring::{compute_scoring, ScoringResult};

const DEFAULT_CUTOFF: f32 = 10.0;


/// Random search algorithm for pose space exploration.
/// This method explores the pose space to find low-energy finding configuration by generating
/// N random poses within the search space, returning the one with the highest score.
pub fn random_search(num_samples: usize, search_space: &Space, ligand: &Ligand, receptor: &Receptor) -> (Pose, ScoringResult) {
    let mut best_sample = Pose::random_pose(search_space, ligand);
    let mut best_score = compute_scoring(ligand, receptor, &best_sample, DEFAULT_CUTOFF);

    for i in 0..num_samples-1 {
        let curr_sample = Pose::random_pose(search_space, ligand);
        let curr_score = compute_scoring(ligand, receptor, &curr_sample, DEFAULT_CUTOFF);

        if curr_score.total < best_score.total {
            best_score = curr_score;
            best_sample = curr_sample;
        }
    }

    (best_sample, best_score)
}



pub fn mc_simulated_annealing(
    n_iterations: usize,
    search_space: &Space,
    ligand: &Ligand,
    receptor: &Receptor,
    step_sizes: &StepSizes,
    temperature: f32) -> (Pose, ScoringResult) {
    let mut curr_pose = Pose::random_pose(search_space, ligand);
    let mut curr_score = compute_scoring(ligand, receptor, &curr_pose, DEFAULT_CUTOFF);
    let mut best_pose = curr_pose.clone();
    let mut best_score = curr_score.clone();
    let mut rng = thread_rng();

    for _ in 1..n_iterations {
        let new_pose = curr_pose.perturb(search_space, step_sizes);
        let new_score = compute_scoring(ligand, receptor, &new_pose, DEFAULT_CUTOFF);
        let delta_e = new_score.total - curr_score.total;

        if delta_e < 0.0 {
            curr_pose = new_pose;
            curr_score = new_score;
        } else {
            let p = (-delta_e/temperature).exp();
            if rng.gen_range(0.0..1.0) < p {
                curr_pose = new_pose;
                curr_score = new_score;
            }
        }

        if curr_score.total < best_score.total {
            best_score = curr_score.clone();
            best_pose = curr_pose.clone();
        }
    }

    (best_pose, best_score)
}