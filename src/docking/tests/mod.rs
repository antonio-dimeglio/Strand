/// Test cases for molecular docking algorithms
///
/// This module contains various test scenarios for validating docking implementations:
/// - Simple molecular systems (water dimer) for basic validation
/// - Real protein-ligand complexes (benzamidine-trypsin) for realistic testing

pub mod water_dimer_test;
pub mod benzamidine_trypsin_test;

// Re-export test runners for easy access
pub use water_dimer_test::run_water_dimer_test;
pub use benzamidine_trypsin_test::run_benzamidine_trypsin_test;
