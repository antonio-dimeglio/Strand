#![feature(core_float_math)]
mod server;
mod docking;
mod core;

use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(version, about = "Strand - Molecular Docking Application", long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Option<Command>,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Run as server for frontend
    Server,
    /// Run docking tests
    Test {
        #[command(subcommand)]
        test_type: TestType,
    },
}

#[derive(Subcommand, Debug)]
enum TestType {
    /// Run water dimer docking test
    WaterDimer,
    /// Run benzamidine-trypsin docking test
    BenzamidineTrypsin,
    /// Run all docking tests
    All,
}

#[tokio::main]
async fn main() {
    let args = Args::parse();

    match args.command {
        Some(Command::Server) => {
            server::run_server().await;
        }
        Some(Command::Test { test_type }) => {
            match test_type {
                TestType::WaterDimer => {
                    docking::tests::run_water_dimer_test();
                }
                TestType::BenzamidineTrypsin => {
                    docking::tests::run_benzamidine_trypsin_test();
                }
                TestType::All => {
                    docking::tests::run_water_dimer_test();
                    println!(); // Add spacing between tests
                    docking::tests::run_benzamidine_trypsin_test();
                }
            }
        }
        None => {
            println!("Strand - Molecular Docking Application");
            println!("\nUsage:");
            println!("  cargo run -- server                    Run as server for frontend");
            println!("  cargo run -- test water-dimer          Run water dimer test");
            println!("  cargo run -- test benzamidine-trypsin  Run benzamidine-trypsin test");
            println!("  cargo run -- test all                  Run all tests");
            println!("\nUse --help for more information");
        }
    }
}
