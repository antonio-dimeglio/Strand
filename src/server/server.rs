use std::net::SocketAddr;
use axum::{Json, Router};
use axum::routing::post;
use axum::http::Method;
use tokio::net::TcpListener;
use tower_http::cors::{CorsLayer, Any};
use crate::server::protein::DockingRequest;

async fn handler(Json(payload): Json<DockingRequest>) -> Json<serde_json::Value> {
    println!("\n=== Docking Request Received ===");
    println!("Receptor:");
    println!("  Name: {}", payload.receptor.name);
    println!("  Type: {}", payload.receptor.molecule_type);
    println!("  Atoms: {}", payload.receptor.atoms.len());
    println!("  Format: {:?}", payload.receptor.file_format);
    println!("  Has original file: {}", payload.receptor.file_content.is_some());

    println!("\nLigand:");
    println!("  Name: {}", payload.ligand.name);
    println!("  Type: {}", payload.ligand.molecule_type);
    println!("  Atoms: {}", payload.ligand.atoms.len());
    println!("  Format: {:?}", payload.ligand.file_format);
    println!("  Has original file: {}", payload.ligand.file_content.is_some());
    println!("================================\n");

    // Return a success response
    Json(serde_json::json!({
        "status": "success",
        "message": "Docking simulation received",
        "receptor": {
            "name": payload.receptor.name,
            "type": payload.receptor.molecule_type,
            "atoms": payload.receptor.atoms.len(),
            "format": payload.receptor.file_format,
        },
        "ligand": {
            "name": payload.ligand.name,
            "type": payload.ligand.molecule_type,
            "atoms": payload.ligand.atoms.len(),
            "format": payload.ligand.file_format,
        }
    }))
}

pub async fn run_server() {
    // Configure CORS to allow requests from any origin
    let cors = CorsLayer::new()
        .allow_origin(Any)
        .allow_methods([Method::GET, Method::POST, Method::OPTIONS])
        .allow_headers(Any);

    let app = Router::new()
        .route("/dock", post(handler))  // Changed from /data to /dock
        .layer(cors);

    let addr = SocketAddr::from(([127, 0, 0, 1], 3000));
    println!("Strand server running on {}", addr);

    let listener = TcpListener::bind(addr).await.unwrap();
    axum::serve(listener, app.into_make_service())
        .await
        .unwrap();
}