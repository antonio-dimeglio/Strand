use std::collections::HashMap;
use nalgebra::Vector3;
use crate::core::atom::Atom;

pub struct SpatialGrid {
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
    cell_size: f32
}


impl SpatialGrid {

    fn compute_cell(pos: Vector3<f32>, grid_size:f32) -> (i32, i32, i32) {
        let cell_x = (pos.x / grid_size).floor() as i32;
        let cell_y = (pos.y / grid_size).floor() as i32;
        let cell_z = (pos.z / grid_size).floor() as i32;

        return (cell_x, cell_y, cell_z)
    }

    pub fn new(atoms: &Vec<Atom>, cell_size: f32) -> SpatialGrid {
        let mut cells:HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();

        for (i, atom) in atoms.iter().enumerate() {
            let pos = atom.coords;
            let cell = Self::compute_cell(pos, cell_size);
            cells.entry(cell)
                .or_insert_with(|| Vec::new())
                .push(i);
        }

        SpatialGrid {
            cells,
            cell_size
        }
    }

    pub fn get_neighbor_cells(&self, position: &Vector3<f32>) -> Vec<(i32, i32, i32)> {
        let cell_x = (position.x / self.cell_size).floor() as i32;
        let cell_y = (position.y / self.cell_size).floor() as i32;
        let cell_z = (position.z / self.cell_size).floor() as i32;

        let mut res = Vec::with_capacity(27);

        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    res.push((cell_x + dx, cell_y + dy, cell_z + dz));
                }
            }
        }

        res
    }

    pub fn query(&self, position: &Vector3<f32>, atoms: &Vec<Atom>, cutoff: f32) -> Vec<usize> {
        let mut nearby_atoms = Vec::new();

        for neighbor in self.get_neighbor_cells(position) {
            if let Some(neighbor_cell) = self.cells.get(&neighbor) {
                for atom_idx in neighbor_cell {
                    let distance = (position - atoms[*atom_idx].coords).norm();
                    if distance <= cutoff {
                        nearby_atoms.push(*atom_idx);
                    }
                }
            }
        }

        return nearby_atoms
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{Atom, AtomType, HydrogenBond};

    fn make_test_atom(x: f32, y: f32, z: f32) -> Atom {
        let atom_type = AtomType {
            vdw_radius: 1.5,
            well_depth: 0.1,
            partial_charge: 0.0,
            hydrogen_bond: HydrogenBond::Neither,
        };
        Atom::new("C", Vector3::new(x, y, z), atom_type)
    }

    #[test]
    fn test_compute_cell() {
        // Test cell computation with grid size 5.0
        let cell1 = SpatialGrid::compute_cell(Vector3::new(0.0, 0.0, 0.0), 5.0);
        assert_eq!(cell1, (0, 0, 0));

        let cell2 = SpatialGrid::compute_cell(Vector3::new(5.0, 5.0, 5.0), 5.0);
        assert_eq!(cell2, (1, 1, 1));

        let cell3 = SpatialGrid::compute_cell(Vector3::new(7.3, 12.8, 3.2), 5.0);
        assert_eq!(cell3, (1, 2, 0));

        // Test negative coordinates
        let cell4 = SpatialGrid::compute_cell(Vector3::new(-2.0, -3.0, -1.0), 5.0);
        assert_eq!(cell4, (-1, -1, -1));
    }

    #[test]
    fn test_single_atom_nearby() {
        // Single atom, query very close to it
        let atoms = vec![make_test_atom(5.0, 5.0, 5.0)];
        let grid = SpatialGrid::new(&atoms, 5.0);

        // Query right next to the atom
        let query_pos = Vector3::new(5.1, 5.1, 5.1);
        let nearby = grid.query(&query_pos, &atoms, 1.0);

        assert_eq!(nearby.len(), 1);
        assert_eq!(nearby[0], 0);
    }

    #[test]
    fn test_single_atom_too_far() {
        // Single atom, query too far away
        let atoms = vec![make_test_atom(5.0, 5.0, 5.0)];
        let grid = SpatialGrid::new(&atoms, 5.0);

        // Query far from the atom
        let query_pos = Vector3::new(15.0, 15.0, 15.0);
        let nearby = grid.query(&query_pos, &atoms, 1.0);

        assert_eq!(nearby.len(), 0);
    }

    #[test]
    fn test_multiple_atoms_different_cells() {
        // Atoms in different cells, far apart
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),    // Cell (0,0,0)
            make_test_atom(20.0, 20.0, 20.0), // Cell (4,4,4) with size 5.0
            make_test_atom(2.0, 2.0, 2.0),    // Cell (0,0,0)
        ];
        let grid = SpatialGrid::new(&atoms, 5.0);

        // Query near atoms 0 and 2
        let query_pos = Vector3::new(1.0, 1.0, 1.0);
        let nearby = grid.query(&query_pos, &atoms, 5.0);

        // Should find atoms 0 and 2, but not 1
        assert_eq!(nearby.len(), 2);
        assert!(nearby.contains(&0));
        assert!(nearby.contains(&2));
        assert!(!nearby.contains(&1));
    }

    #[test]
    fn test_atoms_in_adjacent_cells() {
        // Atoms in adjacent cells should be found
        let atoms = vec![
            make_test_atom(4.9, 4.9, 4.9),  // Near cell boundary
            make_test_atom(5.1, 5.1, 5.1),  // Just across boundary in adjacent cell
        ];
        let grid = SpatialGrid::new(&atoms, 5.0);

        // Query at the boundary
        let query_pos = Vector3::new(5.0, 5.0, 5.0);
        let nearby = grid.query(&query_pos, &atoms, 1.0);

        // Should find both atoms (they're within 1.0 Å)
        assert_eq!(nearby.len(), 2);
    }

    #[test]
    fn test_cutoff_distance_exact() {
        // Test exact cutoff distance
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),
            make_test_atom(3.0, 0.0, 0.0),  // Exactly 3.0 Å away
        ];
        let grid = SpatialGrid::new(&atoms, 5.0);

        let query_pos = Vector3::new(0.0, 0.0, 0.0);

        // With cutoff 3.0, should find both
        let nearby = grid.query(&query_pos, &atoms, 3.0);
        assert_eq!(nearby.len(), 2);

        // With cutoff 2.9, should find only atom 0
        let nearby = grid.query(&query_pos, &atoms, 2.9);
        assert_eq!(nearby.len(), 1);
        assert_eq!(nearby[0], 0);
    }

    #[test]
    fn test_neighbor_cells_count() {
        // Verify we're checking exactly 27 neighbor cells
        let atoms = vec![make_test_atom(0.0, 0.0, 0.0)];
        let grid = SpatialGrid::new(&atoms, 5.0);

        let neighbors = grid.get_neighbor_cells(&Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(neighbors.len(), 27);
    }

    #[test]
    fn test_many_atoms_in_same_cell() {
        // Multiple atoms in the same cell
        let atoms = vec![
            make_test_atom(1.0, 1.0, 1.0),
            make_test_atom(1.5, 1.5, 1.5),
            make_test_atom(2.0, 2.0, 2.0),
            make_test_atom(2.5, 2.5, 2.5),
        ];
        let grid = SpatialGrid::new(&atoms, 5.0);  // All in cell (0,0,0)

        let query_pos = Vector3::new(1.0, 1.0, 1.0);
        let nearby = grid.query(&query_pos, &atoms, 2.0);

        // Should find atoms within 2.0 Å
        // Distance to atom 0: 0.0
        // Distance to atom 1: sqrt(3*0.5^2) ≈ 0.866
        // Distance to atom 2: sqrt(3*1.0^2) ≈ 1.732
        // Distance to atom 3: sqrt(3*1.5^2) ≈ 2.598 (too far)
        assert_eq!(nearby.len(), 3);
        assert!(nearby.contains(&0));
        assert!(nearby.contains(&1));
        assert!(nearby.contains(&2));
    }

    #[test]
    fn test_negative_coordinates() {
        // Test with negative coordinates
        let atoms = vec![
            make_test_atom(-5.0, -5.0, -5.0),
            make_test_atom(-4.0, -4.0, -4.0),
            make_test_atom(5.0, 5.0, 5.0),
        ];
        let grid = SpatialGrid::new(&atoms, 5.0);

        let query_pos = Vector3::new(-4.5, -4.5, -4.5);
        let nearby = grid.query(&query_pos, &atoms, 2.0);

        // Should find atoms 0 and 1, not 2
        assert_eq!(nearby.len(), 2);
        assert!(nearby.contains(&0));
        assert!(nearby.contains(&1));
    }

    #[test]
    fn test_empty_cells() {
        // Grid with sparse atoms (many empty cells)
        let atoms = vec![
            make_test_atom(0.0, 0.0, 0.0),
            make_test_atom(100.0, 100.0, 100.0),
        ];
        let grid = SpatialGrid::new(&atoms, 5.0);

        // Query in empty region
        let query_pos = Vector3::new(50.0, 50.0, 50.0);
        let nearby = grid.query(&query_pos, &atoms, 5.0);

        assert_eq!(nearby.len(), 0);
    }

    #[test]
    fn test_large_receptor() {
        // Simulate larger receptor (100 atoms)
        let mut atoms = Vec::new();
        for i in 0..100 {
            let x = (i % 10) as f32 * 3.0;
            let y = ((i / 10) % 10) as f32 * 3.0;
            let z = (i / 100) as f32 * 3.0;
            atoms.push(make_test_atom(x, y, z));
        }

        let grid = SpatialGrid::new(&atoms, 5.0);

        // Query in the middle
        let query_pos = Vector3::new(15.0, 15.0, 0.0);
        let nearby = grid.query(&query_pos, &atoms, 5.0);

        // Should find some atoms but not all 100
        assert!(nearby.len() > 0);
        assert!(nearby.len() < 100);

        // All nearby atoms should actually be within cutoff
        for &idx in &nearby {
            let distance = (query_pos - atoms[idx].coords).norm();
            assert!(distance <= 5.0, "Atom {} at distance {} > cutoff", idx, distance);
        }
    }

    #[test]
    fn test_correctness_vs_naive() {
        // Compare grid results with naive O(N) search
        let atoms = vec![
            make_test_atom(1.0, 2.0, 3.0),
            make_test_atom(4.0, 5.0, 6.0),
            make_test_atom(7.0, 8.0, 9.0),
            make_test_atom(2.0, 3.0, 4.0),
            make_test_atom(5.0, 6.0, 7.0),
        ];

        let grid = SpatialGrid::new(&atoms, 5.0);
        let query_pos = Vector3::new(3.0, 4.0, 5.0);
        let cutoff = 4.0;

        // Grid search
        let mut grid_result = grid.query(&query_pos, &atoms, cutoff);
        grid_result.sort();

        // Naive search
        let mut naive_result = Vec::new();
        for (i, atom) in atoms.iter().enumerate() {
            let distance = (query_pos - atom.coords).norm();
            if distance <= cutoff {
                naive_result.push(i);
            }
        }
        naive_result.sort();

        // Should get identical results
        assert_eq!(grid_result, naive_result);
    }
}