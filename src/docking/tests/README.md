### 1. Water Dimer Test (`water_dimer_test.rs`)

**System**: Two water molecules (H₂O)
- Total atoms: 6 (3 per molecule)
- Simple hydrogen bonding interaction

**Expected Result**:
- Optimal O-O distance: ~2.8 Å
- Hydrogen bond: O...H distance ~1.8 Å
- Linear O-H...O arrangement

### 2. Benzamidine-Trypsin Test (`benzamidine_trypsin_test.rs`)


**System**:
- **Receptor**: Trypsin binding pocket (simplified, 14 atoms)
  - Key residues: ASP189 (negative charge), SER195 (catalytic)
- **Ligand**: Benzamidine (18 atoms, C₇H₉N₂⁺)
  - Positively charged amidinium group

**Expected Result**:
- Benzamidine binds in S1 specificity pocket
- Salt bridge: Amidinium (NH₂⁺) ↔ Asp189 carboxylate (COO⁻)
- Distance: ~3-4 Å between charged groups
- Known experimental Ki: ~20 µM

**Reference**: PDB ID 3PTB
## Running Tests

To run tests in standalone mode (without server):

```bash
# Run water dimer test
cargo run -- --test water-dimer

# Run benzamidine-trypsin test
cargo run -- --test benzamidine-trypsin

# Run all tests
cargo run -- --test all
```

## Test Data Format

All test cases use the `MoleculeData` format defined in `src/server/protein.rs`:

- **Atoms**: Full 3D coordinates, element types, partial charges
- **Bonds**: Connectivity information
- **Metadata**: Residue names, chain IDs, atom types

## Future Test Cases

Consider adding:
- **Aspirin-COX2**: Non-covalent binding with hydrophobic interactions
- **HIV protease inhibitors**: Symmetric binding site
- **Metal coordination**: Zinc-binding ligands
- **Covalent docking**: Warhead-reactive site
