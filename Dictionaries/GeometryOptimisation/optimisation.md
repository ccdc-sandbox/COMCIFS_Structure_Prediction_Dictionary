# Optimisation Algortihms Dictionary

## Current Fields

| Data Field      | Type | Definition                      | Constraints | Units | Example                                                                                       |
|-----------------|------|---------------------------------|-------------|-------|-----------------------------------------------------------------------------------------------|
| _opt_algorithm  | char | Geometry optimisation algorithm |             |       | - BFGS<br>- L-BFGS<br>- Quasi-Newton<br>- FIRE<br>- Steepest Descent<br>- Conjugate Gradient  |


## Proposed Fields
### General Inputs
| Data Field                   | Type | Definition                                                                                                                                                                                                  | Constraints           | Units  | Example                |
|------------------------------|------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|--------|------------------------|
| _opt_cell                    | char | It can be "fixed" for no cell optimisation, "isotropic" or "anisotropic" for cell relaxation calculations                                                                                                   |                       |        |                        |
| _opt_atoms                   | char | It can be "fixed" for no cell optimisation, "all" for all-atoms geometry optimization, "hydrogens" for optimisation only H atoms, "non-hydrogens" for non-H atoms or a list of atoms for custom relaxation. |                       |        |                        |
| _opt_relax_force_convergence | numb | Convergence criteria for stopping the geometry optimisation. Present in TCOD as `_dft_atom_relax_force_conv`                                                                                                |                       |        |                        |
