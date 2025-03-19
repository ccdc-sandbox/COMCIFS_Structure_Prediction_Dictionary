# Forcefields Core Dictionary

## Proposed Fields

| Data Field                         | Type | Definition | Constraints | Units | Example                                                                                                        |
|------------------------------------|------|------------|-------------|-------|----------------------------------------------------------------------------------------------------------------|
| **1. Forcefield ID**               |      |            |             |       |                                                                                                                |
| _ff_type                           | char |            |             |       | Fixed charge, Polarizable, Distributed multipole based FFs, Hybrid (mixed intra-inter), Reactive (ReaxFF etc.) |
| _ff_name                           | char |            |             |       | GAFF, OPLS, PW                                                                                                 |
| _ff_version                        | char |            |             |       |                                                                                                                |
| _ff_citation_id                    | char |            |             |       |                                                                                                                |
| **2. FF Intramolecular Terms**     |      |            |             |       |                                                                                                                |
| _ff_bond_type                      | char |            |             |       | Harmonic, Morse, Cubic bond stretching potential, ReaxFF Bond Term etc                                         |
| _ff_angle_type                     | char |            |             |       | Harmonic, Cosine, Restricted bending, Linear, ReaxFF Angle Term, etc                                           |
| _ff_improper_type                  | char |            |             |       | Harmonic                                                                                                       |
| _ff_dihedral_type                  | char |            |             |       | Periodic, Ryckaert-Bellemans, Fourier, ReaxFF Torsion Term                                                     |
| _ff_cross_terms                    | char |            |             |       |                                                                                                                |
| _ff_1-4_pairs                      | bool |            |             |       | True/False                                                                                                     |
| _ff_fudgeLJ                        | numb |            |             |       | 0.5                                                                                                            |
| _ff_fudgeQQ                        | numb |            |             |       | 0.5                                                                                                            |
| **3. Molecular Quantum Chemistry** |      |            |             |       |                                                                                                                |
| _mqc_functional                    | char |            |             |       | B3LYP, PBE, PBE0, MP2                                                                                          |
| _mqc_basis_set                     | char |            |             |       | 6-31G(d,p), cc-pVTZ                                                                                            |                                                                                                             | 
| **4. FF Intermolecular Terms**     |      |            |             |       |                                                                                                                |  
| _ff_electrostatic_func             | char |            |             |       | Point-Charge,  Multipoles                                                                                      | 
| _ff_charges_generation             | char |            |             |       | BCC-AM1, Mulliken-AM1, etc.                                                                                    | 
| _ff_multipoles_generation          | char |            |             |       | DMA                                                                                                            | 
| _ff_long_range_electrostatic       | char |            |             |       | Ewald, PME, etc.                                                                                               | 
| _ff_polarization                   | char |            |             |       | Thole, Anharmonic                                                                                              | 
| _ff_vdw_func                       | char |            |             |       | LJ(C6,C12), LJ(epsilon,sigma), Buckingham, ReaxFF Morse-Potential, 14-7 function                               | 
| _ff_long_range_vdw                 | char |            |             |       | Energy/Pressure/EnergyPressure DispersionCorrection, LJ PME                                                    | 
| _ff_combination_rule               | char |            |             |       | Geometric, Arithmetic, Lorentz-Berthelot                                                                       | 
