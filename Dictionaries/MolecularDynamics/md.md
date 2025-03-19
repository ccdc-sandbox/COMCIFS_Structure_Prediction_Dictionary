# Molecular Dynamics Dictionary

## Proposed Fields
### General Inputs
| Data Field                         | Type      | Definition | Constraints | Units            | Example                                    |
|------------------------------------|-----------|------------|-------------|------------------|--------------------------------------------|
| **1. Simulation**                  |           |            |             |                  |                                            |
| _md_integrator                     | char      |            |             |                  | leap-frog algorithm, velocity Verlet, etc. |
| _md_dt                             | numb      |            |             | ps               | 0.001 ps                                   |
| _md_nsteps                         | numb      |            |             |                  | 1000000                                    |
| **2. Neighbor searching**          |           |            |             |                  |                                            |
| _md_cutoff_scheme                  | char      |            |             |                  | Verlet                                     |
| _md_list_update                    | numb      |            |             |                  | 10                                         |
| _md_list_cutoff                    | numb      |            |             | nm               | 1.0 nm                                     |
| **3. Intermolecular Interactions** |           |            |             |                  |                                            |
| _md_coulomb_modifier               | char      |            |             |                  | Potential-shift                            |
| _md_coulomb_cutoff                 | numb      |            |             | nm               | 1.0 nm                                     |
| _md_vdw_modifier                   | char      |            |             |                  | Potential-shift                            |
| _md_vdw_cutoff                     | numb      |            |             | nm               | 1.0 nm                                     |
| **4. Temperature Control**         |           |            |             |                  |                                            |
| _md_thermostat                     | char      |            |             |                  | berendsen, nose-hoover,  andersen          |
| _md_reference_temperature          | numb      |            |             | K                | 300 K                                      |
| _md_temperature_coupling_time      | numb      |            |             | ps               | 10.0 ps                                    |
| _md_temperature_initial_velocities | char      |            |             |                  | Sampled from Maxwell-Boltzmann             |
| _md_annealing_type                 | char      |            |             |                  | None, Single, Periodic                     |
| _md_annealing_time                 | list      |            |             | ps               | 500.0 500.0 500.0                          |
| _md_annealing_temp                 | list      |            |             | K                | 150 200 250                                |
| **5. Pressure Control**            |           |            |             |                  |                                            |
| _md_barostat                       | char      |            |             |                  | Berendsen, Parrinello-Rahman, C-rescale    |
| _md_pressure_coupling_type         | char      |            |             |                  | isotropic, anisotropic                     |
| _md_pressure_coupling_time         | numb      |            |             | ps               | 10.0 ps                                    |
| _md_compressibility                | numb/list |            |             | Pa<sup>-1</sup>  | 1.0                                        |
| _md_reference_pressure             | numb      |            |             | Pa               | 100000 Pa                                  |
| _md_coordinate_scaling             | char      |            |             |                  | "com" or "all-atoms"                       |
| **6. Constraints**                 |           |            |             |                  |                                            |
| _md_constraints                    | char      |            |             |                  | h-bonds, all, etc.                         |
| _md_constraints_algorithm          | char      |            |             |                  | LINCS, SHAKE                               |
