# CSP Core Dictionary

## Current Fields

| Data Field                                        | Type           | Definition                                                                                                                                                                       | Constraints | Units | Example                                                                                                                                                                                                      |
|---------------------------------------------------|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------|-------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| _csp_chemical_name                                | char           | See name_common and name_systematic from Core CIF dictionary.                                                                                                                    |             |       |                                                                                                                                                                                                              |
| _csp_chemical_formula                             | char           | Need to decide which from Core CIF dictionary are appropriate for describing a component and which for describing overall structure (moiety vs sum - also IUPAC and structural). |             |       |                                                                                                                                                                                                              |
| _csp_chemical_identifier                          | char           | InChI / InChI Key / InChI Version.                                                                                                                                               |             |       |                                                                                                                                                                                                              |
| _csp_chemical_connectivity                        | char           | Potentially using existing chemical_conn data items. <br> Might also want to capture SMILES or a MOL file representation.                                                        |             |       |                                                                                                                                                                                                              |
| _csp_structure_generation_space_group_number_list | char/numb/list | Space group selection could be “all” or a subset (list) specifying which spacegroups were used.                                                                                  |             |       |                                                                                                                                                                                                              |
| _csp_structure_generation_method                  | char           | Structure generation method.                                                                                                                                                     |             |       | - Evolutionary algorithm <br> - Random sampling <br> - Simulated annealing <br> - Monte Carlo sampling <br> - Quasi-random sampling <br> - Monte Carlo Parallel tempering <br> - Particle swarm optimisation |
| _csp_software_citation                            | char           | Details of the software used for structure generation.                                                                                                                           |             |       |                                                                                                                                                                                                              |
| _csp_structure_generation_software_version        | char           | Version of software used for structure generation.                                                                                                                               |             |       |                                                                                                                                                                                                              |

## Proposed Fields
### General Inputs
| Data Field                        | Type | Definition                                                                                                        | Constraints           | Units  | Example                |
|-----------------------------------|------|-------------------------------------------------------------------------------------------------------------------|-----------------------|--------|------------------------|
| _csp_density_lower_limit          | numb | Minimum Cell Density                                                                                              | 0.0:                  | kg m-3 | 800                    |
| _csp_density_upper_limit          | numb | Maximum Cell Density                                                                                              | 0.0:                  | kg m-3 | 1400                   |
| _csp_composition_calculation      | char | "fixed" or "variable" composition calculation                                                                     | "fixed" or "variable" |        |                        |
| _csp_group_label                  | numb | Label used to identify a group of atoms                                                                           |                       |        | "MOL1"                 |
| _csp_group_atoms                  | list | List of atoms belonging to a specific group                                                                       |                       |        | "C1 C2 H1 H2"          |
| _csp_maximum_number_of_components | numb | The maximum number of components (atoms or molecules) in the unit cell                                            | 1:                    |        | 4                      |
| _csp_minimum_number_of_components | numb | The minimum number of components (atoms or molecules) in the unit cell                                            | 0:                    |        | 0                      |
| _csp_composition_types            | list | List of atomtypes or groups defining the composition                                                              |                       |        | "Mg O" "Benzene water" |
| _csp_composition_coefficients     | list | List of possible compositions for fixed-composition calculations or extremes for variable-composition simulations |                       |        | "1 1" "2 1"            |

## Evolutionary Algorithms
| Data Field                          | Type | Definition | Constraints | Units | Example |
|-------------------------------------|------|------------|-------------|-------|---------|
| _csp_ea_population_size             | numb |            |             |       |         |
| _csp_ea_initial_population_size     | numb |            |             |       |         |
| _csp_ea_number_of_generations       | numb |            |             |       |         |
| _csp_ea_nextgen_structure_selection | numb |            |             |       |         |
| _csp_ea_parents_structure_selection | numb |            |             |       |         |
| _csp_ea_mutation_fraction           | numb |            |             |       |         |
| _csp_ea_heredity_fraction           | numb |            |             |       |         |
| _csp_ea_permutation_fraction        | numb |            |             |       |         |


