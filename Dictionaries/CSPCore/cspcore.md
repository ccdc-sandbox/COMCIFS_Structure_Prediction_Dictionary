# CSP Core Dictionary

## Input System Fields

| Data Field                   | Type | Definition                                                                                                                                                                       | Constraints | Units | Example |
|------------------------------|------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------|-------|---------|
| `_csp_chemical_name`         | char | See name_common and name_systematic from Core CIF dictionary.                                                                                                                    |             |       |         |
| `_csp_chemical_formula`      | char | Need to decide which from Core CIF dictionary are appropriate for describing a component and which for describing overall structure (moiety vs sum - also IUPAC and structural). |             |       |         |
| `_csp_chemical_identifier`   | char | InChI / InChI Key / InChI Version.                                                                                                                                               |             |       |         |
| `_csp_chemical_connectivity` | char | Potentially using existing chemical_conn data items. <br> Might also want to capture SMILES or a MOL file representation.                                                        |             |       |         |

## Structure Generation Methods

### General Fields

| Data Field                                          | Type           | Definition                                                                                                                                          | Constraints           | Units               | Example                                                                                                                                         |
|-----------------------------------------------------|----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|---------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| `_csp_structure_generation_space_group_number_list` | char/numb/list | Space group selection could be “all” or a subset (list) specifying which spacegroups were used.                                                     |                       |                     |                                                                                                                                                 |
| `_csp_structure_generation_method`                  | char           | Structure generation method.                                                                                                                        |                       |                     | - Evolutionary Algorithm <br>- Particle Swarm Optimisation <br>- Simulated Annealing <br>- Monte Carlo Parallel tempering <br>- Random Sampling |
| `_csp_structure_generation_software_citation`       | char           | Details of the software used for structure generation.                                                                                              |                       |                     |                                                                                                                                                 |
| `_csp_structure_generation_software_version`        | char           | Version of software used for structure generation.                                                                                                  |                       |                     |                                                                                                                                                 |
| `_csp_density_lower_limit`                          | numb           | Minimum Cell Density                                                                                                                                | 0.0:                  | kg m<sup>-3</sup>   | 800                                                                                                                                             |
| `_csp_density_upper_limit`                          | numb           | Maximum Cell Density                                                                                                                                | 0.0:                  | kg m<sup>-3</sup>   | 1400                                                                                                                                            |
| `_csp_composition_calculation`                      | char           | "fixed" or "variable" composition calculation                                                                                                       | "fixed" or "variable" |                     |                                                                                                                                                 |
| `_csp_group_label`                                  | numb           | Label used to identify a group of atoms                                                                                                             |                       |                     | "MOL1"                                                                                                                                          |
| `_csp_group_atoms`                                  | list           | List of atoms belonging to a specific group                                                                                                         |                       |                     | "C1 C2 H1 H2"                                                                                                                                   |
| `_csp_maximum_number_of_components`                 | numb           | The maximum number of components (atoms or molecules) in the unit cell                                                                              | 1:                    |                     | 4                                                                                                                                               |
| `_csp_minimum_number_of_components`                 | numb           | The minimum number of components (atoms or molecules) in the unit cell                                                                              | 0:                    |                     | 0                                                                                                                                               |
| `_csp_composition_types`                            | list           | List of atomtypes or groups defining the composition                                                                                                |                       |                     | "Mg O" "Benzene water"                                                                                                                          |
| `_csp_composition_coefficients`                     | list           | List of possible compositions for fixed-composition calculations or extremes for variable-composition simulations                                   |                       |                     | "1 1" "2 1"                                                                                                                                     |
| `_csp_stopping_criteria`                            | char/list      | List of rules for stopping the generation of new structures.                                                                                        |                       |                     | "Max Structures", "Low-Energy Structures Unchanged"                                                                                             |
| `_csp_stopping_max_structures_evaluated`            | numb           | The maximum total number of unique crystal structures that will be generated and evaluated during the search.                                       | >0                    |                     | 10000                                                                                                                                           |
| `_csp_stopping_iterations_without_improvement`      | numb           | The maximum number of consecutive iterations (generations, MC steps, etc.) where the global minimum (or the lowest few structures) does not change. | >0                    |                     | 50                                                                                                                                              |
| `_csp_stopping_energy_range`                        | numb           | An energy threshold for the selection of low-energy structures to be considering in the convergence criteria.                                       | >0                    | kJ mol<sup>-1</sup> | 5                                                                                                                                               |
| `_csp_stopping_structures_range`                    | numb           | The number of low-energy structures to be considering in the convergence criteria.                                                                  | >0                    |                     | 1000                                                                                                                                            |

### Evolutionary Algorithms

For these fields to be used, the following conditions must be respected:

1. `_csp_structure_generation_method` "Evolutionary Algorithm"

| Data Field                            | Type | Definition                                                                                                                                                     | Constraints | Units | Example |
|:--------------------------------------|:-----|:---------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------|:------|:--------|
| `_csp_ea_population_size`             | numb | The number of candidate structures in each generation.                                                                                                         | >0          |       | 100     |
| `_csp_ea_initial_population_size`     | numb | The number of candidate structures in the first generation.                                                                                                    | >0          |       | 50      |
| `_csp_ea_number_of_generations`       | numb | The maximum number of evolutionary cycles the algorithm will run before termination (unless other stopping criteria are met).                                  | >0          |       | 50      |
| `_csp_ea_nextgen_structure_selection` | numb | The number of individuals that survives in the next generation.                                                                                                | >1          |       | 5       |
| `_csp_ea_parents_structure_fraction`  | numb | The fraction of individuals in the current population that is used to generate structures in the next cycle.                                                   | 0-1         |       | 0.75    |
| `_csp_ea_mutation_fraction`           | numb | The fraction of individuals in the population that will undergo mutation in each generation.                                                                   | 0-1         |       | 0.2     |
| `_csp_ea_heredity_fraction`           | numb | The fraction of individuals in the population that will be generated through heredity (crossover/recombination) operations between two or more parents.        | 0-1         |       | 0.6     |
| `_csp_ea_permutation_fraction`        | numb | The fraction of individuals in the population that will undergo a permutation operation (e.g., swapping atom positions within a structure) in each generation. | 0-1         |       | 0.1     |

### Particle Swarm Optimisation Algorithms

For these fields to be used, the following conditions must be respected:

1. `_csp_structure_generation_method` "Particle Swarm Optimisation"

| Data Field                       | Type | Definition                                                                                                                                        | Constraints | Units | Example |
|:---------------------------------|:-----|:--------------------------------------------------------------------------------------------------------------------------------------------------|:------------|:------|:--------|
| `_csp_pso_population_size`       | numb | The number of candidate crystal structures (particles) in the swarm.                                                                              | >0          |       | 50      |
| `_csp_pso_number_of_generations` | numb | The maximum number of optimization cycles (generations or iterations) the PSO algorithm will run.                                                 | >0          |       | 100     |
| `_csp_pso_inertia_weight`        | numb | A parameter controlling the contribution of the previous velocity of the particle to its current velocity.                                        | 0-1         |       | 0.7     |
| `_csp_pso_max_inertia_weight`    | numb | If the inertia weight changes with each iteration, this parameter specify the maximum value it can have.                                          | 0-1         |       | 0.9     |
| `_csp_pso_min_inertia_weight`    | numb | If the inertia weight changes with each iteration, this parameter specify the minimum value it can have.                                          | 0-1         |       | 0.4     |
| `_csp_pso_cognitive_coefficient` | numb | A parameter (also called self-confidence factor) controlling the influence of the particle's own best position found so far on its movement.      | >=0         |       | 2       |
| `_csp_pso_social_coefficient`    | numb | A parameter (also called swarm confidence factor) controlling the influence of the swarm's best position found so far on the particle's movement. | >=0         |       | 2       |
| `_csp_pso_velocity_clamp_max`    | numb | The maximum allowed velocity for each dimension if velocity clamping is enabled.                                                                  | >0          |       | 0.2     |

### Simulated Annealing

For these fields to be used, the following conditions must be respected:

1. `_csp_structure_generation_method` "Simulated Annealing"

| Data Field                    | Type | Definition                                                                                       | Constraints | Units | Example |
|:------------------------------|:-----|:-------------------------------------------------------------------------------------------------|:------------|:------|:--------|
| `_csp_sa_initial_temperature` | numb | The starting temperature of the simulated annealing process.                                     | >0          | K     | 500     |
| `_csp_sa_cooling_rate`        | numb | The parameter that determine how the temperature is decreased over the course of the simulation. | 0-1         |       | 0.95    |
| `_csp_sa_number_of_steps`     | numb | The number of attempted structure generation and acceptance steps performed at each temperature. | >0          |       | 10      |

### Monte Carlo Parallel Tempering

For these fields to be used, the following conditions must be respected:

1. `_csp_structure_generation_method` "Monte Carlo Parallel tempering"

| Data Field                     | Type | Definition                                                                                                                                                                                                 | Constraints                   | Units | Example         |
|:-------------------------------|:-----|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------|:------|:----------------|
| `_csp_mcpt_number_of_replicas` | numb | The number of independent Monte Carlo simulations (replicas) running in parallel at different temperatures.                                                                                                | >1                            |       | 3               |
| `_csp_mcpt_temperature_range`  | list | The range of temperatures (minimum and maximum) at which the replicas are run. Temperatures are typically spaced such that there is sufficient overlap in the energy distributions for efficient exchange. | `[min_T >= 0, max_T > min_T]` | K     | `[0, 300, 600]` |
| `_csp_mcpt_number_of_steps`    | numb | The number of Monte Carlo steps performed by each replica at its assigned temperature in each parallel tempering cycle.                                                                                    | >0                            |       | 100             |

### Random Search

For these fields to be used, the following conditions must be respected:

1. `_csp_structure_generation_method` "Random Sampling"

| Data Field                             | Type | Definition                                                                                            | Constraints | Units | Example                  |
|:---------------------------------------|:-----|:------------------------------------------------------------------------------------------------------|:------------|:------|:-------------------------|
| `_csp_random_numbers_algorithm`        | char | Specifies the type of random algorithm used.                                                          | N/A         |       | "Simple Random", "Sobol" |
| `_csp_random_search_number_of_samples` | numb | The total number of unique crystal structures to be generated and evaluated during the random search. | >0          |       | 5000                     |

## Structure Ranking Methods (High-level)

### General Fields

| Data Field                       | Type | Definition                                                                    | Constraints | Units | Example                                                                                      |
|:---------------------------------|:-----|:------------------------------------------------------------------------------|:------------|:------|:---------------------------------------------------------------------------------------------|
| `_csp_ranking_method`            | char | The energy or scoring model used to rank structures.                          |             |       | - Force Field<br>- Semi-empirical<br>- DFT<br>- Wavefunction<br>- AI<br>- Other              |
| `_csp_calculation_type`          | char | Indicates how atomic positions are changed.                                   |             |       | - Optimisation<br>- Ensemble Average<br>- Single point                                       |
| `_csp_optimisation_algorithm`    | char | The algorithm used to optimise atomic positions                               |             |       | - BFGS<br>- L-BFGS<br>- Quasi-Newton<br>- FIRE<br>- Steepest Descent<br>- Conjugate Gradient |
| `_csp_ranking_software_citation` | char | Details of the software used for structure generation.                        |             |       |                                                                                              |
| `_csp_ranking_software_version`  | char | Version of software used for structure generation.                            |             |       |                                                                                              |
| `_csp_ranking_stage`             | numb | In case of multi-step approaches, the stage of the ranking method.            | >=0         |       | 0                                                                                            |
| `_csp_ranking_stage_id`          | char | In case of multi-step approaches, the stage identifier of the ranking method. |             |       | FF, PBE, PBE0                                                                                |

### Periodic Density Functional Theory

For these fields to be used, the following conditions must be respected:

1. `_csp_ranking_method` "DFT"

| Data Field                                      | Type | Definition                                                  | Constraints | Units | Example                                                                                    |
|:------------------------------------------------|:-----|:------------------------------------------------------------|:------------|:------|:-------------------------------------------------------------------------------------------|
| `_csp_dft_exchange_correlation_functional_type` | char | Specifies the type of exchange-correlation functional used. |             |       | LDA, GGA, meta-GGA, Hybrid                                                                 |
| `_csp_dft_exchange_correlation_functional_name` | char | Specifies the name of exchange-correlation functional used. |             |       | PBE, PBE0, SCAN                                                                            |
| `_csp_dft_basis_set_type`                       | char | Defines the type of basis functions/pseudopoentials used.   |             |       | "Plane-waves", "PAW", "Norm-conserving", "Ultrasoft"                                       |
| `_csp_dft_dispersion_correction`                | char | The Van der Waals correction used.                          |             |       | - Grimme-D2<br/>- Grimme-D3<br/>- Tkatchenko-Scheffler<br/>- Many-body dipersion<br/>- XDM |

### Forcefields

For these fields to be used, the following conditions must be respected:

1. `_csp_ranking_method` "Force Field"

| Data Field                               | Type | Definition                                                                                                                     | Constraints | Units | Example                                                                          |
|:-----------------------------------------|:-----|:-------------------------------------------------------------------------------------------------------------------------------|:------------|:------|:---------------------------------------------------------------------------------|
| `_csp_ff_name`                           | char | Name of the force field.                                                                                                       |             |       |                                                                                  |
| `_csp_ff_intramolecular_term`            | char | The energy evaluation method for intra-molecular interactions.                                                                 |             |       | "Bonded Parameters", "Isolated Molecule Energy"                                  |
| `_csp_ff_electrostatic_term`             | char | Functional form of electrostatic interactions                                                                                  |             |       | Point-Charge, Multipoles                                                         |
| `_csp_ff_vdw_term`                       | char | Functional form of van der Waals interactions                                                                                  |             |       | LJ(C6,C12), LJ(epsilon,sigma), Buckingham, ReaxFF Morse-Potential, 14-7 function |
| `_csp_ff_parameterization_method`        | char | Briefly describes the primary method used to derive the force field parameters.                                                | N/A         |       | "Fitting to gas-phase QM data", "Transferable parameters based on atom types"    |
| `_csp_ff_qm_parameterization_functional` | char | The exchange-correlation functional used in the gas-phase quantum mechanical calculations when fitting force field parameters. |             |       | "MP2", "CCSD(T)", "B3LYP"                                                        |
| `_csp_ff_qm_parameterization_basis_set`  | char | The basis set used in the gas-phase quantum mechanical calculations when fitting force field parameters.                       |             |       | "aug-cc-pVTZ", "6-31G(d,p)"                                                      |
