# CSP Core Dictionary

## Input System Fields

| Category | Data Field     | Type | Definition                                                                                                                                                                       | Constraints | Units | Example |
|----------|----------------|------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------|-------|---------|
| Chemical | `name`         | char | See name_common and name_systematic from Core CIF dictionary.                                                                                                                    |             |       |         |
| Chemical | `formula`      | char | Need to decide which from Core CIF dictionary are appropriate for describing a component and which for describing overall structure (moiety vs sum - also IUPAC and structural). |             |       |         |
| Chemical | `identifier`   | char | InChI / InChI Key / InChI Version.                                                                                                                                               |             |       |         |
| Chemical | `connectivity` | char | Potentially using existing chemical_conn data items. <br> Might also want to capture SMILES or a MOL file representation.                                                        |             |       |         |

## Structure Generation Methods

### General Fields
Category `_csp_structure_generation_[]`: Category for structure generation methods.

| Category             | Data Field                                | Type           | Definition                                                                                                                                          | Constraints           | Units               | Example                                                                                                                                         |
|----------------------|-------------------------------------------|----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|---------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| Structure Generation | `space_group_number_list`                 | char/numb/list | Space group selection could be “all” or a subset (list) specifying which spacegroups were used.                                                     |                       |                     |                                                                                                                                                 |
| Structure Generation | `method`                                  | char           | Structure generation method.                                                                                                                        |                       |                     | - Evolutionary Algorithm <br>- Particle Swarm Optimisation <br>- Simulated Annealing <br>- Monte Carlo Parallel tempering <br>- Random Sampling |
| Structure Generation | `software_citation`                       | char           | Details of the software used for structure generation.                                                                                              |                       |                     |                                                                                                                                                 |
| Structure Generation | `software_version`                        | char           | Version of software used for structure generation.                                                                                                  |                       |                     |                                                                                                                                                 |
| Structure Generation | `density_lower_limit`                     | numb           | Minimum Cell Density                                                                                                                                | 0.0:                  | kg m<sup>-3</sup>   | 800                                                                                                                                             |
| Structure Generation | `density_upper_limit`                     | numb           | Maximum Cell Density                                                                                                                                | 0.0:                  | kg m<sup>-3</sup>   | 1400                                                                                                                                            |
| Structure Generation | `composition_calculation`                 | char           | "fixed" or "variable" composition calculation                                                                                                       | "fixed" or "variable" |                     |                                                                                                                                                 |
| Structure Generation | `group_label`                             | numb           | Label used to identify a group of atoms                                                                                                             |                       |                     | "MOL1"                                                                                                                                          |
| Structure Generation | `group_atoms`                             | list           | List of atoms belonging to a specific group                                                                                                         |                       |                     | "C1 C2 H1 H2"                                                                                                                                   |
| Structure Generation | `maximum_number_of_components`            | numb           | The maximum number of components (atoms or molecules) in the unit cell                                                                              | 1:                    |                     | 4                                                                                                                                               |
| Structure Generation | `minimum_number_of_components`            | numb           | The minimum number of components (atoms or molecules) in the unit cell                                                                              | 0:                    |                     | 0                                                                                                                                               |
| Structure Generation | `composition_types`                       | list           | List of atomtypes or groups defining the composition                                                                                                |                       |                     | "Mg O" "Benzene water"                                                                                                                          |
| Structure Generation | `composition_coefficients`                | list           | List of possible compositions for fixed-composition calculations or extremes for variable-composition simulations                                   |                       |                     | "1 1" "2 1"                                                                                                                                     |
| Structure Generation | `stopping_criteria`                       | char/list      | List of rules for stopping the generation of new structures.                                                                                        |                       |                     | "Max Structures", "Low-Energy Structures Unchanged"                                                                                             |
| Structure Generation | `stopping_max_structures_evaluated`       | numb           | The maximum total number of unique crystal structures that will be generated and evaluated during the search.                                       | >0                    |                     | 10000                                                                                                                                           |
| Structure Generation | `stopping_iterations_without_improvement` | numb           | The maximum number of consecutive iterations (generations, MC steps, etc.) where the global minimum (or the lowest few structures) does not change. | >0                    |                     | 50                                                                                                                                              |
| Structure Generation | `stopping_energy_range`                   | numb           | An energy threshold for the selection of low-energy structures to be considering in the convergence criteria.                                       | >0                    | kJ mol<sup>-1</sup> | 5                                                                                                                                               |
| Structure Generation | `stopping_structures_range`               | numb           | The number of low-energy structures to be considering in the convergence criteria.                                                                  | >0                    |                     | 1000                                                                                                                                            |

### Evolutionary Algorithms
Category `_csp_ea_[]`: Sub group for CSP Structure Generation methods that use Evolutionary Algorithms. For these fields to be used, the `_csp_structure_generation_method` must be set to "Evolutionary Algorithm".

| Category                | Data Field                    | Type | Definition                                                                                                                                                     | Constraints | Units | Example |
|:------------------------|:------------------------------|:-----|:---------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------|:------|:--------|
| Evolutionary Algorithms | `population_size`             | numb | The number of candidate structures in each generation.                                                                                                         | >0          |       | 100     |
| Evolutionary Algorithms | `initial_population_size`     | numb | The number of candidate structures in the first generation.                                                                                                    | >0          |       | 50      |
| Evolutionary Algorithms | `number_of_generations`       | numb | The maximum number of evolutionary cycles the algorithm will run before termination (unless other stopping criteria are met).                                  | >0          |       | 50      |
| Evolutionary Algorithms | `nextgen_structure_selection` | numb | The number of individuals that survives in the next generation.                                                                                                | >1          |       | 5       |
| Evolutionary Algorithms | `parents_structure_fraction`  | numb | The fraction of individuals in the current population that is used to generate structures in the next cycle.                                                   | 0-1         |       | 0.75    |
| Evolutionary Algorithms | `mutation_fraction`           | numb | The fraction of individuals in the population that will undergo mutation in each generation.                                                                   | 0-1         |       | 0.2     |
| Evolutionary Algorithms | `heredity_fraction`           | numb | The fraction of individuals in the population that will be generated through heredity (crossover/recombination) operations between two or more parents.        | 0-1         |       | 0.6     |
| Evolutionary Algorithms | `permutation_fraction`        | numb | The fraction of individuals in the population that will undergo a permutation operation (e.g., swapping atom positions within a structure) in each generation. | 0-1         |       | 0.1     |

### Particle Swarm Optimisation Algorithms

Category `_csp_pso_[]`: Sub group for CSP Structure Generation methods that use Particle Swarm Optimisation. For these fields to be used, the `_csp_structure_generation_method` must be set to "Particle Swarm Optimisation".

| Category                    | Data Field              | Type | Definition                                                                                                                                        | Constraints | Units | Example |
|:----------------------------|:------------------------|:-----|:--------------------------------------------------------------------------------------------------------------------------------------------------|:------------|:------|:--------|
| Particle Swarm Optimisation | `population_size`       | numb | The number of candidate crystal structures (particles) in the swarm.                                                                              | >0          |       | 50      |
| Particle Swarm Optimisation | `number_of_generations` | numb | The maximum number of optimization cycles (generations or iterations) the PSO algorithm will run.                                                 | >0          |       | 100     |
| Particle Swarm Optimisation | `inertia_weight`        | numb | A parameter controlling the contribution of the previous velocity of the particle to its current velocity.                                        | 0-1         |       | 0.7     |
| Particle Swarm Optimisation | `max_inertia_weight`    | numb | If the inertia weight changes with each iteration, this parameter specify the maximum value it can have.                                          | 0-1         |       | 0.9     |
| Particle Swarm Optimisation | `min_inertia_weight`    | numb | If the inertia weight changes with each iteration, this parameter specify the minimum value it can have.                                          | 0-1         |       | 0.4     |
| Particle Swarm Optimisation | `cognitive_coefficient` | numb | A parameter (also called self-confidence factor) controlling the influence of the particle's own best position found so far on its movement.      | >=0         |       | 2       |
| Particle Swarm Optimisation | `social_coefficient`    | numb | A parameter (also called swarm confidence factor) controlling the influence of the swarm's best position found so far on the particle's movement. | >=0         |       | 2       |
| Particle Swarm Optimisation | `velocity_clamp_max`    | numb | The maximum allowed velocity for each dimension if velocity clamping is enabled.                                                                  | >0          |       | 0.2     |

### Simulated Annealing

Category `_csp_sa_[]`: Sub group for CSP Structure Generation methods that use Simulated Annealing. For these fields to be used, the `_csp_structure_generation_method` must be set to "Simulated Annealing".

| Category            | Data Field            | Type | Definition                                                                                       | Constraints | Units | Example |
|:--------------------|:----------------------|:-----|:-------------------------------------------------------------------------------------------------|:------------|:------|:--------|
| Simulated Annealing | `initial_temperature` | numb | The starting temperature of the simulated annealing process.                                     | >0          | K     | 500     |
| Simulated Annealing | `cooling_rate`        | numb | The parameter that determine how the temperature is decreased over the course of the simulation. | 0-1         |       | 0.95    |
| Simulated Annealing | `number_of_steps`     | numb | The number of attempted structure generation and acceptance steps performed at each temperature. | >0          |       | 10      |

### Monte Carlo Parallel Tempering

Category `_csp_mcpt_[]`: Sub group for CSP Structure Generation methods that use Monte Carlo Parallel tempering. For these fields to be used, the `_csp_structure_generation_method` must be set to "Monte Carlo Parallel tempering".

| Category                       | Data Field           | Type | Definition                                                                                                              | Constraints                   | Units | Example         |
|:-------------------------------|:---------------------|:-----|:------------------------------------------------------------------------------------------------------------------------|:------------------------------|:------|:----------------|
| Monte Carlo Parallel Tempering | `number_of_replicas` | numb | The number of independent Monte Carlo simulations (replicas) running in parallel at different temperatures.             | >1                            |       | 3               |
| Monte Carlo Parallel Tempering | `temperatures_list`  | list | The list of temperatures at which the replicas are run.                                                                 | `[min_T >= 0, max_T > min_T]` | K     | `[0, 300, 600]` |
| Monte Carlo Parallel Tempering | `number_of_steps`    | numb | The number of Monte Carlo steps performed by each replica at its assigned temperature in each parallel tempering cycle. | >0                            |       | 100             |

### Random Search

Category `_csp_random_[]`: Sub group for CSP Structure Generation methods that use Random, Quasi-random algorithms. For these fields to be used, the `_csp_structure_generation_method` should be set to "Random Sampling".


| Category      | Data Field                 | Type | Definition                                                                                            | Constraints | Units | Example                  |
|:--------------|:---------------------------|:-----|:------------------------------------------------------------------------------------------------------|:------------|:------|:-------------------------|
| Random Search | `random_numbers_algorithm` | char | Specifies the type of random algorithm used.                                                          | N/A         |       | "Simple Random", "Sobol" |
| Random Search | `number_of_samples`        | numb | The total number of unique crystal structures to be generated and evaluated during the random search. | >0          |       | 5000                     |

## Structure Ranking Methods (High-level)

### General Fields
Category `_csp_structure_ranking_[]`: Category for structure ranking methods.         

| Category          | Data Field               | Type | Definition                                                                    | Constraints | Units | Example                                                                                      |
|-------------------|:-------------------------|:-----|:------------------------------------------------------------------------------|:------------|:------|:---------------------------------------------------------------------------------------------|
| Structure Ranking | `method`                 | char | The energy or scoring model used to rank structures.                          |             |       | - Force Field<br>- Semi-empirical<br>- DFT<br>- Wavefunction<br>- AI<br>- Other              |
| Structure Ranking | `calculation_type`       | char | Indicates how atomic positions are changed.                                   |             |       | - Optimisation<br>- Ensemble Average<br>- Single point                                       |
| Structure Ranking | `optimisation_algorithm` | char | The algorithm used to optimise atomic positions                               |             |       | - BFGS<br>- L-BFGS<br>- Quasi-Newton<br>- FIRE<br>- Steepest Descent<br>- Conjugate Gradient |
| Structure Ranking | `software_citation`      | char | Details of the software used for structure generation.                        |             |       |                                                                                              |
| Structure Ranking | `software_version`       | char | Version of software used for structure generation.                            |             |       |                                                                                              |
| Structure Ranking | `stage`                  | numb | In case of multi-step approaches, the stage of the ranking method.            | >=0         |       | 0                                                                                            |
| Structure Ranking | `stage_id`               | char | In case of multi-step approaches, the stage identifier of the ranking method. |             |       | FF, PBE, PBE0                                                                                |

### Periodic Density Functional Theory
Category `_csp_dft_[]`: Sub group for CSP Structure Ranking methods that use DFT methods. For these fields to be used, the `_csp_ranking_method` should be set to "DFT".

| Category | Data Field                             | Type | Definition                                                  | Constraints | Units | Example                                                                                    |
|----------|:---------------------------------------|:-----|:------------------------------------------------------------|:------------|:------|:-------------------------------------------------------------------------------------------|
| DFT      | `exchange_correlation_functional_type` | char | Specifies the type of exchange-correlation functional used. |             |       | LDA, GGA, meta-GGA, Hybrid                                                                 |
| DFT      | `exchange_correlation_functional_name` | char | Specifies the name of exchange-correlation functional used. |             |       | PBE, PBE0, SCAN                                                                            |
| DFT      | `basis_set_type`                       | char | Defines the type of basis functions/pseudopoentials used.   |             |       | "Plane-waves", "PAW", "Norm-conserving", "Ultrasoft"                                       |
| DFT      | `dispersion_correction`                | char | The Van der Waals correction used.                          |             |       | - Grimme-D2<br/>- Grimme-D3<br/>- Tkatchenko-Scheffler<br/>- Many-body dipersion<br/>- XDM |

### Forcefields
Category `_csp_ff_[]`: Sub group for CSP Structure Ranking methods that use forcefield or mixed inter/intra molecular methods. For these fields to be used, the `_csp_ranking_method` should be set to "Force Field".

| Category    | Data Field                       | Type | Definition                                                                                                                     | Constraints | Units | Example                                                                          |
|-------------|:---------------------------------|:-----|:-------------------------------------------------------------------------------------------------------------------------------|:------------|:------|:---------------------------------------------------------------------------------|
| Forcefields | `name`                           | char | Name of the force field.                                                                                                       |             |       |                                                                                  |
| Forcefields | `intramolecular_term`            | char | The energy evaluation method for intra-molecular interactions.                                                                 |             |       | "Bonded Parameters", "Isolated Molecule Energy"                                  |
| Forcefields | `electrostatic_term`             | char | Functional form of electrostatic interactions                                                                                  |             |       | Point-Charge, Multipoles                                                         |
| Forcefields | `vdw_term`                       | char | Functional form of van der Waals interactions                                                                                  |             |       | LJ(C6,C12), LJ(epsilon,sigma), Buckingham, ReaxFF Morse-Potential, 14-7 function |
| Forcefields | `parameterization_method`        | char | Briefly describes the primary method used to derive the force field parameters.                                                | N/A         |       | "Fitting to gas-phase QM data", "Transferable parameters based on atom types"    |
| Forcefields | `qm_parameterization_functional` | char | The exchange-correlation functional used in the gas-phase quantum mechanical calculations when fitting force field parameters. |             |       | "MP2", "CCSD(T)", "B3LYP"                                                        |
| Forcefields | `qm_parameterization_basis_set`  | char | The basis set used in the gas-phase quantum mechanical calculations when fitting force field parameters.                       |             |       | "aug-cc-pVTZ", "6-31G(d,p)"                                                      |
