# CSP Core Dictionary

## 0. Introduction

### 0.1 Summary

This is the CSP dictionary for describing predicted crystal structures and the methods, parameters and workflows used to
calculate these.

### 0.2 Table of contents

* **1. Input Chemical System** describing the input atoms or molecules for CSP. 
* **2. Structure Generation Methods** describing the workflow used to generate theoretical crystal structures.
* **3. Structure Ranking Methods** describing the energy evaluation models used to generate and rank the structures.
* **4. Output Structure Properties** describing the properties of each output structure, such as their energy or density.
* **5. Conventions** specifying guidelines to avoid multiple labels for the same term.
* **6. Future Developments** describing what is missing from the current dictionary.


## 1. Input Chemical System
This section specifies the atoms used in inorganic CSP or the input molecules for molecular crystal generation.

| Category       | Data Field                     | Type             | Definition                                                                                                        | Constraints           | Units | Example                                         |
|----------------|--------------------------------|------------------|-------------------------------------------------------------------------------------------------------------------|-----------------------|-------|-------------------------------------------------|
| Input          | `name`                         | char             | See name_common and name_systematic from Core CIF dictionary.                                                     |                       |       |                                                 |
| Input          | `composition_calculation`      | char             | "fixed" or "variable" composition calculation                                                                     | "fixed" or "variable" |       |                                                 |
| Input          | `composition_coefficients`     | list             | List of possible compositions for fixed-composition calculations or extremes for variable-composition simulations |                       |       | "1 1" "2 1"                                     |
| Input          | `maximum_number_of_components` | numb             | The maximum number of components (atoms or molecules) in the unit cell                                            | 1:                    |       | 4                                               |
| Input          | `minimum_number_of_components` | numb             | The minimum number of components (atoms or molecules) in the unit cell                                            | 0:                    |       | 0                                               |
| Input Atom     | `types`                        | list             | List of atomic species defining the composition                                                                   |                       |       | "Mg O" (atomic symbols), "12 8" (Atomic number) |
| Input Molecule | `number`                       | char             | Molecule component index.                                                                                         |                       |       |                                                 |
| Input Molecule | `identifier`                   | char             | Label used to identify the molecule.                                                                              |                       |       |                                                 |
| Input Molecule | `smiles`                       | char             | SMILES of the component.                                                                                          |                       |       |                                                 |
| Input Molecule | `molecule_number`              | char             | Molecule component index for each atom.                                                                           |                       |       |                                                 |
| Input Molecule | `molecule_identifier`          | char (Free text) | Label used to identify the molecule for each atom.                                                                |                       |       |                                                 |
| Input Molecule | `atom_label`                   | char (Free text) | Label of atom in the component.                                                                                   |                       |       |                                                 |

Additional details on atoms in molecule ad their connectivity can be specified through the CIF Chemical dictionary,
available at: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/index.html 

### Input Chemical System (old version)

| Category | Data Field     | Type | Definition                                                                                                                                                                       | Constraints | Units | Example |
|----------|----------------|------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------|-------|---------|
| Chemical | `name`         | char | See name_common and name_systematic from Core CIF dictionary.                                                                                                                    |             |       |         |
| Chemical | `formula`      | char | Need to decide which from Core CIF dictionary are appropriate for describing a component and which for describing overall structure (moiety vs sum - also IUPAC and structural). |             |       |         |
| Chemical | `identifier`   | char | InChI / InChI Key / InChI Version.                                                                                                                                               |             |       |         |
| Chemical | `connectivity` | char | Potentially using existing chemical_conn data items. <br> Might also want to capture SMILES or a MOL file representation.                                                        |             |       |         |

## 2. Structure Generation Methods
This section helps delineate the space search range and specify the parameters used for different methods.

### 2.1 General Fields

Category `_csp_structure_generation_[]`: Category for structure generation methods.

| Category             | Data Field                                         | Type                  | Definition                                                                                                                                          | Constraints | Units               | Example                                                                                                                                                                                                                                    |
|----------------------|----------------------------------------------------|-----------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|-------------|---------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Structure Generation | `space_group_number_list`                          | char/numb/list        | Space group selection could be “all” or a subset (list) specifying which spacegroups were used.                                                     |             |                     |                                                                                                                                                                                                                                            |
| Structure Generation | `method`                                           | char/list             | Structure generation method or list of methods.                                                                                                     |             |                     | - Evolutionary Algorithm (Sec. 2.2) <br>- Particle Swarm Optimisation (Sec. 2.3) <br>- Simulated Annealing (Sec. 2.4) <br>- Monte Carlo Parallel tempering (Sec. 2.5) <br>- Random Sampling (Sec. 2.6)<br>- Analogue Templates <br>- Other |
| Structure Generation | `software`                                         | char                  | Name of the software used for structure generation.                                                                                                 |             |                     |                                                                                                                                                                                                                                            |
| Structure Generation | `software_citation`                                | char                  | Details of the software used for structure generation. Either URL to webpage or DOI of the related publication.                                     |             |                     |                                                                                                                                                                                                                                            |
| Structure Generation | `software_version`                                 | char                  | Version of software used for structure generation.                                                                                                  |             |                     |                                                                                                                                                                                                                                            |
| Structure Generation | `density_lower_limit`                              | numb                  | Minimum Cell Density                                                                                                                                | 0.0:        | kg m<sup>-3</sup>   | 800                                                                                                                                                                                                                                        |
| Structure Generation | `density_upper_limit`                              | numb                  | Maximum Cell Density                                                                                                                                | 0.0:        | kg m<sup>-3</sup>   | 1400                                                                                                                                                                                                                                       |
| Structure Generation | `stopping_criteria`                                | char/list (Free text) | List of rules for stopping the generation of new structures.                                                                                        |             |                     | "Max Structures", "Low-Energy Structures Unchanged"                                                                                                                                                                                        |
| Structure Generation | `stopping_criteria_max_structures_evaluated`       | numb                  | The maximum total number of unique crystal structures that will be generated and evaluated during the search.                                       | >0          |                     | 10000                                                                                                                                                                                                                                      |
| Structure Generation | `stopping_criteria_iterations_without_improvement` | numb                  | The maximum number of consecutive iterations (generations, MC steps, etc.) where the global minimum (or the lowest few structures) does not change. | >0          |                     | 50                                                                                                                                                                                                                                         |
| Structure Generation | `stopping_criteria_energy_range`                   | numb                  | An energy threshold for the selection of low-energy structures to be considered in the convergence criteria.                                        | >0          | kJ mol<sup>-1</sup> | 5                                                                                                                                                                                                                                          |
| Structure Generation | `stopping_criteria_structures_range`               | numb                  | The number of low-energy structures to be considered in the convergence criteria.                                                                   | >0          |                     | 1000                                                                                                                                                                                                                                       |

### 2.2 Evolutionary Algorithms

Category `_csp_evolutionary_algorithm_[]`: Subgroup for CSP Structure Generation methods that use Evolutionary
Algorithms. For these fields
to be used, the `_csp_structure_generation_method` must include "Evolutionary Algorithm".

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

### 2.3 Particle Swarm Optimisation Algorithms

Category `_csp_particle_swarm_optimisation_[]`: Subgroup for CSP Structure Generation methods that use Particle Swarm
Optimisation. For these
fields to be used, the `_csp_structure_generation_method` must include "Particle Swarm Optimisation".

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

### 2.4 Simulated Annealing

Category `_csp_simulated_annealing_[]`: Subgroup for CSP Structure Generation methods that use Simulated Annealing. For
these fields to
be used, the `_csp_structure_generation_method` must include "Simulated Annealing".

| Category            | Data Field            | Type | Definition                                                                                       | Constraints | Units | Example |
|:--------------------|:----------------------|:-----|:-------------------------------------------------------------------------------------------------|:------------|:------|:--------|
| Simulated Annealing | `initial_temperature` | numb | The starting temperature of the simulated annealing process.                                     | >0          | K     | 500     |
| Simulated Annealing | `cooling_rate`        | numb | The parameter that determine how the temperature is decreased over the course of the simulation. | 0-1         |       | 0.95    |
| Simulated Annealing | `number_of_steps`     | numb | The number of attempted structure generation and acceptance steps performed at each temperature. | >0          |       | 10      |

### 2.5 Monte Carlo Parallel Tempering

Category `_csp_monte_carlo_parallel_tempering_[]`: Subgroup for CSP Structure Generation methods that use Monte Carlo
Parallel tempering. For
these fields to be used, the `_csp_structure_generation_method` must be set to "Monte Carlo Parallel Tempering".

| Category                       | Data Field           | Type | Definition                                                                                                              | Constraints                   | Units | Example         |
|:-------------------------------|:---------------------|:-----|:------------------------------------------------------------------------------------------------------------------------|:------------------------------|:------|:----------------|
| Monte Carlo Parallel Tempering | `number_of_replicas` | numb | The number of independent Monte Carlo simulations (replicas) running in parallel at different temperatures.             | >1                            |       | 3               |
| Monte Carlo Parallel Tempering | `temperatures_list`  | list | The list of temperatures at which the replicas are run.                                                                 | `[min_T >= 0, max_T > min_T]` | K     | `[0, 300, 600]` |
| Monte Carlo Parallel Tempering | `number_of_steps`    | numb | The number of Monte Carlo steps performed by each replica at its assigned temperature in each parallel tempering cycle. | >0                            |       | 100             |

### 2.6 Random Search

Category `_csp_random_[]`: Subgroup for CSP Structure Generation methods that use Random, Quasi-random algorithms. For
these fields to be used, the `_csp_structure_generation_method` should be set to "Random Sampling".

| Category      | Data Field                 | Type | Definition                                                                                            | Constraints | Units | Example                       |
|:--------------|:---------------------------|:-----|:------------------------------------------------------------------------------------------------------|:------------|:------|:------------------------------|
| Random Search | `random_numbers_algorithm` | char | Specifies the type of random algorithm used.                                                          | N/A         |       | "Pseudorandom", "Quasirandom" |
| Random Search | `number_of_samples`        | numb | The total number of unique crystal structures to be generated and evaluated during the random search. | >0          |       | 5000                          |

## 3. Structure Ranking Methods (High-level)
Within this section, you can define the workflow used to rank the different crystals and give high-level details of the methods used.

### 3.1 General Fields

Category `_csp_structure_ranking_[]`: Category for structure ranking and optimisation methods.

| Category              | Data Field                | Type             | Definition                                                                                                                                                                                                                | Constraints | Units | Example                                                                                                                                                |
|-----------------------|:--------------------------|:-----------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------|:------|:-------------------------------------------------------------------------------------------------------------------------------------------------------|
| Structure Ranking     | `method`                  | char             | The energy or scoring model used to rank structures.                                                                                                                                                                      |             |       | - pDFT (Sec. 3.2) <br>- Forcefield (Sec. 3.3) <br>- Semi-Empirical (Sec. 3.4) <br>- Wavefunction (Sec. 3.5) <br>- ML Potentials (Sec. 3.6) <br>- Other |
| Structure Ranking     | `calculation_type`        | char             | Indicates how atomic positions are changed.                                                                                                                                                                               |             |       | - Optimisation<br>- Ensemble Average<br>- Single point                                                                                                 |
| Structure Ranking     | `software_citation`       | char             | Details of the software used for structure generation.                                                                                                                                                                    |             |       |                                                                                                                                                        |
| Structure Ranking     | `software_version`        | char             | Version of software used for structure generation.                                                                                                                                                                        |             |       |                                                                                                                                                        |
| Structure Ranking     | `stage`                   | numb             | In case of multi-step approaches, the stage of the ranking method.                                                                                                                                                        | >=0         |       | 0                                                                                                                                                      |
| Structure Ranking     | `stage_id`                | char (Free text) | In case of multi-step approaches, the stage identifier of the ranking method.                                                                                                                                             |             |       | FF, PBE, PBE0                                                                                                                                          |
| Geometry Optimisation | `algorithm`               | char             | Geometry optimisation algorithm                                                                                                                                                                                           |             |       | - BFGS<br>- L-BFGS<br>- Quasi-Newton<br>- FIRE<br>- Steepest Descent<br>- Conjugate Gradient                                                           |
| Geometry Optimisation | `cell`                    | char             | It can be "fixed" for no cell optimisation, "isotropic" or "anisotropic" for cell relaxation calculations.                                                                                                                |             |       |                                                                                                                                                        |
| Geometry Optimisation | `atoms`                   | char             | It can be "fixed" for no atoms' position optimisation, "all" for all-atoms geometry optimization, "hydrogens" for optimisation of only H atoms, "non-hydrogens" for non-H atoms or a list of atoms for custom relaxation. |             |       |                                                                                                                                                        |
| Geometry Optimisation | `relax_force_convergence` | numb             | Convergence criteria for stopping the geometry optimisation. Present in TCOD as `_dft_atom_relax_force_conv`.                                                                                                             |             |       |                                                                                                                                                        |

### 3.2 Periodic Density Functional Theory

Category `_dft_[]`: Subgroup for CSP Structure Ranking methods that use pDFT methods (the *p* of *pDFT* is removed in `_dft` for consistency with the TCOD Dictionary). For these fields to be used,
the `_csp_ranking_method` should be set to "pDFT".

| Category | Data Field                             | Type | Definition                                                  | Constraints | Units | Example                                                                                    |
|:---------|:---------------------------------------|:-----|:------------------------------------------------------------|:------------|:------|:-------------------------------------------------------------------------------------------|
| pDFT     | `exchange_correlation_functional_type` | char | Specifies the type of exchange-correlation functional used. |             |       | LDA, GGA, meta-GGA, Hybrid                                                                 |
| pDFT     | `exchange_correlation_functional_name` | char | Specifies the name of exchange-correlation functional used. |             |       | PBE, PBE0, SCAN                                                                            |
| pDFT     | `basis_set_type`                       | char | Defines the type of basis functions/pseudopoentials used.   |             |       | "Plane-waves", "PAW", "Norm-conserving", "Ultrasoft"                                       |
| pDFT     | `dispersion_correction`                | char | The Van der Waals correction used.                          |             |       | - Grimme-D2<br/>- Grimme-D3<br/>- Tkatchenko-Scheffler<br/>- Many-body dipersion<br/>- XDM |

### 3.3 Forcefields

Category `_forcefield_[]`: Subgroup for CSP Structure Ranking methods that use forcefield or mixed inter/intra
molecular
methods. For these fields to be used, the `_csp_ranking_method` should be set to "Forcefield".

| Category   | Data Field                       | Type | Definition                                                                                                                     | Constraints | Units | Example                                                                          |
|------------|:---------------------------------|:-----|:-------------------------------------------------------------------------------------------------------------------------------|:------------|:------|:---------------------------------------------------------------------------------|
| Forcefield | `name`                           | char | Name of the force field.                                                                                                       |             |       |                                                                                  |
| Forcefield | `intramolecular_term`            | char | The energy evaluation method for intramolecular interactions.                                                                  |             |       | "Bonded Parameters", "Isolated Molecule Energy"                                  |
| Forcefield | `electrostatic_term`             | char | Functional form of electrostatic interactions                                                                                  |             |       | Point-Charge, Multipoles                                                         |
| Forcefield | `vdw_term`                       | char | Functional form of van der Waals interactions                                                                                  |             |       | LJ(C6,C12), LJ(epsilon,sigma), Buckingham, ReaxFF Morse-Potential, 14-7 function |
| Forcefield | `parameterization_method`        | char | Briefly describes the primary method used to derive the force field parameters.                                                | N/A         |       | "Fitting to gas-phase QM data", "Transferable parameters based on atom types"    |
| Forcefield | `qm_parameterization_functional` | char | The exchange-correlation functional used in the gas-phase quantum mechanical calculations when fitting force field parameters. |             |       | "MP2", "CCSD(T)", "B3LYP"                                                        |
| Forcefield | `qm_parameterization_basis_set`  | char | The basis set used in the gas-phase quantum mechanical calculations when fitting force field parameters.                       |             |       | "aug-cc-pVTZ", "6-31G(d,p)"                                                      |

### 3.4 Semi-Empirical

Category `_semiempirical_[]`: Subgroup for CSP Structure Ranking methods that use Semi-Empirical methods. For these
fields to
be used, the `_csp_ranking_method` should be set to "Semi-Empirical".

| Category       | Data Field | Type | Definition                                           | Constraints | Units | Example       |
|:---------------|:-----------|:-----|:-----------------------------------------------------|:------------|:------|:--------------|
| Semi-Empirical | `method`   | char | Specifies the name of the Semi-Empirical method used |             |       | AM1, PM3, PM6 |

### 3.5 Wavefunction

Category `_wavefunction_[]`: Subgroup for CSP Structure Ranking methods that use wavefunction methods. For these
fields to be used, the `_csp_ranking_method` should be set to "Wavefunction".

| Category     | Data Field                        | Type | Definition                             | Constraints | Units | Example     |
|:-------------|:----------------------------------|:-----|:---------------------------------------|:------------|:------|:------------|
| Wavefunction | `exchange_correlation_functional` | char | Specifies the name of functional used. |             |       | HF, MP2, CC |
| Wavefunction | `basis_set_type`                  | char | Defines the type of basis used.        |             |       | GTH, NAO    |

### 3.6 ML Potentials

Category `_ml_potential_[]`: Subgroup for CSP Structure Ranking methods that use machine learning potentials methods.
For
these fields to be used, the `_csp_ranking_method` should be set to "ML Potentials".

| Category     | Data Field | Type | Definition                                                                                                                                  | Constraints | Units | Example             |
|:-------------|:-----------|:-----|:--------------------------------------------------------------------------------------------------------------------------------------------|:------------|:------|:--------------------|
| ML Potential | `method`   | char | Specifies the name of the ML Potential used. In case of ML parameterisation of classical forcefields, refer to the Forcefields dictionaries |             |       | ANI, Mace, GAP/SOAP |

### 3.7 Free Energy

Category `_free_energy_[]`: Subgroup for CSP Structure Ranking methods that use free energy methods.

| Category    | Data Field          | Type | Definition                                                          | Constraints | Units | Example           |
|:------------|:--------------------|:-----|:--------------------------------------------------------------------|:------------|:------|:------------------|
| Free Energy | `method`            | char | Specifies the name of the approach used to calculate free energies. |             |       | HA, QHA           |
| Free Energy | `temperatures_list` | list | The temperatures at which free energies are calculated              | >0          | K     | [ 100, 200, 300 ] |

## 4. Single Output Structure Properties

Describes the structure-specific outputs of CSP methods.
Category `_structure_properties_[]`

| Category              | Data Field                | Type | Definition                                                                                                                                | Constraints | Units               | Example         |
|-----------------------|---------------------------|------|-------------------------------------------------------------------------------------------------------------------------------------------|-------------|---------------------|-----------------|
| Theoretical Structure | `temperature`             | numb | The temperature at which the energy and other properties of the theoretical structure were calculated.                                    | \>=0        | K                   | 298.15          |
| Theoretical Structure | `pressure`                | numb | The pressure at which the energy and other properties of the theoretical structure were calculated.                                       |             | Pa                  | 101325.0        |
| Theoretical Structure | `calculated_density`      | numb | The calculated density of the crystal                                                                                                     |             |                     |                 |
| Theoretical Structure | `total_energy`            | numb | The total energy of the theoretical structure, i.e. energy relative to all of the nuclei and electrons seperated to an infinite distance. |             | kJ mol<sup>-1</sup> | -1500.5         |
| Theoretical Structure | `absolute_lattice_energy` | numb | The absolute lattice energy of the crystal, i.e. energy relative to all the molecules seperated to an infinite distance.                  |             | kJ mol<sup>-1</sup> | -1600.8         |
| Theoretical Structure | `absolute_free_energy`    | numb | The absolute free energy of the crystal.                                                                                                  |             | kJ mol<sup>-1</sup> | -1450.2         |
| Theoretical Structure | `free_energy_correction`  | numb | The correction applied to the lattice energy to obtain the free energy, accounting for vibrational and other thermal effects.             |             | kJ mol<sup>-1</sup> | 50.6            |
| Theoretical Structure | `relative_lattice_energy` | numb | The lattice energy of the theoretical structure relative to the lowest energy structure found in the CSP.                                 | \>=0        | kJ mol<sup>-1</sup> | 0.0, 5.2        |
| Theoretical Structure | `energy_uncertainty`      | numb | An estimate of the uncertainty associated with the calculated energy of the theoretical structure.                                        | \>=0        | kJ mol<sup>-1</sup> | 0.1             |
| Theoretical Structure | `score`                   | numb | To allow for methods that may rank by criteria other than energies (e.g., based on stability or other desired properties).                |             |                     | 1, 0.3333, 0.01 |
| Theoretical Structure | `rank`                    | numb | The rank of the structure when ordered by chosen criteria where 1 is considered to be the most favorable or likely structure.             | \>=1        |                     | 2, 7, 12        |

Details on composition, unit cell, symmetry, and atomic coordinates can be specified through the CIF Core dictionary.

## 5. Conventions

A few guidelines are adopted in the description of specific data fields as highlighted in the table below. Except for
*pDFT*, full names are preferred.

| Category             | Data Field | Suggested Input Item   | Alternatives to avoid                                       |
|----------------------|------------|------------------------|-------------------------------------------------------------|
| Structure Generation | `method`   | Random Sampling        | Quasirandom, Pseudorandom (specified in separate datafield) |
| Structure Generation | `method`   | Evolutionary Algorithm | Genetic Algorithm, EA, GA                                   |
| Structure Ranking    | `method`   | Forcefield             | Force Field, Force-Field, FF                                |
| Structure Ranking    | `method`   | pDFT                   | DFT, Density Functional Theory, periodic-DFT                |
| Structure Ranking    | `method`   | Semi-Empirical         | Semi Empirical                                              |

In addition, the Structure Ranking `method` "ML Potential" refers to methods using *ad hoc* descriptors for neural
network training to directly compute energy and forces.
On the other hand, ML models used to parameterize models constants should be classified in the related method.
For example, forcefield constants parametrised with a ML network should be classified as "Forcefield".

## 6. Future Developments

A few areas relevant to CSP have not been explored yet and might be included in later updates of the dictionary.
In general, new or specific methods can use the "Other" option and specify possible publications describing the
workflow.

A list of missing sections is shown below:

* Initial molecule (or list of molecules and conformers) coordinates and properties.
* ML-based Structure Generation methods.
* Clustering algorithms used to remove duplicates.
* While the TCOD dictionary is available for DFT methods and a draft dictionary for forcefield methods is being
  developed, datafields of other energy evaluation methods are limited to basic identification labels. This includes:
  * Semi-Empirical methods
  * ML Potentials
  * Free energy correction methods
* Output structure properties are limited to the energy or score of the crystal. Other measurable properties (the band
  gap for example) are not currently included.
