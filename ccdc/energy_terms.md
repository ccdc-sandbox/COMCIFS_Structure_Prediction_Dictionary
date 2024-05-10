# Energy terms

Throwaway file for discussing what terms to include in the final dictionary

## CCDC blind test dictionary

* There should probably be a conformer energy term: _ccdc_csp_conformer_energy or something, which would be relative to gas phase minima, given theres so much about conformers in the existing dictionary
* _ccdc_csp_free_energy_correction_method: eg. QHA etc. good, keep. "desired"?
* _ccdc_csp_free_energy_correction_software, _ccdc_csp_free_energy_correction_software_version:  good, Optional
* _ccdc_csp_free_energy_correction_stage: preliminary of final - superfluous imo
* _ccdc_csp_classification_energy_absolute: "The total (absolute) energy of the crystal structure, with respect to gas-phase atoms." I think_ccdc_csp_classification_energy_total would line up with the majority of the community (and TCOD). Must be in kJ/mol not eV as suggested, Hartrees would really be the community standard anyway, but kJ/mol is SI. Methods below Periodic DFT don't have this concept, so Optional not Desirable for me
* _ccdc_csp_classification_energy_lattice_relative and _ccdc_csp_classification_ranking: difficulties with this if the landscape is restarted and a new GM is found, I would move from Required to Desirable
* _ccdc_csp_classification_energy_lattice_absolute: only "Required" for me
* _ccdc_csp_classification_energy_error_absolute, _ccdc_csp_classification_energy_error_relative: Did anyone use this? Looks like GRACE would want to in the future
* _ccdc_csp_free_energy_correction_absolute: "correction" should be made more obvious, and total over absolute if we go with my above recommendation
* _ccdc_csp_free_energy_correction_lattice: good, did many use this? I suspect most who do free energy correction will want to just record "free energy", rather than lattice energy and free energy correction. A seperate _ccdc_csp_free_energy_absolute may be wise
* _ccdc_csp_free_energy_relative,_ccdc_csp_free_energy_ranking: good, but not "Required" unless they actually did free energy calcs

## TCOD

to do
