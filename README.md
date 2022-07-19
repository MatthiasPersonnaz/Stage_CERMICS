# pert_BO_stage_CERMICS
Dépôt GitHub pour les codes implémentant la théorie des perturbations avec l'approximation de Born-Oppenheimer lors de mon stage au CERMICS.


# Contenu dans l'ordre anti-chronologique

* `routine_pb_vp_décomp_H2+_rescale.ipynb` : les codes mettant en oeuvre l'approximation de BO avec théorie des perturbations pour $H_2^+$, avec rescaling
* `optim_hardware_graphes.ipynb` : des codes en programmation dynamique, pour optimisez l'usage du GPU dans un cadre plus général que notre problème, en aparté
* `calculs_sagemath_chgt_variables.ipynb` : calcul formel en SageMath pour le changement de variables des coordonnées
* `routine_pb_vp_décomp_H2+.ipynb` : les codes mettant en oeuvre l'approximation de BO avec théorie des perturbations pour $H_2^+$, sans rescaling
* `routine_pb_vp_décomp_osc_harm.ipynb` : théorie des perturbations pour l'oscillateur harmonique quantique 2D*
* `benchmarking_design_fonctions.ipynb` : un ensemble de tests entre les différentes implémentations des algorithmes pour diverses opérations (multiplications matricielles, résolution de systèmes linéaires, etc) et la comparaison CPU/GPU
