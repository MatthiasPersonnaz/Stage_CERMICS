# Vérification de l'approximation de Born-Oppenheimer dans un cadre simplifié
Dépôt GitHub pour les codes implémentant la théorie des perturbations avec l'approximation de Born-Oppenheimer lors de mon stage au CERMICS.


# Contenu dans l'ordre anti-chronologique

* `vérification_cadre_simplifié_DSE_epsilon.jl` contenant le script mettant en oeuvre le sujet du stage, avec rescaling et un code optimisé et épuré (sans les raviables de contrôle)
* `vérification_cadre_simplifié_DSE_epsilon_heavy.ipynb` et sa version scriptée `.jl` : les codes mettant en oeuvre l'approximation de BO avec théorie des perturbations pour $H_2^+$, avec rescaling et calcul des variables de contrôle (grille, laplaciens, fonctions d'onde non adimensionné(e)s)
* `optim_hardware_graphes.ipynb` : des codes en programmation dynamique, pour optimiser l'usage du GPU dans un cadre plus général que notre problème, en aparté
* `calculs_sagemath_chgt_variables.ipynb` : calcul formel en SageMath pour le changement de variables des coordonnées
* `routine_pb_vp_décomp_H2+.ipynb` : les codes mettant en oeuvre l'approximation de BO avec théorie des perturbations pour $H_2^+$, sans rescaling
* `routine_pb_vp_décomp_osc_harm.ipynb` : théorie des perturbations pour l'oscillateur harmonique quantique 2D*
* `benchmarking_design_fonctions.ipynb` : un ensemble de tests entre les différentes implémentations des algorithmes pour diverses opérations (multiplications matricielles, résolution de systèmes linéaires, etc) et la comparaison CPU/GPU
