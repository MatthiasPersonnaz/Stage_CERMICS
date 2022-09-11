# Stage de recherche au CERMICS

Dépôt GitHub pour les codes implémentant la théorie des perturbations avec l'approximation de Born-Oppenheimer lors de mon stage au CERMICS d'avril à août 2022.
Ce stage avait pour but de tester numériquement dans un cadre simplifié les limites de l’approximation de Born-Oppenheimer utilisée en simulation moléculaire pour découpler la dynamique des électrons de celle des noyaux atomiques.

GitHub repository for the codes implementing the perturbation theory with the Born-Oppenheimer approximation during my internship at CERMICS from April to August 2022.
The aim of this internship was to numerically test in a simplified framework the limits of the Born-Oppenheimer approximation used in molecular simulation to decouple the dynamics of electrons from those of atomic nuclei.


## Contenu dans l'ordre anti-chronologique/Contents in anti-chronological order

| Dossier/Folder | Description/Description |
| ----------- | ----------- |
| `benchmarks_tests_perf` | Contient des codes pour comparer les différentes implémentations de certaines fonctions en Julia pour résoudre différents problèmes d'analyse numérique classiques (inversion de systèmes linéaires, etc) avec la comparaison de performances entre CPG/GPU. Contains codes to compare different implementations of certain functions in Julia to solve different classical numerical analysis problems (inversion of linear systems, etc) with performance comparison between CPG/GPU. |
| `chgt_coord` | Calcul formel pour le changement de coordonnées (reparamétrisation) et la suppression d'un degré de liberté. Computer algebra computations for the change of coordinates (reparametrisation) and the removal of a degree of freedom. |
| `figures_explications_calculs` |  Génération des figures explicatives pour les calculs contenus dans le rapport.  Generation of explanatory figures for the computations contained in the report. |
| `optim_hardware_prog_dyna` | Codes de la partie programmation dynamique, pour optimiser l'usage du GPU dans un cadre plus général que notre problème, ainsi que les figures correspondantes. Codes for the dynamic programming part, to optimise the use of the GPU in a more general framework than our problem, and the corresponding figures.|
| `proper_vérification_H2+` | Codes mettant en oeuvre la théorie des perturbations avec le développement en série entière de la masse et le rescaling + calcul des perturbations par gradients conjugués sur CPU, pour l'ion $H_2^+$. Codes implementing perturbation theory with integer series development of the mass and rescaling + calculation of perturbations by conjugate gradients on CPU, for the $H_2^+$ ion. |
| `routines_parall_H2+_osc-harm` | Codes mettant en oeuvre la théorie des perturbations sur des systèmes simples (oscillateur harmonique quantique 2D) sans rescaling, avec la méthode de la "sum over states", en tant que démonstrateur et tests vis-à-vis de l'orthogonalité des modes propres utilisés. Codes implementing perturbation theory on simple systems (2D harmonic quantum oscillator) without rescaling, using the "sum over states" method, as a demonstrator and test for the orthogonality of the eigen modes used. |
