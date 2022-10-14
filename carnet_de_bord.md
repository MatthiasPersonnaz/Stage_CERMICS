# Carnet de Bord
### Semaine 1: mardi 19 au vendredi 22 avril
- Mardi et mercredi (télétravail)
	* familiarisation avec Julia
	* Création de deux codes pour la simulation 2D d'équations de Schrödinger, schéma en différences finies de Crank-Nicholson avec et sans potentiel extérieur
- Jeudi (présentiel): installation au CERMICS
	* explication du sujet de stage avec les équations de l'hamiltonien
- Vendredi (présentiel): 
	* Lecture du Cohen-Tanoudji sur l'approximation de Born-Oppenheimer, et du poly du CERMICS sur ce sujet pour y parvenir et se familiariser avec les différents types d'équations aux valeurs propres en quantique
	* optimisation du code sur des détails
	* essayé de trouver les DL au propre
	* commencement rédaction rapport au propre avec le début des détails sur Crank-Nicholson
	* essai pour retrouver l'équation découplée comme l'explication mais sans succès
	
### Semaine 2: du lundi 25 au vendredi 29 avril
 - Lundi
	* Lecture de docs: Cohen-Tanoudji, quadrillé le terrain du livre et répérage des chapitres intéressants (cf. références plus bas: théorie des perturbations, etc)
	* en fin de journée repris les équations du sujet de stage et clarification avec le tuteur. Commencement rédaction préliminaire du rapport de stage et repérage des références à inclure (*mostly* problématique, cadre et position du problème), lecture sur la commutativité des opérateurs et la formule de Baker-Campbell-Hausdorff et son utilisation en simulation pour la physique quantique justement pour les opérateurs P et X. Je me suis demandé comment quantifier l'erreur qu'on commet en appliquant la formule aux 2 termes de l'hamiltonien 1D pour la méthode de Trotter utilisée dans le code d'Antoine Levitt: et la meilleure façon d'utiliser le schéma d'Euler (évolution: prendre l'exponentielle de la dérivée ou pas). Optimisation (encore) du code pour le schéma de Crank Nicholson: pour résoudre le laplacien implicite à chaque étape de temps; avec la factorisation de BunchKaufman, il se trouve qu'elle donne des résultats bons (20ms par résolution L\ψ) mais qu'utiliser une matrice creuse donne des résultats aussi bons et moins de code.
 - Mardi (télétravail, tuteur à Jussieu)
	* écriture d'un code pour résoudre le hamiltonien aux v aleurs propres dans le cas d'un état stationnaire en 2D et en 1D. Lecture sur la solution exacte de l'oscillateur harmonique 2D: https://courses.physics.illinois.edu/phys485/fa2015/web/dimensions.pdf 
	* Reste à faire
		+ trouver comment utiliser la fonction csreigvsi avec CUDA pour accélérer les calculs
		+ trouve une démo de la majoration de l'approximation du laplacien à 5 points pour trouver la borne supérieure de la différence
		+ comparer la solution numérique avec les solutions exactes
		+ utiliser Krylov pour les matrices Sparse de l'hamiltonien 2D
 - Mercredi:
	* Explications ?? avec un autre stagiaire
	* Résolution du pb sur la fonction csreigvsi: elle ne fonctionne apparemment que sur le type Float32.
	* Utilisation de Krylovkit pour trouver les premières valeurs propres. ça va plus vite si on ne veut pas toutes les valeurs propres (par exemples les plus petites qui nous intéresseront), mais c'est plus lent pour les avoir toutes
	* Conférence de Laurant Vidal: slides intéressantes pour comprendre le problème
	* Rencontre avec Chiara (étudiante italienne)
 - Jeudi:
	* Optimisation des codes (encore) avec KrylovKit, amélioration de l'utilisation mémoire de la fonction coefs
	* Écriture d'un code pour trouver le plus petit donner le hamiltonien approché (avec les coefs des deux oscillateurs harmoniques découplés et les DS faits automatiquement)
 - Vendredi
	* Optimisation du code pour l'assemblage des matrices avec beaucoup moins de remplissage mémoire
	* Lecture sur la théorie des opérateurs. Je me suis posé la question de savoir pourquoi le spectre d'un certain hamiltonien à potentiel donné était nécessairement discret, et pas continu:  [thread stackoverflow](https://physics.stackexchange.com/questions/65636/discreteness-of-set-of-energy-eigenvalues) en parle. Il s'avère qu'il y a des théorèmes pour ça (la théorie des opérateurs en dimension infinie notamment).
	* Lecture sur la parallélisation en JULIA (https://thomaswiemann.com/assets/teaching/Fall2021-Econ-31720/Econ_31720_discussion_6.pdf et https://julialang.org/blog/2019/07/multithreading/)

### Semaine 3: lundi 02 au vendredi 6 mai
 - Lundi:
	* optimisation du code (encore): résolution du problème de recompilation à chaque exécution en évitant la construction du potentiel par compréhension: j'ai utilisé une meshgrid à la place (2odg plus rapide). Maintenant le code principal pour une grille de 100 par 100 prend environ 0.25 sec et 60Mb d'allocations mémoire.
	* Adaptation du laplacien pour deux masses différentes (deux Δx et Δy différents aussi) mais pas deux N différents.
	* Lecture sur la théorie des perturbations (Cohen-Tanoudji, chapitre 11)
	* écriture des équations et DL poussés au bout. [Intérêt de mettre le terme en E₀(y₀) ?]
- Mardi:
	* Recherche du moyen le plus efficace et le mieux en efficacité / rapidité / simplicité du code pour faire une dérivée approchée
- Mercredi:
	* Finalisation de la fonction calculant la contamination pour la correction à l'ordre 1 d'un état non dégénéré. 
	* En fin de journée, RDV tuteur. On a repris et j'ai mieux compris le sens des équations en physique, ainsi que ce qu'il faut faire précisément avec le hamiltonien selon l'axe Ox, càd trouver numériquement ou analytiquement les états propres sur les deux dimensions, fabriquer un  ψ HBO par produit des deux lois conjointes, et appliquer la théorie des perturbations avec cette solution
	* Séminaire sur les EDP du chercheur de l'EPFL.
	* Reprise de tous les codes pour matcher l'algorithme décrit le jour précédent (ligne d'au dessus). Amélioratio de la précision du modèle avec des différences finies à tous ordres. avec méthode simple de multiplication matricielle
- Jeudi:
	* Continuation du travail précédent. Correction de bugs, et finalisation du code qui trouve maintenant les éléments propres numériquement suivant chaque axe avec les hamiltoniens découplés.
- Vendredi:
	* Calcul (algèbre) des perturbations à l'ordre 22 avec la méthode des perturbations. (ex. 8.2 poly MQ), début du code pour appliquer la méthode des perturbations depuis un produit cartésien des modes propres 1D (avec Itertools), prototypage et recherche doc sur Julia pour faire ça.

### Semaine 4: lundi 09 au vendredi 13 mai
- Lundi: 
	* débbuggage du code pour Born-Oppenheimer. De petites erreurs insidieuses ont été corrigées suite à la reprise du code de la semaine dernière (clarification du sens des équations)
- Mardi: 
	* Lecture de Zur Quantentheorie der Molekln (ainsi que la page wikipédia Approximation de Born-Oppenheimer qui résume bien les idées de l'article)
	* Avancement du rapport sur le sens physique et l'approximation HBO
	* Réflexion sur la meilleur manière de résoudre le pb aux valeurs propres pour contrôler l'erreur: d'abord trouver les valeurs propres, avec Krylov qui est très efficace, surtout pour les premières seulement, puis résoudre le problème restant comme un problème de Dirichlet (avec le laplacien), ce qui a l'avantage de pouvoir donner une majoration de l'erreur. Mais avant ça il faut pouvoir contrôler l'erreur qu'on commet sur le calcul de l'énergie propre en discrétisant !! Donc à creuser
- Mercredi:
	* Lecture sur la convergence, lecture de l'article original de Born et Oppenheimer, un peu de l'article wiki correspondant, j'ai essayé de majorer l'erreur qu'on commet en approchant le hamiltonien 1D mais cela demande de connaître la valeur de phi la solution exacte 1D aux bords donc pas vraiment possible... car on ne compare plus les mêmes solutions (pas les mêmes conditions au bord), dans un cas (celui de l'approximation), on considère que notre vecteur propre psi est nul en dehors d'un segment, dans l'autre cas, ce n'est certainement pas le cas
	* rédaction du rapport
	* codage de l'ordre 2 pour la perturbation
	* Fin de journée: RDV avec le tuteur; points évoqués pour poursuivre:
		- rayon de convergence en λ à étudier
		- norme H¹ **à rajouter**
		- rescaling: l'état devrait se stabiliser
		- trouver le régime où les 2 marchent bien
		- faire tendre la masse vers +∞
		- en 1D: prendre un potentiel de coulomb régularisé en 1/sqrt(x² + α²)
		- estimer l'erreur des vap en fonction des résidus en fct de Hψ et Eψ
		- formule de Kato Temple
		- voir le gap de KT avec le gap numérique de Born-Oppenheimer
		- idée: implémenter des erreurs a posteriori
		- en dim 6: méthodes de tenseurs
		- (mon idée) automatiser le DL ? avec les différents ordres
	* Discussion avec le post doc du fond de la salle, instructive à propos de la recherche.
- Jeudi:
	* Lecture sur Kato-Temple, plusieurs ressources, données par Andrea et les dépendances de celles-ci (les références des articles surtout).
	* Rédaction de la méthode des perturbations sur mon rapport, avec une description de l'algorithme de manière plus systématique que le rapport et un poil plus rigoureux (notamment en prenant un espace hilbertien de dimension finie pour pouvoi contrôler la convergence éventuellement; à ce propos: **penser à rédiger le lien entre espace de dimension finie que l'on pend en discrétsant et l'approximation de la réalité** qui est dans un espace fonctionnel, de dimension infinie, mais qu'on peut approcher cf spectre et contenu des documents pooly CERMICS et livre *Methods of Modern Mathematical Physics*.
	* Commencé la lecture du document de Chiara sur
	* Sont laissés en chantier à l'issue de cette journée:
		- mettre à jour le code de l'ordre 1 et 2, potentiellement l'automatisation
		- faire le calcul de remise à l'échelle et l'idée de regarder quand les deux codes donnent pour l'ordre 1 des résultats similaires, et quand il faut utiliser le code non remis à l'échelle et celui avec
		- regarder les méthodes pour trouver des valeurs propres, les plus adaptées pour les matrices hermitiennes et faire le lien dans le rapport.
		- remettre à jour la rédaction pour le code de la méthode de Crank-Nicholson du début du stage, et décrire les optimisations du code avec les nouvelles connaissances de Julia que j'ai maintenant.
		- mettre la discussion de Chiara dans l'intro de mon rapport, à propos du nombre de variables
		- rajouter les discussion de autoadjoint et l'intérêt en informatique d'utiliser les types hermitian avec les algos pour trouver plus facilement les valeurs propres, ainsi que la preuve argumentée de pourquoi dans le cas d'un potentiel confinant, les valeurs propres sont discrètes
		- dire l'avantage de discrétiser au début de mon algo: déjà on ne peut pas faire autrement car les ordinateurs sont finis et notre espace réel est de dimension infinie, et secondement, et aussi on peut représenter ça dans un espace hermitien $, s'aider de la page avec les limitations dans le cours de Cances
- Vendredi:
	* Finalisation de la thorie des perturbations à tout ordre: écriture de l'algorithme
	
### Semaine 5: lundi 16 mai au vendredi suivant
- Lundi:
	* Matin: rencontre avec Chiara pour lui expliquer le problème, elle avait rien compris
	* Après-midi: codage de la méthode des perturbations à tout ordre
- Mardi:
	* Reprise des API avec CUDA pour résolution de valeurs propres avec CUDA.
	* reste à faire éventuellement:
		- Un diagramme comme pour sklearn i.e. un arbre de décision rapide adapté à ma machine d'après mes tests
	* débuggage du code pour la théorie des perturbations à tout ordre. Il apparaît que avec le développement limité du début en espace, dans la perturbation, les résultats en ordre >= 2 ne sont pas viables
	* écriture rapide d'un code pour tester à tout ordre sans passer par la solution produit
- Mercredi:
	* Pause, travail depuis l'appart, gestion de l'appart Saint-Alexandre, ainsi que la mobilité TUB
- Jeudi et vendredi:
	* rédaction du rapport, réorganisation des parties principalement.
	* continuation et remise au propre du notebook démonstrateur qui calcule la théorie des perturbations en 2D pour démontrer les complications en simulation numérique.
	* Beaucoup de lecture des références aussi ajoutées dans le rapport: les livres Modern Mat
	* réflexions: voir notes papier à propos de la stratégie pour les semaines à venir etc

### Semaine 6: lundi 23 mai au vendredi 27 (2 jours fériés + ponts)
 - Bcp de lecture sur les théorèmes exposés dans la rédaction du rapport.
 - J'ai essayé et partiellement réussi à majorer la norme de la différence entre la solution réelle échantillonnée et la solution numérique.
 - Beaucoup de rédaction du rapport
 - Lecture sur les méthodes de Krylov, le chapitre 6 de Numerical Methods for large eigenvalue problems est une bonne vue d'ensemble
 - Sur la fin de la semaine, lecture de docs sur la résolution de valeurs propres, eigenvalues, eigenvectors, algos et types de matrices, qq pages de lecture de ML, et notamment comment utiliser les fonctions présentes dans http://sites.science.oregonstate.edu/~landaur/nacphy/lapack/eigen.html ?
 - Lecture de notre exam de ML avec l'équivalence du pb lasso / elastic net regression (SD-TSIA 210)
 
### Semaine 7 du 30 mai au 3 juin
 - Lundi : correction de la rédaction de la majoration, principalement
 - Mardi: journée vide et sans motivation/idées de travail, qq corrections au rapport
  - Mercredi: perturbations de W, séries entières: exhibé la condition pour qu'on puisse utiliser la série pour l'énergie avec l'intégrale. A propos de l'algo suggéré par le E Cances: le changement d'échelle nécessite de connaitre la constante de raideur. Qui elle même est estimée : il s'agit de la dérivée seconde en le point minimisant E0. Par conséquent on a un peu un serpent qui se mord la queue.
 - Jeudi: fait l'exercice 7.3.1 du poly de l'X Analyse numérique et optimisation, qui donne la formule générale de l'état fondamental pour un potentiel quadratique. Rédaction de la preuve sur le poly. Question code: mis à jour, tout beau tout propre la routine démonstrateur de la théorie des perturbations avec la somme over states en 2D simple, pour voir si ça marche. Et comme déjà montré, ça ne marche pas super. Ajout des codes de visualisation des termes de contamination. On oberve une converence, même si la meilleure approximation reste à l'ordre 1 à la fois pour l'énergie et le mode.
 - Vendredi: un peu de rédac. Corrigé l'erreur sur le calcul avec le conditionnement dans le rapport: le conditionnement le vérifie pas l'inégalité triangulaire, mais il est multiplicatif --> soulève la question d'estimer le conditionnement 
 - Samedi: Réflexion sur la construction d'un algo global pour donner une estimation de la majoration de l'erreur. Cela nécessite de connaître l'erreur sur l'énergie. Donc de connaître suffisamment précisément l'énergie. Pour le moment je n'ai aucun moyen pour approcher l'énergie dans le domaine continu. Pour utiliser le projecteur du livre, passer par CUDA. On peut connaître le conditionnement  de la matrice H°+λ*W (pour la plus grande valeur propre du moins), et pour la plus petite jsp. Ou bien une factorisation (avec Cholesky puis estimation avec la diagonale cf  https://scicomp.stackexchange.com/questions/32762/cholmod-condition-number-estimate). De +, les tests faits avec CUDA ont montré que 150 points ne sont pas envisageables sur mon GPU (matrice dense 22500x22500), pas assez de mémoire graphique, mais 100 par 100 ok: seulement 10000x10000, et c'est efficace (moins d'une seconde de calcul sans beaucoup d'énergie) contre une vingtaine environ de processeur à 100%). Le produit de matrices creuses d'un certain type n'a pas fonctionné et l'erreur "scalar indexing" ne vient pas de mon ordi: https://github.com/JuliaGPU/CUDA.jl/issues/1113#issuecomment-955759875; https://github.com/Jutho/KrylovKit.jl/issues/15. L'erreur de la cellule qui ne s'arrête pas même conclusion: https://discourse.julialang.org/t/how-to-stop-a-cell-in-vs-code/80010.
 - Dimanche: recherches sur le calcul numérique du résolvent. Trouvé un papier sur le calcul de l'exponentielle: http://www.gipsa-lab.grenoble-inp.fr/~ahmad.hably/Documents/19.pdf (Nineteen Dubious Ways to Compute the Exponential of a Matrix) et https://arxiv.org/pdf/1903.11240.pdf (Eigenvalue and Generalized Eigenvalue Problems: Tutorial) et https://math.stackexchange.com/questions/1039453/resolvent-matrix cette question sur la manière d'inverser le résolvent.
 
### Semaine 8: du mardi 07 juin au vendredi 10 (lundi Pentecôte)
 - Mardi: 
	- normalement conf pour coder efficacement en Julia
	- je me suis rendu compte que ma majoration ne fonctionne pas si bien, car la matrice $H_0+V-E$ est singullière par construction. A voir avec Eric de nouvelles pistes
	- également le calcul de la perturbation avec la projection et automatiser cette tâche: si le livre détaille bien, on n'a que des inverses d'énergies donc facile à inverser. Le pb de l'intégrale rentrera en compte.
	- j'ai résolu la question de la fameuse intégrale dans le bouquin ?
	- j'ai trouvé qu'avec getrf! de https://github.com/JuliaGPU/CUDA.jl/blob/master/test/cusolver, on peut trouver le moyen de calculer $A^{-1}B$ sur GPU (pour les denses)
	- article *How accurate is $inv(A)\times b$ ?*: https://arxiv.org/abs/1201.6035 et ce post de blog http://gregorygundersen.com/blog/2020/12/09/matrix-inversion/
	- trouvé un cours qui explique d'où viennent les formules avec résolvent http://eprints.maths.manchester.ac.uk/1136/1/hht08.pdf et donne une sorte de théorème de Cauchy pour la valeur d'une fonction analytique d'une matrice avec une intégrale de contour
	- The generalized inverse of a singular perturbed matrix: doi:10.1007/bf00944850  
	- Laurent expansion of the inverse of perturbed, singular matrices DOI:10.1016/j.jcp.2015.07.006
	- [Estimates of the determinant of a perturbed identity matrix](https://doi.org/10.1016/j.laa.2018.08.009)
	- Computing $Aα$, $\log(A)$, and related matrix functions by contour integrals DOI:10.1137/070700607
	- A fast contour-integral eigensolver for Non-hermitian matrices
 - Mercredi:
	- rencontre avec Chiara et Eric, il a réexpliqué la technique du projecteur pour la résolution au lieu de la somme over states.
 - Jeudi:
	- après les précisions d'hier j'ai repris les codes avec la méthode originale sujet du stage, i.e. séparation des variables et DSE de V puis perturbation. J'ai 2 problèmes, l'un est un message d'erreur inconnu et ridicule quand je plot les heatmap, l'autre est dû à la grille fixée et finie choisie au préalable, lorsque M tend vers l'infini, le paquet d'onde suivant la direction X (coordonée du noyau) s'écrse jusqu'à une certaine limite puis s'arrête, ce qui est normal car la matrice du laplacien dans le hamiltonien a une limite. Pour contrer ce problème, il faudrait augmenter la taille de la grille. 
	Je me suis aussi rendu compte que fabriquer la projection avec la matrice dense de cette taille est inenvisageable sur CPU car plus de mémoire
 - Vendredi: 
	- codage premier jet des gradients conjugués à l'ordre 1 seulement et avec support gpu. Le résultat est étrange, on voit surtout les blobs 
	Vu qu'on déduit deux vecteurs (un sur chaque sous espace de la somme directe, avec psi0 et son orthogonal, le code est en fait ++ simple moins lourd, mais il faut gérer les vecteurs sur GPU, en FP32 (pour moi) ou FP64 (sur une A100 ou cluster)
	- trouvé le site https://juliateachingctu.github.io/Scientific-Programming-in-Julia/dev/lecture_11/lecture/ comprenant principalement la lecture Julia for programming GPU
	- JULIA Performance et GPU
		* Pour le paquet Krylov: http://juliasmoothoptimizers.github.io/Krylov.jl/dev/gpu/ et https://docs.juliahub.com/Krylov/0fcC3/0.7.0/gpu/
		* https://huijzer.xyz/posts/inference/
		* https://www.juliafordatascience.com/performance-tips/
		* https://gdalle.github.io/JuliaPerf-CERMICS/, https://gdalle.github.io/IntroJulia/
	  	* trouvé le [notebook tuto CUDA pour julia](https://github.com/maleadt/juliacon21-gpu_workshop/blob/main/deep_dive/CUDA.ipynb) (pour la création de kernels, de processus etc)
		* [High-performance computing with Julia and GPU](https://developer.nvidia.com/blog/gpu-computing-julia-programming-language/) explique que quand on utilise le point, julia crée un kernel spécialement optimisé qui parallélise le code.
		* [[++] CUDA programming on Julia part. 1](https://irhum.pubpub.org/pub/julia-cuda/release/4)
		* [An Introduction to GPU Programming in Julia ](https://nextjournal.com/sdanisch/julia-gpu-programming) pour écrire des GPU kernels et les utiliser sur le GPU (exemple avec le calcul d'un ensemble de Julia)
		* [CUDA programming with Julia](https://irhum.github.io/blog/cudajulia/)
		* Un tuto que j'ai utilisé pour faire la fractale ensemble de Julia: [Julia for HPC, GPU Programming](https://enccs.github.io/Julia-for-HPC/GPU/)
		* Voir également cette [page](https://cuda.juliagpu.org/stable/usage/overview/#Kernel-programming-with-@cuda) de la documentation qui explique comment faire des kernels
		* Un post de [blog](https://enccs.se/news/2022/07/julia-for-hpc/) sur Julia pour le HPC
		* https://jenni-westoby.github.io/Julia_GPU_examples/dev/Vector_addition/
		* un article de Nature qui explique en quoi Julia est un super langage: https://www.nature.com/articles/d41586-019-02310-3 https://doi.org/10.1038/d41586-019-02310-3
		* Commandes JULIA LAPACK: `LinearAlgebra.LAPACK.gtsv!`, `LinearAlgebra.LAPACK.potrs!`, `hetri` (pour calculer l'inverse)
### Semaine 9: du lundi 13 juin au vendredi 17
 - Lundi:
	- matin: gestion allemagne, logements, vu pour un passeport, fait le CV europass
	- après-midi: RDV avec Eric, il m'a expliqué comment passer du problème 2D au problème 3D (réduction du problème à 3 corps) et également comment prendre une fonction physiquement acceptable pour les tests
 - Mardi:
	- Continué la rédaction + preuve de la méthode par projection et rédigé
	- sur la fin de journée: lecture sur l'écriture des GPU kernels
 - Mercredi: continué, corrigé des irrégularités et pris la version de l'algo des gradients conjugués  trouvables sur Wikipédia
 - Jeudi:
	- Matin: occupé du workflow TUM (application form et Letter of Acceptance iploadé sur mobility online), occupé de demander à Louis comment éditer les bulletins avec l'école, et les crédits ECTS pour la Deutschlandstipendium, ainsi que écrit des mails à Emannuelle Rivet à Champo pour obtenir les bulletins et Musy pour obtenir d'autres infos par docs officiels, correspondance avec Louis, Neltu, Stefan, GV depuis la veille au soir et Robin Gayet
	- Extraction des nouveles coordonnées en calcul formel
	- Après-midi: rencontre avec Victor 
 - Vendredi:
	- Conf d'Eric au collège de France à 11h puis train
## Semaine 10: du lundi 20 au vendredi 25 juin
- Lundi: paperasse pour Berlin, traduction des bulletins
- Mardi:
	* continué les précédents points
	* trouvé la formule pour l'énergie à l'ordre q.
	* automatisation du calcul.
	* Réfléchi comment optimiser la mémoire sur le GPU (faire les calculs en place autant que possible, etc).
	* Fin d'après-midi puis soirée: rencontre avec Victor pour la deuxième fois
 - Mercredi: télétravail. pas bcp bossé sur le stage, j'ai
	* géré la résidence pour TUB. verdict c'est compliqué. parmi les 3 une est trop loin, une est trop bruyante, et la dernière est pour max un semestre. le Bürgermeister Reuter Stiftung m'a répondu non à l'appart dispo sur le site
	* géré les bulletins de prépas avec Champollion et Mattone. Forwardé à Celine Duchêne à ECL pour qu'elle valide
- Jeudi: au labo. sur la fin de journée, ai été énervé par internet des ponts et le Deutschlandstipendium
- Vendredi: télétravail. En pratique, j'ai géré la paperasse du Deutschlandstipendium à nouveau
## Semaine 11: du lundi 27 juin au 1er juillet
- Lundi:
	* OpenAI me donne l'accès à DALL.E.
	* Reprise de l'algotihme, formalisation codage commencé. Le pb principal reste la gestion de la mémoire entre GPU et CPU et jongler entre eux.
- Mardi.
	* Matin: appelé la licence 3 (Ségolène Serre) pour reconnaissance du diplôme. Je vais m'arranger comme avec Champo, ils sont ok.
	* Am: finalisation de l'algo avec profiling GPU/CPU.
	* Fin am, mail saidi
	* autres jours: vide. en essayant d'optimiser les calculs
## Semaine 12: du lundi 04 juillet au 8 juillet
- Lundi: réflexion pour résoudre le pb de recherche opérationnelle: j'arrive à la conclusion qu'on ne peut PAS le résoudre par une approche d'optimisation en maths pour au moins 2 raisons:
		- On ne peut pas facilement exprimer le coût en transferts mémoire comme une fonction différentiable du paramétrage (par linéaire, ou alors utiliser des ReLU et pas adapté, et donc formulation matricielle impossible)
		- Sans se préoccuper des transferts, exprimer le problème sous forme matricielle comme un pb d'optimisation avec KKT est impossible car contraintes sur le vecteur paramétrage non qualifié (doit prendre des valeurs discrètes donc éventuellement de carré valant 1 mais carré non linéaire comme contrainte égalité, et pas exprimable comme une inégalité simplement).
- Mardi: doc formulation du problème avec un graohe en optimisation combinatoire. Lecture sur les bibliothèques de graphes et parcours associés. Codage de l'algorithme de construction du graphe.
- Mercredi et jeudi: finalisation de l'algo de remplissage du graphe, débugg, apprentissage de la bibliohèque MetaGraphs
- Vendredi: Matin: radio des mains
- Samedi et Dimanche (soirée Lisa) + fini l'algo pour le graphe, finalisé le plot. Tout marche
## Semaine 13: du lundi 11/07 au vendredi 18 juillet
- Lundi: finalisation candidature Stipendium. Tests de l'algo avec des cas concrets pour la démonstration.
- Reste de ala semaine: Rédaction du rapport, fin de l'algorithme, preuve de terminaison, la preuve de correction n'existe pas vraiment, j'explique pourquoi c'est une heuristique. Majoration optimale du nombre de noeuds.
## Semaine 14: ROSCOFF. Arrivée Dimanche.
- Lundi: Discussion avec Paul Cezeaux. CCL: optimiser le parenthésage sur GPU, pour gagner drastiquement en complexité temporelle et mémoire. Pour les gradients conjugués, ne pas assembler la matrice mais utiliser un type de LinearOperators. Mais sur CUDA, cg de IterativeSolvers n'accepte pas les fonctions de  LinearOperators (qi est de plus une licence GPL donc moins libre), et Laurent m'a dit que LinearMaps est mieux et conçu pour s'intégrer bien mieux avec IterativeSolvers dans Julia. Et effectivement, en définissant une telle LinearMap par simple produit vecteur-vecteur, j'arrive à faire fonctionner mes gradients conjugués sur GPU. A retenir pour calculer les produits qui se trouvent dans les formules de projection pour trouver le vecteur propre.
- Reste: discussions avec Eric pour le hamiltonien rescalé
## Semaine 15: Retour Marne
- Lundi: 
- Mardi: Je change pour Arpack.jl pour trouver les vecteurs propres du hamiltonien rescalé, car la convergence avec Krylov ne fonctionne pas. J'ai toujours le bug incompréhensible de l'erreur qui arrive parfois mais pas tout le temps, à un endroit du code (là où on transforme les V0r et V0R en opérateurs d'après les println), où ça devrait normalemenet très bien marcher. Première fois que je tombe sur un pareil mystère.
	* Comparaison entre Arpack et KrylovKit pour les deux opérateurs complets originaux 2D rescalé et non rescalé. Conclusion: dans certains cas Arpack peut être beaucoup plus rapide, je ne sais pas pourquoi. Garder à l'esprit que Arpack utilise une librairie différente externe à Julia
	* Suggestions de choses à faire:
		- faire un algorithme théorique avec la preuve qui trouve la meilleure trajectoire (mais force à résoudre un problème d'optimisation combinatoire à chaque étape.
		- Voir si on peut le faire ave le théorème KKT.
		- voir si pour remplacer krylovkit on peut utiliser Arpack ou les commandes avec CUDA directement (bas niveau) avec les commandes de CUDA.CUSOLVER
		- gérer la figure avec le vecteur propre paramétrisé en u, pourquoi la formule théorique ne fonctionne pas !!!
		- créer la perturbation développée
		- faire une descente de gaadient ou dichotomie pour trouver l'énergie minimale
		- faire l'approche par programmation linéaire pour otimiser le code
- Mercredi: correction des bugs et j'ai enfin obtenu la même énergie pour le hamiltonien 2D rescalé et non rescalé. Le soir: Visio avec Eric: 
	- mettre Psi HBO dans Kato Temple avec lambda 2 - lambda 1 au dénominateur
	- comparer E HBO et fondamental BO  grâce à la simu 2D
- Jeudi: train tôt le matin pour Lyon, j'ai corrigé le point de la todolist qui concernait la décorrélation entre la gaussienne analytique et celle numérique, en fait, le vecteur propre que j'affiche du terme Ĥ⁰u ne correspond à rien. Cela fera partie des variables à virer pour la suite (écrémage du code).
- Vendredi:
	- Victor passe mais je peux pas car déjà à Lyon. Visite chez l'ophtalmo en rentrant à Lyon, je peux rien voir à cause des gouttes, day off, reporté au WE.
	- Les infos pour la chambre à Berlin sont enfin arrivées et dans la semaine j'ai épluché ça entre autres et moments de lecture et détente avec Louis.
- Ce WE: continué le code et la modification.
## Semaine 16: à la maison complètement. Eric a le covid et n'a pas beaucoup suivi de ce fait.
- Lundi et mardi: finalisation; le code marche !!!! pour la première fois, et c'est correct. J'ai fixé l'échelle des u.
- Reste de la semaine: un peu off.
## Semaine 17: dernière semaine. Retour à Marne pour les qq derniers jours.
- Lundi: 
	- J'ai fait la candidature auprès de la bourse AuRA
	- Visio avec Eric à 14h:
		* comment l'énergie de réf se comporte vis-à-vis de la discrétisation
		* comment les approximations se comportent % à la discrétisation
		* --> Pour ces deux dernières questions, il s'agit de faire une campagne sur la discrétisation et de monter pour voir l'impact 
		* peut-on améliorer le dénominateur de KT avec des symétries ? (couplage premiers niveaux)
		* calculer le premier état excité de la référence et regarder si l'approx du omega est bonne ou pas (voir si ω0 ≈ E1-E0)
	- pour moi dans le raport: expliquer qu'on met cette différence dans le dénominateur, et expliquer pourquoi le K ne change pas pour les différentes valeurs de la masse
	- Sur la fin de journée, après la visio, je suis resté dans la bibliothèque et j'ai lu quelques livres en probabilités et machine learning
	- Ne pas oublier de changer le code pour obtenir la dérivée cf tableau 
- Mardi:
	- Codage des recommendations correspondantes faites le jour précédent: supprimé les lignes calculant des variables non utiles (notamment celles non rescalées)
	- Moralité: je me suis rendu compte aussi que les gradients conjugués pour la partie parallèle (perturbation selon la direction de Ψ0 est la partie qui prend la majorité du temps de calcul. Ce qui peut signifier que la matrice n'est pas bien conditionnée (presque nulle car par définition, H0-E0 est identiquement nul sur la direction parallèle, j'ai donc changé en prenant la formule écrite dans sum over states, à confirmer par Eric, car il y avait une grosse erreur de copier-coller dans le code auparavant, à savoir que j'avais appliqué le projecteur orthogonal aux deux accumulateurs (et c'est pas bon, pour le deuxième il fallait hypothétiquement le projecteur parallèle mais ça créait le problème d'identiquement nul suivant la direction qui nous intéressait).
	- Ajouté le vecteur de E2-E1 en sortir pour comparer à ω0
- Mercredi: visio avec Eric à 14h. Conclusions dans l'enregistrement de la conversation. A confirmé que la norme permet de trouver la composante de Psi q sur Psi HBO. Essentiellement, Kato-Temple est trop faible. Explications supplémentaires dans l'enregistrement à compléter
	- u(1-Laplacien)u = H_1^2
	- E-E ref se comporte comme le carré de la norme H1: à vérifier
- Jeudi: Train à 18h.
	- Dans l'après-midi, Eric est venu 10 min constater les résultats. Moralité: on observe peut-être des pentes liées par un facteur 3/2 (ou proche) qui apparaît avec les différentes normes des espaces de Sobolev. La pente de la norme H1 au carré se comporte comme l'erreur sur l'énergie (pas au carré) pour la discrétisation augmentant (théorème pas immédiat, selon ce qu'il a dit). Comme ici on fait varier epsilon à savoir la masse, il n'en est plus très sûr. Egalement, il m'a expliqué "ce qu'on montre c'est que dans l'erreur H1, c'est le terme d'impulsion qui domine et la différence H1 se comporte comme l'énergie"
	- Vu Eric une dernière fois avant de partir pour de bon du labo. Il m'a dit qu'il était très content des codes que j'avais fait et du travail effectué.
	- Retour en train à 18h. Lecture de https://www.daur-rankings.com/blog/diplome-conseil. Voyage plutôt agréable. Avec le soleil couchant, la campagne luit, c'était sympa.
		
	


# Références
+ **TOME I**
	- page 25 (49 du pdf):
		* **Séparation des variables.** Etats stationnaires, décomposition dans une base. Comportement asymptotique d'un état stationnaire de diffusion.
		* Une solution stationnaire correspond à un potentiel indépendant du temps. Un état stationnaire est donc un état d'énergie bien définie.
	- Complément D<sub>II</sub> petit b: observables, recherche des valeurs propres d'un opérateur dans une base. *Lorsque E est de dimension finie, nous avons vu (§ D-1-b) qu’il est toujours possible de former une base avec les vecteurs propres d’un opérateur hermitique. Lorsque E est de dimension infinie, il n’en est plus forcément ainsi. C’est pourquoi il est utile d’introduire une notion nouvelle, celle d’observable.*
	- Complément E<sub>II</sub>: propriétés générales de deux observables vérifiant $[P,Q] = iħ$
	- Complément G<sub>II</sub> page 204 (228 du pdf) **très pertinente**, scindement du problème à deux dimensions avec références
	- Chapitre III: page 216 (240 du pdf): les postulats de la MQ
		* équation de Hamilton-Jacobi
		* Les q<sub>i</sub>(t) et p<sub>i</sub>(t) (i = 1, 2, ..., N) sont appelées variables dynamiques fondamentales. Toutes les grandeurs physiques associées au système (énergie, moment cinétique...) peuvent s’exprimer en fonction de ces variables dynamiques fondamentales.
		* Le P2: Toute grandeur physique mesurable A est décrite par un opérateur A agissant dans E ; cet opérateur est une observable.
	- page 506 (530 du pdf): propriétés générales de l'Hamiltonien quantique
		* Si le potentiel V est minoré alors les valeurs propres de l'hamiltonien quantique (donc l'énergi), sont supérieurs ou égale à V et le spectre est discret (cf aussi complément M3)
	- page 624: référence à l'énergie totale d'un système 
+ **TOME 2**
	- page 1199 (292 du pdf): approximation de Born-Oppenheimer à l'ion H<sub>2</sub><sup>+</sup>. Développement du calcul complet.
	- page 1137 (229 du pdf): perturbations (Complément AXI)
+ **Poly MQ E CANCES**
	- page 39 (51 du pdf): les différents observables courants, leur interprétation, 
	- page 44 (56 du pdf): Time-independent Hamiltonians, avec le théorème d'Ehrenfest
	- page 40 (52 du pdf) Heisenberg Uncertainty principle
	- page 48 (60 du pdf) **Back to the Hermonic Oscillator**
	- page 82 (95 du pdf): potentiel pour deux particules qui interagissent ensemble
	



## Bibliographie consultée pour le moment (littérature)
* Article de Kato sur les fameuses bornes utilisées pour faire notre estimation selon Eric:  On the upper and lower bounds of Matrices eiganvalues (DOI 10.1143/JPSJ.4.334)
* [accès sur Sci-Hub] **The Fourier grid Hamiltonian method for bound state eigenvalues and eigenfunctions, C. Clay Marston and Gabriel G. BalintKurti**
	- *We follow here the philosophy of Dirac's book,9 in that the operators in Eq. (1) act on vectors of an abstract Hilbert space and have not yet been cast into any particular representation.*
* accès [ici](https://arxiv.org/pdf/1804.00434.pdf) **Shortcuts to Adiabaticity Assisted by Counterdiabatic Born-Oppenheimer Dynamics**: décrit le découplage en dynamique pour le hamiltonien prescrit par le sujet de stage
* accès [ici](http://ursula.chem.yale.edu/~batista/classes/v572/fourier_grid_hamiltonian.pdf) The Fourier grid Hamiltonian method for bound state eigenvalues
and eigenfunctions, c. Clay Marston and Gabriel G. Balint-Kurti
	- *"When one has a particular problem to work out in quantum mechanics, one can minimize the labor by using a representation in which the representatives of the most important abstract quantities occurring in the problem are as simple as possible."* P.A.M. Dirac, from The Principles of Quantum Mechanics, 1958.
* [Analyse numérique avancée, cours de l'INRIA:](https://who.rocq.inria.fr/Michel.Kern/Docs/MACS2_AnaNumAv_coursMK.pdf), détaille les méthodes itératives pour résoudre des systèmes linéaires et le conditionnement de la matrice 2D et donc les aspects de convergence
* [Norme et conditionnement d'une matrice](https://www.i2m.univ-amu.fr/perso/thierry.gallouet/licence.d/anum.d/anum-tg2.pdf), cours de niveau L3 très bien fait
* [Numerical solutions of linear eigenvalue problems](https://www.cs.ubc.ca/~greif/Publications/bg2017.pdf) par Jessica Bosch and Chen Greif: PDF épuré et assez didactique
* [Scientific Computing: an introduction survey](https://heath.cs.illinois.edu/scicomp/notes/), cours de l'université de l'Illinois traitant 
    - Scientific Computing
    - Systems of Linear Equations
    - Linear Least Squares
    - Eigenvalue Problems
    - Nonlinear Equations
    - Optimization
    - Interpolation
    - Numerical Integration and Differentiation
    - Initial Value Problems for ODEs
    - Boundary Value Problems for ODEs
    - Partial Differential Equations
    - Fast Fourier Transform
    - Random Numbers and Simulation 
* [Différences finies cours de la Sorbonne](https://www.ljll.math.upmc.fr/~frey/ftp/finite-differences.pdf), TB donne les bornes pour le problème du laplacien 2D, donne le conditionnement et les valeurs propres et l'étude du problème dans le temps mais en 1D seulement et donne le conditionnement en 2D de la matrice du laplacien
* [convergence des différences finies pour la résolution du laplacien 2D](https://giref.ulaval.ca/~fkwok/docs/eig_notes.pdf), complète exactement le cours de Centrale et les valeurs propres de la matrice du laplacien à 5 points en 2D, dans *NUMERICAL METHODS FOR SPECTRAL THEORY* (vient de https://www.ams.org/books/conm/720/conm720-endmatter.pdf)
* Livre *Numerical Methods for large eigenvalues problems*: page 69: Bauer–Fike theorem puis 71/285 section 3.3 avec une série d'inégalités très utiles  pour estimer les valeurs propres dans le chapitre 3: "Perturbation Theory" ainsi que le même Kato-Temple épuré de l'article suivant qui cite ce livre:
* [donne Kato-Temple d'une manière claire et épurée] https://arxiv.org/pdf/1603.06100.pdf
* [The finite difference method](https://www.ljll.math.upmc.fr/frey/cours/UdC/ma691/ma691_ch6.pdf) (avec preuve de la convergence et majoration de l'erreur du laplacien, et citations d'Euler et Laplace) :  ainsi que [Guidelines for the numerical approximation of some models arising in geosciences](https://team.inria.fr/ange/files/2021/09/cours_num_AM_JSM.pdf)
* [A note onf tinite difference methode for solving the eigenvalue problem of second order differential equations](https://www.ams.org/journals/mcom/1962-16-079/S0025-5718-1962-0163434-3/S0025-5718-1962-0163434-3.pdf), par R.M.Osborne
* [Analyse fonctionnelle (cours de l'université d'Evry)](http://www.math-evry.cnrs.fr/_media/members/jmatos/analysefonct1415.pdf)
* [Notes on Quantum Mechanics (Weber state university)](https://physics.weber.edu/schroeder/quantum/QuantumBook.pdf)
* [Explications sur la dérivée fréquentielle](https://scicomp.stackexchange.com/questions/35254/computing-numeric-derivative-via-fft-scipy) et [matecdev](https://www.matecdev.com/posts/julia-fft-derivative.html) et [medium.com](https://medium.com/geekculture/numerical-differentiation-via-the-fft-algorithm-calculating-the-spectral-derivative-a8ab18e1abbe)
* [Notes sur la transformée de Fourier](https://www.iro.umontreal.ca/~mignotte/IFT3205/DemoIFT3205_2_correction.pdf)
* [Julia's Tour on numerical methode](http://www.numerical-tours.com/julia/)
* [NumPy GPU acceleration](https://stsievert.com/blog/2016/07/01/numpy-gpu/)
* [Matrix product complexity on Julia, sujet stack overflow](https://discourse.julialang.org/t/algorithmic-complexity-of-matrix-multiplication-in-julia/57468/8)
* [ASYMMETRY HELPS: EIGENVALUE AND EIGENVECTOR ANALYSES OF ASYMMETRICALLY PERTURBED LOW-RANK MATRICES](https://web.stanford.edu/~chen96/papers/asymmetric_eigens.pdf)
* [4th edition of the Mini-school on mathematics for theoretical chemistry and physics](https://wiki.lct.jussieu.fr/gdrnbody/index.php/MINI-SCHOOL_2022) (organized by GDR NBODY with support from CNRS INC, LCT, LJLL, and ERC EMC2)

### Pour les fractales (testé au début du stage puis réalisé qu'en octobre 2022)
* [Calculer la fractale Julia 'Dragon'](https://personalpages.manchester.ac.uk/staff/yanghong.huang/teaching/MATH4041/DS/juliaset.html)
* http://www.hiddendimension.com/FractalMath/Convergent_Fractals_Main.html
* [une présentation du paquet Fatou.jl](https://crucialflow.com/Fatou.jl/) sur crucialflow.com, qui donne des raccourcis pour tracer plein de figures de systèmes dynamiques


### Reste
* [Différences finies à tous ordres](https://web.njit.edu/~jiang/math712/fornberg.pdf)
* [Sur l'approximation de Born-Oppenheimer](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Quantum_States_of_Atoms_and_Molecules_(Zielinksi_et_al)/10%3A_Theories_of_Electronic_Molecular_Structure/10.01%3A_The_Born-Oppenheimer_Approximation)
* [Approximation de Born-Oppenheimer](https://www.techno-science.net/glossaire-definition/Approximation-de-Born-Oppenheimer.html)
* [Thèse sous la direction de Eric Cances qui parle d'estimation d'erreurs en particulier](https://tel.archives-ouvertes.fr/tel-01689793/document)
* [Introduction to Fourier transform](https://www2.math.upenn.edu/~ccroke/chap4-5.pdf)
* [Solving PDEs with Julia](https://math.stackexchange.com/questions/1657756/is-there-in-discrete-fourier-transform-a-theorem-that-corresponds-to-transforms) (avec des références de bibliothèques Julia)
* [Notes on FFT-based differentiation](https://math.mit.edu/~stevenj/fft-deriv.pdf)
* [La méthode pour dériver avec FFT](https://medium.com/geekculture/numerical-differentiation-via-the-fft-algorithm-calculating-the-spectral-derivative-a8ab18e1abbe) (celui que j'ai donné à Louis) 
* [Solving Schrödinger equation with finite differences method](http://www.scielo.org.mx/pdf/rmfe/v54n2/v54n2a3.pdf)
* [Are there non separable solutions ?](https://physics.stackexchange.com/questions/582447/solutions-of-the-harmonic-oscillator-are-not-always-a-combination-of-separable) 
* [Non separable solutions to the Schrödinger equation sur stackexchange](https://physics.stackexchange.com/questions/255251/non-separable-solutions-of-the-schroedinger-equation )
* [Are there non separable solutions for Schrödinger equation](https://physics.stackexchange.com/questions/190881/analytical-non-separable-solution-for-schrodinger-equation/190896)
* [Solutions of the OHQ3D are not always a combination of separable solutions](https://physics.stackexchange.com/questions/582447/solutions-of-the-harmonic-oscillator-are-not-always-a-combination-of-separable)
* [Why can we assume the solutions are separable](https://physics.stackexchange.com/questions/109750/when-can-we-assume-that-the-wavefunction-is-separable)
* Conseils sur le Machine Learning: https://www.linkedin.com/in/damienbenveniste/recent-activity/shares/
	- [New AI regulations](https://hbr.org/2021/04/new-ai-regulations-are-coming-is-your-organization-ready)
	- [A practical guide to maintaining ML in production](https://eugeneyan.com/writing/practical-guide-to-maintaining-machine-learning/)
	- ai.googleblog.com
* Théorie de Sturm-Liouville
* [Théorie des opérateurs (cours)](http://math.univ-lyon1.fr/~aubrun/enseignement/operateurs/cours.pdf) par Guillaume AUBRUN, université de la Réunion
* [Discreteness of energy eigenvalues](https://physics.stackexchange.com/questions/65636/discreteness-of-set-of-energy-eigenvalues)
* [Create beautiful publication woth plots in Julia](https://nextjournal.com/leandromartinez98/tips-to-create-beautiful-publication-quality-plots-in-julia)
* [A Critical Study of the Finite Difference and Finite Element Methods for the Time Dependent Schrödinger Equation](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.116.5132&rep=rep1&type=pdf) (surtout paragraphe 1.5 page 12 sur les théorèmes importants dont Erehenfest
* [Slides sur le 2D Harmonic Oscillator](https://ocw.nctu.edu.tw/upload/classbfs120904402944051.pdf) qui explique les solutions séparables et le démontre
* [L'approximation de Born-Oppenheimer, Chemistry Lab](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Book%3A_Quantum_States_of_Atoms_and_Molecules_(Zielinksi_et_al)/10%3A_Theories_of_Electronic_Molecular_Structure/10.01%3A_The_Born-Oppenheimer_Approximation)
* [Imaginary time evolution quantum systems](https://physics.stackexchange.com/questions/557225/why-do-we-use-the-imaginary-time-evolution-in-simulations-of-some-quantum-system) explique pourquoi on prend le temps sur l'axe imaginaire
* [Split Operator method avec formule de Campbell](https://www.algorithm-archive.org/contents/split-operator_method/split-operator_method.html)
* Sur le même ton (type tuto python), [solving the Schrödinger equation using the CrankNicholson method](https://artmenlope.github.io/solving-the-2d-schrodinger-equation-using-the-crank-nicolson-method/)
* [The Split-Operator Method (Algorithm's archive, autres sections bien)](https://www.algorithm-archive.org/contents/split-operator_method/split-operator_method.html)
* [The simple Hamrmonic Oscillator](https://faculty.washington.edu/seattle/physics227/reading/reading-22-23.pdf)
* [GPU vs CPU benchmarks with Flux.jl](https://thomaswiemann.com/GPU-vs-CPU-benchmarks-with-Flux.jl)
#### De Towards Data Science
- Data Visualization with Julia and VSCode | by Alan Jones | Towards Data Science
- In-Depth Look At Data Structures In Julia | by Emmett Boudreau | Towards Data Science

