using Logging
Logging.disable_logging(LogLevel(3))

using LinearAlgebra
using SparseArrays
using KrylovKit
using Plots
using Colors
using IterTools
using Polynomials, SpecialPolynomials
using IterativeSolvers
using LinearMaps

@show Threads.nthreads();



function assembly(T,J,N,N²)
    println("assembly"); flush(stdout)
    Λ = spzeros(N²,N²);
    @views for i in 1:N
        Λ[1+(i-1)*N:i*N,1+(i-1)*N:i*N] .= T[:,:]
    end
    
    @views for i in 1:N-1
        Λ[1+(i-1)* N : i    *N, 1+(i)  *N :  (i+1)*N] .= J[:,:] # remplissage blocs etradiag supérieurs
        Λ[1+(i)  * N : (i+1)*N, 1+(i-1)*N :      i*N] .= J[:,:] # remplissage blocs extradiag inférieurs
    end
    return Λ
end

function hamiltonian_2D(δr, δr², δR, δR², N, N², V, m, M)
    println("hamiltonian_2D"); flush(stdout)
    # BLOC EXTRA-DIAGONAL 
    J = -1/(2*M*δR²)*sparse(I,N,N);
    
    # BLOC DIAGONAL
    diag = -(1/(m*δr²)+1/(M*δR²))*ones(Float64,N);  # vecteur de la diagonale RQ: le *2 sur la diagonale du laplacien discret se simplifie avec 1/2m
    extra_diag = 1/(2*m*δr²)*ones(Float64,N-1);     # vecteur de la sur- et sous-diagonale remplies de 1
    T = - SymTridiagonal(diag, extra_diag);    # matrice des blocs diagonaux, qui est tridiagonale symétrique

    # assemblage du laplacien
    Λ = assembly(T,J,N,N²)
    V = Diagonal(reshape(V,N²)); # a priori inutile de convertir en sparse avant

    H = Λ+V;
    H, Λ, V # opérateurs (matrices de taille N²×N²)
end


function hamiltonian_1D(Δd², N, V, m)
    println("hamiltonian_1D"); flush(stdout)
    Λ = -1/(Δd²*2*m)*SymTridiagonal(-2*ones(Float64,N),ones(Float64,N-1));
    V = Diagonal(V); # V_fun.(LinRange(r_min,r_max,N)) pour créer l'argument
    H = Λ + V; 
    H, Λ, V
end



function hamiltonian_2D_rescaled(δr, δr², δu, δu², N, N², V, m, ϵ², K)
    println("hamiltonian_2D_rescaled"); flush(stdout)
    # BLOC EXTRA-DIAGONAL 
    J = -K*ϵ²/(2*δu²)*sparse(I,N,N);
    
    # BLOC DIAGONAL
    diag = -(1/(m*δr²)+K*ϵ²/δu²)*ones(Float64,N); # vecteur de la diagonale RQ: le *2 sur la diagonale du laplacien discret se simplifie avec 1/2m
    extra_diag = 1/(2*m*δr²)*ones(Float64,N-1);     # vecteur de la sur- et sous-diagonale remplies de 1
    T = - SymTridiagonal(diag, extra_diag);    # matrice des blocs diagonaux, qui est tridiagonale symétrique

    # assemblage du laplacien
    Λ = assembly(T,J,N,N²);
    V = Diagonal(reshape(V,N²)); # a priori inutile de convertir en sparse avant

    H = Λ + V;
    H, Λ, V # opérateurs (matrices de taille N²×N²)
end

x = variable(Polynomial{Rational{Int}})
H = [SpecialPolynomials.basis(Hermite, i)(x) for i in 0:3] # /!\ au décalage d'incice
ϕ1Dk = (n,x,m,k) -> (k*m)^(.125)*2^(-n/2)*1/sqrt(factorial(n))*π^(-1/4)*H[n+1]((k*m)^(.25)*x)*exp(-sqrt(k*m)/2*x^2);
ϕ1Du = (n,u,m,k) -> (k*m)^(.125)*2^(-n/2)*1/sqrt(factorial(n))*π^(-1/4)*H[n+1](u)*exp(-1/2*u^2); # u tel que r = ϵu


function get_parameters(r_min, r_max, R_min, R_max, N, m, M)
    println("get_parameters"); flush(stdout)
    δr = (r_max-r_min)/(N-1);
    δR = (R_max-R_min)/(N-1);
    δr² = δr*δr;
    δR² = δR*δR;
    N² = N^2;

    # CONSTRUCTION DE LA MESHGRID
    rs = Vector(r_min:δr:r_max); # sur l'axe donne r ↦ r en direction coordonnée électron
    Rs = Vector(R_min:δR:R_max); # sur l'axe donne R ↦ R en direction coordonnée distance noyau-noyau
    rg = ones(N)' .* rs;         # sur la grille donne   (r,R) ↦ r
    Rg = Rs' .* ones(N);         # sur la grille donne   (r,R) ↦ R
    V = zeros(N,N);              # sur la grille donnera (r,R) ↦ V(r,R) après évaluation ci-dessous

    # CONSTRUCTION DU POTENTIEL ORIGINAL ET DU HAMILTONIEN SUR GRILLE paramétré en R
    V = @. V_nucl_el(rg, Rg) + V_nucl_nucl(Rg); # potentiel d'interaction sur la grille tous phénomènes compris
    Ĥ, Λ, V̂ = hamiltonian_2D(δr, δr², δR, δR², N, N², V, m, M);

    # STRUCTURE DU LAPLACIEN 1D
    LS = SymTridiagonal(-2*ones(Float64,N), ones(Float64,N-1)); 

    # LAPLACIENS SUR AXES
    Λr = -1/(δr²*2*m)*LS;  # laplacien sur l'axe r
    ΛR = -1/(δR²*2*M)*LS;  # laplacien sur l'axe R
    return δr, δR, δr², δR², N², rs, Rs, rg, Rg, V, Ĥ, Λ, V̂, LS, Λr, ΛR
end


function get_lowest_surface_energy(Λr, δR, R_min, rs, N)
    println("get_lowest_surface_energy"); flush(stdout)
    # à changer pour faire une dichotomie ou une biblithèque d'optimisation
    # RECHERCHE DU R₀  minimisant l'énergie de l'état fondamental
    lE₀ = zeros(N);
    Base.Threads.@threads for j in 1:N
        Vx = sparse(Diagonal(Vector(V_nucl_el.(rs,R_min+j*δR)))) # potentiel en x à R=jδR fixé (i.e. à distance noyau-noyau fixé)
        vals, ~, infos = KrylovKit.eigsolve(Λr+Vx, N, 1, :SR, krylovdim=kdim);
        @assert infos.converged > 0;
        lE₀[j]     = infos.converged>=1 ? vals[1] + V_nucl_nucl(R_min + j*δR)  : NaN;
        # on récupère l'énergie  propre du niveau fondamental sur la tranche à R fixé
    end  

    # CALCUL DU R₀ ET DES RAIDEURS
    E₀_at_R₀, ind_R₀ = findmin(lE₀);       # trouver l'énergie de surface minimale
    R₀               = ind_R₀*δR + R_min;  # définir le paramètre nucléaire minimisant l'énergie de surface
    K = 1/(δR)^2 * dot([−1/560 8/315 −1/5 8/5 −205/72 8/5 −1/5 8/315 −1/560], view(lE₀, ind_R₀-4:ind_R₀+4));  # on calcule la dérivée seconde  à l'ordre 8 par rapport à y de E₀ en R₀
    # constante de raideur sur l'axe (Oy) pour le hamiltonien non perturbé
    return lE₀, E₀_at_R₀, ind_R₀, R₀, K
end

function get_rescaling(N)
    println("get_rescaling"); flush(stdout)
    # NOUVEAUX PARAMETRES SUR AXE avec rescaling
    u_min = -3.;     # à fixer de manière à ce que les conditions de Dirichlet soient satisfaites aux bords          
    u_max =  3.;     # idem
    δu = (u_max-u_min)/(N-1);
    δu² = δu^2;
    us = Vector(u_min:δu:u_max);  # sur l'axe donne u ↦ u   
    ug = us' .* ones(N);          # sur la grille donne (r,u) ↦ u
    return u_min, u_max, δu, δu², us, ug
end


function get_parameters_rescaled(u_min, u_max, δu, δu², us, ug, K, M, N, N², LS, V, Rs, rs, rg, Rg, δR, δr, δR², δr², Λ, Λr, ΛR, ind_R₀, R₀)
    println("get_parameters_rescaled"); flush(stdout)
    ϵ = 1/sqrt(sqrt(K*M)); # paramètre de redimensionnement
    ϵ² = ϵ^2;

    Λu = -K*ϵ²/δu²*LS;            # laplacien sur l'axe (Ou)

    # FONCTIONS POTENTIELS HBO NON PERTURBÉS SUR AXES SÉPARÉS  
    V₀rs  = V[:,ind_R₀];          # cf formule 3.19 deuxième ligne du rapport sans la constante .- E₀_at_R₀
    V₀Rs  = .5*K*(Rs.-R₀).^2;     # cf formule 3.19 première ligne du rapport sans la constante .+ E₀_at_R₀
    V₀us  = .5*K*(ϵ*us).^2        # cf formule 3.23 première ligne du rapport sans la constante .+ E₀_at_R₀

    # OPÉRATEURS POTENTIELS HBO NON PERTURBÉS SUR AXES SÉPARÉS 
    V̂⁰rs  = Diagonal(V₀rs);
    V̂⁰Rs  = Diagonal(V₀rs);
    V̂⁰us  = Diagonal(V₀us);

    # HAMILTONIENS NON PERTURBÉS SUR AXES SÉPARÉS
    Ĥ⁰r  = Λr + V̂⁰rs;
    Ĥ⁰R  = ΛR + V̂⁰Rs;

    # OPÉRATEURS RESCALED SUR GRILLE
    V_res = @. V_nucl_el(rg, R₀.+ϵ*ug) + V_nucl_nucl(R₀.+ϵ*ug); # (r,u) ↦ V(r,u) 

    # FONCTION POTENTIEL HBO NON PERTURBÉ SUR GRILLE
    # formule 3.19: somme des deux premières lignes en potentiel:
    # que deux car seulement deux paramétrages (R et u)
    V₀Rg = @. V[:,ind_R₀]*ones(N)'  + .5*K*(Rg.-R₀).^2; # (r,R) ↦ V(r,R₀) + 1/2*(∂²E₀/∂R²)(R₀)(R-R₀)²
    V₀ug = @. V_nucl_el(rg, R₀) + V_nucl_nucl(R₀)  + .5*K*ϵ²*ug.^2;   # (r,u) ↦ V(r,u₀) + Kϵ²/2*(∂²E₀/∂u²)(u₀)(u-u₀)²

    # OPÉRATEUR POTENTIEL NON PERTURBÉ SUR GRILLE
    V̂⁰Rg = Diagonal(reshape(V₀Rg, N²));
    V̂⁰ug = Diagonal(reshape(V₀ug, N²));

    # CONSTRUCTION DU POTENTIEL ET DU HAMILTONIEN NON PERTURBÉS HBO SUR GRILLE
    Ĥ⁰ = Λ + V̂⁰Rg; # hamiltonien HBO non perturbé paramétré en R
    𝔥₀,~,~  = hamiltonian_2D_rescaled(δr, δr², δu, δu², N, N², V₀ug, m, ϵ², K); # hamiltonien HBO non perturbé paramétré en u

    # pour + tard: symétriser un peu la construction de ces variables

    # OPÉRATEUR RESCALÉ ORIGINAL SUR GRILLE
    𝔥,~,~  = hamiltonian_2D_rescaled(δr, δr², δu, δu², N, N², V_res, m, ϵ², K) # 𝔥 : Ψ(r,u) ∈ L²(ℝ^N^2) ↦ -1/2m × ∂²/∂r² -1/2M × ∂²/∂u² + V(r,u) le hamiltonien original complet

    # FONCTION PERTURBATION Vp (développement de Taylor de la perturbation)
    # ordre 1 en espace: (∂V/∂R)(r,R₀)×(R-R₀)
    ∂R_of_V_at_rR₀   = 1/δR*  V[:,ind_R₀-4:ind_R₀+4] * [1/280 −4/105 1/5 −4/5 0. 4/5 −1/5 4/105 −1/280]';       # vecteur, donne r      ↦ ∂V/∂R(r,R₀)
    V₁               = ∂R_of_V_at_rR₀*ones(N)' .* (Rg.-R₀);                                                     # matrice, donne (r, R) ↦ ∂V/∂R(r,R₀)

    # ordre 2 en espace: 1/2×[(∂²V/∂R²)(r,R₀) - d²E₀/dR²(R₀)](R-R₀)² 
    ∂²RR_of_V_at_rR₀ = 1/δR²* V[:,ind_R₀-4:ind_R₀+4] * [−1/560 8/315 −1/5 8/5 −205/72 8/5 −1/5 8/315 −1/560]';  # vecteur, donne r      ↦ ∂²V/∂R²(r,R₀)
    V₂               = .5*(∂²RR_of_V_at_rR₀*ones(N)' .- K) .* (Rg.-R₀).^2;                                      # matrice, donne (r, R) ↦ ∂V/∂R(r,R₀)

    # perturbation totale ordre 1 + ordre 2 en espace (troncature à automatiser éventuellement à tout ordre)
    Vp = V₁ + V₂;

    # OPÉRATEUR PERTURBATION paramétré en R
    Ŵ  = Diagonal(reshape(Vp, N²)); # ψ∈L²(ℝ^N²) ↦ Vp×ψ

    # OPÉRATEURS HAMILTONIENS NON PERTURBÉS RESCALÉS SUR AXES SÉPARÉS
    𝔥u =   K*ϵ²*(-1/2/δu²*LS + 1/2*Diagonal(us.^2));  # 𝔥u : ψ(u) ∈ L²(ℝ^N) ↦ 1/2 × Kϵ²(-∂²/∂u² + u²)ψ pour la solution-produit tensoriel
    # 𝔥u = Λu + V̂⁰us;
    𝔥r = Λr + V̂⁰rs;                                   # 𝔥r : ψ(r) ∈ L²(ℝ^N) ↦ 1/2 ×    (-∂²/∂r² + R²)ψ pour la solution-produit tensoriel

    𝔴  =  @.  V_nucl_el(rg, R₀.+ϵ*ug) + V_nucl_nucl(R₀)  - V_nucl_el(rg, R₀) - V_nucl_nucl(R₀+ϵ*ug) +  1/2*ϵ²*K*ug.^2 # (r,u) ↦ V(r,R₀) - V(r,R₀+ϵu) + Kϵ²/2*u²  (éq. 3.23 du rapport)
    𝔴₁ = ϵ*∂R_of_V_at_rR₀.*ug;
    𝔴₂ = .5*ϵ^2*(∂²RR_of_V_at_rR₀*ones(N)' .- K) .*ug.^2;
    Vp_res = 𝔴₁ + 𝔴₂;
        
    Ŵu = Diagonal(reshape(Vp_res,N²)); # opérateur correspondant à la perturbation paramétrée en u
    return ϵ, ϵ², Λu, V₀rs, V₀rs, V₀us, Ĥ⁰r, Ĥ⁰R, V_res, V₀Rg, V₀ug, V̂⁰Rg, V̂⁰ug, Ĥ⁰, V₁, V₂, Vp, Vp_res, Ŵ, 𝔥, 𝔥u, 𝔥r, 𝔴, Ŵu, 𝔥₀
end



function decompose_hamiltonian_rescaled(r_min, r_max, R_min, R_max, N, m, lM, kdim, Qmax)
    println("decompose_hamiltonian_rescaled"); flush(stdout)
    l = length(lM);
    l_Ψ_pert = zeros(N*N,l);
    l_Ψ_true = zeros(N*N,l);
    l_Ψ_HBO  = zeros(N*N,l);
    l_E_true = zeros(l);
    l_E_pert = zeros(l);
    l_Ψ_err  = zeros(l);
    l_E_err  = zeros(l);

    u_min, u_max, δu, δu², us, ug = get_rescaling(N);

    ############# COMMENCER LA BOUCLE ICI POUR LA MASSE (ce qui précède ne change pas si M change) #############
    ind_M = 1;

    for M in lM
        ### CALCUL DE L'ÉNERGIE DE SURFACE ###
        δr, δR, δr², δR², N², rs, Rs, rg, Rg, V, Ĥ, Λ, V̂, LS, Λr, ΛR = get_parameters(r_min, r_max, R_min, R_max, N,  m, M);
        lE₀, E₀_at_R₀, ind_R₀, R₀, K = get_lowest_surface_energy(Λr, δR, R_min, rs, N);
### CALCUL DES OPÉRATEURS ###
        ϵ, ϵ², Λu, V₀rs, V₀rs, V₀us, Ĥ⁰r, Ĥ⁰R, V_res, V₀Rg, V₀ug, V̂⁰Rg, V̂⁰ug, Ĥ⁰, V₁, V₂, Vp, Vp_res, Ŵ, 𝔥, 𝔥u, 𝔥r, 𝔴, Ŵu, 𝔥₀ = get_parameters_rescaled(u_min, u_max, δu, δu², us, ug, K, M, N, N², LS, V, Rs, rs, rg, Rg, δR, δr, δR², δr², Λ, Λr, ΛR, ind_R₀, R₀);


        ### CALCUL DE LA SOLUTION-PRODUIT HARMONIC-BORN-OPPENHEIMER ###
        lE⁰x, lϕ⁰x, infos_x = KrylovKit.eigsolve(Ĥ⁰r, N, 1, :SR, krylovdim=10); # ou bien avec 𝔥r, car Ĥ⁰r = 𝔥r
        @assert infos_x.converged ≥ 1;
        
        lE⁰u, lϕ⁰u, infos_u = KrylovKit.eigsolve(𝔥u, N, 1, :SR, krylovdim=10);
        @assert infos_u.converged ≥ 1;


        ΨHBO = lϕ⁰x[1] * lϕ⁰u[1]'; # normé car les deux le sont déjà en sortant de Krylov
        ΨHBO = reshape(ΨHBO, N²)
        EHBO = lE⁰x[1] + lE⁰u[1];
        densité_HBO = N^2/(u_max-u_min)/(r_max-r_min);


        l_Ψ_pert[:,ind_M] = ΨHBO;
        l_Ψ_HBO[:,ind_M]  = ΨHBO;
        l_E_pert[ind_M]   = EHBO;


        ### CALCUL DE LA SOLUTION 2D POUR RÉFÉRENCE DU HAMILTONIEN D'INTÉRÊT PARAMÉTRÉ EN u ###
        lE, lϕ, info_2d = KrylovKit.eigsolve(𝔥, N², 1, :SR, krylovdim=40); # KrylovKit.eigsolve plus rapide que Arpack.eigs globalement
        @assert info_2d.converged ≥ 1;


        l_Ψ_true[:,ind_M] = lϕ[1];
        l_E_true[ind_M] = lE[1];
        ### CALCUL DES PERTURBATIONS ###
        Ψ₀ = copy(ΨHBO);
        W =  copy(Ŵu); # W: sparse
        H₀ = copy(𝔥₀); # sparse
        E₀ = EHBO;

        # Qmax   = 3;

        proj = x -> dot(Ψ₀,x)*Ψ₀; # on gagne ~1 prdre de grandeur en temps en utilisant dot au lieu du produit matriciel Ψ₀*(Ψ₀'*x)
        Π_ort  = LinearMap(x -> x - proj(x), N²); # ne pas assembler
        Π_par  = LinearMap(x -> proj(x), N²);

        P_ort  = LinearMap(x -> Π_ort(H₀*Π_ort(x)-E₀*Π_ort(x)), N²); # Π⟂(H₀-E₀)Π⟂
        P_par  = LinearMap(x -> Π_par(H₀*Π_par(x)-E₀*Π_par(x)), N²); # Π∥(H₀-E₀)Π∥

        llΨ    = zeros(Float64, N², Qmax); # création liste des termes d'énergie   en fonction de l'ordre q
        llE    = zeros(Float64, Qmax);      # création liste des termes de vecteurs en fonction de l'ordre q

        Wl     = LinearMap(x -> W*x, N²);
        b      = -Π_ort(W*Ψ₀); # -Π⟂WΨ₀ dans 3.63


        ### GRADIENTS CONJUGUÉS ordre 1 et sauvegarde ###
        llΨ[:,1] = cg(P_ort, b);
        llE[1]    = Ψ₀'*Wl(Ψ₀); # terme d'énergie à l'ordre 1
        
        l_Ψ_pert[:,ind_M] += ϵ*llΨ[:,1];
        l_E_pert[ind_M]   += ϵ*llE[1];



        ### GRADIENTS CONJUGUÉS ordres 2+ et sauvegarde ###
        WlmE₁ = LinearMap(x -> W*x-llE[1]*x, N²);


        R_ort  = -Π_ort*WlmE₁;
        R_par  = -Π_par*WlmE₁;

        acc_b   = zeros(N²);
        acc_ort = zeros(N²);
        acc_par = zeros(N²);



        for q ∈ 2:Qmax
            # calcul énergie ordre q
            llE[q] = Ψ₀'*WlmE₁(llΨ[:,q-1]); # premier terme de 3.65

            for i ∈ 1:(q-2)
                llE[q] = llE[q] - llE[q-i]* Ψ₀'*llΨ[:,i]; # somme du second terme dans 3.65
            end

            # calcul état ordre q
            fill!(acc_b, 0.);
            @views acc_ort[:] = llE[q]*Ψ₀; # dernier terme de la somme de LHS dans 3.66 à i=0
            @views acc_par[:] = llE[q]*Ψ₀; # dernier terme de la somme de LHS dans 3.67 à i=0

            for i ∈ 1:(q-2)
                @. acc_ort[:] = acc_ort + llE[q-i]*llΨ[:,i] # autres termes de la somme dans LHS de 3.66
                @. acc_par[:] = acc_par + llE[q-i]*llΨ[:,i] # autres termes de la somme dans LHS de 3.67
            end
            acc_ort[:] = -Π_ort(WlmE₁(llΨ[:,q-1])) + Π_ort(acc_ort); # LHS de 3.66 complet
            acc_par[:] = -Π_ort(WlmE₁(llΨ[:,q-1])) + Π_ort(acc_par); # LHS de 3.66 complet
            

            acc_ort[:] = cg(P_ort, acc_ort);
            acc_par[:] = cg(P_par, acc_par)
            llΨ[:,q] = Array(acc_ort + acc_ort);
            l_Ψ_pert[:,ind_M] += ϵ^q*llΨ[:,q];
            l_E_pert[ind_M]  += ϵ^q*llE[q];
        end
   
        l_Ψ_err[ind_M] = norm(l_Ψ_pert[:,ind_M] - l_Ψ_true[:,ind_M]);
        l_E_err[ind_M] = abs(l_E_pert[ind_M] - l_E_true[ind_M]);
        ind_M += 1;
    end
    return l_E_pert, l_E_true, l_Ψ_pert, l_Ψ_true, l_Ψ_HBO, l_Ψ_err, l_E_err, # résultats
           δr, δR, δr², δR², N², rs, Rs, rg, Rg, V, Ĥ, Λ, V̂, LS, Λr, ΛR, # paramètres
           u_min, u_max, δu, δu², us, ug # rescaling
end



me = 1; mp = 500;
M=(2*mp^3+mp^2*me)/(2*mp*(me+mp));
m=(2*mp^2*me)/(mp*(2*mp+me)); 
r_min=-5.; r_max=5.; R_min=0.0; R_max=3.5; N=80; ω=1.; kdim=30; # augmenter la dimension Krylov quand N est grand: pour N=200, prendre kdim=50
β=1.5; η=.5; V0=1.5; σ=1.;


function V_nucl_el(r,R)
     return -V0*( exp(-(r-R/2)^2/2/σ^2) + exp(-(r+R/2)^2/2/σ^2) ) # potentiels interaction électron avec les 2 noyaux
end

function V_nucl_nucl(R)
     return + β/sqrt(η^2+R^2) # potentiel interaction des 2 noyaux entre eux
end

Qmax = 2;
lM = [20, 100, 500, 1000, 2000, 3000, 5000];
l_E_pert, l_E_true, l_Ψ_pert, l_Ψ_true, l_Ψ_HBO, l_Ψ_err, l_E_err,
           δr, δR, δr², δR², N², rs, Rs, rg, Rg, V, Ĥ, Λ, V̂, LS, Λr, ΛR, 
           u_min, u_max, δu, δu², us, ug,  = decompose_hamiltonian_rescaled(r_min, r_max, R_min, R_max, N, m, lM, kdim, Qmax);

heatmap(rs, us, reshape(l_Ψ_pert[:,2],N,N)'.^2, xlabel="coordonnée électronique r", ylabel="coordonnée nucléaire R")
heatmap(rs, us, reshape(l_Ψ_true[:,2],N,N)'.^2, xlabel="coordonnée électronique r", ylabel="coordonnée nucléaire R")
heatmap(rs, us, reshape(l_Ψ_HBO[:,2],N,N)'.^2, xlabel="coordonnée électronique r", ylabel="coordonnée nucléaire R")

plot(lM, l_E_err, yaxis=:log, seriestype = :scatter, label="erreur énergie", xlabel="masse M", ylabel="|E - Eₚ|",size=(400,200))
plot(lM, l_Ψ_err, yaxis=:log, seriestype = :scatter, label="résidu", xlabel="masse M", ylabel="|Ψ - Ψₚ|",size=(400,200))