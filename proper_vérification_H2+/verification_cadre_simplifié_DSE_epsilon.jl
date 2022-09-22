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
using DSP # pour la convolution et le fenêtrage
using FFTW # pour les dérivées spectrales
# using ImageFiltering # pour la convolution qui mime le comportement de matlab
# using Optim
# using Roots

@show Threads.nthreads();



function hamiltonian_1D(Δd², N, V, m)
    println("hamiltonian_1D"); flush(stdout)
    Λ = -1/(Δd²*2*m)*SymTridiagonal(-2*ones(Float64,N),ones(Float64,N-1));
    V = Diagonal(V); # V_fun.(LinRange(r_min,r_max,N)) pour créer l'argument
    H = Λ + V; 
    return H, Λ, V
end



#### séparation des deux contributions du laplacien 2D ####
function laplacian_2D_rescaled_dim_elec(N, N²) # partie du laplacien 2D pour la dimension électronique: manque juste le facteur 1/(m*δr²)
    Λ          = spzeros(N²,N²);
    diag       = 2*ones(Float64,N);        # termes uniquement sur la diagonale du bloc diagonal coefficientés de 2
    extra_diag = -ones(Float64,N-1);       # ainsi que sur l'extra diagonale du bloc diagonal coefficientés de -1
    T          = SymTridiagonal(diag, extra_diag);
    @views for i in 1:N
        Λ[1+(i-1)*N:i*N,1+(i-1)*N:i*N] .= T[:,:]
    end
    return Λ # symtridiagonal
end

function laplacian_2D_rescaled_dim_nucl(N, N²) # partie du laplacien 2D pour la dimension nucléaire: manque juste le facteur K*ϵ²/δu²
    Λ = spzeros(N²,N²)
    T = 2*sparse(I,N,N);     # termes uniquement sur la diagonale du bloc diagonal coefficientés de 2
    J = -sparse(I,N,N);      # ainsi que sur la diagonale du du bloc extra-diagonal coefficientés de -1
    @views for i in 1:N
        Λ[1+(i-1)*N:i*N,1+(i-1)*N:i*N] .= T[:,:]
    end
    @views for i in 1:N-1
        Λ[1+(i-1)* N : i    *N, 1+(i)  *N :  (i+1)*N] .= J[:,:]
        Λ[1+(i)  * N : (i+1)*N, 1+(i-1)*N :      i*N] .= J[:,:]
    end
    return Λ # sparse
end
#############################################################



x = variable(Polynomial{Rational{Int}})
H = [SpecialPolynomials.basis(Hermite, i)(x) for i in 0:3] # /!\ au décalage d'incice
ϕ1Dk = (n,x,m,k) -> (k*m)^(.125)*2^(-n/2)*1/sqrt(factorial(n))*π^(-1/4)*H[n+1]((k*m)^(.25)*x)*exp(-sqrt(k*m)/2*x^2);
ϕ1Du = (n,u,m,k) -> (k*m)^(.125)*2^(-n/2)*1/sqrt(factorial(n))*π^(-1/4)*H[n+1](u)*exp(-1/2*u^2); # u tel que r = ϵu



me = 1; mp = 500; Qmax=1;
M=(2*mp^3+mp^2*me)/(2*mp*(me+mp));
m=(2*mp^2*me)/(mp*(2*mp+me)); 
r_min=-5.; r_max=5.; R_min=0.0; R_max=3.5; N=200; ω=1.;
kdim1d=20; kdim2d = 70;
β=1.5; η=.5; V0=1.5; σ=1.;

# INTERACTION ÉLECTRON-NUCLÉON
function V_noy_el(D)
    return -V0*exp(-D^2/2/σ^2) # potentiel d'interaction entre UN NOYAU et UN ÉLECTRON
end

function ∂V_noy_el(D)
    return V0*D/σ^2*exp(-D^2/2/σ^2) # ∂V/∂D (V nucléon-électron)(D) (pour la descente de gradient et trouver E₀)
end

function ∂2V_noy_el(D)
    return D^3*V0*exp(-1/2*D^2/σ^2)/σ^6 - 3*D*V0*exp(-1/2*D^2/σ^2)/σ^4 # ∂²V/∂D² (V nucléon-électron)(D) (pour ∂²E₀/∂R²)
end

# INTERACTION ÉLECTRON-2NOYAUX
function V_nucl_el(r,R)
    return V_noy_el(r-R/2) + V_noy_el(r+R/2) # potentiel d'interaction de l'ÉLECTRON avec les DEUX NOYAUX
end

function ∂R_V_nucl_el(r,R)
    return -.5*∂V_noy_el(r-R/2)   + .5*∂V_noy_el(r+R/2) # ∂V/∂R (V 2noyaux-électrons)(R)
end

function ∂2R_V_nucl_el(r,R)
    return +.25*∂2V_noy_el(r-R/2) + .25*∂2V_noy_el(r+R/2) # ∂²V/∂R² (V 2noyaux-électrons)(R)
end


# INTERACTION NUCLÉON-NUCLÉON
function V_nucl_nucl(R)
    return β/sqrt(η^2+R^2) # potentiel d'interaction des DEUX NOYAUX ENTRE EUX
end

function ∂V_nucl_nucl(R)
    return -β*R/((η^2+R^2)*sqrt(η^2+R^2)) # ∂V/∂R (V nucl-nucl)(R)
end

function ∂2V_nucl_nucl(R)
    return 3*R^2*β/(R^2 + η^2)^2.5 - β/(R^2 + η^2)^1.5 # ∂²V/∂R² (V nucl-nucl)(R)
end




function get_rescaling(N)
    println("calcul du rescaling"); flush(stdout)
    # NOUVEAUX PARAMETRES SUR AXE avec rescaling
    u_min = -2.65;     # à fixer de manière à ce que les conditions de Dirichlet soient satisfaites aux bords          
    u_max = +2.65;     # idem
    δu = (u_max-u_min)/(N-1);
    δu² = δu^2;
    us = Vector(u_min:δu:u_max);  # sur l'axe donne u ↦ u   
    ug = us' .* ones(N);          # sur la grille donne (r,u) ↦ u
    return u_min, u_max, δu, δu², us, ug
end


function get_params_grille(r_max, r_min, R_max, R_min, N)
    println("calcul paramètres de grille"); flush(stdout)
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
    return δr, δR, δr², δR², N², rs, Rs, rg, Rg, V
end


function get_lowest_surface_energy_breadth_first(Λr, δR, R_min, rs, Rs, N, kdim1d)
    println("calcul du minimum de l'énergie de surface par sondage"); flush(stdout)
    # à changer pour faire une dichotomie ou une biblithèque d'optimisation
    # RECHERCHE DU R₀  minimisant l'énergie de l'état fondamental
    lE₀ = zeros(N);
    Base.Threads.@threads for j in 1:N
        Vx = sparse(Diagonal(Vector(V_nucl_el.(rs,Rs[j])))) # potentiel en x à R=jδR fixé (i.e. à distance noyau-noyau fixé)
        vals, ~, infos = KrylovKit.eigsolve(Λr+Vx, N, 1, :SR, krylovdim=kdim1d, issymmetric=true, ishermitian=true);
        @assert infos.converged > 0;
        lE₀[j]     = infos.converged>=1 ? vals[1] + V_nucl_nucl(Rs[j])  : NaN;
        # on récupère l'énergie  propre du niveau fondamental sur la tranche à R fixé
    end  
    # plot de l'énergie pour le rapport
    plot(Rs, lE₀, xlabel="R", ylabel="E₀(R)")
    savefig(raw"./plot_energies.pdf") # dans C:/Users/matth
    
    # CALCUL DU R₀ ET DE LA DÉRIVÉE SECONDE (RAIDEUR) DE L'ÉNERGIE avec les différences finies seulement
    E₀_at_R₀, ind_R₀ = findmin(lE₀);       # trouver l'énergie de surface minimale
    R₀               = ind_R₀*δR + R_min;  # définir le paramètre nucléaire minimisant l'énergie de surface
    diff2_ctr_o8 = [−1/560 8/315 −1/5 8/5 −205/72 8/5 −1/5 8/315 −1/560];
    K = 1/(δR)^2 * dot(diff2_ctr_o8, view(lE₀, ind_R₀-4:ind_R₀+4));  # on calcule la dérivée seconde  à l'ordre 8 par rapport à y de E₀ en R₀
    # constante de raideur sur l'axe (Oy) pour le hamiltonien non perturbé
    return lE₀, E₀_at_R₀, ind_R₀, R₀, K
end



function get_lowest_surface_energy_secant(Λr, δR, R_min, rs, N, R2, R1, ϵ, kdim1d, itermax)
    println("calcul du minimum de l'énergie de surface par l'algorithme des sécantes"); flush(stdout)
    Hr = copy(Λr);
    ∂Hr = Diagonal(zeros(N));
    # à l'avenir transformer le ∂Hr en un simple vecteur, le type sert ici uniquement à l'analogie avec l'objet mathématique
    iter_nb = 0;
    dE_rp = ϵ+1; # initialisation pour rentrer dans la boucle
    dE_lp = 0;
    R_rp = R1;
    R_lp = R2;
    while (iter_nb < itermax && abs(dE_rp) > ϵ)
        # calcul des états
        Hr[diagind(Hr)] = diag(Λr) + Vector(V_nucl_el.(rs, R_rp)) .+ V_nucl_nucl(R_rp);
        vals, vecs, infos = KrylovKit.eigsolve(Hr, N, 1, :SR, krylovdim=kdim1d, issymmetric=true, ishermitian=true);
        @assert infos.converged > 0;
        ψ_rp = Vector(vecs[1]);
        Hr[diagind(Hr)] = diag(Λr) + Vector(V_nucl_el.(rs, R_lp)) .+ V_nucl_nucl(R_rp);
        vals, vecs, infos = KrylovKit.eigsolve(Hr, N, 1, :SR, krylovdim=kdim1d, issymmetric=true, ishermitian=true);
        @assert infos.converged > 0;
        ψ_lp = Vector(vecs[1]);

        # calcul des dérivées
        ∂Hr[diagind(∂Hr)] = ∂R_V_nucl_el.(rs, R_rp) .+ ∂V_nucl_nucl(R_rp); # 1/2*(∂V_noy_el.(rs.+R_rp/2)-∂V_noy_el.(rs.-R_rp/2));
        dE_rp = dot(ψ_rp,∂Hr,ψ_rp); 
        ∂Hr[diagind(∂Hr)] = ∂R_V_nucl_el.(rs, R_lp) .+ ∂V_nucl_nucl(R_lp);
        dE_lp = dot(ψ_lp,∂Hr,ψ_lp);

        # méthode des sécantes pour la dérivée de l'énergie par rapport à R
        R_rp_mem = R_rp;
        step = dE_rp*(R_rp-R_lp)/(dE_rp-dE_lp);
        R_rp -= step;
        R_lp = R_rp_mem;
        iter_nb += 1;
        @show iter_nb, step, dE_rp; # montrer la progression lors de l'appel
    end
    
    ### CALCUL DE LA DÉRIVÉE SECONDE (RAIDEUR) DE L'ÉNERGIE
    # avec les différences finies
    Hr = copy(Λr);
    E_décal = zeros(9);
    ind_décal = Vector(-4:4);
    ψ_at_R₀ = zeros(N);
    for (i,m) in enumerate(ind_décal)
        Hr[diagind(Hr)] = diag(Λr) + Vector(V_nucl_el.(rs, R_rp+i*δR) .+ V_nucl_nucl(R_rp+i*δR)); 
        vals, vecs, infos = KrylovKit.eigsolve(Hr, N, 1, :SR, krylovdim=kdim1d, issymmetric=true, ishermitian=true);
        if i == 0
            ψ_at_R₀ = vecs[1]; # obtenir l'état en R₀ pour le calcul de la dérivée seconde de E(R) en R₀
        end
        @assert infos.converged ≥ 1;
        E_décal[i] = vals[1];
    end
    diff2_ctr_o8 = [−1/560 8/315 −1/5 8/5 −205/72 8/5 −1/5 8/315 −1/560];
    d2EdR2_fd = 1/(δR)^2 * dot(diff2_ctr_o8, view(E_décal, :));

    # avec Ehrenfest modifié (sur GPU ?)
    ∂Hr[diagind(∂Hr)] = ∂R_V_nucl_el.(rs, R_rp) .+ ∂V_nucl_nucl(R_rp);
    Hr[diagind(Hr)] = diag(Λr) + Vector(V_nucl_el.(rs, R_rp) .+ V_nucl_nucl(R_rp)); 
    B = -∂Hr*ψ_at_R₀;
    A = Hr - E_décal[5]*Diagonal(I, N);
    println("gradients conjgués dérivée seconde"); flush(stdout)
    dψdR = cg(A, B);
    println("norme de dψdR: $(norm(dψdR))"); flush(stdout)
    println("gradients conjugués dérivée seconde terminés"); flush(stdout)
    
    ∂2Hr = Diagonal(zeros(N));
    ∂2Hr[diagind(∂2Hr)] = ∂2R_V_nucl_el.(rs, R_rp) .+ ∂2V_nucl_nucl(R_rp); 
    d2EdR2_ehr = dot(ψ_at_R₀, ∂2Hr, ψ_at_R₀) + 2*dot(dψdR, ∂Hr, ψ_at_R₀);


    return dE_rp, R_rp, iter_nb, d2EdR2_fd, d2EdR2_ehr
end




function plot_énergies_dérivées(∂lE₀, N, Λr, rs, Rs, kdim1d, skip, δR)
    println("calcul différences dérivées Ehrenfest et différences finies et graphes"); flush(stdout)
    skip = div(ker_size,2);
    # code pour les dérivées avec Ehrenfest: plot les erreurs
    d∂E_fin_diff = zeros(N-2*skip);
    for ind_t in (skip+1):(N-skip)
        R_t = Rs[ind_t]
        Hr = copy(Λr);
        Hr[diagind(Hr)] = diag(Λr) + Vector(V_nucl_el.(rs, R_t) .+ V_nucl_nucl(R_t)); 
        vals, vecs, infos = KrylovKit.eigsolve(Hr, N, 1, :SR, krylovdim=kdim1d, issymmetric=true, ishermitian=true);
        @assert infos.converged ≥ 1;
        ψ_t = Vector(vecs[1]);

        # plot(rs, ψ_t, title="ψ_t")

        ∂Hr = Diagonal(zeros(N));
        ∂Hr[diagind(∂Hr)] = ∂R_V_nucl_el.(rs, R_t) .+ ∂V_nucl_nucl(R_t)
        # plot(rs, diag(∂Hr), title="diag(∂Hr)")

        dE_t_ehr = dot(ψ_t,∂Hr,ψ_t) # Ehrenfest
        dE_t_fdi = ∂lE₀[ind_t+4]  # finite differences
        d∂E_fin_diff[ind_t-skip] = abs(dE_t_ehr-dE_t_fdi)
    end

    # PLOT DE LA DIFFÉRENCE ENTRE DÉRIVÉES EHRENFEST ET DÉRIVÉES DIFFÉRENCES FINIES
    plot(Rs[skip+1:N-skip], d∂E_fin_diff, yaxis=:log, xlabel="R", ylabel="erreur absolue", minorticks=1,
    title="Approximation de ∂E₀/∂R: Différence entre\ndifférences finies centrées d'ordre 8 et théorème d'Ehrenfest)\nδR≈$(round(δR,digits=4))",
    label="⟨ψ(R)|dH/dR(R)|ψ(R)⟩- ∫E(R)⋅K₈", legend=true, titlefontsize=11,
    ylims=(1e-14,1e-8))
    savefig("./erreur_dérivée.pdf")


    # PLOT DE L'ÉNERGIE ET SES DÉRIVÉES SUR UN MÊME GRAPHE
    # plot de la dérivée de l'énergie avec la méthode fréquentielle
    # y = [sin(2*x)+cos(5*x) for x in Rs]
    # dydx = [2*cos(2*x)-5*sin(5*x) for x in Rs]
    k = fftfreq(N, 1/δR); # cf doc fs est l'inverse de l'espacement des échantillons
    lE₀_fft = real(ifft(2*π*1im*k.*fft(lE₀.*DSP.Windows.tukey(N, α; padding=0, zerophase=false)))); # .*DSP.Windows.hanning(N)
    plot(Rs,lE₀_fft, label="∂E₀(R) par FFT fenêtrée Tukey($α)", ylims=(-2.4,1.0), linestyle=:solid, color="red")
    # plot!(Rs, dydx)

    # plot de la dérivée de l'énergie avec les différences finies
    plot!(Rs[skip+1:N-skip], ∂lE₀[ker_size:N], xlabel="R", label="∂E₀(R) par différences finies centrées o8",
    legend=true, title="Énergie fondamentale du hamiltonien\nélectronique selon R et sa dérivée",
    titlefontsize=11, color="blue")
    # plot de l'énergie
    plot!(Rs, lE₀, xlabel="R", label="E₀(R)", legend=true, color="green")
    savefig(raw"./énergie_fdmtle_hamiltonien.pdf")

    return d∂E_fin_diff
end







# cette fonction reprend le code du début de la fonction complète et fait des plots
function tests_nv_méthodes(r_min, r_max, R_min, R_max, N, m, kdim1d, kdim2d, Qmax, α)
    u_min, u_max, δu, δu², us, ug = get_rescaling(N);
    δr, δR, δr², δR², N², rs, Rs, rg, Rg, V = get_params_grille(r_max, r_min, R_max, R_min, N);

    # vérifier si les maths en virgule flottante sont correctes: https://0.30000000000000004.com/
    @assert length(rs) == length(Rs) == length(us) == N;

    # CONSTRUCTION DU POTENTIEL ORIGINAL ET DU HAMILTONIEN SUR GRILLE paramétré en R
    V[:,:] = @. V_nucl_el(rg, Rg) + V_nucl_nucl(Rg); # potentiel d'interaction sur la grille tous phénomènes compris

    # STRUCTURE DU LAPLACIEN 1D
    LS = SymTridiagonal(-2*ones(Float64,N), ones(Float64,N-1));
    
        # PLOT POTENTIEL ORIGINAL
        heatmap(rs, Rs, V', xlabel="r", ylabel="R", title="Potentiel de H₂⁺ (r,R) ↦ V(r,R) = Vnn(R)+Ven(r+R/2)+Ven(r-R/2)",
        c = cgrad(:deep, rev=true), titlefontsize=11,
        bbox_inches="tight")
        savefig("./potentiel_original.pdf")

    # LAPLACIENS SUR AXES INDEPENDANTS DE M
    Λr = -1/(δr²*2*m)*LS;  # laplacien sur l'axe r

    ### CALCUL DE L'ÉNERGIE DE SURFACE PAR LA MÉTHODE DES SÉCANTES ###
    dE_rp, R_sec, n_iter, d2EdR2_fd, d2EdR2_ehr = get_lowest_surface_energy_secant(Λr, δR, R_min, rs, N, 2/3*R_max+1/3*R_min, 1/3*R_max+2/3*R_min, 1e-12, kdim1d, 30);

    ### CALCUL DE L'ÉNERGIE DE SURFACE PAR PARCOURS EN LARGEUR ###
    lE₀, E₀_at_R₀, ind_R₀, R₀, K = get_lowest_surface_energy_breadth_first(Λr, δR, R_min, rs, Rs, N, kdim1d);
    @show R_sec, R₀

    diff1_ctr_o8 = [1/280, −4/105, 1/5, −4/5, 0., 4/5, −1/5, 4/105, −1/280];
    ∂lE₀ = conv(lE₀, reverse(1/δR*diff1_ctr_o8));
    ker_size = length(diff1_ctr_o8);

    # PLOT DIFFÉRENCE ENTRE DÉRIVÉES D'EHRENFEST ET DIFF FINIES EN LOG
    d∂E_fin_diff = plot_énergies_dérivées(∂lE₀, N, Λr, rs, Rs, kdim1d, ker_size, δR)

    return d∂E_fin_diff, d2EdR2_fd, d2EdR2_ehr, skip, lE₀, ∂lE₀, R₀, ind_R₀, R_sec, u_min, u_max, δu, δu², us, ug, δr, δR, δr², δR², N², rs, Rs, rg, Rg, V, LS, Λr
end


# APPEL AUX TESTS
d∂E_fin_diff, d2EdR2_fd, d2EdR2_ehr, skip,
    lE₀, ∂lE₀, R₀, ind_R₀, R_sec, u_min, u_max, δu, δu², us, ug,
    δr, δR, δr², δR², N², rs, Rs, rg, Rg, V, LS, Λr = tests_nv_méthodes(r_min, r_max, R_min, R_max, N, m, kdim1d, kdim2d, Qmax, .04);



# QUELQUES VÉRIFICATIONS POUR MONTRER QUE LES SÉCANTES ONT BIEN MIEUX FONCTIONNÉ QUE L'ÉCHANTILLONNAGE
Hr = copy(Λr);
Hr[diagind(Hr)] = diag(Λr) + Vector(V_nucl_el.(rs, R_sec) .+ V_nucl_nucl(R_sec)); 
vals, vecs, infos = KrylovKit.eigsolve(Hr, N, 1, :SR, krylovdim=kdim1d, issymmetric=true, ishermitian=true);
@assert infos.converged ≥ 1;
@assert vals[1] < lE₀[ind_R₀];

@show lE₀[ind_R₀]
@assert lE₀[ind_R₀-1] > lE₀[ind_R₀];
@assert lE₀[ind_R₀+1] > lE₀[ind_R₀];

@show dot([1/280, −4/105, 1/5, −4/5, 0., 4/5, −1/5, 4/105, −1/280], lE₀[ind_R₀-4:ind_R₀+4])

@show dot([1/280, −4/105, 1/5, −4/5, 0., 4/5, −1/5, 4/105, −1/280], lE₀[ind_R₀-4-1:ind_R₀+4-1])







function decompose_hamiltonian_rescaled(r_min, r_max, R_min, R_max, N, m, lM, kdim1d, kdim2d, Qmax)
    # INITIALISATION DES LISTES DE VARIABLES-RÉSULTATS
    l = length(lM);
    l_Ψ_pert = zeros(N*N,l);
    l_Ψ_true = zeros(N*N,l);
    l_Ψ_HBO  = zeros(N*N,l);
    l_E_true = zeros(l);
    # l_E_diff = zeros(l); # pour vérifier l'approximation E₁-E₀ ≈ ω₀
    l_E_pert = zeros(l);
    l_Ψ_L2  = zeros(l);
    l_E_err  = zeros(l);
    l_Ψ_H1 = zeros(l);
    λ_approx = zeros(l);
    Kϵ²      = zeros(l);
    résidus_approx  = zeros(l);
    résidus_pert    = zeros(l);
    
############# ICI FIGURENT LES PARAMETRES INCHANGES AVEC LA MASSE ############# à mettre en fonctions à l'avenir
    u_min, u_max, δu, δu², us, ug = get_rescaling(N);
    δr, δR, δr², δR², N², rs, Rs, rg, Rg, V = get_params_grille(r_max, r_min, R_max, R_min, N);

    # vérifier si les maths en virgule flottante sont correctes: https://0.30000000000000004.com/
    @assert length(rs) == length(Rs) == length(us) == N;

    # CONSTRUCTION DU POTENTIEL ORIGINAL ET DU HAMILTONIEN SUR GRILLE paramétré en R
    V[:,:] = @. V_nucl_el(rg, Rg) + V_nucl_nucl(Rg); # potentiel d'interaction sur la grille tous phénomènes compris (non tronqué)



    # STRUCTURE DU LAPLACIEN 1D
    LS = SymTridiagonal(-2*ones(Float64,N), ones(Float64,N-1)); 

    # LAPLACIENS SUR AXES INDEPENDANTS DE M
    Λr = -1/(δr²*2*m)*LS;  # laplacien sur l'axe r

    
    ### CALCUL DE L'ÉNERGIE DE SURFACE ###
    lE₀, E₀_at_R₀, ind_R₀, R₀, K = get_lowest_surface_energy_breadth_first(Λr, δR, R_min, rs, Rs, N, kdim1d);
    # la constante K ne dépend pas de la masse M

    dE_rp, R_test, iter_nb = get_lowest_surface_energy_secant(Λr, δR, R_min, rs, N, 2/3*R_max+1/3*R_min,1/3*R_max+2/3*R_min, 1e-12, kdim1d, 100);

    
    Λ2D_elec = laplacian_2D_rescaled_dim_elec(N,N²);
    Λ2D_nucl = laplacian_2D_rescaled_dim_nucl(N,N²);

    
############# ICI COMMENCE LA BOUCLE POUR LA MASSE (ce qui précède ne change pas si M change) #############
    for (ind_M,M) in enumerate(lM)

        ### CALCUL DES OPÉRATEURS ###
        println("\n## calcul des opérateurs rescalés masse "*string(M)); flush(stdout)
        ϵ = 1/sqrt(sqrt(K*M)); # paramètre de redimensionnement
        ϵ² = ϵ^2;

        # LAPLACIENS SUR AXES
        Λu = -K*ϵ²/δu²*LS;            # laplacien sur l'axe (Ou)
    


        # FONCTIONS POTENTIELS HBO NON PERTURBÉS SUR AXES SÉPARÉS  
        V₀rs  = V[:,ind_R₀];          # cf formule 3.19 deuxième ligne du rapport sans la constante .- E₀_at_R₀
        V₀us  = .5*K*(ϵ*us).^2        # cf formule 3.23 première ligne du rapport sans la constante .+ E₀_at_R₀
    
        # OPÉRATEURS POTENTIELS HBO NON PERTURBÉS SUR AXES SÉPARÉS 
        V̂⁰rs  = Diagonal(V₀rs);
        V̂⁰us  = Diagonal(V₀us);

    
        # OPÉRATEURS RESCALED SUR GRILLE
        V_res = @. V_nucl_el(rg, R₀.+ϵ*ug) + V_nucl_nucl(R₀.+ϵ*ug); # (r,u) ↦ V(r,u) 

    
        # FONCTION POTENTIEL HBO NON PERTURBÉ SUR GRILLE
        # formule 3.19: somme des deux premières lignes en potentiel:
        # que deux car seulement deux paramétrages (R et u)
        V₀Rg  = @. V[:,ind_R₀]*ones(N)'  + .5*K*(Rg.-R₀).^2; # (r,R) ↦ V(r,R₀) + 1/2*(∂²E₀/∂R²)(R₀)(R-R₀)²
        V₀ug  = @. V_nucl_el(rg, R₀) + V_nucl_nucl(R₀)  + .5*K*ϵ²*ug.^2;   # (r,u) ↦ V(r,u₀) + Kϵ²/2*(∂²E₀/∂u²)(u₀)(u-u₀)²
    
        # OPÉRATEUR POTENTIEL NON PERTURBÉ SUR GRILLE
        V̂⁰ug = Diagonal(reshape(V₀ug, N²));
    
        # CONSTRUCTION DU POTENTIEL ET DU HAMILTONIEN NON PERTURBÉS HBO SUR GRILLE
        # création du laplacien 2D sur grille qui factorise les deux cas 𝔥₀ et 𝔥
        Λ𝔥 = K*ϵ²/2/δu²*Λ2D_nucl + 1/(2*m*δr²)*Λ2D_elec;
        
        # OPÉRATEUR RESCALÉ NON PERTURBÉ SUR GRILLE
        𝔥₀ = Λ𝔥 + V̂⁰ug                         # 𝔥 : Ψ(r,u) ∈ L²(ℝ^N^2) ↦ -1/2m × ∂²/∂r² -1/2M × ∂²/∂u² + V(r,u₀) + Kϵ²/2*(∂²E₀/∂u²)(u₀)(u-u₀)² le hamiltonien HBO non perturbé paramétré en u
    
        # OPÉRATEUR RESCALÉ ORIGINAL SUR GRILLE
        𝔥  = Λ𝔥 + Diagonal(reshape(V_res, N²)) # 𝔥 : Ψ(r,u) ∈ L²(ℝ^N^2) ↦ -1/2m × ∂²/∂r² -1/2M × ∂²/∂u² + V(r,u) le hamiltonien original complet

        println("perturbations"); flush(stdout)
        # FONCTION PERTURBATION Vp (développement de Taylor de la perturbation)
        # ordre 1 en espace: (∂V/∂R)(r,R₀)×(R-R₀)
        ∂R_of_V_at_rR₀   = 1/δR*  V[:,ind_R₀-4:ind_R₀+4]*[1/280 −4/105 1/5 −4/5 0. 4/5 −1/5 4/105 −1/280]';       # vecteur, donne r ↦ ∂V/∂R(r,R₀)                                                   
    
        # ordre 2 en espace: 1/2×[(∂²V/∂R²)(r,R₀) - d²E₀/dR²(R₀)](R-R₀)² 
        ∂²RR_of_V_at_rR₀ = 1/δR²* V[:,ind_R₀-4:ind_R₀+4]*[−1/560 8/315 −1/5 8/5 −205/72 8/5 −1/5 8/315 −1/560]';  # vecteur, donne r ↦ ∂²V/∂R²(r,R₀)                                     

    
        # OPÉRATEURS HAMILTONIENS NON PERTURBÉS RESCALÉS SUR AXES SÉPARÉS
        𝔥u = K*ϵ²*(-1/2/δu²*LS + 1/2*Diagonal(us.^2));  # 𝔥u : ψ(u) ∈ L²(ℝ^N) ↦ 1/2 × Kϵ²(-∂²/∂u² + u²)ψ pour la solution-produit tensoriel
        𝔥r = Λr + V̂⁰rs;                                 # 𝔥r : ψ(r) ∈ L²(ℝ^N) ↦ 1/2 ×    (-∂²/∂r² + R²)ψ pour la solution-produit tensoriel
    
        𝔴  =  @.  V_nucl_el(rg, R₀.+ϵ*ug) + V_nucl_nucl(R₀)  - V_nucl_el(rg, R₀) - V_nucl_nucl(R₀+ϵ*ug) +  1/2*ϵ²*K*ug.^2 # (r,u) ↦ V(r,R₀) - V(r,R₀+ϵu) + Kϵ²/2*u²  (éq. 3.23 du rapport)
        𝔴₁ = ϵ*∂R_of_V_at_rR₀.*ug;                              # matrice, donne (r, u) ↦ u × ∂V/∂u(r,u₀)
        𝔴₂ = .5*ϵ^2*(∂²RR_of_V_at_rR₀*ones(N)' .- K) .*ug.^2;   # matrice, donne (r, u) ↦ 1/2 × u² × ∂V/∂u(r,u₀)
        Vp_res = 𝔴₁ + 𝔴₂; # perturbation totale ordre 1 + ordre 2 en espace (troncature à automatiser éventuellement à tout ordre)
            
        Ŵu = Diagonal(reshape(Vp_res,N²)); # opérateur correspondant à la perturbation paramétrée en u
        

        ### CALCUL DE LA SOLUTION-PRODUIT HARMONIC-BORN-OPPENHEIMER ###
        println("calcul solution_produit par Krylov"); flush(stdout)

        lE⁰x, lϕ⁰x, infos_x = KrylovKit.eigsolve(𝔥r, N, 1, :SR, krylovdim=kdim1d); 
        @assert infos_x.converged ≥ 1;
        
        lE⁰u, lϕ⁰u, infos_u = KrylovKit.eigsolve(𝔥u, N, 1, :SR, krylovdim=kdim1d);
        @assert infos_u.converged ≥ 1;


        ΨHBO = lϕ⁰x[1] * lϕ⁰u[1]'; # normé car les deux le sont déjà en sortant de Krylov
        ΨHBO = reshape(ΨHBO, N²)
        EHBO = lE⁰x[1] + lE⁰u[1];
        densité_HBO = N^2/(u_max-u_min)/(r_max-r_min);


        l_Ψ_pert[:,ind_M] = ΨHBO;
        l_Ψ_HBO[:,ind_M]  = ΨHBO;
        l_E_pert[ind_M]   = EHBO;


        println("calcul solution_référence par Krylov"); flush(stdout)
        ### CALCUL DE LA SOLUTION 2D POUR RÉFÉRENCE DU HAMILTONIEN D'INTÉRÊT PARAMÉTRÉ EN u ###
        lE, lϕ, info_2d = KrylovKit.eigsolve(𝔥, N², 1, :SR, krylovdim=kdim2d); # KrylovKit.eigsolve plus rapide que Arpack.eigs globalement
        @assert info_2d.converged ≥ 1;              # mettre 2 pour trouver aussi le second mode propre
        # l_E_diff[ind_M] = lE[2] - lE[1];

        println("## théorie des perturbations"); flush(stdout)
        l_Ψ_true[:,ind_M] = lϕ[1];
        l_E_true[ind_M] = lE[1];

        ### CALCUL DES PERTURBATIONS ###
        Ψ₀ = copy(ΨHBO);
        W =  copy(Ŵu); # W: sparse
        H₀ = copy(𝔥₀); # sparse
        E₀ = EHBO;


        proj = x -> dot(Ψ₀,x)*Ψ₀; # on gagne ~1 ordre de grandeur en temps en utilisant dot au lieu du produit matriciel Ψ₀*(Ψ₀'*x)
        Π_ort  = LinearMap(x -> x - proj(x), N²); # ne pas assembler
        Π_par  = LinearMap(x -> proj(x), N²);

        P_ort  = LinearMap(x -> Π_ort(H₀*Π_ort(x)-E₀*Π_ort(x)), N²); # Π⟂(H₀-E₀)Π⟂
        P_par  = LinearMap(x -> Π_par(H₀*Π_par(x)-E₀*Π_par(x)), N²); # Π∥(H₀-E₀)Π∥

        llΨ    = zeros(Float64, N², Qmax); # création liste des termes d'énergie   en fonction de l'ordre q
        llE    = zeros(Float64, Qmax);      # création liste des termes de vecteurs en fonction de l'ordre q

        Wl     = LinearMap(x -> W*x, N²);
        b      = -Π_ort(W*Ψ₀); # -Π⟂WΨ₀ dans 3.63

        println("gradients conjugués ordre 1"); flush(stdout)
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
            println("calculs termes pour les CG ordre "*string(q)); flush(stdout)
            # calcul énergie ordre q
            llE[q] = Ψ₀'*WlmE₁(llΨ[:,q-1]); # premier terme de 3.65

            for i ∈ 1:(q-2)
                llE[q] = @views llE[q] - llE[q-i]* Ψ₀'*llΨ[:,i]; # somme du second terme dans 3.65
            end

            # calcul état ordre q
            fill!(acc_b, 0.);
            acc_ort[:] = @views llE[q]*Ψ₀; # dernier terme de la somme de LHS dans 3.66 à i=0

            for i ∈ 1:(q-2)
                @. acc_ort[:] = acc_ort + llE[q-i]*llΨ[:,i] # autres termes de la somme dans LHS de 3.66
            end
            acc_ort[:] = -Π_ort(WlmE₁(llΨ[:,q-1])) + Π_ort(acc_ort); # LHS de 3.66 complet
            
            println("calcul gradients conjugués direction orthogonale ordre "*string(q)); flush(stdout)
            acc_ort[:] = cg(P_ort, acc_ort);

            println("calcul coefficients direction parallèle ordre "*string(q)); flush(stdout)
            α = 0.;
            for i ∈ 1:q-1
                α -= @views .5*dot(llΨ[:,i], llΨ[:,q-i]) # coefficient dans la direction parallèle, donnée par la normalisation de Ψ DSE
            end
            
            llΨ[:,q] = @views acc_ort + α*ΨHBO;
            l_Ψ_pert[:,ind_M] += @views ϵ^q*llΨ[:,q];
            l_E_pert[ind_M]   += @views ϵ^q*llE[q];
        end
        
        Ω_norm_H1 = sparse(I, N², N²) +  1/δu²*Λ2D_nucl + 1/δr²*Λ2D_elec; # + car Λ2D_nucl et Λ2D_elec représentent déjà le laplacien
        println("calcul résultats masse ", string(M))
        diff_vectors = l_Ψ_pert[:,ind_M] - l_Ψ_true[:,ind_M];
        l_Ψ_L2[ind_M]  = norm(diff_vectors);
        l_Ψ_H1[ind_M] = sqrt(dot(diff_vectors, Ω_norm_H1, diff_vectors));
        l_E_err[ind_M]  = abs(l_E_pert[ind_M] - l_E_true[ind_M]);
        # calcul inégalité de Kato-Temple:
        λ_approx[ind_M] = dot(l_Ψ_pert[:,ind_M], 𝔥, l_Ψ_pert[:,ind_M]); # numérateur du quotient de rayleigh
        résidus_approx[ind_M]  = norm(𝔥*l_Ψ_pert[:,ind_M] - λ_approx[ind_M]*l_Ψ_pert[:,ind_M]);
        résidus_pert[ind_M]    = norm(𝔥*l_Ψ_pert[:,ind_M] - l_E_pert[ind_M]*l_Ψ_pert[:,ind_M]);
        Kϵ²[ind_M]      = K*ϵ²;
    end
    return l_Ψ_H1, λ_approx, résidus_approx, résidus_pert, Kϵ², l_E_pert, l_E_true, l_Ψ_pert, l_Ψ_true, l_Ψ_HBO, l_Ψ_L2, l_E_err, # résultats
           N², rs, Rs, rg, Rg, V, LS, Λr, # paramètres
           u_min, u_max, δu, δu², us, ug # rescaling
end





lM = [100, 150, 250, 500, 700, 1000, 3000, 5000];
@time l_Ψ_H1, λ_approx, résidus_approx, résidus_pert, Kϵ², l_E_pert, l_E_true, l_Ψ_pert, l_Ψ_true, l_Ψ_HBO, l_Ψ_L2, l_E_err,
            N², rs, Rs, rg, Rg, V, LS, Λr,
            u_min, u_max, δu, δu², us, ug = decompose_hamiltonian_rescaled(r_min, r_max, R_min, R_max, N, m, lM, kdim1d, kdim2d, Qmax);




kato_temple_est = résidus_approx.^2 ./ Kϵ²;
plot(lM, [l_E_err, kato_temple_est, résidus_pert, l_Ψ_L2.^2, l_Ψ_H1.^2],
            xaxis=:log, yaxis=:log, seriestype = :scatter,
            title="Variables de contrôle",
            label=["erreur énergie à la référence: |Eₐ-E|" "Kato-Temple quotient Rayleigh: ||hΨₐ-⟨Ψₐ,h,Ψₐ⟩Ψₐ||²/ω₀ (norme 2)" "résidu ||hΨₐ-EₐΨₐ|| (norme 2)" "erreur état à la référence ||Ψₐ-Ψ||² (norme 2)" "erreur état à la référence ||Ψₐ-Ψ||² (norme H1)"],
            xlabel="masse M", size=(600,400), ylims=(1e-5,1e-1), legend=:bottomleft) 
savefig("./resultats.pdf")

# pour afficher les pentes
abscisses = log10.(lM[1:6]);
ord_E     = log10.(l_E_err[1:6]);
ord_H1_c  = l_Ψ_H1.^2;
ord_H1    = log10.(ord_H1_c[1:6]);
ord_L2_c  = l_Ψ_L2.^2;
ord_L2    = log10.(ord_L2_c[1:6]);

@show (ord_H1[6]-ord_H1[1])/(abscisses[6]-abscisses[1]);
@show (ord_E[6]-ord_E[1])/(abscisses[6]-abscisses[1]);
@show (ord_L2[6]-ord_L2[1])/(abscisses[6]-abscisses[1]);


heatmap(rs, us, reshape(l_Ψ_pert[:,6],N,N)'.^2, xlabel="coordonnée électronique r", ylabel="coordonnée nucléaire R")
heatmap(rs, us, reshape(l_Ψ_true[:,6],N,N)'.^2, xlabel="coordonnée électronique r", ylabel="coordonnée nucléaire R")
heatmap(rs, us, reshape(l_Ψ_HBO[:,6],N,N)'.^2, xlabel="coordonnée électronique r", ylabel="coordonnée nucléaire R")

plot(lM, l_E_err, yaxis=:log, seriestype = :scatter, label="erreur énergie", xlabel="masse M", ylabel="|E - Eₚ|",size=(400,200))
plot(lM, l_Ψ_L2, yaxis=:log, seriestype = :scatter, label="résidu", xlabel="masse M", ylabel="|Ψ - Ψₚ|",size=(400,200))


using CUDA
using CUDA.CUSPARSE
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using KrylovKit
N = 100;
r_cpu = sprand(N*N,N*N,1/N/N);
r_gpu = CuSparseMatrixCSC(r_cpu);
x_cpu = rand(N*N);
x_gpu = cu(x_cpu);

@time KrylovKit.eigsolve(r_cpu, N*N, 1, :SR, krylovdim=20);
CUDA.@time KrylovKit.eigsolve(r_gpu, N*N, 1, :SR, krylovdim=20);

@time      d_cpu = cg(r_cpu, x_cpu);
CUDA.@time d_gpu = cg(r_gpu, x_gpu);