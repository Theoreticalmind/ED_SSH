using LinearAlgebra
using Combinatorics
using Printf

function build_hamiltonian_pbc(L::Int,Np::Int, t1::Float64, t2::Float64, V::Float64,wij_disorder::Vector{Tuple{Tuple{Int,Int},Float64}})


    basis = [Tuple(c) for c in combinations(1:L, Np)]
    basis_size = length(basis)
    state_to_index = Dict(state => i for (i, state) in enumerate(basis))

    H_many_body = zeros(Float64, basis_size, basis_size)

    pbc_sign_factor = (-1)^(Np - 1)

    for i in 1:basis_size
        state = basis[i]

        interaction_energy = 0.0
        for site_j in 1:L
            neighbor = mod1(site_j + 1, L)
            if (site_j in state) && (neighbor in state)
                interaction_energy += V
            end
        end
        H_many_body[i, i] = interaction_energy

        for site_k in state
            if site_k < L
                neighbor = site_k + 1
                if !(neighbor in state)
                    hopping_t = isodd(site_k) ? t1 : t2
                    new_state = Tuple(sort!([s for s in state if s != site_k] ∪ [neighbor]))
                    j = state_to_index[new_state]
                    H_many_body[j, i] += hopping_t
                    H_many_body[i, j] += hopping_t
                end
            else 
                neighbor = 1
                if !(neighbor in state)
                    
                    hopping_t = isodd(L) ? t1 : t2
                    
                    hopping_t_pbc = hopping_t * pbc_sign_factor

                    new_state = Tuple(sort!([s for s in state if s != site_k] ∪ [neighbor]))
                    j = state_to_index[new_state]
                    H_many_body[j, i] += hopping_t_pbc
                    H_many_body[i, j] += hopping_t_pbc
                end
            end
        end
        
    

        #= --- Off-diagonal: Disorder hopping (w_ij) ---
        for (sites, w_ij) in wij_disorder
            site_a, site_b = sites

            # Enforce sublattice symmetry: only hopping between odd <-> even sites
            if isodd(site_a) == isodd(site_b)
                continue
            end

            if (site_a in state) && !(site_b in state)
                new_state_tuple = Tuple(sort!([s for s in state if s != site_a] ∪ [site_b]))
                j = state_to_index[new_state_tuple]
                H_many_body[j, i] += w_ij
                H_many_body[i, j] += w_ij
            elseif !(site_a in state) && (site_b in state)
                new_state_tuple = Tuple(sort!([s for s in state if s != site_b] ∪ [site_a]))
                j = state_to_index[new_state_tuple]
                H_many_body[j, i] += w_ij
                H_many_body[i, j] += w_ij
            end
        end=#
    end
    return Hermitian(H_many_body)
end


function main()
    L = 8
    Np = L ÷ 2
    V = 2.0
    t1 = 1.0
    t2 = 0.8
    #σ = 0.5

    #Random.seed!(1234) 
    wij_disorder = Vector{Tuple{Tuple{Int,Int},Float64}}()
    #=println("Generating random disorder terms with σ = $σ...")
    for site_i in 1:L
        for site_j in 1:L
            if site_i < site_j && isodd(site_i + site_j)
                w_ij = rand(Normal(0, σ / sqrt(L)))
                push!(wij_disorder, ((site_i, site_j), w_ij))
            end
        end
    end=#

    basis = [Tuple(c) for c in combinations(1:L, Np)]
    basis_size = length(basis)

    println("--- Corrected Exact Diagonalization with PBC ---")
    println("System Size L=$L, Particle Number Np=$Np (Half-Filling)")
    println("Interaction V=$V, Hopping t1=$t1, t2=$t2")
    println("Many-body Hilbert space dimension: $basis_size x $basis_size")

    println("\nBuilding the Hamiltonian matrix...")
    H = build_hamiltonian_pbc(L, Np, t1, t2, V,wij_disorder)

    println("Performing exact diagonalization to find eigenvalues...")
    eigenvalues = sort(eigvals(H))
    println("Diagonalization complete.")

    println("\n--- Computed Eigenvalues (Energy Spectrum) ---")
    num_eigenvalues_to_show = min(10, basis_size)
    @printf("Showing the lowest %d eigenvalues:\n", num_eigenvalues_to_show)
    for i in 1:num_eigenvalues_to_show
        @printf("E_%-2d = % .6f\n", i-1, eigenvalues[i])
    end
    #=println("\n--- First 5 generated disorder terms (site_i, site_j), w_ij ---")
    for i in 1:min(5, length(wij_disorder))
        println(wij_disorder[i])
    end=#
end

main()
