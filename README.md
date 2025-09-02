# Exact Diagonalization with Periodic Boundary Conditions (PBC)

# Exact Diagonalization with Periodic Boundary Conditions (PBC)

This project implements exact diagonalization for an interacting 1D fermionic chain with periodic boundary conditions (PBC).
The model includes alternating nearest-neighbor hoppings (`t₁`, `t₂`) and a density-density interaction (`V`).

## Hamiltonian

The Hamiltonian is defined on a 1D chain of length `L` at half-filling (`Nₚ = L/2`):

\[
H = -\sum_{\langle i, j \rangle} t_{ij} (c^\dagger_i c_j + c^\dagger_j c_i) + V \sum_{\langle i, j \rangle} n_i n_j
\]

where:
- \( t_{ij} \) alternates between `t₁` and `t₂` along the chain.
- \( V \) is the nearest-neighbor interaction strength.
- \( c^\dagger_i \), \( c_j \) are fermionic creation and annihilation operators.
- \( n_i = c^\dagger_i c_i \) is the number operator at site \( i \).
- The sums run over all nearest-neighbor pairs \( \langle i, j \rangle \), including the periodic wrap-around.

**Note:**  
Fermionic sign factors are properly included for the periodic boundary (i.e., when hopping from the last to the first site).

---mionic sign factors are included for periodic wrapping.
