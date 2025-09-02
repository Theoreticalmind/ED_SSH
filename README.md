# Exact Diagonalization with Periodic Boundary Conditions (PBC)

This project implements exact diagonalization of an interacting 1D fermionic chain with periodic boundary conditions (PBC).
The model includes nearest-neighbor hopping (t₁, t₂) and density-density interaction (V).

Hamiltonian

The Hamiltonian is defined on a 1D chain of length L at half-filling (Np = L/2):

𝐻
=
−
∑
⟨
𝑖
,
𝑗
⟩
𝑡
𝑖
𝑗
 
(
𝑐
𝑖
†
𝑐
𝑗
+
𝑐
𝑗
†
𝑐
𝑖
)
+
𝑉
∑
⟨
𝑖
,
𝑗
⟩
𝑛
𝑖
𝑛
𝑗
H=−
⟨i,j⟩
∑
	​

t
ij
	​

(c
i
†
	​

c
j
	​

+c
j
†
	​

c
i
	​

)+V
⟨i,j⟩
∑
	​

n
i
	​

n
j
	​


where

𝑡
𝑖
𝑗
t
ij
	​

 alternates between t₁ and t₂,

𝑉
V is the nearest-neighbor interaction strength,

Fermionic sign factors are included for periodic wrapping.
