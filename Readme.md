# Scope and origin of Lcanvas

Lcanvas is a code meant to study the performance scientific codes that operate with sparse data with irregular patterns and the various dimensions that affect them
This code is using: https://www.ict.inaf.it/gitlab/luca.tornatore/nbody-vectorization-toy-lab/ as a starting point.

The main dimensions it aims to explore are:

- The effect of memory Layout of the data of such codes
- The strategies for interleaving memory accesses and computation
- The efficiency of various acceleration strategies (vectorization, GPU-offloading etc)

# A small toy-sized lab for vectorization of Many-body codes pattern

This repository contains the shared work done by (some) participants in the [SPACE CoE](https://www.space-coe.eu/) in experimenting the vectorization achieved when different data layouts (AoS, SoA, SoAos) are adopted.

What we address here is the many-body problem, i.e. the pattern from codes that solve the gravitational interaction among particles that are _not_ “true bodies” but instead represent matter samples in the phase-space. A such, they do not need to deploy highest order $\propto N^2$ integration methods to correclty solve close encounters. At odds, to speed-up the calculation approximate methods can be used while obtaining an adequately small error.

Typically, approximate methods are coupled with a tree structure and a Barnes&Hut algorithm. Hence, the contribution to the force of the “closest” neighbouring particles is accounted in direct summation, as in N-Body codes, while the contribution of “distant” particles is accounted using different approximation (multipole expansion, PM, …).

Here we focus on the loop over the close particles, that resembles the N-body typical loop but for the ffollowing facts:

- the array of neighbours is in general unique per each particles, by definition; then, every particle loops over a different array of partners

- the fact that $n_i$ particles are neighbour of a target particle does not mean that they are close in memory; on the contrary, the very likely pattern is that looping over the neighbours results in memory jumps

- a loop of the kind

  ```c
  for ( int = 0; i < Nneighbours; i++ ) {
      int k = neighbours_list[k];
      force += contributtion( Particles[k] ); }    
  ```

  is difficult to be vectorized, and definitely not automatically vectorizable by the compilers.

The main focus of the toy codes here is to reprodice the pattern aforementioned:

- $N_p$ particles are randomly generated
- A subset of them, of size $N_A$, will be considered “active” - it means that in a real code, in which the timesteps are not all equal, they would need to calculate the force, update velocities, etc.
- Every active particles $i$ is assigned with $N_n^i$ neighbours (where $N_n^i \simeq <N_n>\pm 10\%$) and has to calculate the $\sum_j1/r_{ij}^2$ contribution.

We are experimenting 2 data layouts SoA and SoAos, using $(i)$ a **plain** implementation, $(ii)$ using **builtin vector types** nowadays supported by all the compilers and $(iii)$ **using vector intrinsics**.
Since the memory toll due to the non-cache friendly pattern, we compare 3 different solutions:

**v1**. trying to hide the memory latency by prefetching of “next” neigbours

**v2**. copying all the neighbours into a memory buffer so that the loop is then linear on consequent memory addresses (of course could hardly been done in a real code)

**v3**. evolving the v2 in a more realistic way, having a thread that fetches neighbours particles in a limited buffer while other threads process them

