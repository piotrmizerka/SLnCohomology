# Replication details for [TO FILL](TO FILL)

This code can be used to replicate the results of [arXiv preprint](TO FILL). It provides the necessary computations to prove that the cohomologies $H^2(\text{SL}_3(\mathbb{Z}),\pi_3)$ and $H^3(\text{SL}_4(\mathbb{Z}),\pi_4)$ are non-zero for some orthogonal representations $\pi_3$ and $\pi_4$ not admitting non-trivial invariant vectors.

One can express the rank of the cohomologies as the rank of specific Laplacians built up from the homological data obtained by means of Voronoi cells. We compute the ranks of such Laplacians. We have already precompued these Laplacians since the case $n=4$ required more time and RAM memory to be completed in a reasonable time (probably about 32GB shall suffice). They are saved in [laplacians_computations](./scripts/laplacians_computation) directory. A separate computation of these Laplacians is also possible (provided sufficient computational resources). 

To ensure mathematical rigorousness, all the matrices representing the Laplacians evaluated for particular representations are rationally-valued. For the same reason, we use the package [LinearAlgebraX](https://github.com/scheinerman/LinearAlgebraX.jl) which provides a function to compute the rank of rationally-valued matrices exactly.

For the computations we used Julia in version `1.9.4` but in principle any later version should work.

## Obtaining code
In order to get the replication code, either download it directly from [Zenodo](TO FILL), or issue the following command in the terminal (note that git must be installed in this case)
```bash
git clone https://github.com/piotrmizerka/SLnCohomology
```

## Setting up the environment
Open Julia by running the following command in the terminal in `SLnCohomology` folder
```bash
julia --project=.
```
Next, to set up the proper environment for the replication run in Julia REPL the following
```julia
julia> using Pkg; Pkg.instantiate()
```
Remark: this step needs to be executed only once per installation.

## Running actual replication
In order to replicate the computations for $\text{SL}_n(\mathbb{Z})$ using the precomputed Laplacians, run the following command in the terminal in the `SLnCohomology` folder
```bash
julia --project=. ./scripts/sln_nontrivial_cohomology.jl n
```

The running time of the script will be approximately `3` and `8` minutes on a standard laptop computer for the cases $n=3$ and $n=4$ respectively.

To run the script which recomputes the Laplacians, saving them in [laplacians_computations](./scripts/laplacians_computation) directory, run the following command in the terminal in the `SLnCohomology` folder
```bash
julia --project=. ./scripts/sln_laplacians.jl n
```
Caution: for $n=4$, the running time of the above script may take a few hours and the sufficient amount of RAM memory is recommended (probably about 32GB).
