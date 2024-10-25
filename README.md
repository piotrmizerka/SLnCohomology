# Replication details for [TO FILL](TO FILL)

This code can be used to replicate the results of [arXiv preprint](TO FILL). It provides the necessary computations to prove that the cohomology groups $H^2(\text{SL}_3(\mathbb{Z}),\pi_3)$ and $H^3(\text{SL}_4(\mathbb{Z}),\pi_4)$ are non-zero for some orthogonal representations $\pi_3$ and $\pi_4$ all of whose invariant vectors are trivial.

One can express the rank of the cohomology groups as the rank of specific Laplacians built up from the homological data obtained from Voronoi tesselations of the associated symmetric spaces. We compute the ranks of these Laplacians. We already precomputed the Laplacians since the case $n=4$ required more time and RAM to be completed in a reasonable time (probably about 32GB shall suffice). They are saved in the [laplacians_computations](./scripts/laplacians_computation) directory. A separate computation of these Laplacians is also possible (provided sufficient computational resources). 

To ensure mathematical rigour, the matrices representing the Laplacians evaluated for particular representations are rationally-valued. For the same reason, we use the package [LinearAlgebraX](https://github.com/scheinerman/LinearAlgebraX.jl) which provides a function to compute the rank of rationally-valued matrices exactly.

For the computations we used Julia in version `1.9.4`. On a Windows operating system, please use the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/about).

## Getting Julia 1.9.4

In order to recompute the Laplacians, you need to use Julia in version `1.9.4`. (If you want to use the precomputed Laplacians that we provide in this repository, any later Julia version should work as well.)

To get the correct Julia version, first install [Juliaup](https://github.com/JuliaLang/juliaup), a cross-platform installer for Julia. 

After installing Juliaup, install Julia in version `1.9.4` by running

```bash
juliaup add 1.9.4
```

Run the following command to make version `1.9.4` the default version.

```bash
juliaup default 1.9.4
```

## Obtaining the code
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
using Pkg; Pkg.instantiate()
```
Remark: this step needs to be executed only once per installation.

## Running actual replication
In order to replicate the computations for $\text{SL}_n(\mathbb{Z})$ using the precomputed Laplacians, run the following command in the terminal in the `SLnCohomology` folder, replacing the parameter `n` by `3` or `4`.
```bash
julia --project=. ./scripts/sln_nontrivial_cohomology.jl n
```

The running time of the script will be approximately `3` and `8` minutes on a standard laptop computer for the cases $n=3$ and $n=4$ respectively.

To run the script which recomputes the Laplacians, saving them in [laplacians_computations](./scripts/laplacians_computation) directory, run the following command in the terminal in the `SLnCohomology` folder
```bash
julia --project=. ./scripts/laplacians_computation/sln_laplacians.jl n
```
Caution: for $n=4$, the running time of the above script may take a few hours and the sufficient amount of RAM memory is recommended (probably about 32GB).
