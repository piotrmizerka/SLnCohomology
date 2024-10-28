# Replication details for [TO FILL](TO FILL)

This code can be used to replicate the results of [arXiv preprint](TO FILL). It provides the necessary computations to prove that the cohomology groups $H^2(\text{SL}_3(\mathbb{Z}),\pi_3)$ and $H^3(\text{SL}_4(\mathbb{Z}),\pi_4)$ are non-zero for some orthogonal representations $\pi_3$ and $\pi_4$ all of whose invariant vectors are trivial.

One can express the rank of the cohomology groups as the rank of specific Laplacians built up from the homological data obtained from Voronoi tesselations of the associated symmetric spaces. We compute the ranks of these Laplacians.

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

## Running the actual replication
First compute the Laplacian in the relevant degrees. For this, run the following command in the terminal in the `SLnCohomology` folder, replacing $(n,degree)$ by $(3,3)$ or $(4,6)$.
```bash
julia --project=. ./scripts/laplacians_computation/sln_laplacians.jl n degree
```
For $(4,6)$, this needs around 8 GB of available RAM. You can also replace `degree` by a list of numbers, e.g. writing `2 3 4 5` computes the Laplacians in degrees 2, 3, 4 and 5. But be aware that for `n` equal to 4, degree 4 takes a few hours to compute and some more memory (16 GB of system RAM are enough though).
The Laplacians in the demanded degrees are saved in the [laplacians_computations](./scripts/laplacians_computation) directory. Precomputed versions of the Laplacians in all degrees are also provided on  [Zenodo](TO FILL).

To compute the ranks of these Laplacians for the representations given in the paper, run the following command in the terminal in the `SLnCohomology` folder, replacing the parameter `n` by `3` or `4`.
```bash
julia --project=. ./scripts/sln_nontrivial_cohomology.jl n
```

The running time of the script will be approximately `3` and `8` minutes on a standard laptop computer for the cases $n=3$ and $n=4$ respectively.

There is also a possibility to run all the necessary operations in the [Jupyter notebook](./sln_non_trivial_cohomology.ipynb). The notebook contains procedures for $n=3,4$ and the homology degrees equal to $3$ and $6$ respectively (they correspond to the cohomologies of degrees $2$ and $3$ for $n=3$ and $n=4$).

