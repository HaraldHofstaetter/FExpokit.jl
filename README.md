# FExpokit.jl
Julia interface to [Expokit](http://www.maths.uq.edu.au/expokit/),
"...a software package that provides matrix exponential routines for small dense or very large sparse matrices, real or complex."
The package is called FExpokit.jl  to distinguish it from [Expokit.jl](https://github.com/acroy/Expokit.jl), which is a reimplementation of the Expokit algorithms in Julia;
the letter F indicates that  "Fortran-Expokit" is called.

## Installation
```julia
julia> Pkg.clone("https://github.com/HaraldHofstaetter/FExpokit.jl")
julia> Pkg.build("FExpokit")
```
## Usage

```julia
using FExpokit
```

```julia
w = expv(t, A, v)
```
computes `w = exp(t*A)*v` where `A` is of a type for which  [`Base.LinAlg.A_mul_B!`](https://docs.julialang.org/en/stable/stdlib/linalg/#Base.LinAlg.A_mul_B!) is defined,
`t` is `Real`,
and `v` is of type `Vector{Float64}`or `Vector{Complex128}`.

```julia
w = expv(t, A, v, matrix_times_minus_i=true)
```
computes `w = exp(-1im*t*A)*v`. This is convenient for applications in quantum mechanics.

```julia
w = expv(t, F, v, anorm)
```
Here, `F` is a Julia function which computes `F(x) = A*x`, and `anorm` is an approximation of some norm of `A`.

### Optional parameters

```julia
w = expv(t, A, v,  m=m, tol=tol, symmetric=true, trace=true, anorm=anorm)
```
parameter | description
-------------------|-------------
 `m`               | maximum size of the Krylov basis (default: `m=30`)
 `tol`             | the requested accuracy tolerance on `w`. If `tol=0` or `tol` is too small (`tol<eps`) the default value   `sqrt(eps)` is used.
 `symmetric=true`  | indicates that `A` can be treated as a symmetric matrix. (For complex hermitian matrices use `hermitian=true`).
`matrix_times_minus_i=true` | indicates that `A` is replaced by `-1im*A`. If combined with `hermitian=true` this means that `A` is hermitian (not `-1im*A`, which is skew-hermitian in this case).
 `trace=true`      | indicates that step-by-step info shall be printed.
 `anorm`           | an approximation of some norm of `A`.

### Statistics

```julia
w, stats = expv(t, A, v, statistics=true)
``` 
Here `stats` is a dictionary with the following entries:

mnemonic            |     description
--------------------|---------------------------------------------------
 `:nmult`           | number of matrix-vector multiplications used      
 `:nexph`           | number of Hessenberg matrix exponential evaluated 
 `:nscale`          | number of repeated squaring involved in Pade      
 `:nstep`           | number of integration steps used up to completion 
 `:nreject`         | number of rejected step-sizes                     
 `:ibrkflag`        | set to 1 if "happy breakdown" and 0 otherwise     
 `:mbrkdwn`         | if "happy brkdown", basis-size when it occured    
 `:step_min`        | minimum step-size used during integration         
 `:step_max`        | maximum step-size used during integration         
 `:x_error`         | maximum among all local truncation errors         
 `:s_error`         | global sum of local truncation errors             
 `:tbrkdwn`         | if "happy breakdown", time when it occured        
 `:t_now`           | integration domain successfully covered           
 `:hump`            | `max||exp(sA)||`, `s in [0,t]` (or `[t,0]` if `t<0`)      
 `:scaled_norm_sol` | `||w||/||v||`, scaled norm of the solution w.      



## Examples
To get easy access to the examples, copy them into the home directory:
```julia
julia> cp(joinpath(homedir(), ".julia/v0.4/Expokit/examples/"), joinpath(homedir(), "Expokit_examples"), remove_destination=true)
```
Then 'Expokit_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [Expokit_examples.ipynb](https://github.com/HaraldHofstaetter/FExpokit.jl/blob/master/examples/Expokit_examples.ipynb)

## Technical informations
+ This package contains the file [expokit.f](https://github.com/HaraldHofstaetter/FExpokit.jl/blob/master/deps/src/expokit.f)
  which was (slightly modified) taken from http://www.maths.uq.edu.au/expokit/download.html, see its
  [copyright notice](https://github.com/HaraldHofstaetter/Expokit.jl/blob/master/deps/src/copyright).
  One of my few modifications to this file concerns the external subroutine `matvec` for matrix-vector multiplication, which now has
  a third argument `arg` of type `integer*8`. Correspondingly, the driver routines now also have an additional argument `arg`
  of this type. This allows passing closures via pass-through pointers, see  http://julialang.org/blog/2013/05/callback for a
  description of this technique.
+ In the original Expokit Fotran 77 code there are many `stop` statements to terminate
  the program if an error occurs. This error handling is not advantageous if the code is compiled to a shared library
  mentioned to be called from Julia in an interactive session. The unwanted  program terminations after `stop` are
  avoided by "...subverting the exit mechanism using the obscure long jump mechanism...", see
  http://stackoverflow.com/a/19608029 for this clever idea and
  [fortran_stop_wrapper.c](https://github.com/HaraldHofstaetter/FExpokit.jl/blob/master/deps/src/fortran_stop_wrapper.c)
  for my implementation.
