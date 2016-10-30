# Expokit.jl
Julia interface to [Expokit](http://www.maths.uq.edu.au/expokit/),
"...a software package that provides matrix exponential routines for small dense or very large sparse matrices, real or complex."

## Installation
```julia
julia> Pkg.clone("https://github.com/HaraldHofstaetter/Expokit.jl")
julia> Pkg.build("Expokit")
```

## Examples
To get easy access to the examples, copy them into the home directory:
```julia
julia> cp(joinpath(homedir(), ".julia/v0.4/Expokit/examples/"), joinpath(homedir(), "Expokit_examples"), remove_destination=true)
```
Then 'Expokit_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [Expokit_examples.ipynb](https://github.com/HaraldHofstaetter/Expokit.jl/blob/master/examples/Expokit_examples.ipynb)

## Technical informations
+ This package contains the file [expokit.f](https://github.com/HaraldHofstaetter/Expokit.jl/blob/master/deps/src/expokit.f) 
  which was (slightly modified) taken from http://www.maths.uq.edu.au/expokit/download.html, see its 
  [copyright notice](https://github.com/HaraldHofstaetter/Expokit.jl/blob/master/deps/src/copyright). 
  My only modifications to this file concern the external subroutine `matvec` for matrix-vector multiplication, which now has 
  a third argument `arg` of type `integer*8`. Correspondingly, the driver routines now also have an additional argument `arg` 
  of this type. This allows passing closures via pass-through pointers, see  http://julialang.org/blog/2013/05/callback for a     description of this technique.
+ In the original Expokit Fotran 77 code there are many `stop` statements to terminate
  the program if an error occurs. This is not advantageous if the code is compiled to a shared library to be called
  from Julia in an interactive session. This problem was solved by "...subverting the exit
  mechanism using the obscure long jump mechanism...", see http://stackoverflow.com/a/19608029 for this clever idea and 
  [fortran_stop_wrapper.c](https://github.com/HaraldHofstaetter/Expokit.jl/blob/master/deps/src/fortran_stop_wrapper.c)
  for my implementation.
