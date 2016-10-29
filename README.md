# Expokit.jl
Julia interface to [Expokit](http://www.maths.uq.edu.au/expokit/)

##Installation
```julia
Pkg.clone("https://github.com/HaraldHofstaetter/Expokit.jl")
Pkg.build("Expokit")
```

##Examples
To get easy access to the examples, copy them into the home directory:
```julia
cp(joinpath(homedir(), ".julia/v0.4/Expokit/examples/"), joinpath(homedir(), "Expokit_examples"), remove_destination=true)
```
Then 'Expokit_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [Expokit_examples.ipynb](https://github.com/HaraldHofstaetter/Expokit.jl/blob/master/examples/Expokit_examples.ipynb)
