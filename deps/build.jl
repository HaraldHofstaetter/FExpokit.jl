cd(dirname(@__FILE__))

if (!ispath("lib"))
    run(`mkdir lib`)
end

download("http://www.maths.uq.edu.au/expokit/expokit.tar.gz", "expokit.tar.gz")
run(`tar xzvf expokit.tar.gz`)

run(`mv expokit/data ../examples`)

cd(joinpath(dirname(@__FILE__), "expokit", "fortran" ))
run(`make -f ../../Makefile1`)

run(`mv libexpokit.$(Libdl.dlext) ../../lib`)

