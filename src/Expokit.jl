__precompile__()

module Expokit

export expv

const libexpokit = joinpath(dirname(@__FILE__),
                            "..", "deps", "lib", string("libexpokit.", Libdl.dlext))

function expv(t::Real, matvec::Function, v::Vector{Float64}, anorm::Real; 
              tol::Real=0.0, m::Integer=30, symmetric::Bool=false, trace::Bool=false)
    n = length(v)
    w = zeros(Float64, n)
    lwsp =  max(10, n*(m+1)+n+(m+2)^2+4*(m+2)^2+6+1)
    wsp = zeros(Float64, lwsp)
    liwsp = max(7, m+2)
    iwsp = zeros(Int32, liwsp)
    iflag = zero(Int32)
    matvec_fortran = cfunction(matvec, Void, (Ptr{Float64}, Ptr{Float64}))
    if symmetric
        ccall((:dsexpv_, libexpokit), Void, 
        (Ptr{Int32},   Ptr{Int32},   Ptr{Float64}, Ptr{Float64},   Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32},     Ptr{Int32},   Ptr{Void},      Ptr{Int32},  Ptr{Int32}), 
         &n,           &m,           &t,           v,              w,            &tol, 
         &anorm,       wsp,          &lwsp,        iwsp,           &liwsp,       matvec_fortran, &trace,      &iflag)
    else
        ccall((:dgexpv_, libexpokit), Void, 
        (Ptr{Int32},   Ptr{Int32},   Ptr{Float64}, Ptr{Float64},   Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32},     Ptr{Int32},   Ptr{Void},      Ptr{Int32},  Ptr{Int32}), 
         &n,           &m,           &t,           v,              w,            &tol, 
         &anorm,       wsp,          &lwsp,        iwsp,           &liwsp,       matvec_fortran, &trace,      &iflag)
    end
    if iflag<0
        error("bad input arguments")
    elseif iflag==1
        error("maximum number of steps reached without convergence")
    elseif iflag==2
        error("requested tolerance was too high")
    end
    return w
end       

function _matvec_(v::Ptr{Float64}, w::Ptr{Float64})
    global _A_, _n_
    v1 = pointer_to_array(v, _n_)    
    w1 = pointer_to_array(w, _n_) 
    w1[:] = _A_*v1
    return nothing
end    

function expv(t::Real, A::AbstractArray{Float64,2}, v::Vector{Float64}; 
              tol::Real=0.0, m::Integer=30, symmetric::Bool=isa(A, Hermitian), trace::Bool=false, anorm::Real=-1.0)
    global _A_, _n_
    _n_ = length(v)
    _A_ = A
    if anorm<=0.0
        anorm = norm(A, Inf)
    end
    expv(t, _matvec_, v, anorm, tol=tol, m=m, symmetric=symmetric, trace=trace)
end              

function expv(t::Real, matvec::Function, v::Vector{Complex{Float64}}, anorm::Real; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=false, trace::Bool=false)
    n = length(v)
    w = zeros(Complex{Float64}, n)
    lwsp =  max(10, n*(m+1)+n+(m+2)^2+4*(m+2)^2+6+1)
    wsp = zeros(Complex{Float64}, lwsp)
    liwsp = max(7, m+2)
    iwsp = zeros(Int32, liwsp)
    iflag = zero(Int32)
    matvec_fortran = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}))
    if hermitian
        ccall((:zhexpv_, libexpokit), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32}), 
         &n,           &m,                    &t,           v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec_fortran,
         &trace,       &iflag)
    else
        ccall((:zgexpv_, libexpokit), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32}), 
         &n,           &m,                    &t,           v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec_fortran,
         &trace,       &iflag)
    end
    if iflag<0
        error("bad input arguments")
    elseif iflag==1
        error("maximum number of steps reached without convergence")
    elseif iflag==2
        error("requested tolerance was too high")
    end
    return w
end              


function _matvec_cmplx_(v::Ptr{Complex{Float64}}, w::Ptr{Complex{Float64}})
    global _A_, _n_
    v1 = pointer_to_array(v, _n_)    
    w1 = pointer_to_array(w, _n_) 
    w1[:] = _A_*v1
    return nothing
end    


function expv(t::Real, A::AbstractArray{Complex{Float64},2}, v::Vector{Complex{Float64}}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=isa(A, Hermitian), trace::Bool=false, anorm::Real=-1.0)
    global _A_
    _A_ = A
    if anorm<=0.0
        anorm = norm(A, Inf)
    end
    expv(t, _matvec_cmplx_, v, anorm, tol=tol, m=m, symmetric=symmetric, trace=trace)
end   



end # module Expokit
