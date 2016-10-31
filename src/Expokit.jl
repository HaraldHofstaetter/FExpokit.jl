__precompile__()

module Expokit

export expv


function __init__()
    global const libexpokit = Libdl.dlopen(joinpath(dirname(@__FILE__),
                            "..", "deps", "lib", string("libexpokit.", Libdl.dlext)))
    ccall(Libdl.dlsym(libexpokit, :fortran_stop_handler_init), Void, () )
end


function _expv(t::Real, matvec::Ptr{Void}, v::Vector{Float64}, anorm::Real; 
              tol::Real=0.0, m::Integer=30, symmetric::Bool=false, trace::Bool=false, statistics::Bool=false, arg::Ptr{Void}=convert(Ptr{Void},0))
    n = length(v)
    w = zeros(Float64, n)
    lwsp =  max(10, n*(m+1)+n+(m+2)^2+4*(m+2)^2+6+1)
    wsp = zeros(Float64, lwsp)
    liwsp = max(7, m+2)
    iwsp = zeros(Int32, liwsp)
    iflag = Int32(0) 
    if symmetric
        ccall(Libdl.dlsym(libexpokit, :dsexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},   Ptr{Float64}, Ptr{Float64},   Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32},     Ptr{Int32},   Ptr{Void},    Ptr{Int32},  Ptr{Int32}, Ptr{Void}), 
         &n,           &m,           &t,           v,              w,            &tol, 
         &anorm,       wsp,          &lwsp,        iwsp,           &liwsp,       matvec,       &trace,      &iflag,     arg )
    else
        ccall(Libdl.dlsym(libexpokit, :dgexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},   Ptr{Float64}, Ptr{Float64},   Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32},     Ptr{Int32},   Ptr{Void},    Ptr{Int32},  Ptr{Int32}, Ptr{Void}), 
         &n,           &m,           &t,           v,              w,            &tol, 
         &anorm,       wsp,          &lwsp,        iwsp,           &liwsp,       matvec,       &trace,      &iflag,     arg )
    end
    if ccall(Libdl.dlsym(libexpokit, :has_stopped), Int32, ()) == 1
        error("expokit fortran library stop")
    end    
    if iflag<0
        error("bad input arguments")
    elseif iflag==1
        error("maximum number of steps reached without convergence")
    elseif iflag==2
        error("requested tolerance was too high")
    end
    if statistics
        stat = Dict(
            :nmult => iwsp[1],
            :nexph => iwsp[2],
            :nscale => iwsp[3],
            :nstep => iwsp[4],
            :nreject => iwsp[5],
            :ibrkflag => iwsp[6],
            :mbrkdwn => iwsp[7],
            :step_min => wsp[1],
            :step_max => wsp[2],
            :x_error => wsp[5],
            :s_error => wsp[6],
            :zbrkdwn => wsp[7],
            :hump => wsp[8],
            :scaled_norm_sol => wsp[10],
         )
         return w, stat
    end
    return w
end       


function matvec{T<:AbstractArray{Float64,2}}(v_::Ptr{Float64}, w_::Ptr{Float64}, A_::Ptr{T})
    A = unsafe_pointer_to_objref(A_)::T
    n = size(A,2)
    v = pointer_to_array(v_, n)    
    w = pointer_to_array(w_, n) 
    Base.A_mul_B!(w, A, v) 
    # This *inplace* matrix-vector-multiplication instead of w[:] = A*v
    # results in a noticeable performance improvement
    return nothing
end   


function expv{T<:AbstractArray{Float64,2}}(t::Real, A::T, v::Vector{Float64}; 
              tol::Real=0.0, m::Integer=30, symmetric::Bool=isa(A, Hermitian), trace::Bool=false, anorm::Real=-1.0, statistics::Bool=false)
    if anorm<=0.0
        anorm = norm(A, Inf)
    end
    cmatvec = cfunction(matvec, Void, (Ptr{Float64}, Ptr{Float64}, Ptr{T}))
    arg = pointer_from_objref(A)
    _expv(t, cmatvec, v, anorm, tol=tol, m=m, symmetric=symmetric, trace=trace, statistics=statistics, arg=arg)
end   


function _expv(t::Real, matvec::Ptr{Void}, v::Vector{Complex{Float64}}, anorm::Real; 
               tol::Real=0.0, m::Integer=30, hermitian::Bool=false, trace::Bool=false, statistics::Bool=false, arg::Ptr{Void}=convert(Ptr{Void},0))
    n = length(v)
    w = zeros(Complex{Float64}, n)
    lwsp =  max(10, n*(m+1)+n+(m+2)^2+4*(m+2)^2+6+1)
    wsp = zeros(Complex{Float64}, lwsp)
    liwsp = max(7, m+2)
    iwsp = zeros(Int32, liwsp)
    iflag = zero(Int32) 
    if hermitian
        ccall(Libdl.dlsym(libexpokit, :zhexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32},            Ptr{Void}), 
         &n,           &m,                    &t,           v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec,
         &trace,       &iflag,                arg )
    else
        ccall(Libdl.dlsym(libexpokit, :zgexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32},            Ptr{Void}), 
         &n,           &m,                    &t,           v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec,
         &trace,       &iflag,                arg )
    end
    if ccall(Libdl.dlsym(libexpokit, :has_stopped), Int32, ()) == 1
        error("expokit fortran library stop")
    end    
    if iflag<0
        error("bad input arguments")
    elseif iflag==1
        error("maximum number of steps reached without convergence")
    elseif iflag==2
        error("requested tolerance was too high")
    end
    if statistics
        stat = Dict(
            :nmult => iwsp[1],
            :nexph => iwsp[2],
            :nscale => iwsp[3],
            :nstep => iwsp[4],
            :nreject => iwsp[5],
            :ibrkflag => iwsp[6],
            :mbrkdwn => iwsp[7],
            :step_min => real(wsp[1]),
            :step_max => real(wsp[2]),
            :x_error => real(wsp[5]),
            :s_error => real(wsp[6]),
            :zbrkdwn => real(wsp[7]),
            :hump => real(wsp[8]),
            :scaled_norm_sol => real(wsp[10]),
         )
         return w, stat
    end
    return w
end              


function matvec{T<:AbstractArray{Complex{Float64},2}}(v_::Ptr{Complex{Float64}}, w_::Ptr{Complex{Float64}}, A_::Ptr{T})
    A = unsafe_pointer_to_objref(A_)::T
    n = size(A,2)
    v = pointer_to_array(v_, n)    
    w = pointer_to_array(w_, n) 
    Base.A_mul_B!(w, A, v) 
    # This *inplace* matrix-vector-multiplication instead of w[:] = A*v
    # results in a noticeable performance improvement
    return nothing
end   

function expv{T<:AbstractArray{Complex{Float64},2}}(t::Real, A::T, v::Vector{Complex{Float64}}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=isa(A, Hermitian), trace::Bool=false, anorm::Real=-1.0, statistics::Bool=false)
    if anorm<=0.0
        anorm = norm(A, Inf)
    end
    cmatvec = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{T}))
    arg = pointer_from_objref(A)
    _expv(t, cmatvec, v, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, statistics=statistics, arg=arg)
end  

function expv{T<:AbstractArray{Complex{Float64},2}}(t::Real, A::T, v::Vector{Float64}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=isa(A, Hermitian), trace::Bool=false, anorm::Real=-1.0, statistics::Bool=false)
    if anorm<=0.0
        anorm = norm(A, Inf)
    end
    v1 = v+0im #complexify
    cmatvec = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{T}))
    arg = pointer_from_objref(A)
    _expv(t, cmatvec, v1, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, statistics=statistics, arg=arg)
end  


#Case real matrix A but complex vector v needs its own treatment (of course, complex computation required):

function matvec{T<:AbstractArray{Float64,2}}(v_::Ptr{Complex{Float64}}, w_::Ptr{Complex{Float64}}, A_::Ptr{T})
    A = unsafe_pointer_to_objref(A_)::T
    n = size(A,2)
    v = pointer_to_array(v_, n)    
    w = pointer_to_array(w_, n) 
    Base.A_mul_B!(w, A, v) 
    # This *inplace* matrix-vector-multiplication instead of w[:] = A*v
    # results in a noticeable performance improvement
    return nothing
end 

function expv{T<:AbstractArray{Float64,2}}(t::Real, A::T, v::Vector{Complex{Float64}}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=isa(A, Hermitian), trace::Bool=false, anorm::Real=-1.0, statistics::Bool=false)
    if anorm<=0.0
        anorm = norm(A, Inf)
    end
    cmatvec = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{T}))
    arg = pointer_from_objref(A)
    _expv(t, cmatvec, v, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, statistics=statistics, arg=arg)
end  


# 

function funvec{T<:Union{Float64,Complex{Float64}}}(v_::Ptr{T}, w_::Ptr{T}, nF_::Ptr{Tuple{Int,Function}})
    nF = unsafe_pointer_to_objref(nF_)::Tuple{Int,Function}
    n = nF[1]
    F = nF[2]
    v = pointer_to_array(v_, n)    
    w = pointer_to_array(w_, n) 
    w[:] = F(v)
    return nothing
end 

function expv(t::Real, F::Function, v::Vector{Float64}, anorm::Real;
    tol::Real=0.0, m::Integer=30, symmetric::Bool=false, trace::Bool=false, statistics::Bool=false)
    nF = (length(v), F) 
    cfunvec = cfunction(funvec, Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Tuple{Int,Function}}))
    arg = pointer_from_objref(nF)
    _expv(t, cfunvec, v, anorm, tol=tol, m=m, symmetric=symmetric, trace=trace, statistics=statistics, arg=arg)
end

function expv(t::Real, F::Function, v::Vector{Complex{Float64}}, anorm::Real;
    tol::Real=0.0, m::Integer=30, hermitian::Bool=false, trace::Bool=false, statistics::Bool=false)
    nF = (length(v), F) 
    cfunvec = cfunction(funvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Tuple{Int,Function}}))
    arg = pointer_from_objref(nF)
    _expv(t, cfunvec, v, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, statistics=statistics, arg=arg)
end

#include("acroy.jl")


end # module Expokit
