__precompile__()

module Expokit

export expv


function __init__()
    global const libexpokit = Libdl.dlopen(joinpath(dirname(@__FILE__),
                            "..", "deps", "lib", string("libexpokit.", Libdl.dlext)))
    ccall(Libdl.dlsym(libexpokit, :fortran_stop_handler_init), Void, () )
end


function _expv(t::Real, matvec::Function, v::Vector{Float64}, anorm::Real; 
              tol::Real=0.0, m::Integer=30, symmetric::Bool=false, trace::Bool=false, statistics::Bool=false)
    n = length(v)
    w = zeros(Float64, n)
    lwsp =  max(10, n*(m+1)+n+(m+2)^2+4*(m+2)^2+6+1)
    wsp = zeros(Float64, lwsp)
    liwsp = max(7, m+2)
    iwsp = zeros(Int32, liwsp)
    iflag = Int32(0) 
    matvec_fortran = cfunction(matvec, Void, (Ptr{Float64}, Ptr{Float64}))
    if symmetric
        ccall(Libdl.dlsym(libexpokit, :dsexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},   Ptr{Float64}, Ptr{Float64},   Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32},     Ptr{Int32},   Ptr{Void},      Ptr{Int32},  Ptr{Int32}), 
         &n,           &m,           &t,           v,              w,            &tol, 
         &anorm,       wsp,          &lwsp,        iwsp,           &liwsp,       matvec_fortran, &trace,      &iflag)
    else
        ccall(Libdl.dlsym(libexpokit, :dgexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},   Ptr{Float64}, Ptr{Float64},   Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32},     Ptr{Int32},   Ptr{Void},      Ptr{Int32},  Ptr{Int32}), 
         &n,           &m,           &t,           v,              w,            &tol, 
         &anorm,       wsp,          &lwsp,        iwsp,           &liwsp,       matvec_fortran, &trace,      &iflag)
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

function _matvec_(v::Ptr{Float64}, w::Ptr{Float64})
    global _A_, _n_
    v1 = pointer_to_array(v, _n_)    
    w1 = pointer_to_array(w, _n_) 
    w1[:] = _A_*v1
    return nothing
end   

function _matvec_Av_(v::Ptr{Float64}, w::Ptr{Float64})
    global _Av_, _n_
    v1 = pointer_to_array(v, _n_)    
    w1 = pointer_to_array(w, _n_) 
    w1[:] = _Av_(v1)
    return nothing
end   

function expv(t::Real, A::AbstractArray{Float64,2}, v::Vector{Float64}; 
              tol::Real=0.0, m::Integer=30, symmetric::Bool=isa(A, Hermitian), trace::Bool=false, anorm::Real=-1.0, statistics::Bool=false)
    global _A_, _n_
    _n_ = length(v)
    _A_ = A
    if anorm<=0.0
        anorm = norm(A, Inf)
    end
    _expv(t, _matvec_, v, anorm, tol=tol, m=m, symmetric=symmetric, trace=trace, statistics=statistics)
end   

function expv(t::Real, Av::Function, v::Vector{Float64}, anorm::Real;
              tol::Real=0.0, m::Integer=30, symmetric::Bool=isa(A, Hermitian), trace::Bool=false, statistics::Bool=false)
    global _Av_, _n_
    _n_ = length(v)
    _Av_ = Av
    _expv(t, _matvec_Av_, v, anorm, tol=tol, m=m, symmetric=symmetric, trace=trace, statistics=statistics)
end    



function _expv(t::Real, matvec::Function, v::Vector{Complex{Float64}}, anorm::Real; 
               tol::Real=0.0, m::Integer=30, hermitian::Bool=false, trace::Bool=false, statistics::Bool=false)
    n = length(v)
    w = zeros(Complex{Float64}, n)
    lwsp =  max(10, n*(m+1)+n+(m+2)^2+4*(m+2)^2+6+1)
    wsp = zeros(Complex{Float64}, lwsp)
    liwsp = max(7, m+2)
    iwsp = zeros(Int32, liwsp)
    iflag = zero(Int32) 
    matvec_fortran = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}))
    if hermitian
        ccall(Libdl.dlsym(libexpokit, :zhexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32}), 
         &n,           &m,                    &t,           v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec_fortran,
         &trace,       &iflag)
    else
        ccall(Libdl.dlsym(libexpokit, :zgexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32}), 
         &n,           &m,                    &t,           v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec_fortran,
         &trace,       &iflag)
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


function _matvec_cmplx_(v::Ptr{Complex{Float64}}, w::Ptr{Complex{Float64}})
    global _A_, _n_
    v1 = pointer_to_array(v, _n_)    
    w1 = pointer_to_array(w, _n_) 
    w1[:] = _A_*v1
    return nothing
end    

function _matvec_Av_cmplx_(v::Ptr{Complex{Float64}}, w::Ptr{Complex{Float64}})
    global _Av_, _n_
    v1 = pointer_to_array(v, _n_)    
    w1 = pointer_to_array(w, _n_) 
    w1[:] = _Av_(v1)
    return nothing
end    

function expv(t::Real, A::Union{AbstractArray{Float64,2},AbstractArray{Complex{Float64},2}}, 
              v::Vector{Complex{Float64}}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=isa(A, Hermitian), trace::Bool=false, anorm::Real=-1.0, statistics::Bool=false)
    global _A_, _n_
    _n_ = length(v)
    _A_ = A
    if anorm<=0.0
        anorm = norm(A, Inf)
    end
    _expv(t, _matvec_cmplx_, v, anorm, tol=tol, m=m, symmetric=symmetric, trace=trace, statistics=statistics)
end   

function expv(t::Real, A::AbstractArray{Complex{Float64},2}, 
              v::Vector{Float64}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=isa(A, Hermitian), trace::Bool=false, anorm::Real=-1.0, statistics::Bool=false)
    global _A_, _n_
    _n_ = length(v)
    _A_ = A
    v1 = v+0im #complexify
    if anorm<=0.0
        anorm = norm(A, Inf)
    end
    _expv(t, _matvec_cmplx_, v1, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, statistics=statistics)
end   

function expv(t::Real, Av::Function, v::Vector{Complex{Float64}}, anorm::Real; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=isa(A, Hermitian), trace::Bool=false, statistics::Bool=false)
    global _Av_, _n_
    _n_ = length(v)
    _Av_ = Av
    _expv(t, _matvec_Av_cmplx_, v, anorm, tol=tol, m=m, symmetric=symmetric, trace=trace, statistics=statistics)
end   


#include("acroy.jl")


end # module Expokit
