__precompile__()

module FExpokit

export expv, phiv, expv!, phiv!
export get_lwsp_liwsp_expv, get_lwsp_liwsp_phiv

function __init__()
    global const libexpokit = Libdl.dlopen(joinpath(dirname(@__FILE__),
                            "..", "deps", "lib", string("libexpokit.", Libdl.dlext)))
    ccall(Libdl.dlsym(libexpokit, :fortran_stop_handler_init), Void, () )
end


get_lwsp_liwsp_expv(n::Integer, m::Integer=30) = ( max(10, n*(m+1)+n+(m+2)^2+4*(m+2)^2+6+1), max(7, m+2) )
get_lwsp_liwsp_phiv(n::Integer, m::Integer=30) = ( max(10, n*(m+1)+n+(m+3)^2+4*(m+3)^2+6+1), max(7, m+3) )


function _expv_real!(w::Vector{Float64}, t::Real, matvec::Ptr{Void}, v::Vector{Float64}, anorm::Real; 
              tol::Real=0.0, m::Integer=30, symmetric::Bool=false, trace::Bool=false, statistics::Bool=false, arg::Ptr{Void}=convert(Ptr{Void},0))
    n = length(v)
    m = min(m, n-1)
    lwsp, liwsp = get_lwsp_liwsp_expv(n, m)
    wsp = zeros(Float64, lwsp)
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
    if trace
        Libc.flush_cstdio()
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
            :tbrkdwn => wsp[7],
            :t_now => wsp[8],
            :hump => wsp[9],
            :scaled_norm_sol => wsp[10],
         )
         return w, stat
    end
    return w
end       


function _phiv_real!(w::Vector{Float64}, t::Real, matvec::Ptr{Void}, u::Vector{Float64}, v::Vector{Float64}, anorm::Real; 
                     tol::Real=0.0, m::Integer=30, symmetric::Bool=false, trace::Bool=false, statistics::Bool=false, arg::Ptr{Void}=convert(Ptr{Void},0))
    n = length(v)
    m = min(m, n-1)
    lwsp, liwsp = get_lwsp_liwsp_phiv(n, m)
    wsp = zeros(Float64, lwsp)
    iwsp = zeros(Int32, liwsp)
    iflag = Int32(0)
    if symmetric
        ccall(Libdl.dlsym(libexpokit, :dsphiv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},   Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32}, Ptr{Int32},   Ptr{Void},    Ptr{Int32},  Ptr{Int32}, Ptr{Void}), 
         &n,           &m,           &t,           u,          v,            w,            &tol, 
         &anorm,       wsp,          &lwsp,        iwsp,       &liwsp,       matvec,       &trace,      &iflag,     arg)
    else
        ccall(Libdl.dlsym(libexpokit, :dgphiv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},   Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32},   Ptr{Int32}, Ptr{Void},    Ptr{Int32},  Ptr{Int32}, Ptr{Void}), 
         &n,           &m,           &t,           u,            v,          w,            &tol, 
         &anorm,       wsp,          &lwsp,        iwsp,         &liwsp,     matvec,       &trace,      &iflag,     arg)
    end
    if trace
        Libc.flush_cstdio()
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
            :tbrkdwn => wsp[7],
            :t_now => wsp[8],
         )
         return w, stat
    end
    return w
end     


function _expv_cmplx!(w::Vector{Complex{Float64}},t::Real, matvec::Ptr{Void}, v::Vector{Complex{Float64}}, anorm::Real; 
               tol::Real=0.0, m::Integer=30, hermitian::Bool=false, trace::Bool=false, statistics::Bool=false, 
               arg::Ptr{Void}=convert(Ptr{Void},0), matrix_times_minus_i::Bool=false,
               wsp::Array{Complex{Float64},1}=Complex{Float64}[], iwsp::Array{Int32,1}=Int32[])
    n = length(v)
    m = min(m, n-1)
    lwsp, liwsp = get_lwsp_liwsp_expv(n, m)
    if length(wsp)==0
        wsp = zeros(Complex{Float64}, lwsp)
    elseif length(wsp)<lwsp
        error("wsp too small")
    end
    if length(iwsp)==0
        iwsp = zeros(Int32, liwsp)
    elseif length(iwsp)<liwsp
        error("iwsp too small")
    end
    iflag = zero(Int32) 
    imia = matrix_times_minus_i?1:0
    if hermitian
        ccall(Libdl.dlsym(libexpokit, :zhexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32},            Ptr{Void},    Ptr{Int32}), 
         &n,           &m,                    &t,           v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec,
         &trace,       &iflag,                arg,          &imia)
    else
        ccall(Libdl.dlsym(libexpokit, :zgexpv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32},            Ptr{Void},    Ptr{Int32}), 
         &n,           &m,                    &t,           v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec,
         &trace,       &iflag,                arg,          &imia )
    end
    if trace
        Libc.flush_cstdio()
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
            :tbrkdwn => real(wsp[7]),
            :t_now => real(wsp[8]),
            :hump => real(wsp[9]),
            :scaled_norm_sol => real(wsp[10]),
         )
         return w, stat
    end
    return w
end     


function _phiv_cmplx!(w::Vector{Complex{Float64}}, t::Real, matvec::Ptr{Void}, u::Vector{Complex{Float64}}, v::Vector{Complex{Float64}}, anorm::Real; 
               tol::Real=0.0, m::Integer=30, hermitian::Bool=false, trace::Bool=false, statistics::Bool=false, arg::Ptr{Void}=convert(Ptr{Void},0))
    n = length(v)
    m = min(m, n-1)
    lwsp, liwsp = get_lwsp_liwsp_phiv(n, m)
    wsp = zeros(Complex{Float64}, lwsp)
    iwsp = zeros(Int32, liwsp)
    iflag = zero(Int32) 
    if hermitian
        ccall(Libdl.dlsym(libexpokit, :zhphiv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32},            Ptr{Void}), 
         &n,           &m,                    &t,           u,                     v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec,
         &trace,       &iflag,                arg )
    else
        ccall(Libdl.dlsym(libexpokit, :zgphiv_wrap), Void, 
        (Ptr{Int32},   Ptr{Int32},            Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Float64},
         Ptr{Float64}, Ptr{Complex{Float64}}, Ptr{Int32},   Ptr{Int32},            Ptr{Int32},            Ptr{Void},
         Ptr{Int32},   Ptr{Int32},            Ptr{Void}), 
         &n,           &m,                    &t,           u,                     v,                     w,                     &tol, 
         &anorm,       wsp,                   &lwsp,        iwsp,                  &lwsp,                 matvec,
         &trace,       &iflag,                arg )
    end
    if trace
        Libc.flush_cstdio()
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
            :tbrkdwn => real(wsp[7]),
            :t_now => real(wsp[8]),
         )
         return w, stat
    end
    return w
end





function matvec{T}(v_::Ptr{Float64}, w_::Ptr{Float64}, A_::Ptr{T})
    A = unsafe_pointer_to_objref(A_)::T
    n = size(A,2)
    v = unsafe_wrap(Array, v_, n)
    w = unsafe_wrap(Array, w_, n)
    Base.A_mul_B!(w, A, v) 
    # This *inplace* matrix-vector-multiplication instead of w[:] = A*v
    # results in a noticeable performance improvement
    # and makes it more general
    return nothing
end   

function matvec{T}(v_::Ptr{Complex{Float64}}, w_::Ptr{Complex{Float64}}, A_::Ptr{T})
    A = unsafe_pointer_to_objref(A_)::T
    n = size(A,2)
    v = unsafe_wrap(Array, v_, n)
    w = unsafe_wrap(Array, w_, n)
    Base.A_mul_B!(w, A, v) 
    # This *inplace* matrix-vector-multiplication instead of w[:] = A*v
    # results in a noticeable performance improvement
    # and makes it more general
    return nothing
end   


function expv{T}(t::Real, A::T, v::Vector{Complex{Float64}}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=ishermitian(A), 
              trace::Bool=false, anorm::Real=norm(A, Inf), statistics::Bool=false,
              matrix_times_minus_i::Bool=false,
              wsp::Array{Complex{Float64},1}=Complex{Float64}[], iwsp::Array{Int32,1}=Int32[])
    w = zeros(Complex{Float64}, length(v))
    cmatvec = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{T}))
    arg = pointer_from_objref(A)
    _expv_cmplx!(w, t, cmatvec, v, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, 
                 statistics=statistics, arg=arg, matrix_times_minus_i=matrix_times_minus_i,
                 wsp=wsp, iwsp=iwsp)
end  

function expv!{T}(w::Vector{Complex{Float64}}, t::Real, A::T, v::Vector{Complex{Float64}}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=ishermitian(A), 
              trace::Bool=false, anorm::Real=norm(A, Inf), statistics::Bool=false,
              matrix_times_minus_i::Bool=false,
              wsp::Array{Complex{Float64},1}=Complex{Float64}[], iwsp::Array{Int32,1}=Int32[])
    cmatvec = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{T}))
    arg = pointer_from_objref(A)
    _expv_cmplx!(w, t, cmatvec, v, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, 
                 statistics=statistics, arg=arg, matrix_times_minus_i=matrix_times_minus_i,
                 wsp=wsp, iwsp=iwsp)
end  


function expv{T}(t::Real, A::T, v::Vector{Float64}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=ishermitian(A), trace::Bool=false, anorm::Real=norm(A, Inf), statistics::Bool=false)
    cmatvec = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{T}))
    if eltype(A)==Complex{Float64}
        w = zeros(Complex{Float64}, length(v))
        v1 = v+0im #complexify
        cmatvec = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{T}))
        arg = pointer_from_objref(A)
        return _expv_cmplx!(w, t, cmatvec, v1, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, statistics=statistics, arg=arg)
    else
        w = zeros(Float64, length(v))
        cmatvec = cfunction(matvec, Void, (Ptr{Float64}, Ptr{Float64}, Ptr{T}))
        arg = pointer_from_objref(A)
        return _expv_real!(w, t, cmatvec, v, anorm, tol=tol, m=m, symmetric=hermitian, trace=trace, statistics=statistics, arg=arg)
    end
end   

function expv!{T}(w::Union{Vector{Float64}, Vector{Complex{Float64}}}, t::Real, A::T, v::Vector{Float64}; 
              tol::Real=0.0, m::Integer=30, hermitian::Bool=ishermitian(A), trace::Bool=false, anorm::Real=norm(A, Inf), statistics::Bool=false)
    cmatvec = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{T}))
    if eltype(A)==Complex{Float64}
        v1 = v+0im #complexify
        cmatvec = cfunction(matvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{T}))
        arg = pointer_from_objref(A)
        return _expv_cmplx!(w, t, cmatvec, v1, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, statistics=statistics, arg=arg)
    else
        cmatvec = cfunction(matvec, Void, (Ptr{Float64}, Ptr{Float64}, Ptr{T}))
        arg = pointer_from_objref(A)
        return _expv_real!(w, t, cmatvec, v, anorm, tol=tol, m=m, symmetric=hermitian, trace=trace, statistics=statistics, arg=arg)
    end
end   


#TODO: implement all variants of phiv
#function phiv{T}(t::Real, A::T, u::Vector{Float64}, v::Vector{Float64}; 
#              tol::Real=0.0, m::Integer=30, hermitian::Bool=ishermitian(A), trace::Bool=false, anorm::Real=norm(A, Inf), statistics::Bool=false)
#    cmatvec = cfunction(matvec, Void, (Ptr{Float64}, Ptr{Float64}, Ptr{T}))
#    arg = pointer_from_objref(A)
#    _phiv_real(t, cmatvec, u, v, anorm, tol=tol, m=m, symmetric=hermitian, trace=trace, statistics=statistics, arg=arg)
#end   


function funvec(v_::Ptr{Complex{Float64}}, w_::Ptr{Complex{Float64}}, nF_::Ptr{Tuple{Int,Function,Tuple}})
    nF = unsafe_pointer_to_objref(nF_)::Tuple{Int,Function, Tuple}
    n = nF[1]
    F! = nF[2]
    args = nF[3]
    v = unsafe_wrap(Array, v_, n)
    w = unsafe_wrap(Array, w_, n)
    F!(w, v, args...)
    return nothing
end 

#function expv!(w::Vector{Float64},t::Real, F::Function, v::Vector{Float64}, anorm::Real; 
#               args::Tuple=(), tol::Real=0.0, m::Integer=30, symmetric::Bool=false, trace::Bool=false, statistics::Bool=false)
#    w = zeros(Float64, length(v))
#    nF = (length(v), F, args) 
#    cfunvec = cfunction(funvec, Void, (Ptr{Float64}, Ptr{Float64}, Ptr{Tuple{Int,Function}}))
#    arg = pointer_from_objref(nF)
#    _expv_real!(w, t, cfunvec, v, anorm, tol=tol, m=m, symmetric=symmetric, trace=trace, statistics=statistics, arg=arg)
#end

function expv!(w::Vector{Complex{Float64}},t::Real, F::Function, v::Vector{Complex{Float64}}, anorm::Real;
               args::Tuple=(), tol::Real=0.0, m::Integer=30, hermitian::Bool=false, trace::Bool=false, 
               statistics::Bool=false, matrix_times_minus_i::Bool=false,
               wsp::Array{Complex{Float64},1}=Complex{Float64}[], iwsp::Array{Int32,1}=Int32[])
    nF = (length(v), F, args) 
    cfunvec = cfunction(funvec, Void, (Ptr{Complex{Float64}}, Ptr{Complex{Float64}}, Ptr{Tuple{Int,Function, Tuple}}))
    arg = pointer_from_objref(nF)
    _expv_cmplx!(w, t, cfunvec, v, anorm, tol=tol, m=m, hermitian=hermitian, trace=trace, statistics=statistics, arg=arg,
                matrix_times_minus_i=matrix_times_minus_i, wsp=wsp, iwsp=iwsp)
end



end # module FExpokit
