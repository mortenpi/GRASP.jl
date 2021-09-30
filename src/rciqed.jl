module rciqed

librciqed_so = grasp_path_file = joinpath(@__DIR__, "../deps/libgrasp-rci-qed.so")

function parameter_def(name::Symbol)
    value = Ref{Cint}()
    rv = ccall(
        (:libgrasp_rciqed_parameter_def, librciqed_so), Cint,
        (Cstring, Ref{Cint}),
        name, value,
    )
    rv == 0 || error("Error in libgrasp_rciqed_parameter_def(\"$name\") -> $rv")
    return value[]
end

parameter_def() = NamedTuple(
    name => parameter_def(name)
    for name in (:KEYORB, :NNNP, :NNN1, :NNNW, :NNNWM1, :NNNWM2)
)

const __GRID_INITIALIZED__ = Ref{Bool}(false)

function init(z::Real)
    global __GRID_INITIALIZED__
    ccall(
        (:libgrasp_rciqed_9290_init_pnc, librciqed_so), Cvoid,
        (Cdouble,),
        z,
    )
    __GRID_INITIALIZED__[] = true
    return
end

function grid()
    __GRID_INITIALIZED__[] || error("lib9290 grid not initialized.")
    n, rnt, h, hp = Ref{Cint}(), Ref{Cdouble}(), Ref{Cdouble}(), Ref{Cdouble}()
    r, rp, rpor = Ref{Ptr{Cdouble}}(), Ref{Ptr{Cdouble}}(), Ref{Ptr{Cdouble}}()
    ccall(
        (:libgrasp_rciqed_grid_c, librciqed_so), Cvoid,
        (Ref{Cint}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Ptr{Cdouble}}, Ref{Ptr{Cdouble}}, Ref{Ptr{Cdouble}}),
        n, rnt, h, hp, r, rp, rpor,
    )
    (
        n = n[], rnt = rnt[], h = h[], hp = hp[],
        r = unsafe_wrap(Vector{Cdouble}, r[], n[], own=true),
        rp = unsafe_wrap(Vector{Cdouble}, rp[], n[], own=true),
        rpor = unsafe_wrap(Vector{Cdouble}, rpor[], n[], own=true),
    )
end

function funk(x::Real, n::Integer)
    (x == 0) && (n != 0) && throw(DomainError((x, n), "n must be 0 if x = 0"))
    (n ∉ [0, 1, 3, 5]) && throw(DomainError((x, n), "n must be one of 0, 1, 3, 5"))
    ccall(
        (:libgrasp_rciqed_vp_funk, librciqed_so), Cdouble,
        (Cdouble, Cint),
        x, n,
    )
end

function funl(x::Real, k::Integer)
    (x == 0) && (k != 0) && throw(DomainError((x, k), "k must be 0 if x = 0"))
    (k ∉ [0, 1]) && throw(DomainError((x, k), "k must be one of 0, 1"))
    ccall(
        (:libgrasp_rciqed_vp_funl, librciqed_so), Cdouble,
        (Cdouble, Cint),
        x, k,
    )
end

end
