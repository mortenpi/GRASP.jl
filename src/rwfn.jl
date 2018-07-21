struct OrbitalF90
    npy  :: Cint
    naky :: Cint
    ey   :: Cdouble
    my   :: Cint
    pa   :: Ptr{Cdouble}
    qa   :: Ptr{Cdouble}
    ra   :: Ptr{Cdouble}
end

struct RWFNOrbital
    orbital_n :: Int
    orbital_kappa :: Int
    r :: Vector{Float64}
    p :: Vector{Float64}
    q :: Vector{Float64}
    # ....
    original :: OrbitalF90
    #energy :: Float64
end

function read_rwfn(filename)
    isfile(filename) || error("Unable to open $(filename)")

    orbitals = Ref{Ptr{OrbitalF90}}()
    norbitals = Ref{Cint}()
    ccall( (:rwfnread, libgrasp_so),
        Cvoid, (Cstring, Ref{Cint}, Ref{Ptr{OrbitalF90}}),
        filename, norbitals, orbitals
    )
    orbitals = LCompat.unsafe_wrap(Array{OrbitalF90}, orbitals.x, norbitals.x, own=true)
    map(orbitals) do orb
        RWFNOrbital(
            orb.npy, orb.naky,
            LCompat.unsafe_wrap(Array{Cdouble}, orb.ra, (orb.my,), own=true),
            LCompat.unsafe_wrap(Array{Cdouble}, orb.pa, (orb.my,), own=true),
            LCompat.unsafe_wrap(Array{Cdouble}, orb.qa, (orb.my,), own=true),
            orb
        )
    end
end
