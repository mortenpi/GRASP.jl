struct BlocksInfoF90
    nblocks    :: Cint
    nelectrons :: Cint
    ncsfstotal :: Cint
    norbitals  :: Cint
    nvectotal  :: Cint
    nvecsize   :: Cint
end

struct BlockF90
    blockid       :: Cint
    ncsfs         :: Cint
    nevs          :: Cint
    iatjp         :: Cint
    iaspa         :: Cint
    eav           :: Cdouble
    eigenstates   :: Ptr{Cdouble}
    eigenenergies :: Ptr{Cdouble}
end

struct MixingBlock
    eav :: Float64
    energies :: Vector{Float64}
    states :: Matrix{Float64}
    _block_f90 :: BlockF90
end

struct MixingFile
    nelectrons :: Int
    blocks :: Vector{MixingBlock}
    _blocksinfo_f90 :: BlocksInfoF90
end


function read_rmix(filename)
    isfile(filename) || error("Unable to open $(filename)")

    blockinfo = Ref{BlocksInfoF90}()
    blocks = Ref{Ptr{BlockF90}}()
    status = ccall( (:mixread, libgrasp_so),
        Cint, (Cstring, Ref{BlocksInfoF90}, Ref{Ptr{BlockF90}}),
        filename, blockinfo, blocks
    )
    if status != zero(status)
        error("mixread returned with status $(status)")
    end
    blockinfo = blockinfo.x
    blocks = LCompat.unsafe_wrap(Array{BlockF90}, blocks.x, blockinfo.nblocks, own=true)
    blocks = map(blocks) do block
        energies = LCompat.unsafe_wrap(Array{Cdouble}, block.eigenenergies, (block.nevs,), own=true)
        states = LCompat.unsafe_wrap(Array{Cdouble}, block.eigenstates, (block.ncsfs, block.nevs), own=true)
        MixingBlock(block.eav, energies, states, block)
    end
    MixingFile(blockinfo.nelectrons, blocks, blockinfo)
end
