subroutine mixread(filename_cstr, blockinfo, blocks_ptr) bind(c)
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: iso_c_binding
    use g2k_c_binding
    implicit none

    type, bind(c) :: blocksinfo_t
        integer(c_int) :: nblocks
        integer(c_int) :: nelectrons
        integer(c_int) :: ncsfstotal
        integer(c_int) :: norbitals
        integer(c_int) :: nvectotal
        integer(c_int) :: nvecsize
    end type

    type, bind(c) :: block_t
        integer(c_int) :: blockid
        integer(c_int) :: ncsfs
        integer(c_int) :: nevs
        integer(c_int) :: iatjp, iaspa
        real(c_double) :: eav
        type(c_ptr) :: eigenstates
        type(c_ptr) :: eigenenergies
    end type

    integer, parameter :: dp = real64

    character(kind=c_char), intent(in) :: filename_cstr(1)
    type(blocksinfo_t), intent(out) :: blockinfo
    type(c_ptr), intent(out) :: blocks_ptr
    type(block_t), pointer :: blocks(:)
    type(block_t) :: block

    character(:), allocatable :: filename

    integer :: fhandle, ios
    character(255) :: iom
    character(6) :: g92mix

    integer :: ib, nb, ncfblk, nevblk, iatjp, iaspa

    integer :: i, j, dummy_int
    real(c_double), pointer :: eval(:), evec(:,:)

    filename = from_cstring(filename_cstr)
    open(newunit=fhandle, file=filename, form="unformatted", status="old", iostat=ios, IOMSG=iom)
    if(ios /= 0) then
        print *, "ERROR: Unable to open file:", ios, iom
        stop 1
    endif

    ! Check header
    read(fhandle) g92mix
    if(g92mix /= "G92MIX") then
        print *, "ERROR: Bad header. Not a G92MIX file?"
        stop 1
    endif

    ! Read the global information of the mixing file (number of block etc.).
    ! Corresponds to the following READ in rmixextract:
    !
    !     READ (nfmix) nelec, ncftot, nw, nvectot, nvecsiz, nblock
    !
    read(fhandle) blockinfo%nelectrons, blockinfo%ncsfstotal, blockinfo%norbitals, &
        blockinfo%nvectotal, blockinfo%nvecsize, blockinfo%nblocks

    allocate(blocks(blockinfo%nblocks))
    blocks_ptr = c_loc(blocks)

    do ib = 1, blockinfo%nblocks
        ! ncfblk -- total number of CSFs in the block
        ! nevblk -- number of states in the block
        ! eav -- energy of the state
        read(fhandle) nb, ncfblk, nevblk, blocks(ib)%iatjp, blocks(ib)%iaspa

        blocks(ib)%blockid = nb
        blocks(ib)%ncsfs = ncfblk
        blocks(ib)%nevs = nevblk

        ! The original allocations:
        !   CALL alloc (pnteval, nevblk, 8)
        !   CALL alloc (pntevec, nevblk*ncfblk, 8)
        !   CALL alloc (pntiset, ncfblk, 4)
        allocate(eval(nevblk), evec(ncfblk, nevblk))
        blocks(ib)%eigenstates = c_loc(evec)
        blocks(ib)%eigenenergies = c_loc(eval)

        read(fhandle) (dummy_int, i = 1, nevblk)
        read(fhandle) blocks(ib)%eav, (eval(i), i = 1, nevblk)
        read(fhandle) ((evec(i,j), i = 1, ncfblk), j = 1, nevblk)

        do j = 1, nevblk
            eval(j) = blocks(ib)%eav + eval(j)
        enddo
    enddo

    close(fhandle)

end subroutine mixread
