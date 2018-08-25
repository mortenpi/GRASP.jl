function mixread(filename_cstr, blockinfo, blocks_ptr) bind(c)
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
    integer(c_int) :: mixread

    character(:), allocatable :: filename

    integer :: fhandle, ios
    character(255) :: iom
    character(6) :: g92mix
    integer :: errlineno = 0

    integer :: ib, nb, ncfblk, nevblk, iatjp, iaspa

    integer :: i, j, dummy_int
    real(c_double), pointer :: eval(:), evec(:,:)

    filename = from_cstring(filename_cstr)
    open(newunit=fhandle, file=filename, form="unformatted", status="old", iostat=ios, iomsg=iom)
    if(ios /= 0) then
        print *, "ERROR: Unable to open file:", ios, iom
        mixread = 1
        return
    endif

    ! Check header
    read(fhandle, iostat=ios, iomsg=iom) g92mix
    if(g92mix /= "G92MIX") then
        print *, "ERROR: Bad file header -- not a G92MIX file?"
        mixread = 3
        return
    endif

    ! Read the global information of the mixing file (number of block etc.).
    ! Corresponds to the following READ in rmixextract:
    !
    !     READ (nfmix) nelec, ncftot, nw, nvectot, nvecsiz, nblock
    !
    errlineno=__LINE__; read(fhandle, iostat=ios, iomsg=iom) &
        blockinfo%nelectrons, blockinfo%ncsfstotal, blockinfo%norbitals, &
        blockinfo%nvectotal, blockinfo%nvecsize, blockinfo%nblocks
    if(ios /= 0) go to 999

    allocate(blocks(blockinfo%nblocks))
    blocks_ptr = c_loc(blocks)

    do ib = 1, blockinfo%nblocks
        ! ncfblk -- total number of CSFs in the block
        ! nevblk -- number of states in the block
        ! eav -- energy of the state
        read(fhandle, iostat=ios, iomsg=iom) nb, ncfblk, nevblk, blocks(ib)%iatjp, blocks(ib)%iaspa
        if(ios /= 0) go to 999

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

        errlineno=__LINE__; read(fhandle, iostat=ios, iomsg=iom) (dummy_int, i = 1, nevblk)
        if(ios /= 0) go to 999
        errlineno=__LINE__; read(fhandle, iostat=ios, iomsg=iom) blocks(ib)%eav, (eval(i), i = 1, nevblk)
        if(ios /= 0) go to 999
        errlineno=__LINE__; read(fhandle, iostat=ios, iomsg=iom) ((evec(i,j), i = 1, ncfblk), j = 1, nevblk)
        if(ios /= 0) go to 999

        do j = 1, nevblk
            eval(j) = blocks(ib)%eav + eval(j)
        enddo
    enddo

    close(fhandle)

    mixread = 0 ! no error
    return

    ! Error handling for IO errors (reachable via goto)
    999 continue
    print '(a)', "Terminating mixread() with IO error."
    print '(a,a,":",i0)', " at: ", __FILE__, errlineno
    print '(" while reading: ",a)', filename
    print *, ios, iom
    close(fhandle)
    mixread = 2
    return

end function mixread
