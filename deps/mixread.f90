module g2k_c_binding
    implicit none

contains

    function strlen(str)
        use, intrinsic :: iso_c_binding

        character(kind=c_char), intent(in) :: str(*)
        integer :: strlen

        strlen = 0
        do
            if(str(strlen+1) == c_null_char) then
                exit
            endif
            strlen = strlen + 1
        enddo
    end

    function from_cstring(str)
        use, intrinsic :: iso_c_binding

        character(kind=c_char), intent(in) :: str(*)
        integer :: length
        character(:), allocatable :: from_cstring

        length = strlen(str)
        allocate(character(length) :: from_cstring)
        from_cstring = transfer(str(1:length), from_cstring)
    end

end module g2k_c_binding


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

    integer :: i, j, ivecdum
    real(real64) :: eav
    !integer, allocatable ::
    real(c_double), pointer :: eval(:), evec(:,:)

    ! TESTING
    real(real64) :: rowsum
    real(real64), allocatable :: colsum(:)
    integer(c_int), pointer :: ints(:)

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
    print *, "<pre>blockinfo%nelectrons", blockinfo%nelectrons
    read(fhandle) blockinfo%nelectrons, blockinfo%ncsfstotal, blockinfo%norbitals, &
        blockinfo%nvectotal, blockinfo%nvecsize, blockinfo%nblocks
    print *, "blockinfo%nelectrons", blockinfo%nelectrons
    print *, "blockinfo%ncsfstotal", blockinfo%ncsfstotal
    print *, "blockinfo%norbitals", blockinfo%norbitals
    print *, "blockinfo%nvectotal", blockinfo%nvectotal
    print *, "blockinfo%nvecsize", blockinfo%nvecsize
    print *, "blockinfo%nblocks", blockinfo%nblocks

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

        !CALL alloc (pnteval, nevblk, 8)
        !CALL alloc (pntevec, nevblk*ncfblk, 8)
        !CALL alloc (pntiset, ncfblk, 4)
        allocate(eval(nevblk), evec(ncfblk, nevblk))
        blocks(ib)%eigenstates = c_loc(evec)

        allocate(colsum(nevblk)) ! TEST

        read(fhandle) (ivecdum, i = 1, nevblk)
        read(fhandle) blocks(ib)%eav, (eval(i), i = 1, nevblk)
        read(fhandle) ((evec(i,j), i = 1, ncfblk), j = 1, nevblk)

        print *, blocks(ib)%eav
        print *, size(eval), eval

        print *, "evec=", size(evec), shape(evec)
        do j = 1, nevblk
            colsum(j) = 0.0_dp
        enddo
        do i = 1, ncfblk
            eav = 0.0_dp
            do j = 1, nevblk
                write(*, "(f10.5)", advance="no") evec(i,j)
                eav = eav + abs(evec(i,j))**2
                colsum(j) = colsum(j) + abs(evec(i,j))**2
            enddo
            write(*, "(' -> ', f10.5)", advance="no") eav
            print *
        enddo

        do j = 1, nevblk
            write(*, "(f10.5)", advance="no") colsum(j)
        enddo
        print *

        deallocate(eval)
        deallocate(colsum) ! TEST
    enddo

    close(fhandle)

end subroutine mixread
