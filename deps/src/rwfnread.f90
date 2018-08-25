function rwfnread(filename_cstr, norbitals, orbitals_ptr) bind(c)
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: iso_c_binding
    use g2k_c_binding
    implicit none

    type, bind(c) :: orbital_t
        integer(c_int) :: npy
        integer(c_int) :: naky
        real(c_double) :: ey
        integer(c_int) :: my
        type(c_ptr) :: pa
        type(c_ptr) :: qa
        type(c_ptr) :: ra
    end type

    integer, parameter :: dp = real64

    character(kind=c_char), intent(in) :: filename_cstr(1)
    integer(c_int), intent(out) :: norbitals
    type(c_ptr), intent(out) :: orbitals_ptr
    integer(c_int) :: rwfnread

    character(:), allocatable :: filename

    integer :: fhandle, ios
    character(255) :: iom
    character(6) :: g92rwf
    integer :: errlineno = 0

    integer :: idx, i, npy, naky, my
    real(real64) :: pzy, ey
    real(real64), pointer, dimension(:) :: pa, qa, ra
    type(orbital_t), dimension(1024) :: orbitals
    type(orbital_t), pointer, dimension(:) :: orbitals_allocatable

    filename = from_cstring(filename_cstr)
    open(newunit=fhandle, file=filename, form="unformatted", status="old", iostat=ios, IOMSG=iom)
    if(ios /= 0) then
        print *, "ERROR: Unable to open file:", ios, iom
        rwfnread = 1
        return
    endif

    ! Check header
    read(fhandle, iostat=ios, iomsg=iom) g92rwf
    if(g92rwf /= "G92RWF") then
        print *, "ERROR: Bad header. Not a G92RWF file?"
        rwfnread = 3
        close(fhandle)
        return
    endif

    idx = 0 ! index of the current wavefunction
    do while(.true.)
        ! READ (23,IOSTAT = IOS) npy, naky, ey, my
        ! npy  -- primary quantum number
        ! naky -- kappa quantum number
        ! ey   -- orbital's DC energy (??)
        ! my   -- number of points
        read(fhandle, iostat=ios) npy, naky, ey, my
        if(ios /= 0) then
            !print *, "INFO: No more orbitals to read. Exiting the loop."
            !print *, npy, naky, ey, my, ios
            exit
        endif

        idx = idx + 1
        if(size(orbitals) < idx) then
            print '("ERROR: Hitting a hard-coded limit of ", i5, " orbitals.")', size(orbitals)
            rwfnread = 4
            close(fhandle)
            return
        endif

        orbitals(idx)%npy = npy
        orbitals(idx)%naky = naky
        orbitals(idx)%ey = ey
        orbitals(idx)%my = my

        allocate(pa(orbitals(idx)%my), qa(orbitals(idx)%my), ra(orbitals(idx)%my))

        ! READ (23) PZY,(PA(I),I = 1,MY),(QA(I),I =1 ,MY)
        ! READ (23) (RA(I),I = 1,MY)
        errlineno=__LINE__; read(fhandle, iostat=ios, iomsg=iom) pzy, (pa(i), i = 1, orbitals(idx)%my), (qa(i), i = 1, orbitals(idx)%my)
        if(ios /= 0) go to 999
        errlineno=__LINE__; read(fhandle, iostat=ios, iomsg=iom) (ra(i), i = 1, orbitals(idx)%my)
        if(ios /= 0) go to 999

        orbitals(idx)%ra = c_loc(ra)
        orbitals(idx)%pa = c_loc(pa)
        orbitals(idx)%qa = c_loc(qa)
    end do

    norbitals = idx

    allocate(orbitals_allocatable(norbitals))
    orbitals_allocatable(1:norbitals) = orbitals(1:norbitals)
    orbitals_ptr = c_loc(orbitals_allocatable)

    close(fhandle)
    return

    ! Error handling for IO errors (reachable via goto)
    999 continue
    print '(a)', "Terminating rwfnread() with IO error."
    print '(a,a,":",i0)', " at: ", __FILE__, errlineno
    print '(" while reading: ",a)', filename
    print *, ios, iom
    close(fhandle)
    rwfnread = 2
    return

end function rwfnread
