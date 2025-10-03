program main
    use omp_lib  ! only for wall time
    implicit none

    interface

    subroutine write_g96(g96_unit, time, step, number_of_atoms, Ax, By, Cz, &
        pos, vel)
        implicit none

        integer, intent(in) :: g96_unit
        real(kind=8), intent(in) :: time  ! fs
        integer, intent(in) :: step
        integer, intent(in) :: number_of_atoms
        real(kind=8), intent(in) :: Ax  ! Angstrom
        real(kind=8), intent(in) :: By  ! Angstrom
        real(kind=8), intent(in) :: Cz  ! Angstrom
        real(kind=8), intent(in) :: pos(3, number_of_atoms)  ! Angstrom
        real(kind=8), intent(in), optional :: vel(3, number_of_atoms)  ! atomic

    end subroutine write_g96

    end interface

    !> command line arguments
    integer :: argc

    !> files
    character(len=80) :: file_name
    character(len=80) :: cp2k_cell
    character(len=80) :: cp2k_pos
    character(len=80) :: cp2k_vel
    character(len=80) :: g96

    !> file units
    integer, parameter :: cp2k_cell_unit = 100
    integer, parameter :: cp2k_pos_unit  = 200
    integer, parameter :: cp2k_vel_unit  = 300
    integer, parameter :: g96_unit       = 400

    integer :: stat

    !> file contents
    character(len=80) :: line

    !> cp2k md contents
    integer      :: step
    real(kind=8) :: time    ! fs
    real(kind=8) :: Ax      ! Angstrom
    real(kind=8) :: Ay      ! Angstrom
    real(kind=8) :: Az      ! Angstrom
    real(kind=8) :: Bx      ! Angstrom
    real(kind=8) :: By      ! Angstrom
    real(kind=8) :: Bz      ! Angstrom
    real(kind=8) :: Cx      ! Angstrom
    real(kind=8) :: Cy      ! Angstrom
    real(kind=8) :: Cz      ! Angstrom
    real(kind=8) :: volume  ! Angstrom^3

    integer :: number_of_atoms(2)
    character(len=80) :: title(2)
    character(len=10), allocatable :: symbol(:)
    real(kind=8), allocatable :: pos(:, :)  ! Angstrom
    real(kind=8), allocatable :: vel(:, :)  ! atomic

    logical :: has_vel

    !> summary
    integer :: number_of_frames
    real(kind=8) :: time_start
    real(kind=8) :: time_end
    real(kind=8) :: time_total

    !> loop index
    integer :: i, j, k

    !> parse command line arguments
    argc = iargc()

    if (argc .lt. 2) then
        print "(A)", "Usage:"
        print "(A)", "  $ cp2kmd2g96.x CELL POS [VEL] [G96]"
        stop
    end if

    cp2k_cell = ""
    cp2k_pos  = ""
    cp2k_vel  = ""
    g96       = ""

    has_vel   = .false.

    do i = 1, argc
        call getarg(i, file_name)
        file_name = trim(file_name)

        if (index(file_name, ".cell") .gt. 0) then
            cp2k_cell = file_name
            continue
        else if ((index(file_name, "pos") .gt. 0) &
            .and. (index(file_name, ".xyz") .gt. 0)) then
            cp2k_pos = file_name
            continue
        else if ((index(file_name, "vel") .gt. 0) &
            .and. (index(file_name, ".xyz") .gt. 0)) then
            cp2k_vel = file_name
            has_vel = .true.
            continue
        else if (index(file_name, ".g96") .gt. 0) then
            g96 = file_name
            continue
        end if
    end do

    if (len_trim(g96) .eq. 0) then
        g96 = "traj.g96"
    end if

    !> parse cell, pos, & vel
    open(unit=cp2k_cell_unit, file=trim(cp2k_cell), status="old")
    open(unit=cp2k_pos_unit, file=trim(cp2k_pos), status="old")

    if (has_vel) then
        open(unit=cp2k_vel_unit, file=trim(cp2k_vel), status="old")
    end if

    open(unit=g96_unit, file=trim(g96), status="replace")

    !> skip header line
    read(unit=cp2k_cell_unit, fmt="(A)", iostat=stat) line

    number_of_frames = 0
    time_start = omp_get_wtime()

    do
        !> cell
        read(unit=cp2k_cell_unit, fmt=*, iostat=stat) step, time, Ax, Ay, Az, &
            Bx, By, Bz, Cx, Cy, Cz, volume

        if (stat .ne. 0) then
            exit
        end if

        !> pos
        read(unit=cp2k_pos_unit, fmt=*, iostat=stat) number_of_atoms(1)

        if (stat .ne. 0) then
            exit
        end if

        if (has_vel) then
            read(unit=cp2k_vel_unit, fmt=*, iostat=stat) number_of_atoms(2)

            if (stat .ne. 0) then
                exit
            end if

            if (number_of_atoms(1) .ne. number_of_atoms(2)) then
                print *, "Error: Number of atoms in ", trim(cp2k_pos), &
                    " and ", trim(cp2k_vel), "do not match!"
                stop
            end if
        end if

        if (.not. allocated(symbol)) then
            allocate(symbol(number_of_atoms(1)))
        end if

        if (.not. allocated(pos)) then
            allocate(pos(3, number_of_atoms(1)))
        end if

        if (has_vel) then
            if (.not. allocated(vel)) then
                allocate(vel(3, number_of_atoms(2)))
            end if
        end if

        read(unit=cp2k_pos_unit, fmt=*) title(1)

        if (has_vel) then
            read(unit=cp2k_vel_unit, fmt=*) title(1)
        end if

        do i = 1, number_of_atoms(1)
            read(unit=cp2k_pos_unit, fmt=*) symbol(i), pos(:, i)

            if (has_vel) then
                read(unit=cp2k_vel_unit, fmt=*) symbol(i), vel(:, i)
            end if
        end do

        !> write g96
        if (has_vel) then
            call write_g96(g96_unit, time, step, number_of_atoms(1), &
                Ax, By, Cz, pos, vel)
        else
            call write_g96(g96_unit, time, step, number_of_atoms(1), &
                Ax, By, Cz, pos)
        end if
        
        number_of_frames = number_of_frames + 1
    end do

    time_end = omp_get_wtime()
    time_total = time_end - time_start

    !> deallocate
    deallocate(symbol)
    deallocate(pos)

    if (has_vel) then
        deallocate(vel)
    end if

    !> close
    close(cp2k_cell_unit)
    close(cp2k_pos_unit)

    if (has_vel) then
        close(cp2k_vel_unit)
    end if

    close(g96_unit)

    !> summary
    print *, "Summary"
    print *, "Number of frames: ", number_of_frames
    print *, "Time (s): ", time_total
    print *, "Performance (frame/s): ", dble(number_of_frames) / time_total

end program main

subroutine write_g96(g96_unit, time, step, number_of_atoms, Ax, By, Cz, &
    pos, vel)
    implicit none

    integer, intent(in) :: g96_unit
    real(kind=8), intent(in) :: time  ! fs
    integer, intent(in) :: step
    integer, intent(in) :: number_of_atoms
    real(kind=8), intent(in) :: Ax  ! Angstrom
    real(kind=8), intent(in) :: By  ! Angstrom
    real(kind=8), intent(in) :: Cz  ! Angstrom
    real(kind=8), intent(in) :: pos(3, number_of_atoms)  ! Angstrom
    real(kind=8), intent(in), optional :: vel(3, number_of_atoms)  ! atomic

    integer :: i

    !> write title
    write(g96_unit, fmt="(A)") "TITLE"
    write(g96_unit, fmt="(A)") ""
    write(g96_unit, fmt="(A)") "END"

    !> write time
    write(g96_unit, fmt="(A)") "TIMESTEP"
    !> 1.0 fs = 1.0e-03 ps
    write(g96_unit, fmt="(I15,F15.6)") step, time * 1.0d-03
    write(g96_unit, fmt="(A)") "END"

    !> write position
    write(g96_unit, fmt="(A)") "POSITIONRED"

    do i = 1, number_of_atoms
        !> 1.0 Angstrom = 0.1 nm
        write(g96_unit, fmt="(3F15.9)") pos(:, i) * 0.1d0
    end do

    write(g96_unit, fmt="(A)") "END"

    !> write velocity
    if (present(vel)) then
        write(g96_unit, fmt="(A)") "VELOCITYRED"

        do i = 1, number_of_atoms
            !> 1.0 atomic = 2187.7 nm/ps
            write(g96_unit, fmt="(3F15.9)") vel(:, i) * 2187.69d0
        end do

        write(g96_unit, fmt="(A)") "END"
    end if

    !> write box
    write(g96_unit, fmt="(A)") "BOX"
    !> 1.0 Angstrom -> 0.1 nm
    write(g96_unit, fmt="(3F15.9)") Ax * 0.1d0, By * 0.1d0, Cz * 0.1d0
    write(g96_unit, fmt="(A)") "END"

end subroutine write_g96
