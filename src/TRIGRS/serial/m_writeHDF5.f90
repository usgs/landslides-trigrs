module m_writeHDF5File

    use hdf5

    implicit none

    private

    interface writeHDF5File
        procedure :: writeHDF5File_r3D, writeHDF5File_i2D
    end interface

    public :: writeHDF5File

contains


subroutine writeHDF5File_r3D(nSteps, nRows, nCols, imx1, imax, values1D, values2D, noData, filename)
    !! Takes input from Rex's code, remaps to a 3D buffer, and writes to hdf5.

    integer, intent(in) :: nSteps, nRows, nCols, imx1, imax
    real, intent(in) :: values1D(:)
    real, intent(in) :: values2D(:, :)
    real, intent(in) :: noData
    character(len=*), intent(in) :: filename

    real :: test
    real :: z3(imax)
    real :: buffer(nCols, nRows, nSteps)

    integer :: i, j, k, counter


    do k = 1, nSteps
        z3 = 0.
        do i = 1, imx1
            z3(i) = values1D(i + (k - 1) * imax)
        end do

        counter = 0
        do i = 1, nRows
            do j = 1, nCols
                test = abs(values2D(j, i) - nodata)
                if (test <= 0.1) then
                    buffer(j, i, k) = nodata
                else
                    counter = counter + 1
                    buffer(j, i, k) = z3(counter)
                endif

            enddo
        enddo
    enddo

    call hdf5_r3D(buffer, filename)

end subroutine


subroutine writeHDF5File_i2D(nRows, nCols, values1D, values2D, noData, filename)
    !! Takes input from Rex's code, remaps to a 3D buffer, and writes to hdf5.

    integer, intent(in) :: nRows, nCols
    integer, intent(in) :: values1D(:)
    real, intent(in) :: values2D(:, :)
    real, intent(in) :: noData
    character(len=*), intent(in) :: filename

    real :: test
    integer :: buffer(nCols, nRows)

    integer :: i, j, k, counter

    counter = 0
    do i = 1, nRows
        do j = 1, nCols
            test = abs(values2D(j, i) - nodata)
            if (test <= 0.1) then
                buffer(j, i) = nodata
            else
                counter = counter + 1
                buffer(j, i) = values1D(counter)
            endif

        enddo
    enddo

    call hdf5_i2D(buffer, filename)

end subroutine





subroutine hdf5_r3D(values, filename)

    real, intent(in) :: values(:, :, :)
    character(len=*), intent(in) :: filename

    integer, parameter :: rank = 3

    character(LEN=4), parameter :: dsetname = "data"     ! Dataset name

    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier


    integer(HSIZE_T) :: dims(rank) ! Dataset dimensions
    integer ::   error ! Error flag

    integer :: i
    
    dims = shape(values)

    ! Initialize FORTRAN interface.
    CALL h5open_f(error)

    ! Create a new file using default properties.
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

    ! Create the dataspace.
    CALL h5screate_simple_f(rank, dims, dspace_id, error)

    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, "values", H5T_NATIVE_REAL, dspace_id, dset_id, error)

    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, values, dims, error)

    ! End access to the dataset and release resources used by it.
    CALL h5dclose_f(dset_id, error)

    ! Terminate access to the data space.
    CALL h5sclose_f(dspace_id, error)
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

end subroutine

subroutine hdf5_i2D(values, filename)

    integer, intent(in) :: values(:, :)
    character(len=*), intent(in) :: filename

    integer, parameter :: rank = 2

    character(LEN=4), parameter :: dsetname = "data"     ! Dataset name

    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: dspace_id     ! Dataspace identifier


    integer(HSIZE_T) :: dims(rank) ! Dataset dimensions
    integer ::   error ! Error flag

    integer :: i
    
    dims = shape(values)

    ! Initialize FORTRAN interface.
    CALL h5open_f(error)

    ! Create a new file using default properties.
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

    ! Create the dataspace.
    CALL h5screate_simple_f(rank, dims, dspace_id, error)

    ! Create the dataset with default properties.
    CALL h5dcreate_f(file_id, "values", H5T_NATIVE_INTEGER, dspace_id, dset_id, error)

    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, values, dims, error)

    ! End access to the dataset and release resources used by it.
    CALL h5dclose_f(dset_id, error)

    ! Terminate access to the data space.
    CALL h5sclose_f(dspace_id, error)
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

end subroutine







end module

