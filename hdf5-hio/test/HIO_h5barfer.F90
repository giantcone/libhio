! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the files COPYING and Copyright.html.  COPYING can be found at the root   *
!   of the source code distribution tree; Copyright.html can be found at the  *
!   root level of an installed copy of the electronic HDF5 document set and   *
!   is linked from the top-level documents page.  It can also be found at     *
!   http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
!   access to either file, you may request a copy from help@hdfgroup.org.     *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! G. Cone (HPC-3 LANL) 21 April 2014
! Fortran parallel example.  Copied from Tutorial's example program of
! dataset.f90. Also adding in OpenMP threading capabilities to investigate
! hybrid mode parallel I/O.

!  Note (4 Aug 2015): I'm attempting to modify this code to perform a time-series
!  data writing behavior similar to what would be found in a typical scientific
!  application.  The HDF Group suggested that if I wish to have one 'temporally
!  contiguous' DataSET, that I will need to make use of the 'hyperslab offsets'
!  to append data to the file.  So if I use this approach, I'll need a method 
!  to record those offsets and correlate them with a simulation time-stamp.
!  This could likely be easily done as another dataset that holds the time value
!  and the base offset value of the main dataset.
!
! 2016-09-09: adding in the capability to use the LANL developed 'libhio'
!             library that leverages the DataWarp capability.
!...............................................................................
module hio_compat
! This module will deal with interfacing with the libHIO HDF5 functionality
! provided in C.  I think the best way to deal with this Fortran-C interoperability
! is to present as little as possible to this fortran code, and then have a 
! "go-between" handle code in C that has an explicit Fortran procedure interface
! defined in this module, and returns what is needed by the "vanilla" HDF5 
! program.

use, intrinsic :: iso_c_binding
use mpi
use hdf5
implicit none

private

! define a C compatible enumeration that can be used to specify the HIO
! methods to use as listed in the H5FD_hio_io_t enumeration. 
enum, bind(c)
  enumerator :: H5FD_HIO_BLOCKING = 0          !zero is the default
  enumerator :: H5FD_HIO_NONBLOCKING,    &
                H5FD_HIO_CONTIGUOUS,     &
                H5FD_HIO_STRIDED,        &
                H5FD_HIO_DATASET_SHARED, &
                H5FD_HIO_DATASET_UNIQUE
end enum

! integer,parameter :: HIO_FILE_NAME_SIZE     = 1024
! integer,parameter :: HIO_ELEM_NAME_SIZE     =  256
! integer,parameter :: HIO_CONFIG_FILE_SIZE   =  128
! integer,parameter :: HIO_CONFIG_PREFIX_SIZE =  256


type,bind(c) :: hio_request_t
  type(c_ptr) :: hio_request
end type hio_request_t


type,bind(c) :: hio_settings_t
  integer(c_int)      :: read_blocking, write_blocking, read_io_mode, write_io_mode, dataset_mode
  integer(c_size_t)   :: stride_size
  type(hio_request_t) :: request
  character(kind=c_char),dimension(HIO_FILE_NAME_SIZE)     :: name
  character(kind=c_char),dimension(HIO_ELEM_NAME_SIZE)     :: element_name
  character(kind=c_char),dimension(HIO_CONFIG_FILE_SIZE)   :: config_file
  character(kind=c_char),dimension(HIO_CONFIG_PREFIX_SIZE) :: config_prefix
  integer            :: comm
  integer(c_int64_t) :: setid
  integer(c_int)     :: flags
end type hio_settings_t

interface 
  function h5fd_hio_set_write_io(facct, ) bind(c, name='H5FD_hio_set_write_io')
    integer(hid_t) :: h5fd_hio_set_write_io

end interface

end module hio_compat
     

program paratest

! INTRINSIC FORTRAN MODULES
use, intrinsic :: iso_fortran_env, only: real64, int64, int32
use, intrinsic :: iso_c_binding,   only: c_ptr, c_loc

!use omp_lib
! EXTERNAL LIBRARY MODULES
use hio_compat

implicit none

character(len=12), parameter :: default_fname = "paradump.h5"  ! Default file name
character(len=8), parameter  :: dsetname      = "DBLArray"     ! Dataset name

character(len=64) :: pll_path, filepath, arg_str  ! Path, File path, arg string

integer        :: nargs, nwrt, iwrt  ! # cmd args, # of write instances, index

integer        :: fnamelen        ! File name length
integer(hid_t) :: file_h5id       ! File identifier
integer(hid_t) :: dset_h5id       ! Dataset identifier
integer(hid_t) :: dspc_m_h5id     ! Dataspace identifier in memory
integer(hid_t) :: dspc_f_h5id     ! Dataspace identifier in file
integer(hid_t) :: attr_dspc_h5id  ! Attribute Dataspace identifier in file
integer(hid_t) :: plist_h5id      ! Property list identifier
!integer(hid_t) :: attr_h5id       ! Attribute identifier

integer(int64) :: ndbls, mibibytes
! NOTE: 1 GiB = 2^27 doubles = (2^9)^3  =   512^3 = 134,217,728
!   or  1 GiB = 2^28 singles = (2^14)^2 = 16384^2 = 268,435,456
!       1 MiB = 2^17 doubles = 131,072

real(real64), target, allocatable :: mydata(:)           ! Data to write

!real(real64), target, allocatable :: wrt_times(:) ! An array holding the time
                                                  ! values for each iteration.
                                                  ! Written as Attribute data.   

!integer :: drnk = rank(mydata)    ! Number of Dataspace dimensions
integer :: drnk = 1    ! Number of Dataspace dimensions
!integer :: arnk = rank(wrt_times) ! Number of Dataspace dims for attribute data
integer :: arnk = 1 ! Number of Dataspace dims for attribute data

integer(hssize_t), dimension(1) :: ddims, offset !, & !dspace shape & MPI offsets
!                                  adims            !attr dspace shape

! C-pointer variable to data address in Virtual memory
type(c_ptr) :: addr_mydata

integer :: err_gen, err_h5 ! Error flags

!
! MPI definitions and calls.
!
integer :: err_mpi       ! MPI err_h5 flag
integer :: comm, info, mpi_def_thrd_lvl
integer :: mpi_size, mpi_rank

real(real64) :: t_start, t_stop !, tick_mpi

! set the requested MPI thread-safety level using preprocessor macros
! NOTE: from mpif-common.h the following values for the mpi module 
!       variables are defined:
!
!      parameter (MPI_THREAD_SINGLE=0)
!      parameter (MPI_THREAD_FUNNELED=1)
!      parameter (MPI_THREAD_SERIALIZED=2)
!      parameter (MPI_THREAD_MULTIPLE=3)
#ifndef THREQLVL
#define REQ_THRD_SAFETY MPI_THREAD_SINGLE
#else
#define REQ_THRD_SAFETY THREQLVL
#endif
!...............................................................................
! BEGIN PROGRAM
!...............................................................................
!print *, REQ_THRD_SAFETY
comm = MPI_COMM_WORLD
info  = MPI_INFO_NULL

#ifdef THRDD
call mpi_init_thread(REQ_THRD_SAFETY, mpi_def_thrd_lvl, err_mpi)
#else
call mpi_init(err_mpi)
#endif

call mpi_comm_size(comm, mpi_size, err_mpi)
call mpi_comm_rank(comm, mpi_rank, err_mpi)

! Get command line arguments and ensure this program is being invoked properly.
call invoke_init()
!call mpi_barrier(wcomm, err_mpi)

! calculate the number of doubles to use
ndbls = mibibytes*2_int64**17

! Checking if the right number for ndbls is being calculated:
#ifdef DBG
if (0==mpi_rank) then
  write(*,'(A,I20)') 'Number of doubles to write per MPI rank is ', ndbls
end if
#endif

! Initialize Fortran HDF5 interface
call h5open_f(err_h5)

! Setup file access property list with parallel I/O access.
call h5pcreate_f(H5P_FILE_ACCESS_F, plist_h5id, err_h5)
call h5pset_fapl_mpio_f(plist_h5id, comm, info, err_h5)

!
! Create the file collectively.
!
call h5fcreate_f(filepath, H5F_ACC_TRUNC_F, file_h5id, err_h5, &
                 access_prp = plist_h5id)
call h5pclose_f(plist_h5id, err_h5)

!
! Initialize data buffer with trivial data.
!
allocate(mydata(ndbls), stat=err_gen)

! allocate the space for the attribute data holding timing values for writing
!allocate(wrt_times(nwrt), stat=err_gen)

#ifdef DBG
!Check that this array is allocated properly
if (0==mpi_rank) then
  if (0/=err_gen) then
     write(*,'(A, I8)') 'Allocate return code is ', err_gen
     call mpi_abort(MPI_COMM_WORLD, 255, err_mpi)
  else
     print *, 'Array mydata allocated properly.'
  end if
end if
#endif

! the dataspace description is simply the result of a Fortan 'shape'
ddims = shape(mydata, hssize_t)

#ifdef DEBUG 
if (0==mpi_rank) then
  write(*,'(A,I20,A)') 'mydata contains ', ddims(1), ' elements.'
end if
#endif

! attribute dataspace description also needs to be queried
!adims = shape(wrt_times, hsize_t)

! Fill the allocated memory with simple crap
mydata = real(mpi_rank, real64)

! Create the data space (data layout) for the dataset in memory:
call h5screate_simple_f(drnk, ddims, dspc_m_h5id, err_h5)

! Create the data space (data layout) for the  dataset in the HDF5 file:
call h5screate_simple_f(drnk, int(mpi_size*ddims*nwrt,hssize_t), dspc_f_h5id, err_h5)
!call h5screate_simple_f(drnk, ddims, dspc_f_h5id, err_h5)

!
! Create the DataSET with default properties. I believe this uses the 
! File-based DataSPaCe
!
call h5dcreate_f(file_h5id, dsetname, H5T_NATIVE_DOUBLE, dspc_f_h5id, &
                 dset_h5id, err_h5)


! Create property list for collective dataset write
!
call h5pcreate_f(H5P_DATASET_XFER_F, plist_h5id, err_h5)

! The following defines the "data transfer properties" using MPI-I/O
! "DXPL" = "Data X-fer Property List" is my guess...
!
#ifndef INDY
call h5pset_dxpl_mpio_f(plist_h5id, H5FD_MPIO_COLLECTIVE_F, err_h5)
!
! For independent write use
#else
call h5pset_dxpl_mpio_f(plist_h5id, H5FD_MPIO_INDEPENDENT_F, err_h5)
#endif
!

! Write the dataset collectively. Need to have a C-pointer variable
! that indicates the address of the data to write. I'll also write an 
! 'attribute' that indicates the iteration of data written to file.
!
addr_mydata = c_loc(mydata)

! W R I T E   L O O P
do iwrt = 0, nwrt-1
  ! Ensure the DataSPaCe is closed for this iteration. I'm guessing that the
  ! call to 'h5dget_space_f' is what "reopens" the dataspace?
  call h5sclose_f(dspc_f_h5id, err_h5)

  ! the data offset in the H5 file will be the MPI rank times the shape
  ! time the iteration index
  offset =  (iwrt*mpi_size + mpi_rank)*ddims

#ifdef DBG
  if (mpi_size-1 == mpi_rank) then
    write(*,'(A,I8,A,I20)') 'Offset value for rank ', mpi_rank, ' is ', offset(1)
  end if
#endif

  ! This "returns an ID for a copy of the DataSPaCe for a DataSET" according to
  ! the HDF5 Ref Manual entry for this subroutine. 
  ! This will "recycle" the previous 'dspc_f_h5id' DataSPaCe handle for the 
  ! file view of the dataspace 
  call h5dget_space_f(dset_h5id, dspc_f_h5id, err_h5)
  ! and set the hyperslab DataSPaCe view on the file to write:
  !
  call h5sselect_hyperslab_f(dspc_f_h5id, H5S_SELECT_SET_F, offset, ddims, err_h5)
  
  call mpi_barrier(comm, err_mpi)
  t_start = mpi_wtime()

  ! Write the actual Dataset in parallel to the HDF5 file
  call h5dwrite_f(dset_h5id, H5T_NATIVE_DOUBLE, addr_mydata, err_h5, &
                  file_space_id=dspc_f_h5id, mem_space_id=dspc_m_h5id, &
                  xfer_prp = plist_h5id)

  call mpi_barrier(comm, err_mpi)
  t_stop = mpi_wtime()

  ! Store the time to write this instance to the wrt_times array
!  wrt_times(iwrt) = t_stop - t_start

#ifdef WRTOUT
  if (0 == mpi_rank) then
    write(*,'(A,I4,A,I9,A,F8.3,A)') 'Wallclock time to write instance', iwrt+1, &
                                    ' of size ', mpi_size*mibibytes, & 
                                    ' MiB is: ', (t_stop-t_start)/60, ' min.'
  end if
#endif

  ! Now update the values of the data in memory:
  mydata = mydata + sqrt(mydata)*(-1)**iwrt

end do 

! Now Create the ATTRibute dataspace
!call h5screate_simple_f(arnk, adims, attr_dspc_h5id, err_h5)

! Create an ATTRibute that will be associated with the above DataSET
!call h5acreate_f(file_h5id, 'Write times (seconds)', H5T_NATIVE_INTEGER, &
!                 attr_dspc_h5id, attr_h5id, err_h5)

! Write the attribute to the HDF5 file. I believe this is done in parallel
! This is only needed to be done once to indicate how many iterations of 
! junk was written to the file.
!call h5awrite_f(attr_h5id, H5T_NATIVE_INTEGER, c_loc(wrt_times), err_h5)


!
! Deallocate data buffer.
!
deallocate(mydata)
!deallocate(wrt_times)
!
! Close HDF5 resources.
!

call h5sclose_f(dspc_f_h5id, err_h5)    ! file Dataspace descriptor
call h5sclose_f(dspc_m_h5id, err_h5)    ! memory Dataspace descriptor
!call h5sclose_f(attr_dspc_h5id, err_h5) ! Attribute dataspace descriptor
!call h5aclose_f(attr_h5id, err_h5)      ! Attribute handle
call h5dclose_f(dset_h5id, err_h5)      ! Dataset handle
call h5pclose_f(plist_h5id, err_h5)     ! File property list handle
call h5fclose_f(file_h5id, err_h5)      ! HDF5 file

!
! Close FORTRAN interface
!
call h5close_f(err_h5)       ! end HDF5 interface

call mpi_barrier(comm, err_mpi)
call mpi_finalize(err_mpi) ! end MPI communications
stop

!..............................................................................8
! C O N T A I N E D   P R O C E D U R E S                                      0
!...............................................................................
contains

! invoke_init: grabs command line arguments and ensures saneness of invocation.
subroutine invoke_init()

! Use the environmental variable HDF5_PARAPREFIX to specify the 
! directory of the parallel file-system to use. Otherwise, the program
! will default to /scratch/$USER

  call get_environment_variable('HDF5_PARAPREFIX', pll_path, status=err_gen)

  if (0 /= err_gen) then
    call get_environment_variable('USER', arg_str, status=err_h5)
    pll_path = '/scratch3/' // trim(arg_str) // '/'
    if (0 == mpi_rank) then 
      print *, 'Using ' // trim(pll_path) // ' as default write path.'
    end if
  end if

  nargs = command_argument_count()
  ! ensure program is called with correct arguments
  if ( 3 /= nargs .AND. 0==mpi_rank) then
    call usage()
    if (0 == mpi_rank) print *, 'Using default values for file size and path.'
  !  call mpi_abort(wcomm, 255, err_mpi)
  end if

  ! Get amount of MibiBytes of data to write from the command line:
  call get_command_argument(1, arg_str, status=err_gen)
  if (0 /= err_gen) then
    if (0 == mpi_rank) print *, 'Defaulting to 1GiB/MPI rank'
    mibibytes = 1024
  else
    read(arg_str,'(I8)') mibibytes
  end if

  ! Get the number of iterations of data to write to file (size multiplier)
  call get_command_argument(2, arg_str, status=err_gen)
  if (0 /= err_gen) then
    if (0 == mpi_rank) print *, 'Defaulting to ONE write instance.'
    nwrt = 1
  else
    read(arg_str,'(I8)') nwrt
  end if
 
  ! Get the file name to write data to.
  call get_command_argument(3, arg_str, status=err_gen)
  if (0 == err_gen) then
    filepath = trim(pll_path) // '/' // trim(arg_str)
  else
    filepath = trim(pll_path) // '/' // trim(default_fname)
  end if

  if (0 == mpi_rank) print *, "Using filename = ", filepath
end subroutine invoke_init

  ! usage: informs user on STDOUT the proper method to invoke this program.
subroutine usage()

  character(len=32) :: pname

  ! get the name of the executable:
  call get_command_argument(0, pname, status=err_gen)

  write(*,'(3A)') "Usage: ", trim(pname), "  MibiBytes(int) Iterations(int) Filename(string)"
end subroutine usage

end program
