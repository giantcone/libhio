program hio_example

  use hio  ! this automatically uses the 'mpi' and 'iso_c_binding' Fortran modules
  use, intrinsic :: iso_fortran_env
  implicit none

  type(c_ptr) :: my_hio_context
  integer     :: my_rank, nranks
  
  integer(kind(MPI_COMM_WORLD)) :: mympi_comm
  integer(kind(MPI_ERROR))      :: retval_mpi
  integer(kind(HIO_SUCCESS))    :: retval_hio

  integer        :: argstat, rstat
  integer(int64) :: ndbls, gibibytes
  character(len=32) :: arg_str
  real(c_float),dimension(:),allocatable,target :: crap

  !.......................................................................................
  ! B E G I N   P R O G R A M
  !.......................................................................................

  ! Now grab the number of GibiBytes to write out 
  call get_command_argument(1, arg_str, status=argstat)
  read(arg_str,*,iostat=rstat) gibibytes
  if (0 /= arg_str .or. 0 /= rstat) call usage()
  
  call initialize_hio(stat_hio)

  call data_init()

  call checkpoint(32, 8, 1)

  call finalize_hio()

  call mpi_comm_free(mympi_comm, ierr=retval_mpi)

  call mpi_finalize(ierr=retval_mpi)
  
  stop
  !.....................................................................................................
  ! C O N T A I N E D    P R O C E D U R E S
  !.....................................................................................................
contains

  subroutine initialize_hio()

    call mpi_comm_rank(MPI_COMM_WORLD, my_rank, retval_mpi)
    call mpi_comm_size(MPI_COMM_WORLD, nranks, retval_mpi)

    retval_hio = hio_init_mpi(my_hio_context, mympi_comm, "hio_example.input.dec", "#HIO", &
    "hio_example_context")

    if (HIO_SUCCESS /= retval_hio) then
       ! NOTE: The C hio example code makes use of the 'hio_err_print_all' function.  That C
       ! function makes use of 'variadic arguments' which cannot be easily handled by the Fortran-C
       ! interface capabilities.  I'll need to likely code up some Fortran/C handling functions to
       ! somehow emulate variadic argument capability from Fortran to C.  For now, I'll just use
       ! native Fortran mechanisms.
       write(error_unit,'(A,I5)') 'Error from hio_init_mpi. Rank # ', my_rank
       call mpi_abort(comm_handle, retval_mpi)
    end if

  end subroutine initialize_hio
!......................................................................................................
  subroutine finalize_hio(retval_hio)
    integer(kind(HIO_SUCCESS)),intent(out) :: retval_hio

    retval_hio = hio_fini(my_hio_context)

  end subroutine finalize_hio
!......................................................................................................
  subroutine data_init()

  integer        :: memstat
  integer(int64) :: indx

  ! calculate the number of 32-bit reals for the amount of GibiBytes requested to dump per MPI process
  ndbls = gibibytes*2**28

  allocate(crap(ndbls), stat=memstat)
  if (0/=memstat) then
     write(*,'(A,I5)') 'Problem allocating data on MPI rank ', my_rank
     call mpi_abort(MPI_COMM_WORLD, retval_mpi)
  end if

  crap = real(my_rank, real32)

  end subroutine data_init
  
!......................................................................................................
  subroutine checkpoint(wcount, stride, tstep)

    ! Unlike the C analogue to this example code, in Fortran the hio module does not provide an
    ! interface to the functions that make use of the variadic arguments. So the functions that
    ! print out error messages directly from libhio will not be called here. Later some method to
    ! call the variadic functions from Fortran (or perhaps write limited handle C functions that are
    ! more ammenable to being called from Fortran) will be written.
    
    integer(c_size_t),intent(in)   :: wcount, stride
    integer(c_int),intent(in)      :: tstep

    type(c_ptr)                    :: element_ptr, set_ptr, data_ptr
    type(c_ptr),dimension(3)       :: req_ptr
    integer(c_size_t),dimension(3) :: bytes_xfrd
    integer(type(HIO_SCP_NOT_NOW)) :: hint_enum
    integer(type(HIO_FLAG_CREAT))  :: flag_val
    integer(c_int)                 :: wrt_offset
    integer(c_size_t)              :: bytes_wrtn, sz_c_size, sz_c_float
    integer(c_intptr_t)            :: addr_offset
    real(c_float)                  :: whatev


    ! type sizes 
    ! NOTE: a 'sizeof' equivalent in Fortran would be bit_size(var)/8
    whatev     = 3.14159
    wrt_offset = 0
    sz_c_int   = bit_size(wrt_offset)/8
    sz_c_size  = bit_size(wcount)/8
    sz_c_float = bit_size(whatev)/8


    hint_enum = hio_dataset_should_checkpoint(my_hio_context, "restart")

    if (HIO_SCP_NOT_NOW == hint_enum) return

    !/* Open the dataset associated with this timestep. The timestep can be used
    ! * later to ensure ordering when automatically rolling back to a prior    
    ! * timestep on data failure. Specify that files in this dataset have unique
    ! * offset spaces. */
    
    ret_enum = hio_dataset_alloc(my_hio_context, my_hio_set, "restart", tstep, flag_val, HIO_SET_ELEMENT_UNIQUE)

    ! Set the stage mode on the dataset

    retval_hio = hio_config_set_value(set_ptr, "datawarp_stage_mode", "lazy")

    retval_hio = hio_config_set_value(set_ptr, "stripe_width", "4096")

    retval_hio = hio_dataset_open(set_ptr)

    ! FILE mode bitmask
    flag_val = ior(HIO_FLAG_CREAT, HIO_FLAG_WRITE)

    retval_hio = hio_element_open(set_ptr, element_ptr, "data", flag_val)
    
    bytes_wrtn = hio_element_write(element_ptr, wrt_offset, 0, wcount, 1, sz_c_size)

    wrt_offset = wrt_offset + bytes_wrtn

    ! for use with the libHIO C function calls, I'll need to set the C pointer to the Fortran array that
    ! is used for holding the data. In the C example, there is pointer math that I believe is being done.
    ! I think I should be able to perform the same type of pointer math here as well.
    data_ptr = c_loc(crap)

    bytes_wrtn = hio_element_write_strided(element_ptr, wrt_offset, 0, data_ptr, wcount, sz_c_float, stride)

    wrt_offset = wrt_offset + bytes_wrtn

    addr_offset = int(sz_z_float, c_intptr_t)

    bytes_wrtn = hio_element_write_strided(element_ptr, wrt_offset, 0, data_ptr+addr_offset, wcount, &
                    sz_c_int, stride)

    ! Force flush all data locally (the following close operation should do this too.)
    retval_hio = hio_element_flush(element_ptr, HIO_FLUSH_MODE_LOCAL)

    !/* close the file and force all data to be written to the backing store */!
    retval_hio = hio_element_close(c_loc(element_ptr))

    ! collective call finalizing the dataset
    retval_hio = hio_dataset_close(set_ptr)

    !/* Dataset performance variables coule be retrieved between close and free */
    retval_hio = hio_dataset_free(c_loc(set_ptr))
    
    wrt_offset = my_rank*sz_c_size

    retval_hio = hio_element_write_nb (element_ptr, req_ptr, wrt_offset, 0, wcount, 1, sz_c_size);

    retval_hio = hio_element_write_strided_nb(element_ptr, req_ptr, wrt_offset, 0, wcount, 1, bit_size(sz_c_size)/8)
