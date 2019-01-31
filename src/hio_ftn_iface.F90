!................................................................................
! Copyright (c) 2018 Triad National Security, LLC.
! $COPYRIGHT$
!  
! Additional copyrights may follow
!
! L I B H I O    F O R T R A N    I N T E R F A C E    M O D U L E
!
! This module will deal with interfacing with the libHIO HDF5 functionality
! provided in C.  I think the best way to deal with this Fortran-C interoper-
! ability is to present as little as possible to this fortran code, and then have
! a go-between handle code in C that has an explicit Fortran procedure interface
! defined in this module.
!
! For those C functions that take a character type argument, there will need to
! be intermediate handle functions that properly modify the Fortran character
! string into one that C and handle with a End-Of-Line null character.
!................................................................................
module hio

use, intrinsic :: iso_c_binding
#ifdef HIO_USE_MPI
use mpi
#endif
implicit none

!private

!#define HIO_FILE_NAME_SIZE 1024
!#define HIO_ELEM_NAME_SIZE 256
!#define HIO_CONFIG_FILE_SIZE 128
!#define HIO_CONFIG_PREFIX_SIZE 256

! A note about the "opaque data structures" used in 'hio.h':
!
!   typedef struct hio_STUFF *hio_STUFF_t;
!
! In Fortran, these type of structure pointer are treated as C pointers to an
! "opaque stuctured datatype" using the type(c_ptr)[,value] shown in the argument list
! of some of the interface definitions below.


! The HIO Enumerated types are defined here for convenience and to keep in line
! with those defined in the C header file
!
! 'hio_return_t' enumerated type
!
! HINT: To define a type in Fortran that makes use of this enumeration,
! do the following:
!
!     integer(kind(HIO_SUCCESS)) :: VARIABLE_OF_ENUM_TYPE

enum, bind(c)
  enumerator :: &
  !
  HIO_SUCCESS  =  0, &
  ! Generic HIO error
  HIO_ERROR = -1, &
  ! Permissions error or operation not permitted 
  HIO_ERR_PERM = -2, &
  ! Short read or write
  HIO_ERR_TRUNCATE = -3, &          
  ! Out of resources 
  HIO_ERR_OUT_OF_RESOURCE = -4, &
  ! Not found 
  HIO_ERR_NOT_FOUND       = -5, &
  ! Not available 
  HIO_ERR_NOT_AVAILABLE = -6, &
  ! Bad parameter 
  HIO_ERR_BAD_PARAM     = -7, &
  ! Dataset id already exists
  HIO_ERR_EXISTS        = -8, &
  ! Temporary IO error. Try the IO again later.
  HIO_ERR_IO_TEMPORARY  = -Z'00010001', &
  ! Permanent IO error. IO to the current data root is no longer available. 
  HIO_ERR_IO_PERMANENT  = -Z'00010002'
end enum

! 'hio_flags_t' enumerated type
enum, bind(c)
  enumerator :: &
  !/** Open the dataset read-only Can not be combined with
  !   * HIO_FLAG_WRITE at this time. */
  HIO_FLAG_READ     = 1,&
  !/** Open the dataset write-only. Can not be combined with
  !  * HIO_FLAG_READ at this time. */
  HIO_FLAG_WRITE    = 2,&
  !/** Create a new element. If the element exists and HIO_FLAG_TRUNC
  ! * is not specified an error is returned. */
  HIO_FLAG_CREAT    = 64,&
  !/** Remove all existing data associated with the dataset. */
  HIO_FLAG_TRUNC    = 512,&
  !/** Append to an existing dataset. Not supported at this time. */
  HIO_FLAG_APPEND   = 1024
end enum

! 'hio_flush_mode_t' enumerated type
enum, bind(c)
  enumerator :: &
  ! /** Locally flush data. This mode ensures that the user buffers can
  !  * be reused by the application. It does not ensure the data has
  !  * been written out to the backing store. This is the equivalent of
  !  * the fflush() call from streaming IO. */
  HIO_FLUSH_MODE_LOCAL = 0, &
  ! /** Ensure all data has been written out to the backing store. On
  !  * posix this is the equivalent of fsync(). */
  HIO_FLUSH_MODE_COMPLETE = 1
end enum

! 'hio_config_type_t' enumerated type
enum, bind(c)
  enumerator :: &
  !/** C99/C++ boolean. Valid values: true, false, 0, 1 */
  HIO_CONFIG_TYPE_BOOL,&
  !/** NULL-terminated string */
  HIO_CONFIG_TYPE_STRING,&
  !/** 32-bit signed integer */
  HIO_CONFIG_TYPE_INT32,&
  !/** 32-bit unsigned integer */
  HIO_CONFIG_TYPE_UINT32,&
  !/** 64-bit signed integer */
  HIO_CONFIG_TYPE_INT64,&
  !/** 64-bit unsigned integer */
  HIO_CONFIG_TYPE_UINT64,&
  !/** IEEE-754 floating point number */
  HIO_CONFIG_TYPE_FLOAT,&
  !/** IEEE-754 double-precision floating point number */
  HIO_CONFIG_TYPE_DOUBLE
end enum

! 'hio_dataset_mode_t' enumerated type
enum, bind(c)
  enumerator :: &
  ! /** Element(s) in the dataset have unique offset spaces. This mode
  !  * is equivalent to an N-N IO pattern where N is the number of
  !  * ranks in the communicator used to create the hio context. */
  HIO_SET_ELEMENT_UNIQUE, &
  ! /** Element(s) in the dataset have shared offset spaces. This mode
  !  * is equivalent to an N-1 IO pattern where N is the number of
  !  * ranks in the communicator used to create the hio context. */
  HIO_SET_ELEMENT_SHARED
end enum

! 'hio_recommendation_t' enumerated type
enum, bind(c)
  enumerator :: &
  !/** Do not attempt to checkpoint at this time */
  HIO_SCP_NOT_NOW, &
  !/** Checkpoint strongly recommended */
  HIO_SCP_MUST_CHECKPOINT
end enum

! 'hio_unlink_mode_t' enumerated type
enum, bind(c)
  enumerator :: &
  !/** Unlink dataset id only in the current active data root */
  HIO_UNLINK_MODE_CURRENT, &
  !/** Unlink first matching dataset id instance */
  HIO_UNLINK_MODE_FIRST,   &
  !/** Unlink all matching dataset id instances */
  HIO_UNLINK_MODE_ALL
end enum

!................................................................................
! LIBHIO PROCEDURE (C FUNCTION) INTERFACES
!................................................................................
interface
   
  integer(kind(HIO_SUCCESS)) &
       function hio_init_single (new_context, config_file, config_file_prefix, &     
       context_name) bind(c)
    import :: c_ptr, c_char  
    type(c_ptr),value                   :: new_context
    character(kind=c_char),dimension(*) :: config_file
    character(kind=c_char),dimension(*) :: config_file_prefix
    character(kind=c_char),dimension(*) :: context_name
  end function hio_init_single
  

#if !defined(MPI_VERSION)

   integer(kind(HIO_SUCCESS)) &
        function hio_init_mpi (new_context, comm, config_file, config_file_prefix, &
        name) bind(c)
    import :: c_ptr, c_char
    type(c_ptr)                         :: new_context
    type(c_ptr)                         :: comm
    character(kind=c_char),dimension(*) :: config_file, config_file_prefix, &
                                           name
  end function hio_init_mpi


#else
! This following interface may 

   integer(kind(HIO_SUCCESS)) &
        function hio_init_mpi (new_context, comm, config_file, config_file_prefix, &
        name) bind(c)
    import :: c_ptr, c_char, mpi_comm  
    type(c_ptr)                         :: new_context
    integer(kind(MPI_COMM_WORLD))       :: comm
    character(kind=c_char),dimension(*) :: config_file, config_file_prefix, &
                                           name
  end function hio_init_mpi

#endif


   integer(kind(HIO_SUCCESS)) &
        function hio_fini(context) bind(c)
   import :: c_ptr
   type(c_ptr)                :: context
 end function hio_fini



   integer(kind(HIO_SUCCESS)) &
        function hio_err_get_last (context, error_string) bind(c)
     import :: c_ptr
     type(c_ptr)                :: context
     type(c_ptr)                :: error_string ! 'char **error_string'
   end function hio_err_get_last



! The following three C functions, dealing with printing stuff, make use of what are
! known as 'variadic functions' that can accept a variable number of arguments.
! I may need to create an interface that is abstract. It may also be the case that these
! functions will need to be unavailable to calling from a Fortran code.

  !  integer(kind(HIO_SUCCESS)) &
  !       function hio_err_print_last(context, output, frmt, va_args) bind(c)
  !    import :: c_ptr, c_char
  !    type(c_ptr)                         :: output
  !    character(kind=c_char),dimension(*) :: frmt
  !    character
  ! end function hio_err_print_last

  !  integer(kind(HIO_SUCCESS)) &
  !  function hio_err_print_all(context, output, frmt, va_args) bind(c)
  !    import :: c_ptr, c_char
  !    type(c_ptr)                         :: output
  !    character(kind=c_char),dimension(*) :: frmt
  !    character
  !  end function hio_err_print_all

  !  integer(kind(HIO_SUCCESS)) &
  !  function hio_print_vars(object, type_rx, name_rx, output, frmt, va_args) bind(c)
  !    import :: c_ptr, c_char
  !    type(c_ptr), value                  :: object
  !    character(kind=c_char),dimension(*) :: type_rx, name_rx, format
  !    type(c_ptr)                         :: output

!END VA_ARGS C FUNCTIONS

   integer(kind(HIO_SUCCESS)) &
   function hio_dataset_alloc(context, set_out, name, set_id, flags, mode) &
       bind(c)
     import :: c_ptr, c_char, c_int, c_int64_t, HIO_SET_ELEMENT_UNIQUE
     type(c_ptr)                           :: set_out
     character(kind=c_char),dimension(*)   :: name
     integer(kind=c_int64_t)               :: set_id
     integer(kind=c_int)                   :: flags
     integer(kind(HIO_SET_ELEMENT_UNIQUE)) :: mode
   end function hio_dataset_alloc

   integer(kind(HIO_SUCCESS)) &
   function hio_dataset_free(dset) bind(c)
     import :: c_ptr     
     type(c_ptr)                :: dset     
   end function hio_dataset_free

   integer(kind(HIO_SUCCESS)) &
   function hio_dataset_get_id(dset, set_id) bind(c)
     import :: c_ptr, c_int64_t
     type(c_ptr),value          :: dset
     integer(kind=c_int64_t)    :: set_id
   end function hio_dataset_get_id

   integer(kind(HIO_SUCCESS)) &
   function hio_dataset_unlink(context, name, set_id, mode) bind(c)
     import :: c_ptr, c_char, c_int64_t, HIO_UNLINK_MODE_ALL
     type(c_ptr),value                        :: context
     character(kind=c_char),dimension(*)      :: name
     integer(kind=c_int64_t), value           :: set_id
     integer(kind(HIO_UNLINK_MODE_ALL)),value :: mode
   end function hio_dataset_unlink

   integer(kind(HIO_SUCCESS)) &
   function hio_element_open(dset, elem_out, name, flags) bind(c)
     import :: c_ptr, c_char, c_int
     type(c_ptr),value                   :: dset
     type(c_ptr)                         :: elem_out
     character(kind=c_char),dimension(*) :: name
     integer(kind=c_int),value           :: flags
   end function hio_element_open

   integer(kind(HIO_SUCCESS)) &
   function hio_element_size(element, e_size) bind(c)
     import :: c_ptr, c_int64_t
     type(c_ptr),value          :: element
     integer(kind=c_int64_t)    :: e_size
   end function hio_element_size

   integer(kind(HIO_SUCCESS)) &
   function hio_dataset_construct(dataset, destination, flags) bind(c)
     import :: c_ptr, c_char, c_int
     type(c_ptr),value                   :: dataset
     character(kind=c_char),dimension(*) :: destination
     integer(kind=c_int),value           :: flags
   end function hio_dataset_construct

   integer(kind(HIO_SUCCESS)) &
   function hio_element_close(element) bind(c)
     import :: c_ptr
     type(c_ptr) :: element
   end function hio_element_close

   integer(kind(HIO_SUCCESS)) &
   function hio_complete(element)
     import :: c_ptr
     type(c_ptr),value :: element
   end function hio_complete

   integer(kind(HIO_SUCCESS)) &
   function hio_request_test(requests, nrequests, bytes_trns, complete) bind(c)
     import :: c_ptr, c_int, c_size_t, c_bool
     type(c_ptr)               :: requests
     integer(kind=c_int),value :: nrequests
     integer(kind=c_size_t)    :: bytes_trns
     logical(kind=c_bool)      :: complete
   end function hio_request_test

   integer(kind(HIO_SUCCESS)) &
   function hio_request_wait(requests, nrequests, bytes_trns) bind(c)
     import :: c_ptr, c_int, c_size_t
     type(c_ptr)               :: requests
     integer(kind=c_int),value :: nrequests
     integer(kind=c_size_t)    :: bytes_trns
   end function hio_request_wait

   integer(kind(HIO_SUCCESS)) &
   function hio_dataset_should_checkpoint(context, name) bind(c)
     import :: c_ptr, c_char
     type(c_ptr),value                   :: context
     character(kind=c_char),dimension(*) :: name
   end function hio_dataset_should_checkpoint

   integer(kind(HIO_SUCCESS)) &
   function hio_config_set_value(object, variable, strval) bind(c)
     import :: c_ptr, c_char
     type(c_ptr),value                   :: object
     character(kind=c_char),dimension(*) :: variable
     character(kind=c_char),dimension(*) :: strval
   end function hio_config_set_value

   integer(kind(HIO_SUCCESS)) &
   function hio_config_get_value(object, variable, strval) bind(c)
     import :: c_ptr, c_char
     type(c_ptr),value                   :: object
     character(kind=c_char),dimension(*) :: variable
     type(c_ptr)                         :: strval
   end function hio_config_get_value
   
   
   integer(kind(HIO_SUCCESS)) &
   function hio_config_get_count(object, count) bind(c)
     import :: c_ptr, c_int
     type(c_ptr),value   :: object
     integer(kind=c_int) :: count
   end function hio_config_get_count
   
   integer(kind(HIO_SUCCESS)) &
   function hio_config_get_info(object, index, name, conftype, readonly)
     import :: c_ptr, c_int, c_bool, HIO_CONFIG_TYPE_BOOL
     type(c_ptr),value                   :: object
     integer(kind=c_int)                 :: index
     type(c_ptr)                         :: name
     integer(kind(HIO_CONFIG_TYPE_BOOL)) :: conftype
     logical(kind=c_bool)                :: readonly
   end function hio_config_get_info

   integer(kind(HIO_SUCCESS)) &
   function hio_perf_get_count(object, count) bind(c)
     import :: c_ptr, c_int
     type(c_ptr),value   :: object
     integer(kind=c_int) :: count
   end function hio_perf_get_count
   
   integer(kind(HIO_SUCCESS)) &
   function hio_perf_get_info(object, index, name, conftype) bind(c)
     import :: c_ptr, c_int, HIO_CONFIG_TYPE_BOOL
     type(c_ptr),value :: object, name
     integer(kind=c_int) :: index
     integer(kind(HIO_CONFIG_TYPE_BOOL)) :: conftype
   end function hio_perf_get_info

   integer(kind(HIO_SUCCESS)) &
   function hio_perf_get_value(object, variable, value, value_len) bind(c)
     import :: c_ptr, c_char, c_size_t
     type(c_ptr),value                   :: object
     character(kind=c_char),dimension(*) :: variable
     type(c_ptr)                         :: value
     integer(kind=c_size_t)              :: value_len
   end function hio_perf_get_value

   integer(kind(HIO_SUCCESS)) &
   function hio_object_get_name(object, name_out) bind(c)
     import :: c_ptr
     type(c_ptr),value :: object
     type(c_ptr)       :: name_out
   end function hio_object_get_name
     
end interface

! Abstract interfaces to various HIO operations where the calling API
! have similar signatures.

abstract interface
   integer(kind(HIO_SUCCESS)) &
   function hio_dataset_opcl(dset) bind(c)
     import :: c_ptr
     type(c_ptr),value          :: dset
   end function hio_dataset_opcl
   
   integer(c_size_t) &
   function hio_element_rw(element, offset, reserved0, ptr, count, size) bind(c)
   import :: c_int64_t, c_long, c_ptr, c_size_t
     type(c_ptr), value            :: element
     integer(kind=c_long),value    :: offset
     integer(kind=c_int64_t),value :: reserved0
     type(c_ptr)                   :: ptr
     integer(kind=c_size_t),value  :: count
     integer(kind=c_size_t),value  :: size
   end function hio_element_rw

   integer(kind(HIO_SUCCESS)) &
   function hio_element_rw_nb(element, request, offset, reserved0, ptr, count, size) bind(c)
     import :: c_ptr, c_long, c_int64_t, c_size_t
     type(c_ptr), value            :: element
     type(c_ptr)                   :: request
     integer(kind=c_long),value    :: offset
     integer(kind=c_int64_t),value :: reserved0
     type(c_ptr)                   :: ptr
     integer(kind=c_size_t),value  :: count
     integer(kind=c_size_t),value  :: size
   end function hio_element_rw_nb

   integer(c_size_t) &
   function hio_element_rw_strided(element, offset, reserved0, ptr, count, size, stride) bind(c)
     import :: c_ptr, c_long, c_int64_t, c_size_t
     type(c_ptr), value            :: element
     integer(kind=c_long),value    :: offset
     integer(kind=c_int64_t),value :: reserved0
     type(c_ptr)                   :: ptr
     integer(kind=c_size_t),value  :: count
     integer(kind=c_size_t),value  :: size
     integer(kind=c_size_t),value  :: stride
   end function hio_element_rw_strided

   integer(kind(HIO_SUCCESS)) &
   function hio_element_rw_strided_nb(element, request, offset, reserved0, ptr, count, size, stride) bind(c)
     import :: c_ptr, c_long, c_int64_t, c_size_t
     type(c_ptr), value            :: element
     type(c_ptr)                   :: request
     integer(kind=c_long),value    :: offset
     integer(kind=c_int64_t),value :: reserved0
     type(c_ptr)                   :: ptr
     integer(kind=c_size_t),value  :: count
     integer(kind=c_size_t),value  :: size
     integer(kind=c_size_t),value  :: stride
   end function hio_element_rw_strided_nb

   integer(kind(HIO_SUCCESS)) &
   function hio_flush(item, mode) bind(c)
     import :: c_ptr, HIO_FLUSH_MODE_LOCAL
     type(c_ptr),value                         :: item
     integer(kind(HIO_FLUSH_MODE_LOCAL)),value :: mode
   end function hio_flush
end interface

procedure(hio_dataset_opcl)          :: hio_dataset_open, hio_dataset_close
procedure(hio_element_rw)            :: hio_element_read, hio_element_write
procedure(hio_element_rw_nb)         :: hio_element_read_nb, hio_element_write_nb
procedure(hio_element_rw_strided)    :: hio_element_read_strided, hio_element_write_strided
procedure(hio_element_rw_strided_nb) :: hio_element_read_strided_nb, hio_element_write_strided_nb
procedure(hio_flush)                 :: hio_element_flush, hio_dataset_flush

end module hio
