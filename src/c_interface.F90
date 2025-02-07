
module w90_library_c
!Fortran 2018: Assumed-length character dummy argument ‘keyword’ at (1) of procedure ‘cset_option_int’ with BIND(C) attribute
  use iso_c_binding
  use w90_library
  implicit none

  public

  type, bind(c) :: w90_data
    type(c_ptr) :: caddr
  end type

contains

  subroutine w90_create(w90_obj) bind(c)
    !! return a c-pointer to a instance of the wannier90 library data structure
    type(lib_common_type), pointer :: common_data
    type(w90_data) :: w90_obj
    if (c_associated(w90_obj%caddr)) return
    allocate (common_data)
    w90_obj%caddr = c_loc(common_data)
  end subroutine

  subroutine w90_delete(w90_obj) bind(c)
    !! deallocates/clears a c-pointer to a instance of the wannier90 library data structure
    implicit none
    type(w90_data) :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    if (.not. c_associated(w90_obj%caddr)) return
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    deallocate (w90_fptr)
    w90_obj%caddr = C_NULL_PTR
  end subroutine

  subroutine w90_disentangle_c(w90_obj, ierr) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_disentangle(w90_fptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_get_centres_c(w90_obj, centres) bind(c)
    ! returns the centres of calulated mlwfs
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: centres
    type(lib_common_type), pointer :: w90_fptr
    real(kind=8), pointer :: fcentres(:, :)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(centres, fcentres, [3, w90_fptr%num_wann])
    call w90_get_centres(w90_fptr, fcentres)
  end subroutine

  subroutine w90_get_gkpb_c(w90_obj, gkpb) bind(c)
    ! return the g-offset of adjacent k-points in finite difference scheme
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: gkpb
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: nfptr(:, :, :)
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(gkpb, nfptr, [3, w90_fptr%num_kpts, w90_fptr%kmesh_info%nntot])
    call w90_get_gkpb(w90_fptr, nfptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_get_nn_c(w90_obj, n) bind(c)
    ! return the number of adjacent k-points in finite difference scheme
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: n
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: ndat
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(n, ndat)
    call w90_get_nn(w90_fptr, ndat, istdout, istderr, ierr)
  end subroutine

  subroutine w90_get_nnkp_c(w90_obj, nnkp) bind(c)
    ! return the indexing of adjacent k-points in finite difference scheme
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: nnkp
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int), pointer :: nfptr(:, :)
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(nnkp, nfptr, [w90_fptr%num_kpts, w90_fptr%kmesh_info%nntot])
    call w90_get_nnkp(w90_fptr, nfptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_get_spreads_c(w90_obj, spreads) bind(c)
    ! returns the spreads of calulated mlwfs
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: spreads
    type(lib_common_type), pointer :: w90_fptr
    real(kind=8), pointer :: fspreads(:)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(spreads, fspreads, [w90_fptr%num_wann])
    call w90_get_spreads(w90_fptr, fspreads)
  end subroutine

  subroutine w90_input_setopt_c(w90_obj, seedname, ierr) bind(c)
    ! specify parameters through the library interface
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: seedname
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_input_setopt(w90_fptr, seedname, istdout, istderr, ierr)
  end subroutine

  subroutine w90_input_reader_c(w90_obj, ierr) bind(c)
    ! read (optional) parameters from .win file
    implicit none
    type(w90_data), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istderr, istdout, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_input_reader(w90_fptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_project_overlap_c(w90_obj, ierr) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_project_overlap(w90_fptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_set_eigval_c(w90_obj, eigval_cptr) bind(c)
    ! copy a pointer to eigenvalue data
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: eigval_cptr
    real(8), pointer :: eigval_fptr(:, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(eigval_cptr, eigval_fptr, [w90_fptr%num_bands, w90_fptr%num_kpts])
    call w90_set_eigval(w90_fptr, eigval_fptr)
  end subroutine

  subroutine w90_set_m_local_c(w90_obj, m_cptr) bind(c)
    ! copy a pointer to m-matrix data
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: m_cptr
    complex(8), pointer :: fptr(:, :, :, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(m_cptr, fptr, [w90_fptr%num_bands, w90_fptr%num_bands, w90_fptr%kmesh_info%nntot, w90_fptr%num_kpts])
    call w90_set_m_local(w90_fptr, fptr)
  end subroutine

  subroutine w90_set_u_matrix_c(w90_obj, a_cptr) bind(c)
    ! copy pointer to u-matrix
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: a_cptr
    complex(8), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(a_cptr, a_fptr, [w90_fptr%num_wann, w90_fptr%num_wann, w90_fptr%num_kpts]) ! these are reversed wrt c
    call w90_set_u_matrix(w90_fptr, a_fptr)
  end subroutine

  subroutine w90_set_u_opt_c(w90_obj, a_cptr) bind(c)
    ! copy pointer to u-matrix (also used for initial projections)
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value :: a_cptr
    complex(kind=8), pointer :: a_fptr(:, :, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(a_cptr, a_fptr, [w90_fptr%num_bands, w90_fptr%num_wann, w90_fptr%num_kpts]) ! these are reversed wrt c
    call w90_set_u_opt(w90_fptr, a_fptr)
  end subroutine

  subroutine w90_wannierise_c(w90_obj, ierr) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(lib_common_type), pointer :: w90_fptr
    integer(kind=c_int) :: istdout, istderr, ierr
    call w90_get_fortran_stderr(istderr)
    call w90_get_fortran_stdout(istdout)
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_wannierise(w90_fptr, istdout, istderr, ierr)
  end subroutine

  subroutine w90_set_option_double_c(w90_obj, keyword, cdble) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    real(kind=c_double), value  :: cdble
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, cdble)
  end subroutine

  subroutine w90_set_option_double1d_c(w90_obj, keyword, arg_cptr, x) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    type(c_ptr), value  :: arg_cptr
    integer(kind=c_int), value :: x
    type(lib_common_type), pointer :: w90_fptr
    real(kind=8), pointer :: fptr(:)
    call c_f_pointer(arg_cptr, fptr, [x])
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, fptr)
  end subroutine

  subroutine w90_set_option_double2d_c(w90_obj, keyword, arg_cptr, x, y) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    type(c_ptr), value  :: arg_cptr
    integer(kind=c_int), value :: x, y
    type(lib_common_type), pointer :: w90_fptr
    real(kind=8), pointer :: fptr(:, :)
    call c_f_pointer(arg_cptr, fptr, [y, x]) ! these are reversed wrt c
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, fptr)
  end subroutine

  subroutine w90_set_option_int_c(w90_obj, keyword, cint) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    integer(kind=c_int), value  :: cint
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, cint)
  end subroutine

  subroutine w90_set_option_int1d_c(w90_obj, keyword, arg_cptr, x) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value ::  arg_cptr
    character(*, kind=c_char) :: keyword
    integer(kind=c_int), value  :: x
    integer, pointer :: fptr(:)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(arg_cptr, fptr, [x])
    call w90_set_option(w90_fptr, keyword, fptr)
  end subroutine

  subroutine w90_set_option_int2d_c(w90_obj, keyword, arg_cptr, x, y) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    type(c_ptr), value ::  arg_cptr
    character(*, kind=c_char) :: keyword
    integer(kind=c_int), value  :: x, y
    integer(kind=c_int), pointer :: fptr(:, :)
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call c_f_pointer(arg_cptr, fptr, [y, x]) ! these are reversed wrt c
    call w90_set_option(w90_fptr, keyword, fptr)
  end subroutine

  subroutine w90_set_option_logical_c(w90_obj, keyword, bool) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    logical(kind=c_bool), value  :: bool
    logical :: fbool
    type(lib_common_type), pointer :: w90_fptr
    fbool = bool
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, fbool)
  end subroutine

  subroutine w90_set_option_text_c(w90_obj, keyword, text) bind(c)
    implicit none
    type(w90_data), value :: w90_obj
    character(*, kind=c_char) :: keyword
    character(*, kind=c_char) :: text
    type(lib_common_type), pointer :: w90_fptr
    call c_f_pointer(w90_obj%caddr, w90_fptr)
    call w90_set_option(w90_fptr, keyword, text)
  end subroutine

end module w90_library_c
