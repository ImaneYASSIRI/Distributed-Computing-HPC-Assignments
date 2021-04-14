module bind_c_mod_dqc4xnth

use mod_dqc4xnth, only: solve_1d_linearconv_pyccel
use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function bind_c_solve_1d_linearconv_pyccel(n0_u, u, n0_un, un, nt, nx, &
      dt, dx, c) bind(c) result(Out_0001)

  implicit none

  integer(C_INT64_T), value :: n0_u
  real(C_DOUBLE), intent(inout) :: u(0:n0_u-1)
  integer(C_INT64_T), value :: n0_un
  real(C_DOUBLE), intent(inout) :: un(0:n0_un-1)
  integer(C_INT64_T), value :: nt
  integer(C_INT64_T), value :: nx
  real(C_DOUBLE), value :: dt
  real(C_DOUBLE), value :: dx
  integer(C_INT64_T), value :: c
  integer(C_INT64_T) :: Out_0001

  Out_0001 = solve_1d_linearconv_pyccel(u, un, nt, nx, dt, dx, c)

end function bind_c_solve_1d_linearconv_pyccel
!........................................

end module bind_c_mod_dqc4xnth
