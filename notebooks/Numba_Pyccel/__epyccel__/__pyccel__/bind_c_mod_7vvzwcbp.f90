module bind_c_mod_7vvzwcbp

use mod_7vvzwcbp, only: solve_1d_diff_pyccel
use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function bind_c_solve_1d_diff_pyccel(n0_u, u, n0_un, un, nt, nx, dt, dx, &
      nu) bind(c) result(Out_0001)

  implicit none

  integer(C_INT64_T), value :: n0_u
  real(C_DOUBLE), intent(inout) :: u(0:n0_u-1)
  integer(C_INT64_T), value :: n0_un
  real(C_DOUBLE), intent(inout) :: un(0:n0_un-1)
  integer(C_INT64_T), value :: nt
  integer(C_INT64_T), value :: nx
  real(C_DOUBLE), value :: dt
  real(C_DOUBLE), value :: dx
  real(C_DOUBLE), value :: nu
  integer(C_INT64_T) :: Out_0001

  Out_0001 = solve_1d_diff_pyccel(u, un, nt, nx, dt, dx, nu)

end function bind_c_solve_1d_diff_pyccel
!........................................

end module bind_c_mod_7vvzwcbp
