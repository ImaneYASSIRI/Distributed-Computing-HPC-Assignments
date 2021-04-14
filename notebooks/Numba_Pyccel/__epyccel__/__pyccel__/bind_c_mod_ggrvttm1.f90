module bind_c_mod_ggrvttm1

use mod_ggrvttm1, only: solve_2d_diff_pyccel
use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function bind_c_solve_2d_diff_pyccel(n0_u, n1_u, u, n0_un, n1_un, un, nt &
      , dt, dx, dy, nu) bind(c) result(Out_0001)

  implicit none

  integer(C_INT64_T), value :: n0_u
  integer(C_INT64_T), value :: n1_u
  real(C_DOUBLE), intent(inout) :: u(0:n1_u-1,0:n0_u-1)
  integer(C_INT64_T), value :: n0_un
  integer(C_INT64_T), value :: n1_un
  real(C_DOUBLE), intent(inout) :: un(0:n1_un-1,0:n0_un-1)
  integer(C_INT64_T), value :: nt
  real(C_DOUBLE), value :: dt
  real(C_DOUBLE), value :: dx
  real(C_DOUBLE), value :: dy
  integer(C_INT64_T), value :: nu
  integer(C_INT64_T) :: Out_0001

  Out_0001 = solve_2d_diff_pyccel(u, un, nt, dt, dx, dy, nu)

end function bind_c_solve_2d_diff_pyccel
!........................................

end module bind_c_mod_ggrvttm1
