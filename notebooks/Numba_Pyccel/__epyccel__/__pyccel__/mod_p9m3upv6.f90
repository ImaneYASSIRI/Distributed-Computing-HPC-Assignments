module mod_p9m3upv6

use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function solve_1d_diff_pyccel(u, un, nt, nx, dt, dx, nu) result(Out_0001 &
      )

  implicit none

  integer(C_INT64_T) :: Out_0001
  real(C_DOUBLE), intent(inout) :: u(0:)
  real(C_DOUBLE), intent(inout) :: un(0:)
  integer(C_INT64_T), value :: nt
  integer(C_INT64_T), value :: nx
  real(C_DOUBLE), value :: dt
  real(C_DOUBLE), value :: dx
  real(C_DOUBLE), value :: nu
  integer(C_INT64_T) :: n
  integer(C_INT64_T) :: i

  !fill the update of u
  do n = 0_C_INT64_T, nt - 1_C_INT64_T, 1_C_INT64_T
    do i = 0_C_INT64_T, nx - 1_C_INT64_T, 1_C_INT64_T
      un(i) = u(i)
    end do
    do i = 1_C_INT64_T, nx - 1_C_INT64_T, 1_C_INT64_T
      u(i) = un(i) - nu * dt / dx ** 2_C_INT64_T * (u(i) - 2_C_INT64_T * &
      un(i) + un(i - 1_C_INT64_T))
    end do
  end do
  Out_0001 = 0_C_INT64_T
  return

end function solve_1d_diff_pyccel
!........................................

end module mod_p9m3upv6
