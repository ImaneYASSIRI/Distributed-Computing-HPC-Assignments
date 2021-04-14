module mod_dqc4xnth

use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function solve_1d_linearconv_pyccel(u, un, nt, nx, dt, dx, c) result( &
      Out_0001)

  implicit none

  integer(C_INT64_T) :: Out_0001
  real(C_DOUBLE), intent(inout) :: u(0:)
  real(C_DOUBLE), intent(inout) :: un(0:)
  integer(C_INT64_T), value :: nt
  integer(C_INT64_T), value :: nx
  real(C_DOUBLE), value :: dt
  real(C_DOUBLE), value :: dx
  integer(C_INT64_T), value :: c
  integer(C_INT64_T) :: n
  integer(C_INT64_T) :: i

  do n = 0_C_INT64_T, nt - 1_C_INT64_T, 1_C_INT64_T
    do i = 0_C_INT64_T, nx - 1_C_INT64_T, 1_C_INT64_T
      un(i) = u(i)
    end do
    do i = 1_C_INT64_T, nx - 1_C_INT64_T, 1_C_INT64_T
      u(i) = un(i) - c * dt / dx * (un(i) - un(i - 1_C_INT64_T))
    end do
  end do
  Out_0001 = 0_C_INT64_T
  return

end function solve_1d_linearconv_pyccel
!........................................

end module mod_dqc4xnth
