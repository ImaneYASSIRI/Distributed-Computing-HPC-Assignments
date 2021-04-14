module mod_cehjx3r8

use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function solve_2d_linearconv_pyccel(u, un, nt, dt, dx, dy, c) result( &
      Out_0001)

  implicit none

  integer(C_INT64_T) :: Out_0001
  real(C_DOUBLE), intent(inout) :: u(0:,0:)
  real(C_DOUBLE), intent(inout) :: un(0:,0:)
  integer(C_INT64_T), value :: nt
  real(C_DOUBLE), value :: dt
  real(C_DOUBLE), value :: dx
  real(C_DOUBLE), value :: dy
  integer(C_INT64_T), value :: c
  integer(C_INT64_T) :: row
  integer(C_INT64_T) :: col
  integer(C_INT64_T) :: n
  integer(C_INT64_T) :: i
  integer(C_INT64_T) :: j

  row = size(u, 2_C_INT64_T - 0_C_INT64_T, C_INT64_T)
  col = size(u, 2_C_INT64_T - 1_C_INT64_T, C_INT64_T)
  !fill the update of u and v
  do n = 0_C_INT64_T, nt - 1_C_INT64_T, 1_C_INT64_T
    do i = 0_C_INT64_T, row - 1_C_INT64_T, 1_C_INT64_T
      do j = 0_C_INT64_T, col - 1_C_INT64_T, 1_C_INT64_T
        un(j, i) = u(j, i)
      end do
    end do
    do i = 1_C_INT64_T, row - 1_C_INT64_T, 1_C_INT64_T
      do j = 1_C_INT64_T, col - 1_C_INT64_T, 1_C_INT64_T
        u(j, i) = un(j, i) - c * dt / dx * (un(j, i) - un(j, i - &
      1_C_INT64_T)) - c * dt / dy * (un(j, i) - un(j - 1_C_INT64_T, i))
      end do
    end do
  end do
  Out_0001 = 0_C_INT64_T
  return

end function solve_2d_linearconv_pyccel
!........................................

end module mod_cehjx3r8
