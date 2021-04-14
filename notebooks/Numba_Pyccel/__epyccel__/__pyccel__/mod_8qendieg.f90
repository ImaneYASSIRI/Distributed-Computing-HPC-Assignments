module mod_8qendieg

use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function solve_2d_diff_pyccel(u, un, nt, dt, dx, dy, nu) result(Out_0001 &
      )

  implicit none

  integer(C_INT32_T) :: Out_0001
  real(C_DOUBLE), intent(inout) :: u(0:,0:)
  real(C_DOUBLE), intent(inout) :: un(0:,0:)
  integer(C_INT32_T), value :: nt
  real(C_DOUBLE), value :: dt
  real(C_DOUBLE), value :: dx
  real(C_DOUBLE), value :: dy
  real(C_DOUBLE), value :: nu
  integer(C_INT32_T) :: row
  integer(C_INT32_T) :: col
  integer(C_INT32_T) :: n
  integer(C_INT32_T) :: i
  integer(C_INT32_T) :: j

  row = size(u, 2, C_INT32_T)
  col = size(u, 1, C_INT32_T)
  !#Assign initial conditions
  !set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
  u(Int(0.5_C_DOUBLE / dx, C_INT32_T):Int(1_C_INT32_T / dx + 1_C_INT32_T &
      , C_INT32_T) - 1_C_INT32_T, Int(0.5_C_DOUBLE / dy, C_INT32_T):Int &
      (1_C_INT32_T / dy + 1_C_INT32_T, C_INT32_T) - 1_C_INT32_T) = &
      2_C_INT32_T
  !fill the update of u and v
  do n = 0_C_INT32_T, nt-1_C_INT32_T, 1_C_INT32_T
    do i = 0_C_INT32_T, row-1_C_INT32_T, 1_C_INT32_T
      do j = 0_C_INT32_T, col-1_C_INT32_T, 1_C_INT32_T
        un(j, i) = u(j, i)
      end do
    end do
    do i = 1_C_INT32_T, row - 1_C_INT32_T-1_C_INT32_T, 1_C_INT32_T
      do j = 1_C_INT32_T, col - 1_C_INT32_T-1_C_INT32_T, 1_C_INT32_T
        u(j, i) = un(j, i) + nu * dt / dx ** 2_C_INT32_T * (un(j, i + &
      1_C_INT32_T) - 2_C_INT32_T * un(j, i) + un(j, i - 1_C_INT32_T)) + &
      nu * dt / dy ** 2_C_INT32_T * (un(j + 1_C_INT32_T, i) - &
      2_C_INT32_T * un(j, i) + un(j - 1_C_INT32_T, i))
      end do
    end do
  end do
  Out_0001 = 0_C_INT32_T
  return

end function solve_2d_diff_pyccel
!........................................

end module mod_8qendieg
