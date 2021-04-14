module mod_fseu68qo

use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function solve_2d_poisson_pyccel(p, pd, b, nx, ny, nt, dx, dy) result( &
      Out_0001)

  implicit none

  integer(C_INT32_T) :: Out_0001
  real(C_DOUBLE), intent(inout) :: p(0:,0:)
  real(C_DOUBLE), intent(inout) :: pd(0:,0:)
  real(C_DOUBLE), intent(inout) :: b(0:,0:)
  integer(C_INT32_T), value :: nx
  integer(C_INT32_T), value :: ny
  integer(C_INT32_T), value :: nt
  real(C_DOUBLE), value :: dx
  real(C_DOUBLE), value :: dy
  integer(C_INT32_T) :: row
  integer(C_INT32_T) :: col
  integer(C_INT32_T) :: it
  integer(C_INT32_T) :: i
  integer(C_INT32_T) :: j

  row = size(p, 2, C_INT32_T)
  col = size(p, 1, C_INT32_T)
  !Source
  b(Int(Real(nx, C_DOUBLE) / 4.0_C_DOUBLE, C_INT32_T), Int(Real(ny, &
      C_DOUBLE) / 4.0_C_DOUBLE, C_INT32_T)) = 100_C_INT32_T
  b(Int(Real(3_C_INT32_T * nx, C_DOUBLE) / 4.0_C_DOUBLE, C_INT32_T), Int &
      (Real(3_C_INT32_T * ny, C_DOUBLE) / 4.0_C_DOUBLE, C_INT32_T)) = &
      -100_C_INT32_T
  do it = 0_C_INT32_T, nt-1_C_INT32_T, 1_C_INT32_T
    do i = 0_C_INT32_T, nx-1_C_INT32_T, 1_C_INT32_T
      pd(:, i) = p(:, i)
    end do
    do j = 2_C_INT32_T, row-1_C_INT32_T, 1_C_INT32_T
      do i = 2_C_INT32_T, col-1_C_INT32_T, 1_C_INT32_T
        p(i - 1_C_INT32_T, j - 1_C_INT32_T) = ((pd(i, j - 1_C_INT32_T) + &
      pd(i - 2_C_INT32_T, j - 1_C_INT32_T)) * dy ** 2_C_INT32_T + (pd(i &
      - 1_C_INT32_T, j) + pd(i - 1_C_INT32_T, j - 2_C_INT32_T)) * dx ** &
      2_C_INT32_T - b(i - 1_C_INT32_T, j - 1_C_INT32_T) * dx ** &
      2_C_INT32_T * dy ** 2_C_INT32_T) / (2_C_INT32_T * (dx ** &
      2_C_INT32_T + dy ** 2_C_INT32_T))
      end do
    end do
    p(:, 0_C_INT32_T) = 0_C_INT32_T
    p(:, ny - 1_C_INT32_T) = 0_C_INT32_T
    p(0_C_INT32_T, :) = 0_C_INT32_T
    p(nx - 1_C_INT32_T, :) = 0_C_INT32_T
  end do
  Out_0001 = 0_C_INT32_T
  return

end function solve_2d_poisson_pyccel
!........................................

end module mod_fseu68qo
