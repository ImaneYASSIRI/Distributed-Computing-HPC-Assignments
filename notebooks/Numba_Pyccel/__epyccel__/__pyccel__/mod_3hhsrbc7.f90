module mod_3hhsrbc7

use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function solve_2d_nonlinearconv_pyccel(u, un, v, vn, nt, dt, dx, dy, c) &
      result(Out_0001)

  implicit none

  integer(C_INT64_T) :: Out_0001
  real(C_DOUBLE), intent(inout) :: u(0:,0:)
  real(C_DOUBLE), intent(inout) :: un(0:,0:)
  real(C_DOUBLE), intent(inout) :: v(0:,0:)
  real(C_DOUBLE), intent(inout) :: vn(0:,0:)
  integer(C_INT64_T), value :: nt
  real(C_DOUBLE), value :: dt
  real(C_DOUBLE), value :: dx
  real(C_DOUBLE), value :: dy
  integer(C_INT64_T), value :: c
  integer(C_INT64_T) :: row
  integer(C_INT64_T) :: col
  real(C_DOUBLE), allocatable :: u0(:,:)
  real(C_DOUBLE), allocatable :: v0(:,:)
  integer(C_INT64_T) :: i
  integer(C_INT64_T) :: j
  integer(C_INT64_T) :: n

  !##Assign initial conditions
  !#set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
  u(Int(0.5_C_DOUBLE / dx, C_INT64_T):Int(1_C_INT64_T / dx + 1_C_INT64_T &
      , C_INT64_T) - 1_C_INT64_T, Int(0.5_C_DOUBLE / dy, C_INT64_T):Int &
      (1_C_INT64_T / dy + 1_C_INT64_T, C_INT64_T) - 1_C_INT64_T) = &
      2_C_INT64_T
  !#set hat function I.C. : v(.5<=x<=1 && .5<=y<=1 ) is 2
  v(Int(0.5_C_DOUBLE / dx, C_INT64_T):Int(1_C_INT64_T / dx + 1_C_INT64_T &
      , C_INT64_T) - 1_C_INT64_T, Int(0.5_C_DOUBLE / dy, C_INT64_T):Int &
      (1_C_INT64_T / dy + 1_C_INT64_T, C_INT64_T) - 1_C_INT64_T) = &
      2_C_INT64_T
  row = size(u, 2_C_INT64_T - 0_C_INT64_T, C_INT64_T)
  col = size(u, 2_C_INT64_T - 1_C_INT64_T, C_INT64_T)
  allocate(u0(0:col - 1_C_INT64_T, 0:row - 1_C_INT64_T))
  u0 = 1.0_C_DOUBLE
  allocate(v0(0:col - 1_C_INT64_T, 0:row - 1_C_INT64_T))
  v0 = 1.0_C_DOUBLE
  do i = 0_C_INT64_T, row - 1_C_INT64_T, 1_C_INT64_T
    do j = 0_C_INT64_T, col - 1_C_INT64_T, 1_C_INT64_T
      u0(j, i) = u(j, i)
      v0(j, i) = v(j, i)
    end do
  end do
  !fill the update of u and v
  do n = 0_C_INT64_T, nt - 1_C_INT64_T, 1_C_INT64_T
    do i = 0_C_INT64_T, row - 1_C_INT64_T, 1_C_INT64_T
      do j = 0_C_INT64_T, col - 1_C_INT64_T, 1_C_INT64_T
        un(j, i) = u(j, i)
        vn(j, i) = v(j, i)
      end do
    end do
    do i = 1_C_INT64_T, row - 1_C_INT64_T, 1_C_INT64_T
      do j = 1_C_INT64_T, col - 1_C_INT64_T, 1_C_INT64_T
        u(j, i) = un(j, i) - u0(j, i) * dt / dx * (un(j, i) - un(j, i - &
      1_C_INT64_T)) - vn(j, i) * dt / dy * (un(j, i) - un(j - &
      1_C_INT64_T, i))
        v(j, i) = vn(j, i) - v0(j, i) * dt / dx * (vn(j, i) - vn(j, i - &
      1_C_INT64_T)) - vn(j, i) * dt / dy * (vn(j, i) - vn(j - &
      1_C_INT64_T, i))
      end do
    end do
  end do
  Out_0001 = 0_C_INT64_T
  return

end function solve_2d_nonlinearconv_pyccel
!........................................

end module mod_3hhsrbc7
