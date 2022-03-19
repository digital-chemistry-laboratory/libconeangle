module test_integration
  use stdlib_sorting, only: sort
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use coneangle_env, only: wp
  use coneangle_main, only: cone_angle
  implicit none(type, external)
  private

  public :: collect_suite_integration

contains

!> Collect all exported unit tests
  subroutine collect_suite_integration(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("PdCO", test_PdCO), &
                new_unittest("Pdbpy", test_Pdbpy), &
                new_unittest("PdPMe3", test_PdPMe3), &
                new_unittest("PdCO_close", test_PdCO_close, should_fail=.true.), &
                new_unittest("PdCO_mismatch", test_PdCO_mismatch, should_fail=.true.), &
                new_unittest("PdCO_bounds", test_PdCO_bounds, should_fail=.true.) &
                ]

  end subroutine collect_suite_integration

  subroutine test_PdCO(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(3, 3), radii(3), alpha, axis(3), ref_axis(3)
    integer :: tangent_atoms(3), ref_tangent_atoms(3), stat, i
    character(:), allocatable :: errmsg

    ! Run cone angle calculation
    coordinates = reshape([0._wp, 0._wp, -0.52_wp, 0._wp, 0._wp, 1.76_wp, 0._wp, 0._wp, 2.86_wp], [3, 3])
    radii = [2.1_wp, 1.7_wp, 1.52_wp]
    call cone_angle(coordinates, radii, 1, alpha, axis, tangent_atoms, stat, errmsg)

    ! Check cone angle
    call check(error, alpha, 96.423_wp, thr=0.001_wp)
    if (allocated(error)) return

    ! Check cone indices
    ref_tangent_atoms = [2, 0, 0]
    do i = 1, 3
      call check(error, tangent_atoms(i), ref_tangent_atoms(i))
      if (allocated(error)) return
    end do

    ! Check cone axis
    ref_axis = [0._wp, 0._wp, 1._wp]
    do i = 1, 3
      call check(error, axis(i), ref_axis(i), thr=0.0001_wp)
      if (allocated(error)) return
    end do

  end subroutine test_PdCO

  subroutine test_PdCO_close(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(3, 3), radii(3), alpha, axis(3)
    integer :: tangent_atoms(3), stat
    character(:), allocatable :: errmsg

    ! Run cone angle calculation
    coordinates = reshape([0._wp, 0._wp, 0.1_wp, 0._wp, 0._wp, 1.76_wp, 0._wp, 0._wp, 2.86_wp], [3, 3])
    radii = [2.1_wp, 1.7_wp, 1.52_wp]
    call cone_angle(coordinates, radii, 1, alpha, axis, tangent_atoms, stat, errmsg)

    ! Check stat calculation failed.
    call check(error, stat, 0)
    if (allocated(error)) return

  end subroutine test_PdCO_close

  subroutine test_PdCO_mismatch(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(3, 3), radii(2), alpha, axis(3)
    integer :: tangent_atoms(3), stat
    character(:), allocatable :: errmsg

    ! Run cone angle calculation
    coordinates = reshape([0._wp, 0._wp, -0.52_wp, 0._wp, 0._wp, 1.76_wp, 0._wp, 0._wp, 2.86_wp], [3, 3])
    radii = [2.1_wp, 1.7_wp]
    call cone_angle(coordinates, radii, 1, alpha, axis, tangent_atoms, stat, errmsg)

    ! Check stat calculation failed.
    call check(error, stat, 0)
    if (allocated(error)) return

  end subroutine test_PdCO_mismatch

  subroutine test_PdCO_bounds(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(3, 3), radii(3), alpha, axis(3)
    integer :: tangent_atoms(3), stat
    character(:), allocatable :: errmsg

    ! Run cone angle calculation
    coordinates = reshape([0._wp, 0._wp, -0.52_wp, 0._wp, 0._wp, 1.76_wp, 0._wp, 0._wp, 2.86_wp], [3, 3])
    radii = [2.1_wp, 1.7_wp, 1.52_wp]
    call cone_angle(coordinates, radii, 0, alpha, axis, tangent_atoms, stat, errmsg)

    ! Check stat calculation failed.
    call check(error, stat, 0)
    if (allocated(error)) return

  end subroutine test_PdCO_bounds

  subroutine test_Pdbpy(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(3, 21), radii(21), alpha, axis(3), ref_axis(3)
    integer :: tangent_atoms(3), ref_tangent_atoms(3), stat, i
    character(:), allocatable :: errmsg

    ! Run cone angle calculation
    coordinates = reshape( &
                  [-1.899000e-03_wp, -1.846004e+00_wp, -2.007000e-03_wp, -1.343609e+00_wp, &
                   -8.359200e-02_wp, 4.197960e-01_wp, 1.345137e+00_wp, -8.591400e-02_wp, &
                   -4.193200e-01_wp, -7.399130e-01_wp, 1.062540e+00_wp, 8.798000e-03_wp, &
                   -1.473982e+00_wp, 2.183787e+00_wp, -3.961990e-01_wp, -2.865881e+00_wp, &
                   2.130461e+00_wp, -3.962120e-01_wp, -3.487476e+00_wp, 9.490540e-01_wp, &
                   1.374000e-02_wp, -2.690562e+00_wp, -1.192190e-01_wp, 4.151720e-01_wp, &
                   -9.542920e-01_wp, 3.073408e+00_wp, -7.377720e-01_wp, -3.451265e+00_wp, &
                   2.985788e+00_wp, -7.207940e-01_wp, -4.569174e+00_wp, 8.579540e-01_wp, &
                   3.831000e-02_wp, -3.135734e+00_wp, -1.043925e+00_wp, 7.700380e-01_wp, &
                   2.691916e+00_wp, -1.229550e-01_wp, -4.128390e-01_wp, 3.489163e+00_wp, &
                   9.451910e-01_wp, -1.169700e-02_wp, 2.868283e+00_wp, 2.127761e+00_wp, &
                   3.960000e-01_wp, 1.476591e+00_wp, 2.182376e+00_wp, 3.941950e-01_wp, &
                   7.424670e-01_wp, 1.061320e+00_wp, -1.040300e-02_wp, 3.136172e+00_wp, &
                   -1.048666e+00_wp, -7.660710e-01_wp, 4.570699e+00_wp, 8.531320e-01_wp, &
                   -3.487800e-02_wp, 3.453586e+00_wp, 2.983188e+00_wp, 7.202120e-01_wp, &
                   9.563590e-01_wp, 3.072406e+00_wp, 7.338950e-01_wp], [3, 21])
    radii = [2.1_wp, 1.55_wp, 1.55_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.1_wp, 1.1_wp, 1.1_wp, &
             1.1_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp]
    call cone_angle(coordinates, radii, 1, alpha, axis, tangent_atoms, stat, errmsg)

    ! Check cone angle
    call check(error, alpha, 190.799_wp, thr=0.001_wp)
    if (allocated(error)) return

    ! Check cone indices
    ref_tangent_atoms = [12, 18, 0]
    call sort(ref_tangent_atoms)
    call sort(tangent_atoms)
    do i = 1, 3
      call check(error, tangent_atoms(i), ref_tangent_atoms(i))
      if (allocated(error)) return
    end do

    ! Check cone axis
    ref_axis = [0.00200079_wp, 0.99998474_wp, 0.00514848_wp]
    do i = 1, 3
      call check(error, axis(i), ref_axis(i), thr=0.0001_wp)
      if (allocated(error)) return
    end do
  end subroutine test_Pdbpy

  subroutine test_PdPMe3(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(3, 14), radii(14), alpha, axis(3), ref_axis(3)
    integer :: tangent_atoms(3), ref_tangent_atoms(3), stat, i
    character(:), allocatable :: errmsg

    ! Run cone angle calculation
    coordinates = reshape( &
                  [-7.433970e-01_wp, -9.920000e-04_wp, -1.518770e+00_wp, 1.005896e+00_wp, &
                   2.102000e-03_wp, -1.636180e-01_wp, 2.689042e+00_wp, -1.142000e-03_wp, &
                   -9.464240e-01_wp, 2.797104e+00_wp, -8.888730e-01_wp, -1.576886e+00_wp, &
                   3.486761e+00_wp, 3.481000e-03_wp, -1.921520e-01_wp, 2.795646e+00_wp, &
                   8.800330e-01_wp, -1.586258e+00_wp, 1.164898e+00_wp, 1.431506e+00_wp, &
                   1.010324e+00_wp, 1.222204e+00_wp, 2.365419e+00_wp, 4.431760e-01_wp, &
                   2.057780e+00_wp, 1.338752e+00_wp, 1.642300e+00_wp, 2.773700e-01_wp, &
                   1.479634e+00_wp, 1.648431e+00_wp, 1.165320e+00_wp, -1.429939e+00_wp, &
                   1.006886e+00_wp, 1.222850e+00_wp, -2.362601e+00_wp, 4.376840e-01_wp, &
                   2.776990e-01_wp, -1.479572e+00_wp, 1.644733e+00_wp, 2.058112e+00_wp, &
                   -1.338384e+00_wp, 1.639153e+00_wp], [3, 14])
    radii = [2.1_wp, 1.8_wp, 1.7_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.7_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.7_wp, 1.1_wp, 1.1_wp, 1.1_wp]
    call cone_angle(coordinates, radii, 1, alpha, axis, tangent_atoms, stat, errmsg)

    ! Check cone angle
    call check(error, alpha, 117.110_wp, thr=0.001_wp)
    if (allocated(error)) return

    ! Check cone indices
    ref_tangent_atoms = [6, 10, 13]
    call sort(ref_tangent_atoms)
    call sort(tangent_atoms)
    do i = 1, 3
      call check(error, tangent_atoms(i), ref_tangent_atoms(i))
      if (allocated(error)) return
    end do

    ! Check cone axis
    ref_axis = [7.89524317e-01_wp, 5.22544038e-05_wp, 6.13719276e-01_wp]
    do i = 1, 3
      call check(error, axis(i), ref_axis(i), thr=0.0001_wp)
      if (allocated(error)) return
    end do

  end subroutine test_PdPMe3
end module test_integration
