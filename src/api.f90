!> C API
module coneangle_api
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use coneangle_main, only: cone_angle
  implicit none(type, external)

contains
  subroutine cone_angle_c(n_atoms, coordinates, radii, index_metal, alpha, axis, tangent_atoms, stat) bind(c, name="cone_angle")
    !! Calculate cone angle, cone axis and tangent atoms
    !> Number of atoms
    integer(c_int), value, intent(in) :: n_atoms
    !> Coordinates (Å)
    real(c_double), intent(in) :: coordinates(3, n_atoms)
    !> vdW radii (Å)
    real(c_double), intent(in) :: radii(n_atoms)
    !> Index of metal atom
    integer(c_int), value, intent(in) :: index_metal
    !> Cone angle (degrees)
    real(c_double), intent(out) :: alpha
    !> Cone axis (Å)
    real(c_double), intent(out) :: axis(3)
    !> Indices of atoms tangent to cone
    integer(c_int), intent(out) :: tangent_atoms(3)
    !> Return code
    integer(c_int), intent(out) :: stat
    !> Error message
    character(:), allocatable :: errmsg

    call cone_angle(coordinates, radii, index_metal, alpha, axis, tangent_atoms, stat, errmsg)
  end subroutine cone_angle_c
end module coneangle_api
