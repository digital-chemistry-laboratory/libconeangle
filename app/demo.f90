program demo
  use coneangle_main, only: cone_angle

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp) :: coordinates(3, 3), radii(3), alpha, axis(3)
  integer :: tangent_atoms(3), stat
  character(:), allocatable :: errmsg

  coordinates = reshape([0._dp, 0._dp, -0.52_dp, 0._dp, 0._dp, 1.76_dp, 0._dp, 0._dp, 2.86_dp], [3, 3])
  radii = [2.1_dp, 1.7_dp, 1.52_dp]
  call cone_angle(coordinates, radii, 1, alpha, axis, tangent_atoms, stat, errmsg)
  write (*, *) "Cone angle:", alpha
  write (*, *) "Cone axis:", axis
  write (*, *) "Tangent atoms:", tangent_atoms
end program demo
