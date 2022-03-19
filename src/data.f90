!> Data
module conenagle_data
  use coneangle_env, only: wp
  implicit none(type, external)

  private
  public :: RAD_TO_DEG, PI

  !> Pi
  real(wp), parameter :: PI = acos(-1.0_wp)
  !> Radians to degrees conversion factor
  real(wp), parameter :: RAD_TO_DEG = 180._wp/PI
end module conenagle_data
