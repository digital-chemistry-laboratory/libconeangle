!> Environment and accuracy
module coneangle_env
  use, intrinsic :: iso_fortran_env, only: int64, stderr => error_unit
  implicit none(type, external)

  private
  public :: int64, stderr, wp

  !> Double precision real numbers
  integer, parameter :: dp = selected_real_kind(15, 307), wp = dp
end module coneangle_env
