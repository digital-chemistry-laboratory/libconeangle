!> Utility and math functions
module coneangle_util
  use coneangle_env, only: wp
  implicit none(type, external)

  private
  public :: binomial, cross_product, get_distances, solve_quadratic, vector_angle

contains
  pure function binomial(n, k)
    !! Returns binomial coefficient \( \binom{n}{k} \)
    !> n
    integer, intent(in) :: n
    !> k
    integer, intent(in) :: k
    !> Binomial coefficient
    integer :: binomial

    binomial = nint(exp(log_gamma(n + 1.0_wp) - log_gamma(n - k + 1.0_wp) - log_gamma(k + 1.0_wp)))
  end function binomial

  pure function cross_product(a, b)
    !! Returns cross product between two vectors  \( a \times b \)
    !> First vector
    real(wp), intent(in) :: a(3)
    !> Second vector
    real(wp), intent(in) :: b(3)
    !> Cross product
    real(wp) :: cross_product(3)

    cross_product(1) = a(2)*b(3) - b(2)*a(3)
    cross_product(2) = a(3)*b(1) - b(3)*a(1)
    cross_product(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product

  pure function get_distances(a, b) result(distances)
    !! Calculate distance between point and group of points
    !> Point
    real(wp), intent(in) :: a(:)
    !> Group of points
    real(wp), intent(in) :: b(:, :)
    !> Distances
    real(wp) :: distances(size(b, 2))

    distances = norm2(spread(a, 2, size(b, 2)) - b, 1)
  end function get_distances

  pure function solve_quadratic(a, b, c) result(roots)
    !! Return roots for quadratic equation  \( ax^2 + bx + c = 0 \)
    !> a*x**2
    real(wp), intent(in) :: a
    !> b*x
    real(wp), intent(in) :: b
    !> c
    real(wp), intent(in) :: c
    !> Roots
    complex(wp) :: roots(2)

    associate (delta => sqrt(cmplx(b**2 - 4*a*c, kind=wp)))
      roots(1) = (-b + delta)/(2*a)
      roots(2) = (-b - delta)/(2*a)
    end associate
  end function solve_quadratic

  pure function vector_angle(a, b) result(angle)
  !! Calculate angle between vectors
    !> First vector
    real(wp), intent(in) :: a(3)
    !> Second vector
    real(wp), intent(in) :: b(3)
    !> Angle (radians)
    real(wp) :: angle

    angle = atan2(norm2(cross_product(a, b)), dot_product(a, b))
  end function vector_angle
end module coneangle_util
