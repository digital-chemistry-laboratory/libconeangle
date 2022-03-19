!> Subroutines and functions to search for cones
module coneangle_cones
  use stdlib_array, only: trueloc
  use stdlib_error, only: check
  use stdlib_math, only: is_close
  use stdlib_sorting, only: sort_index
  use coneangle_env, only: int64, stderr, wp
  use conenagle_data, only: PI
  use coneangle_util, only: binomial, cross_product, solve_quadratic, vector_angle
  implicit none(type, external)

  private
  public :: search_1_cones, search_2_cones, search_3_cones, test_inside, prune_cones

contains
  subroutine search_1_cones(alphas, axes, alpha, axis, indices)
    !! Search for cones tangent to one atom
    !> Cone angles
    real(wp), intent(in) :: alphas(:)
    !> Cone axes
    real(wp), intent(in) :: axes(:, :)
    !> Angle of found cone
    real(wp), intent(out) :: alpha
    !> Angle of found cone
    real(wp), intent(out) :: axis(3)
    !> Indices of tangent atoms of found cone
    integer, intent(out) :: indices(3)

    integer :: idx

    ! Take out largest cone
    idx = maxloc(alphas, 1)
    alpha = alphas(idx)
    axis = axes(:, idx)
    indices = [idx, 0, 0]
  end subroutine search_1_cones

  subroutine search_2_cones(alphas, axes, alpha, axis, indices)
    !! Search for cones tangent to two atoms
    !> Cone angles
    real(wp), intent(in) :: alphas(:)
    !> Cone axes
    real(wp), intent(in) :: axes(:, :)
    !> Angle of found cone
    real(wp), intent(out) :: alpha
    !> Axis of found cone
    real(wp), intent(out) :: axis(3)
    !> Indices of tangent atoms of found cone
    integer, intent(out) :: indices(3)

    real(wp) :: a_ij, b_ij, c_ij, beta_i, beta_j, beta_ij, alpha_ij, axis_ij(3)
    integer :: i, j

    alpha = 0
    do i = 1, size(alphas)
      do j = i + 1, size(alphas)
        ! Calculate angle
        beta_i = alphas(i)
        beta_j = alphas(j)
        beta_ij = acos(dot_product(axes(:, i), axes(:, j)))
        alpha_ij = (beta_ij + beta_i + beta_j)/2
        if (alpha_ij < alpha) cycle

        ! Calculate axis vector
        a_ij = (1/sin(beta_ij))*sin(0.5*(beta_ij + beta_i - beta_j))
        b_ij = (1/sin(beta_ij))*sin(0.5*(beta_ij - beta_i + beta_j))
        c_ij = 0
        axis_ij = a_ij*axes(:, i) + b_ij*axes(:, j) + c_ij
        axis_ij = axis_ij/spread(norm2(axis_ij, 1), 1, size(axis_ij))
        alpha = alpha_ij
        axis = axis_ij
        indices = [i, j, 0]
      end do
    end do

  end subroutine search_2_cones

  subroutine search_3_cones(alphas, axes, coordinates, alpha, axis, indices, stat, errmsg)
    !! Search for cones tangent to three atoms
    !> Cone angles
    real(wp), intent(in) :: alphas(:)
    !> Cone axes
    real(wp), intent(in) :: axes(:, :)
    !> Coordinates of ligand
    real(wp), intent(in) :: coordinates(:, :)
    !> Angle of found cone
    real(wp), intent(inout) :: alpha
    !> Axis of found cone
    real(wp), intent(out) :: axis(3)
    !> Indices of tangent atoms of found cone
    integer, intent(out) :: indices(3)
    !> Return code
    integer, intent(out) :: stat
    !> Error message
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, j, k, l, sign_p, n_combinations, i_all, index_min
    integer(int64) :: sort_idx(4)
    real(wp) :: beta_i, beta_j, beta_k, beta_ij, gamma, &
                m_i(3), m_j(3), m_k(3), u(3), v(3), N(3, 3), P(3, 3), A, B, C, D, &
                p_2, p_1, p_0, cos_roots(4), angles(4), D_tests(4), &
                a_ij, b_ij, c_ij, physical_angles(2), p_(3), axis_(3), upper_bound, &
                lower_bound
    real(wp), allocatable :: axes_all(:, :), alphas_all(:)
    integer, allocatable :: indices_all(:, :)
    logical, allocatable :: keep_all(:)
    complex(wp) :: roots(2)

    n_combinations = binomial(size(alphas), 3)
    allocate (axes_all(3, n_combinations*2))
    allocate (alphas_all(n_combinations*2))
    allocate (indices_all(3, n_combinations*2))
    ! Get three atom cones
    i_all = 1
    do i = 1, size(alphas)
      do j = i + 1, size(alphas)
        do k = j + 1, size(alphas)
          ! Set up angles and axes
          beta_i = alphas(i)
          beta_j = alphas(j)
          beta_k = alphas(k)

          m_i = axes(:, i)
          m_j = axes(:, j)
          m_k = axes(:, k)

          ! Set up angles between atom vectors
          beta_ij = acos(dot_product(m_i, m_j))

          ! Set up matrices
          u = [cos(beta_i), cos(beta_j), cos(beta_k)]
          v = [sin(beta_i), sin(beta_j), sin(beta_k)]
          N(1, :) = cross_product(m_j, m_k)
          N(2, :) = cross_product(m_k, m_i)
          N(3, :) = cross_product(m_i, m_j)
          P = matmul(N, transpose(N))
          gamma = dot_product(m_i, cross_product(m_j, m_k))

          ! Set up coefficients of quadratic equation
          associate (A_temp => matmul(matmul(u, P), reshape(u, [3, 1])), &
                     B_temp => matmul(matmul(v, P), reshape(v, [3, 1])), &
                     C_temp => matmul(matmul(u, P), reshape(v, [3, 1])))
            A = A_temp(1)
            B = B_temp(1)
            C = C_temp(1)
          end associate
          D = gamma**2

          ! Solve quadratic equation
          p_2 = (A - B)**2 + 4*C**2
          p_1 = 2*(A - B)*(A + B - 2*D)
          p_0 = (A + B - 2*D)**2 - 4*C**2

          roots = solve_quadratic(p_2, p_1, p_0)
          if (.not. any(is_close(abs(aimag(roots)), 0._wp))) then
            errmsg = "Complex roots encountered."
            stat = 1
            return
          end if

          cos_roots([1, 3]) = acos(real(roots))
          cos_roots([2, 4]) = 2*PI - cos_roots([1, 3])

          ! Calculate physical angles
          angles = cos_roots/2
          D_tests = abs(A*cos(angles)**2 + B*sin(angles)**2 + 2*C*sin(angles)*cos(angles) - D)
          call sort_index(D_tests, sort_idx)
          physical_angles = angles(sort_idx(:2))

          ! Calculate axis vectors
          associate (pa => physical_angles)
            do l = 1, size(pa)
              a_ij = (cos(pa(l) - beta_i) - cos(pa(l) - beta_j)*cos(beta_ij))/sin(beta_ij)**2
              b_ij = (cos(pa(l) - beta_j) - cos(pa(l) - beta_i)*cos(beta_ij))/sin(beta_ij)**2
              c_ij = sqrt(1 - a_ij**2 - b_ij**2 - 2*a_ij*b_ij*cos(beta_ij))

              p_ = matmul(transpose(N), u*cos(pa(l)) + v*sin(pa(l)))
              sign_p = int(sign(1._wp, gamma)*sign(1._wp, dot_product(p_, cross_product(m_i, m_j))))
              if (sign_p /= int(sign(1._wp, c_ij))) c_ij = -c_ij
              axis_ = a_ij*m_i + b_ij*m_j + c_ij/sin(beta_ij)*cross_product(m_i, m_j)

              alphas_all(i_all) = pa(l)
              axes_all(:, i_all) = axis_
              indices_all(:, i_all) = [i, j, k]
              i_all = i_all + 1
            end do
          end associate
        end do
      end do
    end do
    ! Get upper and lower bound to apex angle
    upper_bound = get_upper_bound(alphas, axes, coordinates)
    lower_bound = alpha

    ! Remove cones from which are outside the bounds
    keep_all = (is_close(alphas_all, upper_bound) .or. alphas_all < upper_bound) &
               .and. (is_close(alphas_all, upper_bound) .or. alphas_all > lower_bound)

    ! Remove cones which don't contain all atoms
    associate (temp_indices => trueloc(keep_all))
      do concurrent(i=1:size(temp_indices))
        j = temp_indices(i)
        keep_all(j) = all(test_inside(alphas_all(j), axes_all(:, j), alphas, axes))
      end do
    end associate

    index_min = minloc(alphas_all, dim=1, mask=keep_all)
    alpha = alphas_all(index_min)
    axis = axes_all(:, index_min)
    indices = indices_all(:, index_min)
    stat = 0
  end subroutine search_3_cones

  pure function get_upper_bound(alphas, axes, coordinates) result(upper_bound)
    !! Calculate upper bound for cone angle tangent to three atoms
    !> Cone angles
    real(wp), intent(in) :: alphas(:)
    !> Coordinates of ligand
    real(wp), intent(in) :: axes(:, :)
    !> Upper bound
    real(wp), intent(in) :: coordinates(:, :)
    !> Cone axes
    real(wp) :: upper_bound

    real(wp) :: centroid_vector(3), vertex_angle, angle_sums(size(alphas))
    integer :: i

    centroid_vector = sum(coordinates, 2)/size(coordinates, 2)
    centroid_vector = centroid_vector/norm2(centroid_vector)

    do concurrent(i=1:size(alphas))
      vertex_angle = acos(dot_product(centroid_vector, axes(:, i)))
      angle_sums(i) = alphas(i) + vertex_angle
    end do
    upper_bound = maxval(angle_sums)
  end function get_upper_bound

  pure function prune_cones(alphas, axes) result(is_inside)
    !! Returns mask of cones that are inside other cones
    !> Cone angles
    real(wp), intent(in) :: alphas(:)
    !> Cone axes
    real(wp), intent(in) :: axes(:, :)
    !> Mask of cones inside another cone
    logical :: is_inside(size(alphas))

    logical :: is_inside_atom(size(alphas))
    integer :: i

    ! Loop over atoms and check whether they are inside
    is_inside = .false.
    do i = 1, size(alphas)
      if (is_inside(i)) cycle
      is_inside_atom = test_inside(alphas(i), axes(:, i), alphas, axes)
      is_inside_atom(i) = .false.
      is_inside = is_inside .or. is_inside_atom
    end do
  end function prune_cones

  pure function test_inside(alpha, axis, alphas, axes) result(is_inside)
    !! Returns mask of cones inside test cone
    !> Cone angle of test cone
    real(wp), intent(in) :: alpha
    !> Cone axis of test cone
    real(wp), intent(in) :: axis(3)
    !> Cone angles
    real(wp), intent(in) :: alphas(:)
    !> Cone axes
    real(wp), intent(in) :: axes(:, :)
    !> Mask of cones inside test cone
    logical :: is_inside(size(alphas))

    integer :: i
    real(wp) :: angle

    do concurrent(i=1:size(alphas))
      if (alpha == alphas(i)) then
        is_inside(i) = .true.
      else
        angle = vector_angle(axis, axes(:, i))
        is_inside(i) = (alpha > (angle + alphas(i))) .or. is_close(alpha, (angle + alphas(i)))
      end if
    end do
  end function test_inside
end module coneangle_cones
