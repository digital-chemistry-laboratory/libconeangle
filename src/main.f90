!> Main cone angle function
module coneangle_main
  use stdlib_array, only: trueloc
  use coneangle_env, only: stderr, wp
  use coneangle_cones, only: prune_cones, search_1_cones, search_2_cones, search_3_cones, test_inside
  use conenagle_data, only: RAD_TO_DEG
  use coneangle_util, only: get_distances
  implicit none(type, external)

  private
  public :: cone_angle

contains
  subroutine cone_angle(coordinates, radii, index_metal, alpha, axis, tangent_atoms, stat, errmsg)
    !! Calculate cone angle, cone axis and tangent atoms
    !> Coordinates (Å)
    real(wp), intent(in) :: coordinates(:, :)
    !> vdW radii (Å)
    real(wp), intent(in) :: radii(:)
    !> Index of metal atom
    integer, intent(in) :: index_metal
    !> Cone angle (degrees)
    real(wp), intent(out) :: alpha
    !> Cone axis (Å)
    real(wp), intent(out) :: axis(3)
    !> Indices of atoms tangent to cone
    integer, intent(out) :: tangent_atoms(3)
    !> Return code
    integer, intent(out) :: stat
    !> Error message
    character(:), allocatable, intent(out) :: errmsg

    real(wp) :: distances(size(radii)), coordinates_centered(3, size(radii)), &
                axes(3, size(radii)), alphas(size(radii)), center(3)
    integer, allocatable :: indices(:)
    integer :: n_atoms, indices_ligand(size(radii) - 1), atom_indices(size(radii)), i
    logical :: mask(size(radii)), mask_ligand(size(radii))
    logical, allocatable :: is_inside(:)

    ! Validate input
    if (size(coordinates, dim=2) /= size(radii)) then
      errmsg = "Mismatch in dimension between coordinates and radii."
      stat = 1
      return
    end if

    if (index_metal < lbound(coordinates, dim=2) .or. index_metal > ubound(coordinates, dim=2)) then
      errmsg = "Metal index out of bounds."
      stat = 1
      return
    end if

    ! Set up masks
    n_atoms = size(radii)
    mask = .true.
    mask(index_metal) = .false.
    mask_ligand = mask

    indices = trueloc(mask)
    indices_ligand = indices

    ! Check whether atoms are within vdW radius of metal atom
    distances = get_distances(coordinates(:, index_metal), coordinates) - radii

    if (any(distances(indices) < 0)) then
      errmsg = "Atoms within vdW radius of metal atom."
      stat = 1
      return
    end if

    ! Center coordinate system on metal atom
    center = coordinates(:, index_metal)
    coordinates_centered = coordinates - spread(center, 2, size(coordinates, 2))

    ! Calculate cone axes and angles
    distances = norm2(coordinates_centered, 1)
    axes(:, indices) = coordinates_centered(:, indices)/spread(distances(indices), 1, 3)
    alphas(indices) = asin(radii(indices)/distances(indices))

    ! Search over one atom cones
    call search_1_cones(alphas(indices), axes(:, indices), alpha, axis, tangent_atoms)
    is_inside = test_inside(alpha, axis, alphas(indices), axes(:, indices))

    if (.not. all(is_inside)) then
      ! Prune out atoms inside other atoms cone
      is_inside = prune_cones(alphas(indices), axes(:, indices))
      mask(indices) = mask(indices) .and. .not. is_inside
      indices = trueloc(mask)

      ! Search over two atom cones
      call search_2_cones(alphas(indices), axes(:, indices), alpha, axis, tangent_atoms)
      is_inside = test_inside(alpha, axis, alphas(indices), axes(:, indices))
    end if

    if (.not. all(is_inside)) then
      ! Search over three atom cones
      call search_3_cones(alphas(indices), axes(:, indices), &
                          coordinates_centered(:, indices_ligand), alpha, axis, tangent_atoms, stat, errmsg)
      if (stat /= 0) then
        return
      end if
    end if

    ! Map cone indices onto original atom indices
    atom_indices = [(i, i=1, size(radii))]
    associate (temp_indices => atom_indices(indices))
      where (tangent_atoms /= 0)
        tangent_atoms = temp_indices(tangent_atoms)
      end where
    end associate

    ! Convert cone angle to degrees
    alpha = 2*alpha*RAD_TO_DEG

    stat = 0
  end subroutine cone_angle
end module coneangle_main
