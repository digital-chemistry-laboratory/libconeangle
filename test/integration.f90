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
                new_unittest("PdCO_bounds", test_PdCO_bounds, should_fail=.true.), &
                new_unittest("No_cone", test_no_cone, should_fail=.true.) &
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
    call check(error, stat, 0, message=errmsg)
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
    call check(error, stat, 0, message=errmsg)
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
    call check(error, stat, 0, message=errmsg)
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

  subroutine test_no_cone(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(3, 114), radii(114), alpha, axis(3)
    integer :: tangent_atoms(3), stat
    character(:), allocatable :: errmsg

    ! Run cone angle calculation
    coordinates = reshape( &
                  [-0.150731_wp, -0.19074103_wp, -1.92057203_wp, &
                   3.37788306_wp, 0.05130647_wp, 3.61481396_wp, &
                   2.54925444_wp, 0.10907458_wp, 2.4732705_wp, &
                   1.16592121_wp, 0.19499813_wp, 2.7295883_wp, &
                   0.67750873_wp, 0.32027115_wp, 4.03553064_wp, &
                   -0.66459397_wp, 0.33472346_wp, 4.29526449_wp, &
                   -1.5896386_wp, 0.22951037_wp, 3.26335973_wp, &
                   -2.97600312_wp, 0.16907618_wp, 3.47545234_wp, &
                   -3.37273754_wp, 0.14349526_wp, 4.82745361_wp, &
                   -1.12082686_wp, 0.16072819_wp, 1.92806203_wp, &
                   0.25612043_wp, 0.13107103_wp, 1.67287536_wp, &
                   1.02958224_wp, 0.02476611_wp, 0.01490941_wp, &
                   1.82059756_wp, 1.7253701_wp, -0.29829352_wp, &
                   0.70905762_wp, 2.73064994_wp, 0.07168378_wp, &
                   1.17203462_wp, 4.18917778_wp, -0.1847218_wp, &
                   2.40741965_wp, 4.4864593_wp, 0.70070731_wp, &
                   3.54375578_wp, 3.49571492_wp, 0.34372851_wp, &
                   3.93048362_wp, 3.66565495_wp, -1.14642971_wp, &
                   2.69574266_wp, 3.37243964_wp, -2.03329984_wp, &
                   1.55037404_wp, 4.35567902_wp, -1.67841579_wp, &
                   2.20574022_wp, 1.92384277_wp, -1.78094054_wp, &
                   3.05713071_wp, 2.04198933_wp, 0.57964817_wp, &
                   2.15728409_wp, -1.50864448_wp, -0.07241542_wp, &
                   1.99383629_wp, -2.27757743_wp, -1.40472633_wp, &
                   2.86192227_wp, -3.56211624_wp, -1.41764689_wp, &
                   4.35051981_wp, -3.16593263_wp, -1.25560966_wp, &
                   4.53490985_wp, -2.40802958_wp, 0.08272307_wp, &
                   4.11213979_wp, -3.32596808_wp, 1.25833128_wp, &
                   2.6248743_wp, -3.72959249_wp, 1.09122431_wp, &
                   2.43708516_wp, -4.48368875_wp, -0.24816519_wp, &
                   1.7486867_wp, -2.45094363_wp, 1.08450343_wp, &
                   3.64796946_wp, -1.13454067_wp, 0.07974954_wp, &
                   -2.15873843_wp, 0.02557711_wp, 0.87494921_wp, &
                   -2.83375913_wp, 1.15424363_wp, 0.36371904_wp, &
                   -3.79170515_wp, 0.96990112_wp, -0.61710415_wp, &
                   -4.13561567_wp, -0.28988183_wp, -1.07699306_wp, &
                   -3.55863471_wp, -1.3919975_wp, -0.47177744_wp, &
                   -2.58828911_wp, -1.26278612_wp, 0.50348897_wp, &
                   -2.1100334_wp, -2.48589077_wp, 1.26782191_wp, &
                   -3.1332671_wp, -2.87618502_wp, 2.33469829_wp, &
                   -1.78019071_wp, -3.675817_wp, 0.36281612_wp, &
                   -5.05101947_wp, -0.44670698_wp, -2.27345372_wp, &
                   -4.23171871_wp, -0.79345761_wp, -3.5188832_wp, &
                   -6.16277239_wp, -1.46827493_wp, -2.02460123_wp, &
                   -2.6389261_wp, 2.53465503_wp, 0.97093586_wp, &
                   -2.6731693_wp, 3.66659474_wp, -0.06385225_wp, &
                   -3.67816057_wp, 2.79319683_wp, 2.06269304_wp, &
                   3.09870866_wp, -0.79978635_wp, 4.23676883_wp, &
                   3.30075044_wp, 0.98089069_wp, 4.17942435_wp, &
                   4.41250643_wp, -0.07494945_wp, 3.29490855_wp, &
                   1.36612294_wp, 0.39439745_wp, 4.85692446_wp, &
                   -1.00413428_wp, 0.41590767_wp, 5.31160633_wp, &
                   -4.4574898_wp, 0.04228008_wp, 4.87101366_wp, &
                   -3.08700057_wp, 1.07437316_wp, 5.31820892_wp, &
                   -2.91961721_wp, -0.70934337_wp, 5.33392628_wp, &
                   0.45320601_wp, 2.60280398_wp, 1.12970462_wp, &
                   -0.19085692_wp, 2.5009366_wp, -0.50473306_wp, &
                   0.35492371_wp, 4.88682367_wp, 0.07034344_wp, &
                   2.14421904_wp, 4.3943906_wp, 1.76081388_wp, &
                   2.74277135_wp, 5.51883252_wp, 0.54920794_wp, &
                   4.42373568_wp, 3.69813931_wp, 0.97994855_wp, &
                   4.75063997_wp, 2.98584776_wp, -1.40467617_wp, &
                   4.29620677_wp, 4.68241516_wp, -1.32997566_wp, &
                   2.96839395_wp, 3.48856684_wp, -3.09739297_wp, &
                   0.67652035_wp, 4.15885389_wp, -2.31075326_wp, &
                   1.85887259_wp, 5.38767303_wp, -1.88008983_wp, &
                   2.96534664_wp, 1.19348819_wp, -2.07803376_wp, &
                   1.33404653_wp, 1.74391149_wp, -2.4132579_wp, &
                   2.80000551_wp, 1.90433112_wp, 1.63257294_wp, &
                   3.86967954_wp, 1.34944608_wp, 0.35361523_wp, &
                   2.30473569_wp, -1.63346618_wp, -2.23048312_wp, &
                   0.943587_wp, -2.53527988_wp, -1.56287883_wp, &
                   2.72089862_wp, -4.09188083_wp, -2.37664043_wp, &
                   4.66676824_wp, -2.53530469_wp, -2.09457136_wp, &
                   4.98553074_wp, -4.058925_wp, -1.27717134_wp, &
                   5.59425957_wp, -2.11901399_wp, 0.20154492_wp, &
                   4.25663636_wp, -2.80691628_wp, 2.21277006_wp, &
                   4.74609151_wp, -4.2194651_wp, 1.28838594_wp, &
                   2.32101659_wp, -4.38137367_wp, 1.9295423_wp, &
                   1.39110083_wp, -4.78849902_wp, -0.36828869_wp, &
                   3.03412685_wp, -5.40277855_wp, -0.25327918_wp, &
                   1.84793775_wp, -1.92904579_wp, 2.04037241_wp, &
                   0.69355273_wp, -2.71517113_wp, 0.97742744_wp, &
                   3.94364514_wp, -0.47412877_wp, -0.74156635_wp, &
                   3.79104303_wp, -0.58608794_wp, 1.01295901_wp, &
                   -4.2742406_wp, 1.82643368_wp, -1.05108578_wp, &
                   -3.84581055_wp, -2.37458531_wp, -0.79712727_wp, &
                   -1.18368241_wp, -2.19257929_wp, 1.78763211_wp, &
                   -3.3536679_wp, -2.03497956_wp, 2.99220055_wp, &
                   -2.77928716_wp, -3.70502614_wp, 2.94884454_wp, &
                   -4.07120364_wp, -3.18222374_wp, 1.87169859_wp, &
                   -1.32736176_wp, -4.49395165_wp, 0.92444772_wp, &
                   -1.0952821_wp, -3.39356505_wp, -0.43550308_wp, &
                   -2.6789532_wp, -4.06780733_wp, -0.11110905_wp, &
                   -5.52924589_wp, 0.53389222_wp, -2.44222512_wp, &
                   -4.86089039_wp, -0.88510392_wp, -4.40472287_wp, &
                   -3.70765326_wp, -1.73971975_wp, -3.3807132_wp, &
                   -3.47439401_wp, -0.03299876_wp, -3.71557572_wp, &
                   -6.72153597_wp, -1.23264259_wp, -1.11797902_wp, &
                   -5.75277225_wp, -2.47091635_wp, -1.90793774_wp, &
                   -6.86768348_wp, -1.50183613_wp, -2.85638946_wp, &
                   -1.64495071_wp, 2.53371097_wp, 1.44508932_wp, &
                   -3.68359058_wp, 3.8164286_wp, -0.4422523_wp, &
                   -2.03981694_wp, 3.4530384_wp, -0.92276893_wp, &
                   -2.34994421_wp, 4.61433395_wp, 0.36862794_wp, &
                   -4.68247165_wp, 2.80296196_wp, 1.63957451_wp, &
                   -3.51501765_wp, 3.75298046_wp, 2.55428919_wp, &
                   -3.65502807_wp, 2.01491302_wp, 2.8234341_wp, &
                   -1.09122209_wp, 1.25105272_wp, -2.2802411_wp, &
                   0.94507581_wp, -0.26624979_wp, -3.32603572_wp, &
                   -1.11868769_wp, -1.6330022_wp, -2.21693533_wp, &
                   -1.73355125_wp, -2.57286054_wp, -2.41694099_wp, &
                   -1.68224675_wp, 2.19788714_wp, -2.51393017_wp, &
                   2.00766997_wp, -0.2091706_wp, -3.73779609_wp], &
                  [3, 114])
    radii = [1.97_wp, 1.7_wp, 1.52_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.52_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.8_wp, 1.7_wp, 1.7_wp, &
             1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, &
             1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, &
             1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, &
             1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, &
             1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, &
             1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, &
             1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.1_wp, 1.7_wp, 1.7_wp, 1.7_wp, 1.52_wp, &
             1.52_wp, 1.52_wp]

    call cone_angle(coordinates, radii, 1, alpha, axis, tangent_atoms, stat, errmsg)

    ! Check stat calculation failed.
    call check(error, stat, 0, message=errmsg)
    if (allocated(error)) return

  end subroutine test_no_cone
end module test_integration
