module test_suite1
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  integer, parameter :: qp = selected_real_kind(33)
  integer, parameter :: dp = selected_real_kind(15)
  
  private

  public :: collect_suite1

contains

!> Collect all exported unit tests
subroutine collect_suite1(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("shear modulus", test_shear_modulus), &
    new_unittest("beta star binary", test_beta_star_binary), &
    new_unittest("beta star ternary", test_beta_star_ternary) &
    ]

end subroutine collect_suite1


subroutine test_shear_modulus(error)
  use KWN_model_functions, only: calculate_shear_modulus
  type(error_type), allocatable, intent(out) :: error
  
  real(qp) :: val
  real(qp) :: expected
  
  val = calculate_shear_modulus(800.0_qp)

  expected = 19311670851.89_qp
  call check(error, val, expected, rel=.true., thr=0.001_qp)
    
end subroutine test_shear_modulus



subroutine test_beta_star_binary(error)
  use KWN_model_functions, only: calculate_beta_star
  use KWN_data_types, only: tKwnpowerlawMicrostructure

  type(error_type), allocatable, intent(out) :: error

  real(qp) :: val
  real(qp) :: expected

  real(qp) :: radius_crit = 1d-9
  real(qp) :: lattice_param = 3

  type(tKwnpowerlawMicrostructure) :: dst

  integer :: en = 1

   allocate(dst%c_matrix                 (2,1), source=0.0_qp)
   allocate(dst%diffusion_coefficient    (2,1), source=0.0_qp)

  dst%diffusion_coefficient(1,:) = 2
  dst%diffusion_coefficient(2,:) = 2
  dst%c_matrix(1,:) = 1
  dst%c_matrix(2,:) = 0
  
  val = calculate_beta_star(radius_crit, lattice_param,  en, dst)
  
  expected = 0.310280E-18
  call check(error, val, expected, rel=.true., thr=0.001_qp)

end subroutine test_beta_star_binary


subroutine test_beta_star_ternary(error)
  use KWN_model_functions, only: calculate_beta_star
  use KWN_data_types, only: tKwnpowerlawMicrostructure


  type(error_type), allocatable, intent(out) :: error

  real(qp) :: val
  real(qp) :: expected
  
  real(qp) :: radius_crit = 1d-9
  real(qp) :: lattice_param = 3

  type(tKwnpowerlawMicrostructure) :: dst

  integer :: en = 1

  allocate(dst%c_matrix                 (2,1), source=0.0_qp)
  allocate(dst%diffusion_coefficient    (2,1), source=0.0_qp)

  dst%diffusion_coefficient(1,:) = 2
  dst%diffusion_coefficient(2,:) = 2
  dst%c_matrix(1,:) = 1
  dst%c_matrix(2,:) = 1
  
  val = calculate_beta_star(radius_crit, lattice_param,  en, dst)
  
  expected = 0.155140E-18
  call check(error, val, expected, rel=.true., thr=0.001_qp)

end subroutine test_beta_star_ternary



end module test_suite1
