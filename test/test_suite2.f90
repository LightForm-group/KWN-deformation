module test_suite2
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  integer, parameter :: qp = selected_real_kind(33)
  integer, parameter :: dp = selected_real_kind(15)
  
  private

  public :: collect_suite2

contains

!> Collect all exported unit tests
subroutine collect_suite2(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("shear modulus", test_shear_modulus, should_fail=.true.), &
    new_unittest("beta star", test_beta_star, should_fail=.true.) &
    ]

end subroutine collect_suite2


subroutine test_shear_modulus(error)
  use KWN_model_functions, only: calculate_shear_modulus
  type(error_type), allocatable, intent(out) :: error
  
  real(qp) :: val
  real(qp) :: expected
  
  val = calculate_shear_modulus(800.0_qp)

  expected = 3.3
  !expected = 19311670851.89_qp
  call check(error, val, expected, rel=.true., thr=0.001_qp)
    
end subroutine test_shear_modulus


subroutine test_beta_star(error)
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
  
  val = calculate_beta_star(radius_crit, lattice_param, en, dst)
  
  expected = 33.0
  call check(error, val, expected, rel=.true.)

end subroutine test_beta_star

end module test_suite2