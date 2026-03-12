module declarations

  use type_defs
  implicit none
  
!  public :: mass_number, proton_number, neutron_number
!  public :: Nmax1_2N_file, Nmax12_2N_file
!  public :: Nmax1_3N_file, Nmax12_3N_file, Nmax123_3N_file
!  public :: Nmax1, Nmax12, Nmax123
!  public :: bcs_n, bcs_p
!  public :: n_sp_levels
!  public :: hbar_omega
!  public :: epsilon
!  public :: lev_p,lev_n
!  public :: Up_HFB, Vp_HFB, Un_HFB, Vn_HFB
!  public :: rho_p, rho_n, kappa_p, kappa_n
!  public :: T1b_p, T1b_n
  
  
  integer :: mass_number, proton_number, neutron_number
  integer :: Nmax1_2N_file, Nmax12_2N_file
  integer :: Nmax1_3N_file, Nmax12_3N_file, Nmax123_3N_file
  integer :: Nmax1, Nmax12, Nmax123
  logical :: bcs_n, bcs_p, hf_p_diag, hf_n_diag

  integer :: n_sp_levels,n_sp_levels_3N,j2_lev_max, j2_lev3_max, n_sp_levels_calc, j2_lev_max_calc
  real(kind=8) :: epsilon
  real(kind=8) :: hbar_omega

  integer :: iter_max, iter_counter,iter_bcs,max_iter
  real(kind=8) :: E_kin, E_prot, E_neut, E_pn, E_pair, E_HFB, E_prot_2b, E_neut_2b, E_pn_2b
  real(kind=8) :: fermi_energy_p, fermi_energy_n

  type(level_type),dimension(:),allocatable, save :: lev_p,lev_n,lev,levhf_p,levhf_n
  type(level3b_type),dimension(:),allocatable,save :: lev3

  real(kind=8), allocatable, save :: h_p(:,:), h_n(:,:), delta_p(:,:), delta_n(:,:)
  real(kind=8), allocatable, save :: e_p(:), e_n(:)
  
  integer, allocatable, save :: lev3ord(:)
  INTEGER, ALLOCATABLE, SAVE :: lpoint(:,:,:,:,:,:)
  integer (kind=8), allocatable, save :: idim3j(:)
  integer (kind=8), allocatable, save :: idim3jt(:,:)

  real(kind=8), allocatable, SAVE :: Up_HFB(:,:), Vp_HFB(:,:), Un_HFB(:,:), Vn_HFB(:,:)
  real(kind=8), allocatable, SAVE :: rho_p(:,:), rho_n(:,:), kappa_p(:,:), kappa_n(:,:)
  real(kind=8), allocatable, SAVE :: T1b_p(:,:), T1b_n(:,:)

  real(kind=8), ALLOCATABLE, SAVE :: Vpp(:,:,:,:,:)
  real(kind=8), ALLOCATABLE, SAVE :: Vnn(:,:,:,:,:)
  real(kind=8), ALLOCATABLE, SAVE :: Vpn(:,:,:,:,:)

  integer, allocatable, SAVE :: lp1(:), lp2(:)

  real(kind=8), ALLOCATABLE, SAVE :: Vpp_ar(:)
  real(kind=8), ALLOCATABLE, SAVE :: Vnn_ar(:)
  real(kind=8), ALLOCATABLE, SAVE :: Vpn_ar(:)
        
  integer(kind=8), ALLOCATABLE, SAVE :: Vpp_pair(:)
  integer(kind=8), ALLOCATABLE, SAVE :: Vnn_pair(:)
  integer(kind=8), ALLOCATABLE, SAVE :: Vpn_pair(:)

  real(kind=4), ALLOCATABLE, SAVE :: V3B_ar(:)
  integer(kind=8), ALLOCATABLE, SAVE :: V3B_pair(:)        

  REAL(kind=4), ALLOCATABLE, SAVE :: V3B(:,:)
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: V3BNO2(:,:,:,:,:,:,:)
  integer(kind=8) :: ipozjt1,ipozt1,ipoz,ipozjt


end module declarations