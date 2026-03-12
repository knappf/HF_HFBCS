module type_defs
  implicit none

TYPE level_type
!      SEQUENCE
    integer           :: index       ! the number of level
    integer           :: par         ! parity
	integer           :: Nosc        ! the quantum number N=2*n+l
	integer           :: nr          ! the quantum number n
	integer           :: l           ! the quantum number l
    integer           :: j2          ! the quantum number 2*j
    real(kind=8)      :: e_ho        ! s.p.-energy e_i = hbar*omega*(N+3/2)
	real(kind=8)      :: e_hf        ! s.p. HF energy
    real(kind=8)      :: qe          ! quasi-s.p. HFB energy
    real(kind=8)      :: u           ! Bogolyubov coefficient u_i 
    real(kind=8)      :: v           ! Bogolyubov coefficient v_i 
END TYPE level_type

TYPE level3b_type
!      SEQUENCE
        INTEGER           :: index        ! the number of level
        INTEGER           :: i            ! the i-th state
        INTEGER           :: j            ! the j-th state
        INTEGER           :: k            ! the k-th state
        INTEGER           :: Jab          ! the angular momentum [ji x jj]Jab
        INTEGER           :: JJ           ! 2-times the total angular momentum 
        INTEGER           :: Tab          ! the isospin [ti x tj]Tab
        INTEGER           :: TT           ! 2-times the total isospin
        INTEGER           :: ord          ! position in J ordered basis
END TYPE level3b_type

end module type_defs