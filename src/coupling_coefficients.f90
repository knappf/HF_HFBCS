module coupling_coeffs

use fwigxjpf
implicit none

contains

real(kind=8) function sixj_real(a,b,c,d,e,f) 

! 6j coefficient, double precision arguments

real(kind=8):: a,b,e,d,c,f, value

value = fwig6jj(int(2*a),int(2*b),int(2*c),int(2*d),int(2*e),int(2*f))
sixj_real = value

end function sixj_real
! -------------------------------------------------------------------------
real(kind=8) function  delta1(i,j)
  
    integer :: i,j
  if(i == j) then
    delta1=1.d0
        else
    delta1=0.d0
  endif
end function delta1
! -------------------------------------------------------------------------


end module coupling_coeffs