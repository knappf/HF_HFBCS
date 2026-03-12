subroutine read_input
   
  use declarations
  !implicit none
  ! Subroutine to read input parameters for HF-BCS calculations

  write(*,*)"Reading input parameters..."

  open(unit=1, file='input_HFBCS.dat', status='old', action='read', form='formatted')

  read(1,*) mass_number, proton_number
  read(1,*) Nmax1_2N_file,Nmax12_2N_file
  read(1,*) Nmax1_3N_file,Nmax12_3N_file,Nmax123_3N_file
  read(1,*) Nmax1,Nmax12,Nmax123
  read(1,*) hbar_omega
  read(1,*) epsilon

  close(1)

  neutron_number = mass_number - proton_number

  write(*,*) "Input parameters read:"
  write(*,*) "Mass number (A): ", mass_number
  write(*,*) "Proton number (Z): ", proton_number
  write(*,*) "Model space 2N_file: ", Nmax1_2N_file, Nmax12_2N_file
  write(*,*) "Model space 3N_file: ", Nmax1_3N_file, Nmax12_3N_file, Nmax123_3N_file
  write(*,*) "Model space used in calculation: ", Nmax1, Nmax12, Nmax123
  write(*,*) "Harmonic oscillator frequency (hbar_omega): ", hbar_omega
  write(*,*) "Convergence threshold (epsilon): ", epsilon

  if(mod(proton_number,2) == 1) then
    write(*,*) 'Error! Wrong input Z - Z must be even.'
    stop
  endif

  if(mod(neutron_number,2) == 1) then
    write(*,*) 'Error! Wrong input N - N must be even.'
    stop
  endif

! temporarily set BCS flags to false
  bcs_n=.true.
  bcs_p=.true.

  if(proton_number.eq.2.or.proton_number.eq.8.or.proton_number.eq.20.or.proton_number.eq.28.or.proton_number.eq.50.or.proton_number.eq.82.or.proton_number.eq.126) bcs_p = .false.
  if(neutron_number.eq.2.or.neutron_number.eq.8.or.neutron_number.eq.20.or.neutron_number.eq.28.or.neutron_number.eq.50.or.neutron_number.eq.82.or.neutron_number.eq.126) bcs_n = .false.

return
end subroutine read_input