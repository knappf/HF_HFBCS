module HFBCS_solver

  use declarations
  use v3body_no2b_bin
  use v2body
  implicit none

contains
    ! Here would go the subroutines and functions for solving the HF-BCS equations

subroutine HFBCS_iteration
   
   real(kind=8), allocatable :: rho0_p(:,:), rho0_n(:,:), kappa0_p(:,:), kappa0_n(:,:)
   real(kind=8) :: Ediff,Ediff_bcs,E_HFB_energy,dens_average_factor
   integer :: iter_hf_initial


! alllocate densities
    allocate(rho_p(n_sp_levels_calc,n_sp_levels_calc), rho_n(n_sp_levels_calc,n_sp_levels_calc))
    allocate(rho0_p(n_sp_levels_calc,n_sp_levels_calc), rho0_n(n_sp_levels_calc,n_sp_levels_calc))
    rho_n=0.d0
    rho_p=0.d0

!  testing neutron BCSs

    dens_average_factor=0.9d0
!    bcs_n = .true.
!    bcs_p = .true.
    hf_n_diag = .true.
    hf_p_diag = .true.

    

    if (bcs_p) then
        write(*,*) 'BCS for protons is ON '
        allocate(kappa_p(n_sp_levels_calc,n_sp_levels_calc))
        allocate(kappa0_p(n_sp_levels_calc,n_sp_levels_calc))
        kappa_p=0.d0
    endif

    if (bcs_n) then
        write(*,*) 'BCS for neutrons is ON '
        allocate(kappa_n(n_sp_levels_calc,n_sp_levels_calc))
        allocate(kappa0_n(n_sp_levels_calc,n_sp_levels_calc))
        kappa_n=0.d0
    endif

!  initial densities
    call densities_update
!  for BCS iteration, calculation of Fermi energies and initial steps must be added here

!    call HFB_energy

!  allocation of h_p, h_n matrices

    allocate(h_p(n_sp_levels_calc,n_sp_levels_calc),h_n(n_sp_levels_calc,n_sp_levels_calc))
    allocate(e_p(n_sp_levels_calc),e_n(n_sp_levels_calc))
    if (bcs_p .or. bcs_n) allocate(delta_p(n_sp_levels_calc,n_sp_levels_calc),delta_n(n_sp_levels_calc,n_sp_levels_calc))

!   call h_fields(h_p,h_n)
!   call d_fields   will be here

!   HF-BCS iterations

    max_iter=500
    iter_counter=0
    iter_hf_initial=4
    Ediff=10000.d0
    Ediff_bcs=10000.d0

    do while ((Ediff > epsilon .and. iter_counter < max_iter).and.(.not.(bcs_p.or.bcs_n)).or.(Ediff_bcs > epsilon .and. (bcs_p.or.bcs_n) .and. iter_counter <= iter_hf_initial))
        
        iter_counter=iter_counter+1
        ! Update fields
        call HFB_energy
        Ediff=dabs(E_HFB_energy - E_HFB)
        E_HFB_energy=E_HFB
        call h_fields(h_p,h_n)

        if ((bcs_p .or. bcs_n ).and. (iter_counter == iter_hf_initial)) then
            iter_bcs=0
!            Ediff_bcs=1d9
            if (bcs_n)hf_n_diag = .false.
            if (bcs_p)hf_p_diag = .false.
            do while (Ediff_bcs > epsilon .and. iter_bcs < max_iter)
                iter_counter=iter_counter+1
                iter_bcs=iter_bcs+1
                rho0_n=rho_n
                rho0_p=rho_p
                kappa0_n=kappa_n
                kappa0_p=kappa_p
                call h_fields(h_p,h_n)
                call d_fields(delta_p,delta_n)
                call diag_hf(h_p,e_p,h_n,e_n)
                call hf_levels(h_p,h_n,e_p,e_n,iter_counter)
                if (bcs_n) then
                    call shift_delta_h(delta_n,h_n,fermi_energy_n,iter_bcs)
                    hf_n_diag = .true.
                    hf_p_diag = .false.
                    call diag_hf(h_p,e_p,h_n,e_n)
                    call hfbcs_levels(h_p,h_n,e_p,e_n,delta_p,delta_n,iter_counter)
                    call densities_update
                    rho_n=dens_average_factor*rho0_n+(1.d0-dens_average_factor)*rho_n
                    rho_p=dens_average_factor*rho0_p+(1.d0-dens_average_factor)*rho_p
                    kappa_n=dens_average_factor*kappa0_n+(1.d0-dens_average_factor)*kappa_n                
                    call HFB_energy
                    Ediff_bcs=dabs(E_HFB_energy - E_HFB)
                    E_HFB_energy=E_HFB
                    hf_n_diag = .false.
                    hf_p_diag = .true.
                    continue
                endif

                if (bcs_p) then
                    call shift_delta_h(delta_p,h_p,fermi_energy_p,iter_bcs)
                    hf_p_diag = .true.
                    hf_n_diag = .false.
                    call diag_hf(h_p,e_p,h_n,e_n)
                    call hfbcs_levels(h_p,h_n,e_p,e_n,delta_p,delta_n,iter_counter)
                    call densities_update
                    rho_n=dens_average_factor*rho0_n+(1.d0-dens_average_factor)*rho_n
                    rho_p=dens_average_factor*rho0_p+(1.d0-dens_average_factor)*rho_p
                    kappa_p=dens_average_factor*kappa0_p+(1.d0-dens_average_factor)*kappa_p                
                    call HFB_energy
                    Ediff_bcs=dabs(E_HFB_energy - E_HFB)
                    E_HFB_energy=E_HFB
                    hf_p_diag = .false.
                    hf_n_diag = .true.
                    continue
                endif
            enddo
            continue
        endif 

        if ((bcs_p.or.bcs_n).and.iter_counter < iter_hf_initial) then
            call diag_hf(h_p,e_p,h_n,e_n)
            call hf_levels(h_p,h_n,e_p,e_n,iter_counter)
            call densities_update
            continue
        endif

        if (.not.(bcs_p.or.bcs_n).and.iter_counter < max_iter) then
            call diag_hf(h_p,e_p,h_n,e_n)
            call hf_levels(h_p,h_n,e_p,e_n,iter_counter)
            call densities_update
            continue
        endif


    end do
  
    call write_summary

return
end subroutine HFBCS_iteration

!--------------------------------------------------------------------------
subroutine densities_update
    integer :: id
    id = n_sp_levels_calc

    call dgemm('N','T',id,id,id,1.d0,Vp_HFB,max(1,id),Vp_HFB,max(1,id),0.d0,rho_p,max(1,id))
    if (bcs_p) call dgemm('N','T',id,id,id,1.d0,Vp_HFB,max(1,id),Up_HFB,max(1,id),0.d0,kappa_p,max(1,id))

    call dgemm('N','T',id,id,id,1.d0,Vn_HFB,max(1,id),Vn_HFB,max(1,id),0.d0,rho_n,max(1,id))
    if (bcs_n) call dgemm('N','T',id,id,id,1.d0,Vn_HFB,max(1,id),Un_HFB,max(1,id),0.d0,kappa_n,max(1,id))
return
end subroutine densities_update

!--------------------------------------------------------------------------

subroutine HFB_energy

    real(kind=8), allocatable :: Vpp_gen(:,:,:,:),Vnn_gen(:,:,:,:)
    integer :: i,j,k,l,m,n,Jp
    real(kind=8) :: valp,valn

    E_kin=0.d0
    E_prot=0.d0
    E_neut=0.d0
    E_pn=0.d0
    E_pair=0.d0
    E_HFB=0.d0

! Kinetic energy contribution to the total energy

    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            E_kin=E_kin+T1b_p(i,j)*rho_p(j,i)*dble(lev_p(j)%j2+1)+T1b_n(i,j)*rho_n(j,i)*dble(lev_n(j)%j2+1)
        enddo
    enddo

!    Two-body interaction energy contribution to the total energy
!    note that lev_p and lev_n are identical, so only one of them is used in the following loops

    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            if(lev_p(i)%j2 == lev_p(j)%j2) then
                do k=1,n_sp_levels_calc
                    do l=1,n_sp_levels_calc
                        if(lev_p(k)%j2 == lev_p(l)%j2) then
                            do Jp=0,j2_lev_max_calc
                                E_prot=E_prot+0.5d0*Vpp_me(i,k,j,l,Jp)*rho_p(l,k)*rho_p(j,i)*dble(2*Jp+1)
                                E_neut=E_neut+0.5d0*Vnn_me(i,k,j,l,Jp)*rho_n(l,k)*rho_n(j,i)*dble(2*Jp+1)
                                E_pn=E_pn+1.d0*Vpn_me(i,k,j,l,Jp)*rho_p(j,i)*rho_n(l,k)*dble(2*Jp+1)
                            enddo
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo


    E_prot_2b=E_prot
    E_neut_2b=E_neut
    E_pn_2b=E_pn

  
!    Radek's version of 3-body NO2B contribution to the energy
!    note that lev_p and lev_n are identical, so only one of them is used in the following loops

!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP& PRIVATE(i,j,l,m,n,Jp)
!$OMP DO REDUCTION(+:E_prot,E_neut,E_pn)
    do k=1,n_sp_levels_calc
        do n=1,n_sp_levels_calc
            if ((lev_p(k)%j2 == lev_p(n)%j2).and.(lev_p(k)%l == lev_p(n)%l)) then
                do i=1,n_sp_levels_calc
                    do l=1,n_sp_levels_calc
                        if ((lev_p(i)%j2 == lev_p(l)%j2).and.(lev_p(i)%l == lev_p(l)%l)) then

                            do j=1,n_sp_levels_calc
                                do Jp=abs(lev_p(i)%j2-lev_p(j)%j2)/2,(lev_p(i)%j2+lev_p(j)%j2)/2
                                    do m=1,n_sp_levels_calc
                                        if ((lev_p(j)%j2 == lev_p(m)%j2).and.(lev_p(j)%l == lev_p(m)%l)) then

!                                            E_prot=E_prot+(1.d0/6.d0)*rho_p(lp1(l),lp1(i))*rho_p(lp1(m),lp1(j))*rho_p(lp1(n),lp1(k))*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)
!                                            E_neut=E_neut+(1.d0/6.d0)*rho_n(lp1(l),lp1(i))*rho_n(lp1(m),lp1(j))*rho_n(lp1(n),lp1(k))*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)
!                                            E_pn=E_pn +(1.d0/3.d0)*rho_p(lp1(l),lp1(i))*rho_p(lp1(m),lp1(j))*rho_n(lp1(n),lp1(k))*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1) &
!&    +(1.d0/6.d0)*rho_p(lp1(l),lp1(i))*rho_p(lp1(m),lp1(j))*rho_n(lp1(n),lp1(k))*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3) &
!&    +(1.d0/6.d0)*(rho_p(lp1(l),lp1(i))*rho_n(lp1(m),lp1(j))*rho_n(lp1(n),lp1(k)))*((3.d0/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,0,0,1)+(dsqrt(3.d0)/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,1,0,1) &
!&    +(dsqrt(3.d0)/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,0,1,1)+ (1.d0/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1)+ V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3))

                                            E_prot=E_prot+(1.d0/6.d0)*rho_p(l,i)*rho_p(m,j)*rho_p(n,k)*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)
                                            E_neut=E_neut+(1.d0/6.d0)*rho_n(l,i)*rho_n(m,j)*rho_n(n,k)*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)
                                            E_pn=E_pn +(1.d0/3.d0)*rho_p(l,i)*rho_p(m,j)*rho_n(n,k)*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1) &
&    +(1.d0/6.d0)*rho_p(l,i)*rho_p(m,j)*rho_n(n,k)*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3) &
&    +(1.d0/6.d0)*(rho_p(l,i)*rho_n(m,j)*rho_n(n,k))*((3.d0/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,0,0,1)+(dsqrt(3.d0)/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,1,0,1) &
&    +(dsqrt(3.d0)/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,0,1,1)+ (1.d0/2.d0)*V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1)+ V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3))


                                        endif
                                    enddo
                                enddo
                            enddo
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

!      Calculation of the pairing energy

    if(bcs_p.or.bcs_n) then
    if (.not.allocated(Vpp_gen)) allocate(Vpp_gen(n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc),Vnn_gen(n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc))
       Vpp_gen=0.d0
       Vnn_gen=0.d0
    
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP& PRIVATE(j,k,l,valp,valn,m,n)

    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            if ((lev_p(i)%j2 == lev_p(j)%j2).and.(lev_p(i)%l == lev_p(j)%l)) then
                do k=1,n_sp_levels_calc
                    do l=1,n_sp_levels_calc
                        if ((lev_p(k)%j2 == lev_p(l)%j2).and.(lev_p(k)%l == lev_p(l)%l)) then
                            valp=0.d0
                            valn=0.d0

                            do m=1,n_sp_levels_calc
                                do n=1,n_sp_levels_calc
                                    if ((lev_p(m)%j2 == lev_p(n)%j2).and.(lev_p(m)%l == lev_p(n)%l)) then
                                        valp=valp+rho_p(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)+rho_n(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)/3.d0 &
     &                                  +rho_n(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,1)*2.d0/3.d0
                                        valn=valn+rho_n(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)+rho_p(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)/3.d0 &
     &                                  +rho_p(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,1)*2.d0/3.d0
                                    endif 
                                enddo
                            enddo

                            Vpp_gen(i,j,k,l)=Vpp_me(i,j,k,l,0)+valp
                            Vnn_gen(i,j,k,l)=Vnn_me(i,j,k,l,0)+valn
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo

!$OMP END PARALLEL



!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP& PRIVATE(j,k,l)
!$OMP DO REDUCTION(+:E_pair)

    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            do k=1,n_sp_levels_calc
                do l=1,n_sp_levels_calc
                    if(bcs_p) E_pair=E_pair+0.25d0*Vpp_gen(i,j,k,l)*kappa_p(i,j)*kappa_p(l,k)*dsqrt(dble((lev_p(i)%j2+1)*(lev_p(k)%j2+1)))
                    if(bcs_n) E_pair=E_pair+0.25d0*Vnn_gen(i,j,k,l)*kappa_n(i,j)*kappa_n(l,k)*dsqrt(dble((lev_n(i)%j2+1)*(lev_n(k)%j2+1)))
                enddo
            enddo
        enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    end if


!       deallocate(Vpp_gen,Vnn_gen)


    E_HFB = E_kin + E_prot + E_neut + E_pn + E_pair

    write(*,*)'Iteration # ', iter_counter,'     E_HF(BCS)  = ',E_HFB,'MeV'

!       write(*,*) 'HFB energy                 =',E_HFB,'  MeV'
!       write(*,*) 'Kinetic energy             =',E_kin,'  MeV'
!       write(*,*) 'Proton interaction energy  =',E_prot,'  MeV'
!       write(*,*) 'Neutron interaction energy =',E_neut,'  MeV'
!       write(*,*) 'P-N interaction energy     =',E_pn,'  MeV'
!       write(*,*) 'Pairing energy             =',E_pair,'  MeV'

       
    return
end subroutine HFB_energy
!--------------------------------------------------------------------------
subroutine h_fields(h_p,h_n)
   
    real(kind=8) :: h_p(n_sp_levels_calc,n_sp_levels_calc),h_n(n_sp_levels_calc,n_sp_levels_calc)
    real(kind=8) :: h2(n_sp_levels_calc,n_sp_levels_calc),h3(n_sp_levels_calc,n_sp_levels_calc)
    integer :: i,j,k,l,m,n,Jp
    real(kind=8) :: val, val2,v3b1,v3b2,v3b3,v3b4,v3b5

!  calculation of proton field h_p
    h2=0.d0
    h3=0.d0
    h_p=0.d0

    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            val=0.d0
            do k=1,n_sp_levels_calc
                do l=1,n_sp_levels_calc
                    if((lev_p(i)%j2 == lev_p(j)%j2).and.(lev_p(k)%j2 == lev_p(l)%j2)) then
                        do Jp=0,j2_lev_max_calc
                            val=val+Vpp_me(i,k,j,l,Jp)*rho_p(l,k)*dble(2*Jp+1)/dble(lev_p(i)%j2+1)+Vpn_me(i,k,j,l,Jp)*rho_n(l,k)*dble(2*Jp+1)/dble(lev_p(i)%j2+1)
                        enddo
                    endif
                enddo
            enddo
            h2(i,j)=T1b_p(i,j)+val
        enddo
    enddo


!      formulae according to Radek Folprecht
    do i=1,n_sp_levels_calc
        do l=1,n_sp_levels_calc
            if ((lev_p(i)%j2 == lev_p(l)%j2).and.(lev_p(i)%l == lev_p(l)%l)) then
                val2=0.d0
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP& PRIVATE(j,m,k,n,Jp,v3b1,v3b2,v3b3,v3b4,v3b5)
!$OMP DO REDUCTION(+:val2)
         
                do j=1,n_sp_levels_calc
                    do m=1,n_sp_levels_calc
                        if ((lev_p(j)%j2 == lev_p(m)%j2).and.(lev_p(j)%l == lev_p(m)%l)) then
                            do k=1,n_sp_levels_calc
                                do n=1,n_sp_levels_calc
                                    if ((lev_p(k)%j2 == lev_p(n)%j2).and.(lev_p(k)%l == lev_p(n)%l)) then
                                        do Jp=abs(lev_p(i)%j2-lev_p(j)%j2)/2,(lev_p(i)%j2+lev_p(j)%j2)/2
                                            v3b1=V3BNO2_me(i,j,k,l,m,n,Jp,0,0,1)
                                            v3b2=V3BNO2_me(i,j,k,l,m,n,Jp,1,0,1)
                                            v3b3=V3BNO2_me(i,j,k,l,m,n,Jp,0,1,1)
                                            v3b4=V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1)
                                            v3b5=V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)
                                            val2=val2+1.d0/dble(lev_p(i)%j2+1)*((0.25d0*(v3b1+1.d0/dsqrt(3.d0)*v3b2+1.d0/dsqrt(3.d0)*v3b3+1.d0/3.d0*v3b4+2.d0/3.d0*v3b5) &         
     &                                      *rho_n(m,j)*rho_n(n,k)) &
     &                                      +(1.d0/3.d0*(2.d0*v3b4+v3b5)*rho_p(m,j)*rho_n(n,k)) &
     &                                      +0.5d0*v3b5*rho_p(m,j)*rho_p(n,k))
                                        enddo  
                                    endif 
                                enddo   ! k
                            enddo  ! m 
                        endif          
                    enddo    ! m
                enddo   ! j
!$OMP END DO
!$OMP END PARALLEL
                 h3(i,l)=val2
            endif
        enddo  !  l
    enddo  !  i 

    h_p=h2+h3

!  calculation of neutron field h_n
    h2=0.d0
    h3=0.d0
    h_n=0.d0

    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            val=0.d0
            do k=1,n_sp_levels_calc
                do l=1,n_sp_levels_calc
                    if((lev_n(i)%j2 == lev_n(j)%j2).and.(lev_n(k)%j2 == lev_n(l)%j2)) then
                        do Jp=0,j2_lev_max_calc
                            val=val+Vnn_me(i,k,j,l,Jp)*rho_n(l,k)*dble(2*Jp+1)/dble(lev_n(i)%j2+1)+Vpn_me(k,i,l,j,Jp)*rho_p(l,k)*dble(2*Jp+1)/dble(lev_n(i)%j2+1)
                        enddo
                    endif
                enddo
            enddo
            h2(i,j)=T1b_n(i,j)+val
        enddo
    enddo


!      formulae according to Radek Folprecht
    do i=1,n_sp_levels_calc
        do l=1,n_sp_levels_calc
            if ((lev_n(i)%j2 == lev_n(l)%j2).and.(lev_n(i)%l == lev_n(l)%l)) then
                val2=0.d0
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP& PRIVATE(j,m,k,n,Jp,v3b1,v3b2,v3b3,v3b4,v3b5)
!$OMP DO REDUCTION(+:val2)
         
                do j=1,n_sp_levels_calc
                    do m=1,n_sp_levels_calc
                        if ((lev_n(j)%j2 == lev_n(m)%j2).and.(lev_n(j)%l == lev_n(m)%l)) then
                            do k=1,n_sp_levels_calc
                                do n=1,n_sp_levels_calc
                                    if ((lev_n(k)%j2 == lev_n(n)%j2).and.(lev_n(k)%l == lev_n(n)%l)) then
                                        do Jp=abs(lev_n(i)%j2-lev_n(j)%j2)/2,(lev_n(i)%j2+lev_n(j)%j2)/2
                                            v3b1=V3BNO2_me(i,j,k,l,m,n,Jp,0,0,1)
                                            v3b2=V3BNO2_me(i,j,k,l,m,n,Jp,1,0,1)
                                            v3b3=V3BNO2_me(i,j,k,l,m,n,Jp,0,1,1)
                                            v3b4=V3BNO2_me(i,j,k,l,m,n,Jp,1,1,1)
                                            v3b5=V3BNO2_me(i,j,k,l,m,n,Jp,1,1,3)
                                            val2=val2+1.d0/dble(lev_n(i)%j2+1)*((0.25d0*(v3b1+1.d0/dsqrt(3.d0)*v3b2+1.d0/dsqrt(3.d0)*v3b3+1.d0/3.d0*v3b4+2.d0/3.d0*v3b5) &
                                            *rho_p(m,j)*rho_p(n,k)) &
     &                                      +(1.d0/3.d0*(2.d0*v3b4+v3b5)*rho_n(m,j)*rho_p(n,k)) &
     &                                      +0.5d0*v3b5*rho_n(m,j)*rho_n(n,k))
                                        enddo  
                                    endif 
                                enddo   ! k
                            enddo  ! m 
                        endif          
                    enddo    ! m
                enddo   ! j
!$OMP END DO
!$OMP END PARALLEL
                 h3(i,l)=val2
            endif
        enddo  !  l
    enddo  !  i 

    h_n=h2+h3


    return
end subroutine h_fields
!--------------------------------------------------------------------------

subroutine d_fields(delta_p,delta_n)
    
    real(kind=8) :: delta_p(n_sp_levels_calc,n_sp_levels_calc),delta_n(n_sp_levels_calc,n_sp_levels_calc)
    real(kind=8), allocatable :: Vpp_gen(:,:,:,:),Vnn_gen(:,:,:,:)
    real(kind=8) ::valp,valn, val
    integer :: i,j,k,l,m,n



!    if(bcs_p.or.bcs_n) then
!        if (.not.allocated(Vpp_gen)) allocate(Vpp_gen(n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc),Vnn_gen(n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc))
!        Vpp_gen=0.d0
!        Vnn_gen=0.d0
!    endif

    if(bcs_p) then
        if (.not.allocated(Vpp_gen)) allocate(Vpp_gen(n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc))
        Vpp_gen=0.d0
    

!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP& PRIVATE(j,k,l,valp,valn,m,n)

    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            if ((lev_p(i)%j2 == lev_p(j)%j2).and.(lev_p(i)%l == lev_p(j)%l)) then
                do k=1,n_sp_levels_calc
                    do l=1,n_sp_levels_calc
                        if ((lev_p(k)%j2 == lev_p(l)%j2).and.(lev_p(k)%l == lev_p(l)%l)) then
                            valp=0.d0
!                            valn=0.d0

                                do m=1,n_sp_levels_calc
                                    do n=1,n_sp_levels_calc
                                        if ((lev_n(m)%j2 == lev_n(n)%j2).and.(lev_n(m)%l == lev_n(n)%l)) then
                                            valp=valp+rho_p(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)+rho_n(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)/3.d0+rho_n(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,1)*2.d0/3.d0
!                                           valn=valn+rho_n(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)+rho_p(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)/3.d0+rho_p(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,1)*2.d0/3.d0
                                        endif
                                    enddo
                                enddo
                                Vpp_gen(i,j,k,l)=Vpp_me(i,j,k,l,0)+valp
!                                Vnn_gen(i,j,k,l)=Vnn_me(i,j,k,l,0)+valn
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
!$OMP END PARALLEL

    delta_p=0.d0
       
!$OMP PARALLEL DEFAULT(SHARED)&
!$OMP& PRIVATE(j,k,l,val)       
    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            if((lev_p(i)%j2 == lev_p(j)%j2).and.(lev_p(i)%l == lev_p(j)%l)) then
                val=0.d0
                    do k=1,n_sp_levels_calc
                        do l=1,n_sp_levels_calc
                            if((lev_p(k)%j2 == lev_p(l)%j2).and.(lev_p(k)%l == lev_p(l)%l)) val=val+Vpp_gen(k,l,i,j)*kappa_p(l,k)*dsqrt(dble(lev_p(k)%j2+1))  
                        enddo
                    enddo
                delta_p(i,j)=0.5d0*val/dsqrt(dble(lev_p(i)%j2+1))
            endif
        enddo
    enddo
!$OMP END PARALLEL      
    endif

    if(bcs_n) then
        if (.not.allocated(Vnn_gen)) allocate(Vnn_gen(n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc))
        Vnn_gen=0.d0
    

!$OMP PARALLEL DEFAULT(SHARED)& 
!$OMP& PRIVATE(j,k,l,valp,valn,m,n)

    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            if ((lev_n(i)%j2 == lev_n(j)%j2).and.(lev_n(i)%l == lev_n(j)%l)) then
                do k=1,n_sp_levels_calc
                    do l=1,n_sp_levels_calc
                        if ((lev_n(k)%j2 == lev_n(l)%j2).and.(lev_n(k)%l == lev_n(l)%l)) then
!                            valp=0.d0
                            valn=0.d0

                                do m=1,n_sp_levels_calc
                                    do n=1,n_sp_levels_calc
                                        if ((lev_n(m)%j2 == lev_n(n)%j2).and.(lev_n(m)%l == lev_n(n)%l)) then
!                                            valp=valp+rho_p(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)+rho_n(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)/3.d0+rho_n(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,1)*2.d0/3.d0
                                            valn=valn+rho_n(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)+rho_p(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,3)/3.d0+rho_p(n,m)*V3BNO2_me(i,j,m,k,l,n,0,1,1,1)*2.d0/3.d0
                                        endif
                                    enddo
                                enddo
!                                Vpp_gen(i,j,k,l)=Vpp_me(i,j,k,l,0)+valp
                                Vnn_gen(i,j,k,l)=Vnn_me(i,j,k,l,0)+valn
                        endif
                    enddo
                enddo
            endif
        enddo
    enddo
!$OMP END PARALLEL

    delta_n=0.d0
       
!$OMP PARALLEL DEFAULT(SHARED)& 
!$OMP& PRIVATE(j,k,l,val)       
    do i=1,n_sp_levels_calc
        do j=1,n_sp_levels_calc
            if((lev_n(i)%j2 == lev_n(j)%j2).and.(lev_n(i)%l == lev_n(j)%l)) then
                val=0.d0
                    do k=1,n_sp_levels_calc
                        do l=1,n_sp_levels_calc
                            if((lev_n(k)%j2 == lev_n(l)%j2).and.(lev_n(k)%l == lev_n(l)%l)) val=val+Vnn_gen(k,l,i,j)*kappa_n(l,k)*dsqrt(dble(lev_n(k)%j2+1))  
                        enddo
                    enddo
                delta_n(i,j)=0.5d0*val/dsqrt(dble(lev_n(i)%j2+1))
            endif
        enddo
    enddo
!$OMP END PARALLEL
    endif

!       deallocate(Vpp_gen,Vnn_gen)

    return
end subroutine d_fields
!--------------------------------------------------------------------------
subroutine hf_levels(h_p,h_n,e_p,e_n,iter)
    real(kind=8) :: h_p(n_sp_levels_calc,n_sp_levels_calc),h_n(n_sp_levels_calc,n_sp_levels_calc)
    real(kind=8) :: h_vu(n_sp_levels_calc,n_sp_levels_calc)
    real(kind=8) :: e_p(n_sp_levels_calc),e_n(n_sp_levels_calc)
    integer :: jlp(n_sp_levels_calc),jln(n_sp_levels_calc)
    real(kind=8) :: d_max,e_filled,e_empty
    integer :: i,j,j_pointer,ll,jj,nn,mm,lhf,jhf,iter

!   proton levels
    do i=1,n_sp_levels_calc
        d_max=0.d0
        do j=1,n_sp_levels_calc
            if(dabs(h_p(j,i)) > d_max) then
                d_max=dabs(h_p(j,i))
                j_pointer=j
            endif
        enddo
        jlp(i)=1000*lev_p(j_pointer)%l+lev_p(j_pointer)%j2
    enddo
    

    do i=1,n_sp_levels_calc
        lhf=jlp(i)/1000
        jhf=mod(jlp(i),1000)
        levhf_p(i)%par=(-1)**lhf
        levhf_p(i)%l=lhf
        levhf_p(i)%j2=jhf
        levhf_p(i)%e_hf=e_p(i)
        levhf_p(i)%qe=dabs(e_p(i)-fermi_energy_p)
    enddo

    do ll=0,(j2_lev_max_calc+1)/2
        do jj=1,j2_lev_max_calc,2
            nn=0
            do i=1,j2_lev_max_calc
                if(lev_p(i)%j2==jj.and.lev_p(i)%l==ll) then
                    nn=nn+1
                    mm=0
                    do j=1,n_sp_levels_calc
                        if(levhf_p(j)%j2==jj.and.levhf_p(j)%l==ll) then
                            mm=mm+1
                            if(mm==nn) levhf_p(j)%u=lev_p(i)%u
                            if(mm==nn) levhf_p(j)%v=lev_p(i)%v
                        endif
                    enddo
                endif
            enddo
        enddo
    enddo

    Vp_HFB=0.d0
    Up_HFB=0.d0
    h_vu=0.d0    

    do i=1,n_sp_levels_calc
        Up_HFB(i,i)=levhf_p(i)%u
        Vp_HFB(i,i)=levhf_p(i)%v
    enddo

    call dgemm('N','N',n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,1.d0,h_p,n_sp_levels_calc,Vp_HFB,n_sp_levels_calc,0.d0,h_vu,n_sp_levels_calc)

    Vp_HFB=h_vu
    h_vu=0.d0

    call dgemm('N','N',n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,1.d0,h_p,n_sp_levels_calc,Up_HFB,n_sp_levels_calc,0.d0,h_vu,n_sp_levels_calc)
       
    Up_HFB=h_vu
    h_vu=0.d0
       
    if(.not.bcs_p) then
        e_filled=0.d0
        e_empty=0.d0
        
        do i=1,n_sp_levels_calc
            if(dabs(levhf_p(i)%v-1.d0)< 1.d-7) e_filled=levhf_p(i)%e_hf
            if(dabs(levhf_p(n_sp_levels_calc-i+1)%u-1.d0) < 1.d-7) e_empty=levhf_p(n_sp_levels_calc-i+1)%e_hf
        enddo
        fermi_energy_p=0.5d0*(e_empty+e_filled)
!        fermi_energy_p=e_filled
    endif

    if(bcs_p.and.iter <=2) then
        do i=1,n_sp_levels_calc
            if(dabs(levhf_p(i)%v) > 1.d-7) fermi_energy_p=levhf_p(i)%e_hf
        enddo
    endif


!   neutron levels

    do i=1,n_sp_levels_calc
        d_max=0.d0
        do j=1,n_sp_levels_calc
            if(dabs(h_n(j,i)) > d_max) then
                d_max=dabs(h_n(j,i))
                j_pointer=j
            endif
        enddo
        jlp(i)=1000*lev_n(j_pointer)%l+lev_n(j_pointer)%j2
    enddo
    
    if (.not.allocated(levhf_n)) allocate(levhf_n(n_sp_levels_calc))

    do i=1,n_sp_levels_calc
        lhf=jlp(i)/1000
        jhf=mod(jlp(i),1000)
        levhf_n(i)%par=(-1)**lhf
        levhf_n(i)%l=lhf
        levhf_n(i)%j2=jhf
        levhf_n(i)%e_hf=e_n(i)
        levhf_n(i)%qe=dabs(e_n(i)-fermi_energy_n)
    enddo

    do ll=0,(j2_lev_max_calc+1)/2
        do jj=1,j2_lev_max_calc,2
            nn=0
            do i=1,j2_lev_max_calc
                if(lev_n(i)%j2==jj.and.lev_n(i)%l==ll) then
                    nn=nn+1
                    mm=0
                    do j=1,n_sp_levels_calc
                        if(levhf_n(j)%j2==jj.and.levhf_n(j)%l==ll) then
                            mm=mm+1
                            if(mm==nn) levhf_n(j)%u=lev_n(i)%u
                            if(mm==nn) levhf_n(j)%v=lev_n(i)%v
                        endif
                    enddo
                endif
            enddo
        enddo
    enddo

    Vn_HFB=0.d0
    Un_HFB=0.d0
    h_vu=0.d0

    do i=1,n_sp_levels_calc
        Un_HFB(i,i)=levhf_n(i)%u
        Vn_HFB(i,i)=levhf_n(i)%v
    enddo

    call dgemm('N','N',n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,1.d0,h_n,n_sp_levels_calc,Vn_HFB,n_sp_levels_calc,0.d0,h_vu,n_sp_levels_calc)

    Vn_HFB=h_vu
    h_vu=0.d0

    call dgemm('N','N',n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,1.d0,h_n,n_sp_levels_calc,Un_HFB,n_sp_levels_calc,0.d0,h_vu,n_sp_levels_calc)
       
    Un_HFB=h_vu
    h_vu=0.d0

    if(.not.bcs_n) then
        e_filled=0.d0
        e_empty=0.d0
        
        do i=1,n_sp_levels_calc
            if(dabs(levhf_n(i)%v-1.d0)< 1.d-7) e_filled=levhf_n(i)%e_hf
            if(dabs(levhf_n(n_sp_levels_calc-i+1)%u-1.d0) < 1.d-7) e_empty=levhf_n(n_sp_levels_calc-i+1)%e_hf
        enddo
        fermi_energy_n=0.5d0*(e_empty+e_filled)
!        fermi_energy_n=e_filled
    endif

    if(bcs_n.and.iter <=2) then
        do i=1,n_sp_levels_calc
            if(dabs(levhf_n(i)%v) > 1.d-7) fermi_energy_n=levhf_n(i)%e_hf
        enddo
    endif



    open(1,file='HF_levels.dat',status='unknown',form='formatted')
    write(1,*)'Iteration # ', iter
    write(1,*)
    write(1,*)'proton HF levels'
    write(1,*)'   i    l    2j          e        qe         v'
    write(1,*)'-----------------------------------------------'
    do i = 1, n_sp_levels_calc
        write(1,'(3i5,5x,s3f10.5)') i ,levhf_p(i)%l,levhf_p(i)%j2, levhf_p(i)%e_hf,levhf_p(i)%qe,levhf_p(i)%v
    end do
    write(1,*)
    write(1,*)'neutron HF levels' 
    write(1,*)'   i    l    2j          e        qe         v'
    write(1,*)'-----------------------------------------------'
    do i = 1, n_sp_levels_calc
        write(1,'(3i5,5x,s3f10.5)') i ,levhf_n(i)%l,levhf_n(i)%j2, levhf_n(i)%e_hf,levhf_n(i)%qe,levhf_n(i)%v
    end do
    close(1)

 

    return
end subroutine hf_levels
!--------------------------------------------------------------------------
subroutine hfbcs_levels(h_p,h_n,e_p,e_n,delta_p,delta_n,iter)

    real(kind=8) :: h_p(n_sp_levels_calc,n_sp_levels_calc),h_n(n_sp_levels_calc,n_sp_levels_calc)
    real(kind=8) :: h_vu(n_sp_levels_calc,n_sp_levels_calc)
    real(kind=8) :: delta_p(n_sp_levels_calc,n_sp_levels_calc),delta_n(n_sp_levels_calc,n_sp_levels_calc)
    real(kind=8) :: e_p(n_sp_levels_calc),e_n(n_sp_levels_calc)
    real(kind=8) :: u_sum,v_sum,d_max,h_elem,denom,delta_fermi
    integer:: i,j,i_it_pair,j_pointer,ll,jj,iter
    integer :: jlp(n_sp_levels_calc),jln(n_sp_levels_calc)

    if(bcs_n) then
        i_it_pair=0
        v_sum=1.d9
        do while(dabs(v_sum-dble(neutron_number)) >= 1.d-7.and.i_it_pair <= 100000)
            i_it_pair=i_it_pair+1
                do i=1,n_sp_levels_calc
                    d_max=0.d0
                    do j=1,n_sp_levels_calc
                        if(dabs(h_n(j,i)).gt.d_max) then
                            d_max=dabs(h_n(j,i))
                            j_pointer=j
                        endif
                    enddo
            jln(i)=1000*lev_n(j_pointer)%l+lev_n(j_pointer)%j2
        enddo

        do i=1,n_sp_levels_calc
            ll=jln(i)/1000
            jj=mod(jln(i),1000)
            h_elem=e_n(i)
            denom=dsqrt(h_elem**2.d0+delta_n(i,i)**2.d0)
            u_sum=0.5d0*(1.d0+h_elem/denom) !0.5d0*(1.d0+h_elem/dsqrt(w2(nklm)))
            v_sum=0.5d0*(1.d0-h_elem/denom) !0.5d0*(1.d0-h_elem/dsqrt(w2(nklm)))
            levhf_n(i)%par=(-1)**ll
            levhf_n(i)%l=ll
            levhf_n(i)%j2=jj
            levhf_n(i)%qe=denom !dsqrt(w2(nklm))
            levhf_n(i)%e_hf=h_elem+fermi_energy_n
            levhf_n(i)%u=dsqrt(u_sum)*(-1)**levhf_n(i)%l
            levhf_n(i)%v=dsqrt(v_sum)
        enddo

        Vn_HFB=0.d0
        Un_HFB=0.d0
        do i=1,n_sp_levels_calc
            Un_HFB(i,i)=levhf_n(i)%u
            Vn_HFB(i,i)=levhf_n(i)%v!*(-1)**((levhf_n(i)%j2-1)/2)
        enddo

        call dgemm('N','N',n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,1.d0,h_n,n_sp_levels_calc,Vn_HFB,n_sp_levels_calc,0.d0,h_vu,n_sp_levels_calc)
        Vn_HFB=h_vu
        call dgemm('N','N',n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,1.d0,h_n,n_sp_levels_calc,Un_HFB,n_sp_levels_calc,0.d0,h_vu,n_sp_levels_calc)
        Un_HFB=h_vu

        v_sum=0.d0
        do j=1,n_sp_levels_calc
            v_sum=v_sum+levhf_n(j)%v**2.d0*dble(levhf_n(j)%j2+1)
        enddo

        delta_fermi=-0.001d0*(v_sum-dble(neutron_number))
        fermi_energy_n=fermi_energy_n+delta_fermi
        do j=1,n_sp_levels_calc
            e_n(j)=e_n(j)-delta_fermi
        enddo
        enddo
    endif

    if(bcs_p) then
        i_it_pair=0
        v_sum=1.d9
        do while(dabs(v_sum-dble(proton_number)) >= 1.d-7.and.i_it_pair <= 100000)
            i_it_pair=i_it_pair+1
                do i=1,n_sp_levels_calc
                    d_max=0.d0
                    do j=1,n_sp_levels_calc
                        if(dabs(h_p(j,i)).gt.d_max) then
                            d_max=dabs(h_p(j,i))
                            j_pointer=j
                        endif
                    enddo
            jlp(i)=1000*lev_p(j_pointer)%l+lev_p(j_pointer)%j2
        enddo

        do i=1,n_sp_levels_calc
            ll=jlp(i)/1000
            jj=mod(jlp(i),1000)
            h_elem=e_p(i)
            denom=dsqrt(h_elem**2.d0+delta_p(i,i)**2.d0)
            u_sum=0.5d0*(1.d0+h_elem/denom) !0.5d0*(1.d0+h_elem/dsqrt(w2(nklm)))
            v_sum=0.5d0*(1.d0-h_elem/denom) !0.5d0*(1.d0-h_elem/dsqrt(w2(nklm)))
            levhf_p(i)%par=(-1)**ll
            levhf_p(i)%l=ll
            levhf_p(i)%j2=jj
            levhf_p(i)%qe=denom !dsqrt(w2(nklm))
            levhf_p(i)%e_hf=h_elem+fermi_energy_p
            levhf_p(i)%u=dsqrt(u_sum)*(-1)**levhf_p(i)%l
            levhf_p(i)%v=dsqrt(v_sum)
        enddo

        Vp_HFB=0.d0
        Up_HFB=0.d0
        do i=1,n_sp_levels_calc
            Up_HFB(i,i)=levhf_p(i)%u
            Vp_HFB(i,i)=levhf_p(i)%v!*(-1)**((levhf_p(i)%j2-1)/2)
        enddo

        call dgemm('N','N',n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,1.d0,h_p,n_sp_levels_calc,Vp_HFB,n_sp_levels_calc,0.d0,h_vu,n_sp_levels_calc)
        Vp_HFB=h_vu
        call dgemm('N','N',n_sp_levels_calc,n_sp_levels_calc,n_sp_levels_calc,1.d0,h_p,n_sp_levels_calc,Up_HFB,n_sp_levels_calc,0.d0,h_vu,n_sp_levels_calc)
        Up_HFB=h_vu

        v_sum=0.d0
        do j=1,n_sp_levels_calc
            v_sum=v_sum+levhf_p(j)%v**2.d0*dble(levhf_p(j)%j2+1)
        enddo

        delta_fermi=-0.001d0*(v_sum-dble(proton_number))
        fermi_energy_p=fermi_energy_p+delta_fermi
        do j=1,n_sp_levels_calc
            e_p(j)=e_p(j)-delta_fermi
        enddo
        enddo
    endif


    open(1,file='HFBCS_levels.dat',status='unknown',form='formatted')
    write(1,*)'Iteration # ', iter
    write(1,*)
    write(1,*)'proton HF levels'
    write(1,*)'   i    l    2j          e        qe         v'
    write(1,*)'-----------------------------------------------'
    do i = 1, n_sp_levels_calc
        write(1,'(3i5,5x,s3f10.5)') i ,levhf_p(i)%l,levhf_p(i)%j2, levhf_p(i)%e_hf,levhf_p(i)%qe,levhf_p(i)%v
    end do
    write(1,*)
    write(1,*)'neutron HF levels' 
    write(1,*)'   i    l    2j          e        qe         v'
    write(1,*)'-----------------------------------------------'
    do i = 1, n_sp_levels_calc
        write(1,'(3i5,5x,s3f10.5)') i ,levhf_n(i)%l,levhf_n(i)%j2, levhf_n(i)%e_hf,levhf_n(i)%qe,levhf_n(i)%v
    end do
    close(1)

    return
end subroutine hfbcs_levels
!--------------------------------------------------------------------------
subroutine diag_hf(h_p,e_p,h_n,e_n)

    real(kind=8) :: h_p(n_sp_levels_calc,n_sp_levels_calc),h_n(n_sp_levels_calc,n_sp_levels_calc)
    real(kind=8) :: e_p(n_sp_levels_calc),e_n(n_sp_levels_calc)
    integer :: info, lwork
    real(kind=8), allocatable :: work(:)
    
    lwork=3*n_sp_levels_calc-1
    if (.not.allocated(work)) allocate(work(lwork))
 
    ! LAPACK routine to diagonalize symmetric matrix
    if (hf_p_diag) then
        call dsyev('V','U', n_sp_levels_calc, h_p, n_sp_levels_calc, e_p, work, lwork, info)
        if(info /= 0) then
            write(*,*) 'Diagonalization  of h_p failed!!!'
            stop
        endif
    endif

    if (hf_n_diag) then
        call dsyev('V','U', n_sp_levels_calc, h_n, n_sp_levels_calc, e_n, work, lwork, info)
        if(info /= 0) then
            write(*,*) 'Diagonalization  of h_n failed!!!'
            stop
        endif
    endif

    return
end subroutine diag_hf
!--------------------------------------------------------------------------
subroutine shift_delta_h(delta,h,fermi_energy,iter_bcs)
    real(kind=8) :: delta(n_sp_levels_calc,n_sp_levels_calc),h(n_sp_levels_calc,n_sp_levels_calc)
    real(kind=8) :: fermi_energy
    integer :: iter_bcs
    integer :: i,j

    do i=1,n_sp_levels_calc
        h(i,i)=h(i,i)-fermi_energy
        if (iter_bcs==1) then 
            do j=1,n_sp_levels_calc
                if(lev_n(i)%j2 == lev_n(j)%j2.and.lev_n(i)%l == lev_n(j)%l) delta(i,j)=delta(i,j)+0.1d0
            enddo
        endif
    enddo

    return

end subroutine shift_delta_h

subroutine write_summary
   
    integer :: i
    
    open(unit=10,file='HFBCS_summary.dat',status='unknown', form='formatted')
    write(10,*) 'HF-BCS calculation summary'
    write(10,*) 'Proton number: ', proton_number
    write(10,*) 'Neutron number: ', neutron_number
    write(10,*) 'hbar_omega: ', hbar_omega
    write(10,*) 'Model space used in calculation: ', Nmax1, Nmax12, Nmax123
    write(10,*) 'Number of single-particle levels: ', n_sp_levels_calc
    write(10,*) 'Total energy: ', E_HFB
    write(10,*) 'Pairing energy: ', E_pair
    write(10,*) 
    
    write(10,*) 'Number of iterations: ', iter_counter
    write(10,*)
    write(10,*)'proton HF levels'
    write(10,*)'   i    l    2j          e        qe         v'
    write(10,*)'-----------------------------------------------'
    do i = 1, n_sp_levels_calc
        write(10,'(3i5,5x,s3f10.5)') i ,levhf_p(i)%l,levhf_p(i)%j2, levhf_p(i)%e_hf,levhf_p(i)%qe,levhf_p(i)%v
    end do
    write(10,*)
    write(10,*)'neutron HF levels' 
    write(10,*)'   i    l    2j          e        qe         v'
    write(10,*)'-----------------------------------------------'
    do i = 1, n_sp_levels_calc
        write(10,'(3i5,5x,s3f10.5)') i ,levhf_n(i)%l,levhf_n(i)%j2, levhf_n(i)%e_hf,levhf_n(i)%qe,levhf_n(i)%v
    end do

    close(10)
    return
end subroutine write_summary

end module HFBCS_solver