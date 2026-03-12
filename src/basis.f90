module basis
    use declarations
    implicit none
   
contains
!--------------------------------------------------------------------------
!  Subroutine to initialize the basis levels
subroutine initialize_basis

integer :: i,j,k
integer :: N,l,jj,ii,nn,im
integer :: Ni,Nj,Nk,J12,itab,itot,it

    n_sp_levels = dim_HO(Nmax1_2N_file)
    n_sp_levels_calc = dim_HO(Nmax1)
    allocate(lev_p(n_sp_levels),lev_n(n_sp_levels))
    im = 0

    do N = 0,Nmax1_2N_file
        do l = N,0,-1
            if(mod(N,2).eq.mod(l,2)) then
                nn = (N - l)/2
                    do jj = 2*l+1, 2*l-1, -2
                        if(jj > 0) then
                            im = im + 1
                            lev_p(im)%index = im
                            lev_p(im)%par = (-1)**l
                            lev_p(im)%Nosc = N
                            lev_p(im)%nr = nn
                            lev_p(im)%l = l
                            lev_p(im)%j2 = jj
                            lev_p(im)%e_ho = hbar_omega * (dble(N)+1.5d0)
                            lev_p(im)%e_hf = lev_p(im)%e_ho
                            lev_p(im)%qe = lev_p(im)%e_ho
                        endif
                    end do
            endif
        end do
    end do

    lev_n = lev_p

    call occupations

!  Myiagi's ordering of sp. levels for reading 2b matrix elements
    allocate(lev(n_sp_levels))
    im = 0

    do N = 0,Nmax1_2N_file
        do l = 0,N
            if(mod(N,2) == mod(l,2)) then
                nn = (N - l)/2
                    do jj = 2*l-1, 2*l+1, 2
                        if(jj > 0) then
                            im = im + 1
                            lev(im)%index = im
                            lev(im)%par = (-1)**l
                            lev(im)%Nosc = N
                            lev(im)%nr = nn
                            lev(im)%l = l
                            lev(im)%j2 = jj
                            lev(im)%e_ho = hbar_omega * (dble(N)+1.5d0)
                            lev(im)%e_hf = lev(im)%e_ho
                            lev(im)%qe = lev(im)%e_ho
                        endif
                    end do
            endif
        end do
    end do

    allocate(lp1(n_sp_levels),lp2(n_sp_levels))
    lp1=0
    lp2=0

    do i=1,n_sp_levels
      do j=1,n_sp_levels
         if(lev_p(i)%nr == lev(j)%nr.and.(lev_p(i)%l == lev(j)%l.and.lev_p(i)%j2 == lev(j)%j2)) then
          lp1(j)=i
          lp2(i)=j
         endif
      enddo
    enddo

    j2_lev_max=0
    do i=1,n_sp_levels
        if(lev_p(i)%j2 > j2_lev_max) j2_lev_max=lev_p(i)%j2
    enddo

    j2_lev_max_calc=0
    do i=1,n_sp_levels_calc
        if(lev_p(i)%j2 > j2_lev_max_calc) j2_lev_max_calc=lev_p(i)%j2
    enddo

    n_sp_levels_3N=dim_HO(Nmax1_3N_file)
    im=0
    do i=1,n_sp_levels_3N
        Ni=lev(i)%Nosc
        if(Ni <= Nmax1_3N_file) then
            do j=1,i 
                Nj=lev(j)%Nosc
                    if(Ni+Nj <= Nmax12_3N_file) then
                        do k=1,n_sp_levels_3N
                            Nk=lev(k)%Nosc
                            if(Ni+Nj+Nk <=Nmax123_3N_file) then
                                do J12=abs(lev(i)%j2-lev(j)%j2)/2,(lev(i)%j2+lev(j)%j2)/2
                                    do itab=0,1
                                        do itot=abs(2*itab-1),2*itab+1,2 !1,3,2
                                            im=im+1
                                        enddo
                                    enddo
                               enddo

                            endif
                        enddo
                    endif
            enddo
        endif
    enddo

    j2_lev3_max=0
    do i=1,n_sp_levels_3N
        if(lev(i)%j2 > j2_lev3_max) j2_lev3_max=lev(i)%j2
    enddo

    allocate(lev3(im))
    allocate(lpoint(n_sp_levels_3N,n_sp_levels_3N,n_sp_levels_3N,0:j2_lev3_max,0:1,3))
    lpoint=0
    im=0

    do i=1,n_sp_levels_3N
        Ni=lev(i)%Nosc
        if(Ni <= Nmax1_3N_file) then
            do j=1,i
                Nj=lev(j)%Nosc
                if(Ni+Nj <= Nmax12_3N_file) then
                    do k=1,n_sp_levels_3N
                        Nk=lev(k)%Nosc
                        if(Ni+Nj+Nk <= Nmax123_3N_file) then
                            do J12=abs(lev(i)%j2-lev(j)%j2)/2,(lev(i)%j2+lev(j)%j2)/2
                                do itab=0,1
                                    do itot=abs(2*itab-1),2*itab+1,2 !1,3,2
                                        im=im+1
                                        lev3(im)%index=im
                                        lev3(im)%i=i
                                        lev3(im)%j=j
                                        lev3(im)%k=k
                                        lev3(im)%Jab=J12
                                        lev3(im)%Tab=itab
                                        lev3(im)%TT=itot
                                        lpoint(i,j,k,J12,itab,itot)=im
                                    enddo
                                enddo
                            enddo
                        endif
                    enddo
                endif
            enddo
        endif
    enddo


    allocate(idim3j(0:j2_lev3_max))
    allocate(idim3jt(0:j2_lev3_max,1:3))
    idim3j=0
    idim3jt=0

!  ordering to J blocks 
    allocate(lev3ord(im))
    jj=0
    do it=1,3,2
        do j=0,j2_lev3_max
            do ii=1,im 
                if (lev3(ii)%Jab.eq.j.and.lev3(ii)%TT.eq.it) then
                    jj=jj+1
                    lev3ord(jj)=ii
                    lev3(ii)%ord=jj
                endif
            enddo
        enddo  
    enddo

return
end subroutine initialize_basis

!--------------------------------------------------------------------------
!  Subroutine to set initial occupations in HF(BCS)

subroutine occupations

!use declarations

integer :: id,ipp,inn,i
logical :: occ_p,occ_n

    id=n_sp_levels_calc

    allocate(Up_HFB(id,id),Un_HFB(id,id),Vp_HFB(id,id),Vn_HFB(id,id))
    Up_HFB=0.d0
    Un_HFB=0.d0
    Vp_HFB=0.d0
    Vn_HFB=0.d0

!   Here the initial occupations are determined                          
    ipp=0
    inn=0
    occ_p=.true.
    occ_n=.true.
       
    do i=1,id
        ipp=ipp+lev_p(i)%j2+1
        inn=inn+lev_n(i)%j2+1

        if(.not.occ_p) then
         lev_p(i)%v=0.d0
         lev_p(i)%u=1.d0*(-1)**lev_p(i)%l ! note the phase according to Suhonen's book
         Vp_HFB(i,i)=0.d0
         Up_HFB(i,i)=1.d0
        endif
        
        if(occ_p) then
         if(ipp <= proton_number) then
          lev_p(i)%v=1.d0
          lev_p(i)%u=0.d0
          Vp_HFB(i,i)=1.d0
          Up_HFB(i,i)=0.d0
         endif

         if(ipp > proton_number) then
          occ_p=.false.
          lev_p(i)%v=(dble(lev_p(i)%j2+1-ipp+proton_number)/dble(lev_p(i)%j2+1))**0.5d0 
          lev_p(i)%u=dsqrt(dble(ipp-proton_number)/dble(lev_p(i)%j2+1))*(-1)**lev_p(i)%l ! note the phase according to Suhonen's book
          Vp_HFB(i,i)=lev_p(i)%v
          Up_HFB(i,i)=lev_p(i)%u
         endif

        endif

        if(.not.occ_n) then
         lev_n(i)%v=0.d0
         lev_n(i)%u=1.d0*(-1)**lev_n(i)%l ! note the phase according to Suhonen's book
         Vn_HFB(i,i)=0.d0
         Un_HFB(i,i)=1.d0
        endif
        if(occ_n) then
         if(inn <= neutron_number) then
          lev_n(i)%v=1.d0
          lev_n(i)%u=0.d0
          Vn_HFB(i,i)=1.d0
          Un_HFB(i,i)=0.d0
         endif
         if(inn > neutron_number) then
          occ_n=.false.
          lev_n(i)%v=(dble(lev_n(i)%j2+1-inn+neutron_number)/dble(lev_n(i)%j2+1))**0.5d0
          lev_n(i)%u=dsqrt(dble(inn-neutron_number)/dble(lev_n(i)%j2+1))*(-1)**lev_n(i)%l ! note the phase according to Suhonen's book
          Vn_HFB(i,i)=lev_n(i)%v
          Un_HFB(i,i)=lev_n(i)%u
         endif
        endif
    enddo

    allocate(levhf_p(n_sp_levels_calc),levhf_n(n_sp_levels_calc))
    do i=1,n_sp_levels_calc
        levhf_p(i) =lev_p(i)
        levhf_n(i) =lev_n(i)
    enddo

return
end subroutine occupations

!--------------------------------------------------------------------------
!  Function to calculate the dimension of the HO basis up to Nmax
function dim_HO(Nmax) result(dim)

integer, intent(in) :: Nmax
integer :: dim
! number of j-shells in a harmonic oscillator basis up to Nmax
dim = (Nmax + 1)*(Nmax + 2)/2

end function dim_HO 

end module basis