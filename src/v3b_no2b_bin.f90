module v3body_no2b_bin
! this module should be completely rewritten in the future, but for the moment it contains only the function to read the 3-body matrix elements from a binary file
! several functions adopted and modified from NuHamil code by T.Myiagi

use declarations

contains

real function V3BNO2_me(i_lev,j_lev,k_lev,l_lev,m_lev,n_lev,Jij,Tij,Tlm,Ttot)

integer :: i,j,k,l,m,n,Jij,Tij,Tlm,Ttot,ii,jj,ipoz,ipozj
integer :: i_lev,j_lev,k_lev,l_lev,m_lev,n_lev 
integer :: ibra,iket,phase_bra,phase_ket

V3BNO2_me=0.0

! we use different order ordering of sp levels, so we need to find the corresponding indices in the Myiagi's ordering

i=lp2(i_lev)
j=lp2(j_lev)
k=lp2(k_lev)
l=lp2(l_lev)
m=lp2(m_lev)
n=lp2(n_lev)

if (i >= j ) then
  ibra=lpoint(i,j,k,Jij,Tij,Ttot)
  phase_bra=1
  else 
  ibra=lpoint(j,i,k,Jij,Tij,Ttot)
  phase_bra=(-1)**(Jij+Tij-(lev(i)%j2+lev(j)%j2)/2)
endif

if (l>=m) then
  iket=lpoint(l,m,n,Jij,Tlm,Ttot)
  phase_ket=1
  else
  iket=lpoint(m,l,n,Jij,Tlm,Ttot)
  phase_ket=(-1)**(Jij+Tlm-(lev(l)%j2+lev(m)%j2)/2)
endif



if(ibra.ne.0.and.iket.ne.0) then
 V3BNO2_me=V3Bno(max(ibra,iket),min(ibra,iket))*float(phase_bra*phase_ket)
endif 

return 

end function V3BNO2_me
!--------------------------------------------------------------------------
real function V3Bno(i,j)
implicit none
integer :: i,j,ii,jj,it
integer(kind=8) :: iv3b,arrl
integer(kind=8) :: low,high,middle,rang
integer(kind=8) :: ipoz,ipozjt,pair_ij
real :: diff,var

V3Bno=0.0

arrl=size(V3B_ar,DIM=1,KIND=8)
low=1
high=arrl
rang=high-low
middle=(low+high)/2

!This is only for testing
!diff=0.0

if (lev3(i)%jab.ne.lev3(j)%jab) return

ipoz=0
ipozjt=0
!do it=1,lev3(j)%TT-2,2
!do jj=0,lev3(j)%jab-1
! ipoz=ipoz+(idim3jt(jj,it)*(idim3jt(jj,it)+1)/2)
! ipozjt=ipozjt+idim3jt(jj,it)
!enddo
!enddo


if (lev3(j)%TT==1) then  
 do jj=0,lev3(j)%jab-1
   ipoz=ipoz+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
   ipozjt=ipozjt+idim3jt(jj,1)
 enddo
endif 

if (lev3(j)%TT==3) then  

!         do jj=0,jmax
!          ipoz=ipoz+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
!          ipozjt=ipozjt+idim3jt(jj,1)
!         enddo
  ipoz=ipozt1
  ipozjt=ipozjt1
  do jj=0,lev3(j)%jab-1
    ipoz=ipoz+(idim3jt(jj,3)*(idim3jt(jj,3)+1)/2)
    ipozjt=ipozjt+idim3jt(jj,3)
  enddo
endif 

ii=lev3(i)%ord-ipozjt
jj=lev3(j)%ord-ipozjt

if ((ii.ne.0).and.(jj.ne.0)) then
  if (i >= j ) then
!     diff=V3B_ar2((int(ii,8)*(int(ii,8)-1))/2+int(jj,8)+ipoz)
     pair_ij=(int(ii,8)*(int(ii,8)-1))/2+int(jj,8)+ipoz

     do while((V3B_pair(middle).ne.pair_ij).and.(rang.gt.0))
        if(pair_ij>V3B_pair(middle)) then
           low=middle+1
        else
           high=middle-1
        end if
        
        rang=high-low
        middle=(high+low)/2
     end do
     if(V3B_pair(middle)==pair_ij) then
       V3Bno=V3B_ar(middle)
     end if

  endif
endif 

!if (abs(diff-V3Bno)>1.0e-6) then
!  write(*,*)'Diff at',i,j,V3Bno,diff
!endif 


if (i<j) then 
  write(*,*)'i<j in V3B !'
  stop
endif 

return
end function V3Bno
!---------------------------------------------------------------------------

subroutine read_3b_no2b_bin(e1max,e2max,e3max,d1max,d2max,d3max)

integer(kind=8) :: n3b_el

integer :: e1max, e2max, e3max, lmax
integer :: d1max, d2max, d3max
integer :: nlev
integer(kind=8) :: cnt,total_cnt,nelms,delms
integer(kind=8) :: pair_ij,tmp4,tmp6,iv3b,iiv3b,jjv3b
integer :: i
integer :: i1, n1, l1, j1, e1, ch12
integer :: i2, n2, l2, j2, e2
integer :: i3, n3, l3, j3, e3, ch3
integer :: i4, n4, l4, j4, e4
integer :: i5, n5, l5, j5, e5
integer :: i6, n6, l6, j6, e6
integer :: ii1, ii2, ii3, ii4, ii5, ii6
integer :: jj, iii1, iii2
integer :: T12, T45, T, P123, P456, J
integer :: ch, ch_t, ch_no2b, bra, ket, iphase
integer :: buffer_size=100000000
integer :: runit = 22
real(kind=4), allocatable :: v(:), v_tmp1(:) 
integer(kind=8), allocatable :: v_tmp2(:)
real(kind=4) :: v3me

nelms=count_3bme_no2b(e1max, e2max, e3max)
delms=count_3bme_no2b(d1max, d2max, d3max)

write(*,*) "Calculated number of 3b No2b MEs: ", nelms

nlev=(e1max+1)*(e1max+2)/2  ! number of single-particle levels
jmax=0
do i=1,nlev
  if (lev(i)%j2 > jmax) jmax=lev(i)%j2
enddo

idim_v3b=0
do T=1,3,2
 do i=0,jmax
   idim_v3b=idim_v3b+idim3jt(i,T)*(idim3jt(i,T)+1)/2
  enddo
enddo

ipozt1=0

ipozjt1=0
  do jj=0,jmax
    ipozt1=ipozt1+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
    ipozjt1=ipozjt1+idim3jt(jj,1)
  enddo

write(*,*)'Loading 3B m.e.'
!write(*,*)'V3B2 array size (GB)', idim_v3b*4.d0/(1024.d0*1024.d0*1024.d0)
!allocate(V3B_ar2(idim_v3b))
!V3B_ar2=0.0                                      
write(*,*)'V3B array size (GB)', nelms*4.d0/(1024.d0*1024.d0*1024.d0)

allocate(v_tmp1(nelms))
allocate(v_tmp2(nelms))
v_tmp1=0.0
v_tmp2=0
iv3b=1

!!allocate(v(nelms))
allocate(v(buffer_size))
v=0.0

open(runit, file='3belem_NO2.stream.bin', action='read', iostat=io,form='unformatted',access='stream')

    lmax = e1max
    cnt = 0
    total_cnt = 0
    do i1 = 1, nlev
      n1 = lev(i1)%nr
      l1 = lev(i1)%l
      j1 = lev(i1)%j2
      e1 = lev(i1)%Nosc
      if(e1 > e1max) cycle
      if(l1 > lmax) cycle
      do i2 = 1, i1
        n2 = lev(i2)%nr
        l2 = lev(i2)%l
        j2 = lev(i2)%j2
        e2 = lev(i2)%Nosc
        if(e1 + e2 > e2max) cycle
        do i3 = 1, nlev
          n3 = lev(i3)%nr
          l3 = lev(i3)%l
          j3 = lev(i3)%j2
          e3 = lev(i3)%Nosc
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle

          P123 = (-1) ** (l1+l2+l3)
          do i4 = 1, i1
            n4 = lev(i4)%nr
            l4 = lev(i4)%l
            j4 = lev(i4)%j2
            e4 = lev(i4)%Nosc

            do i5 = 1, i4
              n5 = lev(i5)%nr
              l5 = lev(i5)%l
              j5 = lev(i5)%j2
              e5 = lev(i5)%Nosc
              if(e4 + e5 > e2max) cycle

              do i6 = 1, nlev
                n6 = lev(i6)%nr
                l6 = lev(i6)%l
                j6 = lev(i6)%j2
                e6 = lev(i6)%Nosc
                if(j3 /= j6) cycle
                if(l3 /= l6) cycle
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2), (j4+j5))/2
                  do T12 = 0, 1
                    do T45 = 0, 1
                      do T = max(abs(2*T12-1),abs(2*T45-1)),&
                            &min(   (2*T12+1),   (2*T45+1)), 2
                        if(cnt==buffer_size) cnt=0
                        if(cnt==0) then
!#ifdef single_precision_three_body_file
                          v(:) = 0.0
!#elif half_precision_three_body_file
!                          v(:) = 0
!#else
!                          v(:) = 0.d0
!#endif
                          if(nelms - total_cnt >= buffer_size) then
                            read(runit) v(:)
                          else
                            read(runit)v(:nelms-total_cnt)
                          end if
                        end if
                        cnt = cnt + 1
                        total_cnt = total_cnt + 1
!!! changes here
                        if(l1 > d1max) cycle
                        if(l2 > d1max) cycle
                        if(l3 > d1max) cycle

                        if(l4 > d1max) cycle
                        if(l5 > d1max) cycle
                        if(l6 > d1max) cycle

                        if(e1 > d1max) cycle
                        if(e2 > d1max) cycle
                        if(e3 > d1max) cycle

                        if(e4 > d1max) cycle
                        if(e5 > d1max) cycle
                        if(e6 > d1max) cycle

                        if(e1 + e2 > d2max) cycle
                        if(e2 + e3 > d2max) cycle
                        if(e3 + e1 > d2max) cycle

                        if(e4 + e5 > d2max) cycle
                        if(e5 + e6 > d2max) cycle
                        if(e6 + e4 > d2max) cycle

!                        if(e1 + e2 + e3 > e3max) cycle
!                        if(e4 + e5 + e6 > e3max) cycle
!                        for the moment only d3max is implemented
                        
                        if(e1 + e2 + e3 > d3max) cycle
                        if(e4 + e5 + e6 > d3max) cycle

                        
!                        ii1 = sps%nlj2idx(n1,l1,j1)
!                        ii2 = sps%nlj2idx(n2,l2,j2)
!                        ii3 = sps%nlj2idx(n3,l3,j3)
!                        ii4 = sps%nlj2idx(n4,l4,j4)
!                        ii5 = sps%nlj2idx(n5,l5,j5)
!                        ii6 = sps%nlj2idx(n6,l6,j6)
                         ii1=i1
                         ii2=i2  
                         ii3=i3
                         ii4=i4
                         ii5=i5
                         ii6=i6

!                         write(999,'(10i5,f10.5)')ii1,ii2,ii3,T12,ii4,ii5,ii6,T45,J,T,v(cnt)

                        if(ii1==ii2 .and. mod(J+T12,2)==0) then
!#ifdef half_precision_three_body_file
!                          if(v(cnt) /= 0 .and. v(cnt) /= -32768) then
!                            write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
!                            write(*,"(10i4,i16,i8)") ii1,ii2,ii3,J,T12,ii4,ii5,ii6,J,T45,cnt,v(cnt)
!                          end if
!#else
                          if(abs(v(cnt)) > 1.d-6) then
                            write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                            write(*,"(10i4,i16,f12.6)") ii1,ii2,ii3,J,T12,ii4,ii5,ii6,J,T45,cnt,v(cnt)
                          end if
!#endif
                          cycle
                        end if

                        if(ii4==ii5 .and. mod(J+T45,2)==0) then
!#ifdef half_precision_three_body_file
!                          if(v(cnt) /= 0 .and. v(cnt) /= -32768) then
!                            write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
!                            write(*,"(10i4,i16,i8)") ii1,ii2,ii3,J,T12,ii4,ii5,ii6,J,T45,cnt,v(cnt)
!                          end if
!#else
                          if(abs(v(cnt)) > 1.d-6) then
                            write(*,*) "Warning: something wrong, this three-body matrix element has to be zero."
                            write(*,"(10i4,i16,f12.6)") ii1,ii2,ii3,J,T12,ii4,ii5,ii6,J,T45,cnt,v(cnt)
                          end if
!#endif
                          cycle
                        end if

!#ifdef half_precision_three_body_file
!                        tmp = v(cnt)
!                        v123 = tmp
!                        call thr%SetThBMENO2B(ii1,ii2,ii3,T12,ii4,ii5,ii6,T45,J,T,v123)
!#else
!!                        call thr%SetThBMENO2B(ii1,ii2,ii3,T12,ii4,ii5,ii6,T45,J,T,v(cnt))
!#endif

         iii1=lpoint(i1,i2,i3,J,T12,T)
         iii2=lpoint(i4,i5,i6,J,T45,T)

        if (iii1 == 0.or.iii2 == 0) then
          write(*,*) 'Error in 3-body matrix elements pointer'
          stop
        endif

         ipoz=0
         ipozjt=0

         if (T==1) then  
         do jj=0,J-1
          ipoz=ipoz+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
          ipozjt=ipozjt+idim3jt(jj,1)
         enddo
         endif 



         if (T==3) then  

!         do jj=0,jmax
!          ipoz=ipoz+(idim3jt(jj,1)*(idim3jt(jj,1)+1)/2)
!          ipozjt=ipozjt+idim3jt(jj,1)
!         enddo
          ipoz=ipozt1
          ipozjt=ipozjt1
         do jj=0,J-1
          ipoz=ipoz+(idim3jt(jj,3)*(idim3jt(jj,3)+1)/2)
          ipozjt=ipozjt+idim3jt(jj,3)
         enddo
         endif 

 
         ii1j=lev3(iii1)%ord-ipozjt
         ii2j=lev3(iii2)%ord-ipozjt

         if((iii1.ne.0).and.(iii2.ne.0)) then
         if (ii1j >= ii2j ) then
!            V3B_ar2((int(ii1j,8)*(int(ii1j,8)-1))/2+int(ii2j,8)+ipoz)=v(cnt)
            if (abs(v(cnt)) > 1.d-9) then
             pair_ij=(int(ii1j,8)*(int(ii1j,8)-1))/2+int(ii2j,8)+ipoz
             v_tmp1(iv3b)=v(cnt)
             v_tmp2(iv3b)=pair_ij
             iv3b = iv3b + 1
            endif               
!         else
!            V3B_ar2((int(ii2j,8)*(int(ii2j,8)-1))/2+int(ii1j,8)+ipoz)=v(cnt)
!            if (abs(v(cnt)) > 1.d-9) then
!             pair_ij=(int(ii2j,8)*(int(ii2j,8)-1))/2+int(ii1j,8)+ipoz
!             v_tmp1(iv3b)=v(cnt)
!             v_tmp2(iv3b)=pair_ij
!             iv3b = iv3b + 1
!            endif
         endif
!          V3B(ii1,ii2)=xv3*V3b_NNN
!          V3B(ii2,ii1)=xv3*V3b_NNN
         endif

                      end do
                    end do
                  end do
                end do


              end do
            end do
          end do

        end do
      end do
    end do
    close(runit)
    deallocate(v)

    write(*,*) "Compressed 3b No2b MEs:", iv3b 
    
!!! sorting array
! fist keep only the non-zero values
    iv3b=iv3b-1
    
    allocate(V3B_ar(0:iv3b))
     V3B_ar=0
     do iiv3b=1,iv3b
       V3B_ar(iiv3b)=v_tmp1(iiv3b)
    enddo
    deallocate(v_tmp1)
 
    allocate(V3B_pair(0:iv3b))
     V3B_pair=0
     do iiv3b=1,iv3b
       V3B_pair(iiv3b)=v_tmp2(iiv3b)
    enddo   
    deallocate(v_tmp2)

! then we order the arrays in growing order of pair_ij
    
    call sort2(iv3b,V3B_pair,V3B_ar)   

! finally we set index 0 to 0 for possible anomaly in v3body.f  

   V3B_ar(0)=0.0
   V3B_pair(0)=0
   write(*,*)'Compressed V3B array size (GB)', iv3b*4.d0/(1024.d0*1024.d0*1024.d0)
 
 ! write(*,*) "Loaded number of 3b No2b MEs: ", total_cnt
  if (total_cnt /= nelms) then
    write(*,*) "WARNING: number of 3b No2b MEs read from file is not equal to calculated number of 3b No2b MEs"
    stop
  end if

end subroutine read_3b_no2b_bin


function count_3bme_no2b(e_cut, e2max, e3max) result(r)
      
    integer, intent(in) :: e_cut, e2max, e3max
    integer(8) :: r
    integer :: nlev
    integer :: i1, l1, j1, e1
    integer :: i2, l2, j2, e2
    integer :: i3, l3, j3, e3
    integer :: i4, l4, j4, e4
    integer :: i5, l5, j5, e5
    integer :: i6, l6, j6, e6
    integer :: T12, T45, T, P123, P456, J

    r = 0

    nlev=(e_cut+1)*(e_cut+2)/2  ! number of single-particle levels
    do i1 = 1, nlev
      l1 = lev(i1)%l
      j1 = lev(i1)%j2
      e1 = lev(i1)%Nosc
      if(e1>e_cut) cycle
      do i2 = 1, i1
        l2 = lev(i2)%l
        j2 = lev(i2)%j2  
        e2 = lev(i2)%Nosc
        if(e1 + e2 > e2max) cycle
        do i3 = 1, nlev
          l3 = lev(i3)%l
          j3 = lev(i3)%j2
          e3 = lev(i3)%Nosc
          if(e1 + e3 > e2max) cycle
          if(e2 + e3 > e2max) cycle
          if(e1 + e2 + e3 > e3max) cycle
          P123 = (-1) ** (l1+l2+l3)
          do i4 = 1, i1
            l4 = lev(i4)%l
            j4 = lev(i4)%j2
            e4 = lev(i4)%Nosc

            do i5 = 1, i4
              l5 = lev(i5)%l
              j5 = lev(i5)%j2
              e5 = lev(i5)%Nosc
              if(e4 + e5 > e2max) cycle
              if(mod(l1+l2+l4+l5,2)==1) cycle

              do i6 = 1, nlev
                l6 = lev(i6)%l
                j6 = lev(i6)%j2
                e6 = lev(i6)%Nosc
                if(j3 /= j6) cycle
                if(l3 /= l6) cycle
                if(e4 + e6 > e2max) cycle
                if(e5 + e6 > e2max) cycle
                if(e4 + e5 + e6 > e3max) cycle

                P456 = (-1) ** (l4+l5+l6)

                if(P123 /= P456) cycle

                do J = max(abs(j1-j2), abs(j4-j5))/2, min((j1+j2), (j4+j5))/2
                  do T12 = 0, 1
                    do T45 = 0, 1
                      do T = max(abs(2*T12-1),abs(2*T45-1)),&
                            &min(   (2*T12+1),   (2*T45+1)), 2

                        r = r + 1

                      end do
                    end do
                  end do
                end do

              end do
            end do
          end do


        end do
      end do
    end do
    !write(*,*) "Number of MEs: ", r
  end function count_3bme_no2b

     SUBROUTINE sort2(n,arr,brr)
     INTEGER (kind=8) M,NSTACK
     PARAMETER (M=7,NSTACK=500)
     INTEGER (kind=8) n,i,ir,j,jstack,k,l,istack(NSTACK)
     INTEGER (kind=8) a, tempa
     REAL (kind=4) b, tempb
     INTEGER (kind=8) arr(n)
     REAL (kind=4) brr(n)
!Sorts an array arr(1:n) into ascending order using Quicksort, while making the corresponding rearrangement of the array brr(1:n).
     jstack=0
     l=1
     ir=n

1    if(ir-l.lt.M) then
        do j=l+1,ir
           a=arr(j)
           b=brr(j)
           do i=j-1,l,-1
              if(arr(i).le.a) goto 2
              arr(i+1)=arr(i)
              brr(i+1)=brr(i)
           enddo
           i=l-1
2          arr(i+1)=a
           brr(i+1)=b
        enddo
        if(jstack.eq.0) return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
     else
        k=(l+ir)/2
        tempa=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=tempa
        tempb=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=tempb
        if(arr(l).gt.arr(ir)) then
           tempa=arr(l)
           arr(l)=arr(ir)
           arr(ir)=tempa
           tempb=brr(l)
           brr(l)=brr(ir)
           brr(ir)=tempb
        endif
        if (arr(l+1).gt.arr(ir)) then
           tempa=arr(l+1)
           arr(l+1)=arr(ir)
           arr(ir)=tempa
           tempb=brr(l+1)
           brr(l+1)=brr(ir)
           brr(ir)=tempb
        endif
        if(arr(l).gt.arr(l+1)) then
           tempa=arr(l)
           arr(l)=arr(l+1)
           arr(l+1)=tempa
           tempb=brr(l)
           brr(l)=brr(l+1)
           brr(l+1)=tempb
        endif
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       continue
           i=i+1
        if(arr(i).lt.a) goto 3
4       continue
           j=j-1
        if(arr(j).gt.a) goto 4
        if(j.lt.i) goto 5
        tempa=arr(i)
        arr(i)=arr(j)
        arr(j)=tempa
        tempb=brr(i)
        brr(i)=brr(j)
        brr(j)=tempb
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK) pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l) then
           istack(jstack)=ir
           istack(jstack-1)=i
           ir=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        endif
     endif   
     goto 1
     END SUBROUTINE sort2
  
end module v3body_no2b_bin
