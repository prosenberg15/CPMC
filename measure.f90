!This subroutine measure the quanity. It is combined with FW and BK measure.

!---------------------------
!The free projection measure
!---------------------------
subroutine freep_measure(i,j)
use param
use lattice_param
use meas_param
use adET_param
use mpi_serial_param
use mc_loop_param
use one_meas_param
implicit none
integer,intent(IN)::i,j  !i is the sample,j is the i_local
integer::sitei,ib,jb

 call measure()

 kin_l(i,j)=dble(K_one)
 v_l(i,j)=dble(V_one)
 e_l(i,j)=dble(E_one)
 var_l(i,j)=var_one
 nu_l(i,j)=dble(nu_one)
 nd_l(i,j)=dble(nd_one)
 do sitei=1,Nbravais,1
    sisj_l(i,sitei,j)=dble(S_one(sitei))
 end do
 do sitei=1,2*Nsite,1
    cicj_l(i,sitei,j)=c_one(sitei)
 end do
 sig(j)=sig(j)+Sig_one
 absig(j)=absig(j)+AbSig_one

 if(j.lt.(max_ad*meastep)) then 
   ET=dble(E_one) !adjust the ET
   if(rank.eq.0) then
      write(*,*) "ADJUST ET FPMC:",ET
      write(*,*) ""
      write(*,*) ""
      write(*,*) ""
   end if
 end if
 !if(rank.eq.0) write(*,*) E_one
 !write(*,*) E_one;pause
end subroutine freep_measure


!--------------------------
!The cpmc and rcpmc measure
!--------------------------
subroutine add_measure(sample,i_release)
use param
use lattice_param
use model_param
use meas_param
use one_meas_param
use method_param
use mpi_serial_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer,intent(IN)::sample,i_release
integer::sitei,sitej,i_beta,ipair,ib,jb,alpha,beta,idir
real(kind=8)::t0,t2


!measure ground state correlations
  call measure()

!measure mixed estimator of energy
  call measure_energy() 

!measure dynamical properties
  if(Nbeta.gt.0)then

   call green_particles()
   call green_holes()
   !if two-body dynamics is desired
   if(I_twob.eq.1)then  
     call response_functions()
   endif
   !if the old method is desired
   if(I_onebf.eq.1)then
     call measure_one_body_dynamics()
   endif
  
  endif 


!accumulate estimations
  do idir=1,Dimen
    ZetaN_l(sample,idir,i_release)=ZetaN_l(sample,idir,i_release)+ZetaN_one(idir)
  enddo
  kin_l(sample,i_release)=kin_l(sample,i_release)+dble(K_one)
  v_l(sample,i_release)=v_l(sample,i_release)+dble(V_one)
  e_l(sample,i_release)=e_l(sample,i_release)+dble(E_one)
  eE_l(sample,i_release)=eE_l(sample,i_release)+dble(EE_one)
  var_l(sample,i_release)=var_l(sample,i_release)+var_one 
  nu_l(sample,i_release)=nu_l(sample,i_release)+dble(nu_one)
  nd_l(sample,i_release)=nd_l(sample,i_release)+dble(nd_one)
  if(I_obdm.eq.1)then
    do sitei=1,2*Nsite
      do sitej=1,2*Nsite
        obdm_l(sample,sitej,sitei,i_release)=obdm_l(sample,sitej,sitei,i_release)+Obdm_one(sitej,sitei)
      enddo
    enddo
  endif

!periodic systems
  if(ipinn.eq.0)then
    do sitei=1,Nbravais,1
     sisj_l(sample,sitei,i_release)=sisj_l(sample,sitei,i_release)+dble(S_one(sitei))
    end do
    do beta=1,Nbands
      do alpha=1,Nbands
        do sitei=1,Nbravais
          ninj_l(sample,sitei,alpha,beta,i_release)=ninj_l(sample,sitei,alpha,beta,i_release)+dble(N_one(sitei,alpha,beta))
        enddo
      enddo
    enddo
    do beta=1,Nbands
      do alpha=1,Nbands
        do sitei=1,Nbravais
          szsz_l(sample,sitei,alpha,beta,i_release)=szsz_l(sample,sitei,alpha,beta,i_release)+dble(Sz_one(sitei,alpha,beta))
        enddo
      enddo
    enddo
  elseif(ipinn.eq.1)then
!pinned systems
    do beta=1,2
      do alpha=1,2
        do sitei=1,Nsite
          nofr_l(sample,sitei,alpha,beta,i_release)=nofr_l(sample,sitei,alpha,beta,i_release)+Nr_one(sitei,alpha,beta)
        enddo
      enddo
    enddo
  else
    write(*,*)'Problem with ipinn '
    stop
  endif

!pairing correlations
  do ipair=1,Npair,1
   do sitei=1,Nbravais,1
     didj_l(sample,sitei,ipair,i_release)=didj_l(sample,sitei,ipair,i_release)+D_one(sitei,ipair)
   end do
  end do
 
  if(I_onebf.eq.1)then
    do i_beta=1,Nbeta,1
      do sitei=1,2*Nsite,1
        cicj_t_l(sample,sitei,i_beta,i_release)=cicj_t_l(sample,sitei,i_beta,i_release)+GP_one(sitei,i_beta)
        cicjh_t_l(sample,sitei,i_beta,i_release)=cicjh_t_l(sample,sitei,i_beta,i_release)+GH_one(sitei,i_beta)
      enddo
    enddo
  endif

!dynamical correlations
  if(Nbeta.gt.0)then
    if(I_twob.eq.1)then
      do i_beta=0,Nbeta,1
        do sitei=1,Nbravais,1
          nupnup_t_l(sample,sitei,i_beta,i_release)=nupnup_t_l(sample,sitei,i_beta,i_release)+nupnup_one(sitei,i_beta)
          ndnnup_t_l(sample,sitei,i_beta,i_release)=ndnnup_t_l(sample,sitei,i_beta,i_release)+ndnnup_one(sitei,i_beta)
        enddo
      enddo
    endif

    do i_beta=0,Nbeta,1
     GreenP_t_l(sample,i_beta,i_release)=GreenP_t_l(sample,i_beta,i_release)+GreenP_one(i_beta)
    enddo
    do i_beta=0,Nbeta,1
     GreenH_t_l(sample,i_beta,i_release)=GreenH_t_l(sample,i_beta,i_release)+GreenH_one(i_beta)
    enddo
  endif
  
  sig(i_release)=sig(i_release)+Sig_one
  absig(i_release)=absig(i_release)+AbSig_one
end subroutine add_measure




!-----------------------------------------------------------------
!The measurement give out the K V E S C sig the result is only one
!measurement, it should work together with other add subroutines
!-----------------------------------------------------------------
subroutine measure()
use param
use lattice_param
use phiT_param
use model_param
use mc_loop_param
use mpi_serial_param
use one_meas_param
use phi_param
use method_param
implicit none
#ifdef MPI
include "mpif.h"
#endif

!denominator and the numerator of the measurement
complex(kind=8)::Kinm,obdmm(2*Nsite,2*Nsite),nofrm(Nsite,2,2)
complex(kind=8)::nu_m,nd_m
complex(kind=8)::Vinm,sisjm(Nbravais),ninjm(Nbravais,Nbands,Nbands),didjm(Nbravais,Npair)
complex(kind=8)::szszm(Nbravais,Nbands,Nbands)
complex(kind=8)::ZetaN_m(Dimen)
complex(kind=8)::denominator  !Determinate the sum of weight

!The local files: <phL|O_local|phr>
complex(kind=8)::ZetaN_local(Dimen,Dtot)
complex(kind=8)::Kinm_local(Dtot)
complex(kind=8),allocatable::nofr_local(:,:,:,:)
complex(kind=8)::nu_local(Dtot),nd_local(Dtot)
complex(kind=8)::Vinm_local(Dtot)
complex(kind=8),allocatable::sisj_local(:,:)
complex(kind=8),allocatable::ninj_local(:,:,:,:)
complex(kind=8),allocatable::szsz_local(:,:,:,:)
complex(kind=8),allocatable::didj_local(:,:,:)

!For the tmp wavefunction file
complex(kind=8),allocatable::phL(:,:,:)
complex(kind=8)::coe(Dtot)
complex(kind=8)::phR(2*Nsite,Ntot)
complex(kind=8),allocatable::ovp_local(:,:,:)
complex(kind=8)::imp_local(Dtot),tot_local
complex(kind=8),allocatable::Amat_local(:,:,:),Bmat_local(:,:,:)
complex(kind=8),allocatable::Amat_pure(:,:)

!For the weight meas
complex(kind=8)::w_meas,tot_meas

!The cycle parameter or tmp parameter
integer::i,j,k,sitei,sitej,sitek,ipair,ib,jb,alpha,beta,idir
integer::i_deb,j_deb
logical::i_debug
complex(kind=8)::tmp,dummy


!Allocation of local estimators
if(ipinn.eq.0)then
  allocate(sisj_local(Nbravais,Dtot))
  allocate(ninj_local(Nbravais,Nbands,Nbands,Dtot))
  allocate(szsz_local(Nbravais,Nbands,Nbands,Dtot))
endif
if(ipinn.eq.1)allocate(nofr_local(Nsite,2,2,Dtot))
if(Npair.gt.0)allocate(didj_local(Nbravais,Npair,Dtot))
if(I_wavefun.eq.1)then
  allocate(phL(2*Nsite,Ntot,Dtot))
  allocate(Amat_local(2*Nsite,2*Nsite,Dtot))
elseif(I_wavefun.eq.2)then
  allocate(Amat_pure(2*Nsite,2*Nsite))
endif
allocate(ovp_local(Ntot,Ntot,Dtot))


!Set zero
ZetaN_m=zero
Kinm=zero
Vinm=zero
nu_m=zero
nd_m=zero
obdmm=zero
if(ipinn.eq.0)then
  sisjm=zero
  ninjm=zero
  szszm=zero
endif
if(ipinn.eq.1)nofrm=zero
if(Npair.gt.0)didjm=zero
denominator=zero


do i=1,Nwalkers,1

   if(I_wavefun.eq.1)then
     call     get_meas_array(i,phL    ,coe,phR,ovp_local,imp_local,tot_local,w_meas,tot_meas)
   elseif(I_wavefun.eq.2)then
     call bcs_get_meas_array(i,phR,coe,w_meas,tot_meas)
   endif

   !To get the Amat = < Psi_{trial}|(ci^+)(cj)| phR >
   if(I_wavefun.eq.1)then
     Amat_local=zero
     do k=1,Dtot,1
       call cal_Amat_withovlpinv_dc(phL(1,1,k),phR(1,1),ovp_local(1,1,k),Amat_local(1,1,k))
     enddo
   elseif(I_wavefun.eq.2)then
     call bcs_green_pure(i,phR,Amat_pure,.true.)
   endif

   !Kinetic energy

   if(I_wavefun.eq.1)then
     call get_klocal(Kinm_local,Amat_local)
   elseif(I_wavefun.eq.2)then
     call bcs_get_klocal(Kinm_local,Amat_pure)
   endif
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,Kinm_local,Kinm)

   !Potential energy

   if(I_wavefun.eq.1)then
     call get_vlocal(Vinm_local,Amat_local)
   elseif(I_wavefun.eq.2)then
     call bcs_get_vlocal(Vinm_local,Amat_pure)
   endif
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,Vinm_local,Vinm)


   !Particles localization

   ZetaN_local=zero
   do k=1,Dtot,1
      if(I_wavefun.eq.1)then
        call get_zetaNlocal(ZetaN_local(1,k),phL(1,1,k),phR(1,1))
      elseif(I_wavefun.eq.2)then     !we do not compute the Resta-Sorella for BCS
!        call bcs_get_zetaNlocal(ZetaN_local(1,k),phR(1,1)) !BCS case mixed estimator
        ZetaN_local=Zero
      endif
      do idir=1,Dimen
        ZetaN_local(idir,k)=ZetaN_local(idir,k)/imp_local(k)
      enddo
   enddo
   do idir=1,Dimen
        call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,          &
     &                     ZetaN_local(idir,1:Dtot),ZetaN_m(idir))
   enddo




   !number of system

   if(I_wavefun.eq.1)then
     call get_nlocal(nu_local,nd_local,Amat_local)
   elseif(I_wavefun.eq.2)then
     call bcs_get_nlocal(nu_local,nd_local,Amat_pure)
   endif
   !Get the numerator nud_m=w_meas*rx_i*<phL|Nup,Ndn|phR>/d_tmpi
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,nu_local,nu_m)
   call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,nd_local,nd_m)


   !The OBDM
     if(I_obdm.eq.1)then
       do sitei=1,2*Nsite
         do sitej=1,2*Nsite
           if(I_wavefun.eq.1)then
             call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,         &
     &                          Amat_local(sitej,sitei,1:Dtot),obdmm(sitej,sitei))
           elseif(I_wavefun.eq.2)then
             call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,         &
     &                          Amat_pure(sitej,sitei),obdmm(sitej,sitei))
           endif
         enddo
       enddo
     endif
  
   !Periodic systems
   if(ipinn.eq.0)then

!The spin correlation
     if(I_wavefun.eq.1)then
       call measure_sisj(Amat_local(1,1,1),sisj_local(1,1))
     elseif(I_wavefun.eq.2)then
       call bcs_measure_sisj(Amat_pure,sisj_local)
     endif
   
     do sitei=1,Nbravais,1
      call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,sisj_local(sitei,1:Dtot),sisjm(sitei))
     enddo

!The density-density correlation
     if(I_wavefun.eq.1)then
       call measure_ninj(Amat_local(1,1,1),ninj_local(1,1,1,1))
     elseif(I_wavefun.eq.2)then
       call bcs_measure_ninj(Amat_pure,ninj_local(1,1,1,1))
     endif

     do beta=1,Nbands
       do alpha=1,Nbands
         do sitei=1,Nbravais
           call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,                      &
     &                       ninj_local(sitei,alpha,beta,1:Dtot),ninjm(sitei,alpha,beta))
         enddo
       enddo
     enddo

!the Sz-Sz correlation
     if(I_wavefun.eq.1)then
       call measure_szsz(Amat_local(1,1,1),szsz_local(1,1,1,1))
     elseif(I_wavefun.eq.2)then
       call bcs_measure_szsz(Amat_pure,szsz_local(1,1,1,1))
     endif

     do beta=1,Nbands
       do alpha=1,Nbands
         do sitei=1,Nbravais
           call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,              &
     &                       szsz_local(sitei,alpha,beta,1:Dtot),szszm(sitei,alpha,beta))
         enddo
       enddo
     enddo

   endif

    
   !Pinned systems
   if(ipinn.eq.1)then
!The spin resolved and band resolved local Green function
     if(I_wavefun.eq.1)then
       call measure_nofr(Amat_local(1,1,1),nofr_local(1,1,1,1))
     elseif(I_wavefun.eq.2)then
       call bcs_measure_nofr(Amat_pure,nofr_local(1,1,1,1))
     endif

     do beta=1,2
       do alpha=1,2
         do sitei=1,Nsite
           call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,           &
     &                       nofr_local(sitei,alpha,beta,1:Dtot),nofrm(sitei,alpha,beta))
         enddo
       enddo
     enddo
   endif
     

   if(Npair.gt.0)then

   !The pairing correlation functions
     if(I_wavefun.eq.1)then
       call measure_pairing(Amat_local(1,1,1),didj_local(1,1,1))
     elseif(I_wavefun.eq.2)then
       call bcs_measure_pairing(Amat_pure,didj_local(1,1,1))
     endif

     do ipair=1,Npair
      do sitei=1,Nbravais,1
        call add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,didj_local(sitei,ipair,1:Dtot),didjm(sitei,ipair))
      enddo
     enddo
   endif

   call add_denominator(i,tot_meas,w_meas,denominator)  
end do


!Get the measurement value to physics_one=numerator/denominator
#ifdef MPI
call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(denominator,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
denominator=tmp
do idir=1,Dimen
  call MPI_ALLREDUCE(ZetaN_m(idir),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
  ZetaN_m(idir)=tmp
enddo
call MPI_ALLREDUCE(Kinm,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
Kinm=tmp
call MPI_ALLREDUCE(Vinm,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
Vinm=tmp
call MPI_ALLREDUCE(nu_m,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
nu_m=tmp
call MPI_ALLREDUCE(nd_m,tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
nd_m=tmp

if(I_obdm.eq.1)then
  do sitei=1,2*Nsite
    do sitej=1,2*Nsite
      call MPI_ALLREDUCE(obdmm(sitej,sitei),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM, &
     &                        MPI_COMM_WORLD,IERR)
      obdmm(sitej,sitei)=tmp
    enddo
  enddo
endif


if(ipinn.eq.0)then

  do sitei=1,Nbravais,1
    call MPI_ALLREDUCE(sisjm(sitei),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
    sisjm(sitei)=tmp
  enddo

  do beta=1,Nbands
    do alpha=1,Nbands
      do sitei=1,Nbravais
        call MPI_ALLREDUCE(ninjm(sitei,alpha,beta),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
        ninjm(sitei,alpha,beta)=tmp
      enddo
    enddo
  enddo

  do beta=1,Nbands
    do alpha=1,Nbands
      do sitei=1,Nbravais
        call MPI_ALLREDUCE(szszm(sitei,alpha,beta),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
        szszm(sitei,alpha,beta)=tmp
      enddo
    enddo
  enddo

endif


if(ipinn.eq.1)then
  do beta=1,2
     do alpha=1,2
       do sitei=1,Nsite
         call MPI_ALLREDUCE(nofrm(sitei,alpha,beta),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM, &
     &                      MPI_COMM_WORLD,IERR)
         nofrm(sitei,alpha,beta)=tmp
       enddo
     enddo
  enddo
endif

do ipair=1,Npair,1
  do sitei=1,Nbravais,1
    call MPI_ALLREDUCE(didjm(sitei,ipair),tmp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
    didjm(sitei,ipair)=tmp
  end do
end do
#endif


do idir=1,Dimen
  ZetaN_one(idir)=ZetaN_m(idir)/denominator
enddo


K_one=Kinm/denominator
V_one=Vinm/denominator
E_one=(Kinm+Vinm)/denominator
nu_one=nu_m/denominator
nd_one=nd_m/denominator

if(I_obdm.eq.1)then
  do sitei=1,2*Nsite
    do sitej=1,2*Nsite
      Obdm_one(sitej,sitei)=obdmm(sitej,sitei)/denominator
    enddo
  enddo
endif

if(ipinn.eq.0)then

  do sitei=1,Nbravais,1
    S_one(sitei)=sisjm(sitei)/denominator
  end do

  do beta=1,Nbands
     do alpha=1,Nbands
       do sitei=1,Nbravais
         N_one(sitei,alpha,beta)=ninjm(sitei,alpha,beta)/denominator
       enddo
     enddo
  enddo

  do beta=1,Nbands
     do alpha=1,Nbands
       do sitei=1,Nbravais
         Sz_one(sitei,alpha,beta)=szszm(sitei,alpha,beta)/denominator
       enddo
     enddo
  enddo
  
endif
if(ipinn.eq.1)then
  do beta=1,2
     do alpha=1,2
       do sitei=1,Nsite
         Nr_one(sitei,alpha,beta)=nofrm(sitei,alpha,beta)/denominator
       enddo
     enddo
  enddo
endif


do ipair=1,Npair,1
  do sitei=1,Nbravais,1
    D_one(sitei,ipair)=didjm(sitei,ipair)/denominator
  end do
end do


Sig_one=denominator
AbSig_one=abs(denominator)


deallocate(ovp_local)
if(allocated(Amat_local))deallocate(Amat_local)
if(allocated(Amat_pure)) deallocate(Amat_pure)
if(allocated(phL))       deallocate(phL)
if(ipinn.eq.0)deallocate(sisj_local,ninj_local,szsz_local)
if(ipinn.eq.1)deallocate(nofr_local)
if(Npair.gt.0)deallocate(didj_local)


end subroutine measure



subroutine bcs_get_meas_array(i,phR,coe,w_meas,tot_meas)
use param
use lattice_param
use model_param
use project_param
use phiT_param
use phi_param
use method_param
use caldet_module
use mpi_serial_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(OUT)::coe(Dtot),phR(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::w_meas,tot_meas
complex(kind=8)::ovp_tmp(Nspin(1),Nspin(1)),tmp,phtmp(2*Nsite,Ntot)
integer::j,k,sitei,sitej

!Get the phR
!exp_halfK = exp( - dt* T/2 )
if(.not.back_pro) then
  call k_to_ph_dc(exp_halfK,phi(1,1,i),phR)
else if(back_pro) then
  call k_to_ph_dc(exp_halfK,phi0(1,1,i),phR)
else
  write(*,*) "Something is wrong with back_pro:",back_pro
end if
coe(:)=one

!Get the tot_meas
if(.not.back_pro) then
  tot_meas=zero
  do k=1,Dtot,1
    if(Nzeta.eq.0)then
      call bcs_over_lap_dc(Fpairing,phR(1,1),ovp_tmp(1:Nspin(1),1:Nspin(1)))
    else
      call unp_bcs_over_lap_dc(Fpairing,Dunpaired,phR(1,1),ovp_tmp(1:Nspin(1),1:Nspin(1)))
    endif
    call inverse_d(ovp_tmp(1:Nspin(1),1:Nspin(1)),Nspin(1),tmp)  
    tot_meas=tot_meas+conjg(coe_multi(k))*tmp
  enddo
else if(back_pro) then
  call k_to_ph_dc(exp_halfK,phi(1,1,i),phtmp)
  tot_meas=zero
  do k=1,Dtot,1
     if(Nzeta.eq.0)then
       call bcs_over_lap_dc(Fpairing,phtmp(1,1),ovp_tmp(1:Nspin(1),1:Nspin(1)))
     else
       call unp_bcs_over_lap_dc(Fpairing,Dunpaired,phtmp(1,1),ovp_tmp(1:Nspin(1),1:Nspin(1)))
     endif
     call inverse_d(ovp_tmp(1:Nspin(1),1:Nspin(1)),Nspin(1),tmp)
     tot_meas=tot_meas+conjg(coe_multi(k))*tmp
  end do
else
  write(*,*) "Something is wrong with back_pro:",back_pro
end if


!Get the w_meas
if(crn.LT.0.d0) then  !CPMC
  w_meas=weight(i)*tot_meas/tot_imp(i)
else
  w_meas=weight(i)*one
end if

end subroutine bcs_get_meas_array






!------------------------------------------------------------
!We write the wave function after exp_mhalf_K, get the w_meas
!tot_meas, we also get the phL and coe, phR, the overlap 
!information between left and right wave functions.
!------------------------------------------------------------
subroutine get_meas_array(i,phL,coe,phR,ovp_local,imp_local,tot_local,w_meas,tot_meas)
use param
use lattice_param
use model_param
use project_param
use phiT_param
use phi_param
use method_param
use caldet_module
use mpi_serial_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(OUT)::phL(2*Nsite,Ntot,Dtot),coe(Dtot),phR(2*Nsite,Ntot)
complex(kind=8),intent(OUT)::ovp_local(Ntot,Ntot,Dtot)
complex(kind=8),intent(OUT)::imp_local(Dtot),tot_local,w_meas,tot_meas
complex(kind=8)::ovp_tmp(Ntot,Ntot),tmp,phtmp(2*Nsite,Ntot)

integer::j,k,sitei,sitej

!test
! phR=zero
! phL=zero
!end test


!Get the phR
if(.not.back_pro) then
  !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phi(1,1,i),2*Nsite,zero,phR,2*Nsite)
  call k_to_ph_dc(exp_halfK,phi(1,1,i),phR)
else if(back_pro) then
  !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phi0(1,1,i),2*Nsite,zero,phR,2*Nsite)
  call k_to_ph_dc(exp_halfK,phi0(1,1,i),phR)
else
  write(*,*) "Something is wrong with back_pro:",back_pro
end if


!Get the phL, coe
if(.not.back_pro) then
  !call zcopy(2*Nsite*Ntot*Dtot,phiT(1,1,1),1,phL(1,1,1),1)
  call copy_wf_T_dc(phiT(1,1,1),phL(1,1,1))
  call zcopy(Dtot,coe_multi(1),1,coe(1),1)
else if(back_pro) then
  call back_prog_phL(i,phL,coe)
else
  write(*,*) "Something is wrong with back_pro:",back_pro
end if




!Get the ovp_local,imp_local,tot_local
tot_local=zero
do k=1,Dtot,1
   !call deter_overlap(2*Nsite,Ntot,phL(1,1,k),phR(1,1),ovp_local(1,1,k))
   call over_lap_dc(phL(1,1,k),phR(1,1),ovp_local(1,1,k))

   !call caldet_dc(ovp_local(1:Ntot,1:Ntot,k),imp_local(k))
   !call inverse_dc(ovp_local(1:Ntot,1:Ntot,k))
   call inverse_d_dc(ovp_local(1:Ntot,1:Ntot,k),imp_local(k))
  ! if(rank.eq.0)write(*,*)'INITIAL OVERLAP ',imp_local(k)
   tot_local=tot_local+conjg(coe(k))*imp_local(k)
end do



!Get the tot_meas
if(.not.back_pro) then
  tot_meas=tot_local
else if(back_pro) then
  !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phi(1,1,i),2*Nsite,zero,phtmp,2*Nsite)
  call k_to_ph_dc(exp_halfK,phi(1,1,i),phtmp)
  tot_meas=zero
  do k=1,Dtot,1
     !call deter_overlap(2*Nsite,Ntot,phiT(1,1,k),phtmp(1,1),ovp_tmp(1,1))
     call over_lap_dc(phiT(1,1,k),phtmp(1,1),ovp_tmp(1,1))
     !call caldet(Ntot,ovp_tmp(1:Ntot,1:Ntot),tmp)
     call caldet_dc(ovp_tmp(1:Ntot,1:Ntot),tmp)
     tot_meas=tot_meas+conjg(coe_multi(k))*tmp
  end do
else
  write(*,*) "Something is wrong with back_pro:",back_pro
end if


!Get the w_meas
if(crn.LT.0.d0) then
  w_meas=weight(i)*tot_meas/tot_imp(i)
  !write(*,*)'tot_meas/tot_imp(i) ',tot_meas/tot_imp(i),'BBBBBBBBBBBBBBBBBBBBB'
else 
  w_meas=weight(i)*one
end if

end subroutine get_meas_array



subroutine bcs_get_zetaNlocal(ZetaN_local,ph_right)
use param
use lattice_param
use model_param
use phiT_param
implicit none
complex(kind=8),intent(OUT)::ZetaN_local(Dimen)
complex(kind=8),intent(IN)::ph_right(2*Nsite,Ntot)
integer::i,ib,isigma,k,ip,idir
real(kind=8)::rvec(Dimen)
complex(kind=8)::dummy
complex(kind=8)::phase(Dimen)
complex(kind=8)::ph(2*Nsite,Ntot)


do idir=1,Dimen

  phase(idir)=Xi*2.d0*Pi/Nl(idir)
  call copy_wf_dc(ph_right,ph)

  k=0
  do isigma=1,2
    do ib=1,Nbands
      do i=1,Nbravais
        k=k+1
        if(ib.eq.1)then
          rvec(idir)=coor(i,idir)
        elseif(ib.eq.2)then
          if(idir.eq.1)then
            rvec(idir)=coor(i,idir)+0.5d0
          else
            rvec(idir)=coor(i,idir)
          endif
        elseif(ib.eq.3)then
          if(idir.eq.1)then
            rvec(idir)=coor(i,idir)
          else
            rvec(idir)=coor(i,idir)+0.5d0
          endif
        endif
        do ip=1,Ntot
         ph(k,ip)=ph(k,ip)*exp(phase(idir)*rvec(idir))
        enddo
      enddo
    enddo
  enddo

  if(Nzeta.eq.0)then
    call bcs_imp_fun_dc(FPairing(1,1),ph,dummy)
  else
    call unp_bcs_imp_fun_dc(FPairing(1,1),DUnpaired(1,1),ph,dummy)
  endif

  ZetaN_local(idir)=dummy

enddo


end subroutine bcs_get_zetaNlocal


subroutine get_zetaNlocal(ZetaN_local,ph_left,ph_right)
use param
use lattice_param
use model_param

implicit none
complex(kind=8),intent(OUT)::ZetaN_local(Dimen)
complex(kind=8),intent(IN)::ph_left(2*Nsite,Ntot),ph_right(2*Nsite,Ntot)

integer::i,ib,isigma,k,ip,idir
real(kind=8)::rvec(Dimen)
complex(kind=8)::dummy
complex(kind=8)::phase(Dimen)
complex(kind=8)::ph(2*Nsite,Ntot),ovp(Ntot,Ntot)

do idir=1,Dimen

  !write(101,*)'idir = ',idir

  phase(idir)=Xi*2.d0*Pi/Nl(idir)

  !write(101,*)'phase(idir) = ',phase(idir)

  call copy_wf_dc(ph_right,ph)

  k=0
  do isigma=1,2
    do ib=1,Nbands
      do i=1,Nbravais
        k=k+1
        if(ib.eq.1)then
          rvec(idir)=coor(i,idir)
        elseif(ib.eq.2)then
          if(idir.eq.1)then
            rvec(idir)=coor(i,idir)+0.5d0
          else
            rvec(idir)=coor(i,idir)
          endif
        elseif(ib.eq.3)then
          if(idir.eq.1)then
            rvec(idir)=coor(i,idir)
          else
            rvec(idir)=coor(i,idir)+0.5d0
          endif
        endif
        do ip=1,Ntot
         ph(k,ip)=ph(k,ip)*exp(phase(idir)*rvec(idir))    
        enddo   
      enddo
    enddo
  enddo

  call over_lap_dc(ph_left(1,1),ph(1,1),ovp(1,1))
  call caldet_dc(ovp,dummy)  

  ZetaN_local(idir)=dummy

  !write(101,*)'dummy ',dummy
  !flush(101)

enddo

!stop 'deb'

end subroutine get_zetaNlocal


subroutine bcs_get_klocal(Kinm_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(OUT)::Kinm_local(Dtot)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite)
integer::k,sitei,sitej
Kinm_local=zero

k=1
do sitei=1,Nsite,1
        do sitej=1,Nsite,1
           Kinm_local(k)=Kinm_local(k)+Hzero(sitei,sitej)*Amat(sitei,sitej)
        enddo
enddo
do sitei=Nsite+1,2*Nsite,1
        do sitej=Nsite+1,2*Nsite,1
           Kinm_local(k)=Kinm_local(k)+Hzero(sitei,sitej)*Amat(sitei,sitej)
        enddo
enddo

end subroutine bcs_get_klocal


!--------------------------------------------
!Get the kinm_local(k)=<phR|K|phL>/ <phR|phL>
!the kinm_local is should from 1~Dtot
!--------------------------------------------
subroutine get_klocal(Kinm_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(OUT)::Kinm_local(Dtot)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
integer::k,sitei,sitej
Kinm_local=zero
if(dtype.EQ.'c') then
   do k=1,Dtot,1
     do sitei=1,2*Nsite,1
        do sitej=1,2*Nsite,1
           Kinm_local(k)=Kinm_local(k)+Hzero(sitei,sitej)*Amat(sitei,sitej,k)
        enddo
     enddo
   end do
else if(dtype.EQ.'d') then
   do k=1,Dtot,1
     do sitei=1,Nsite,1
        do sitej=1,Nsite,1
           Kinm_local(k)=Kinm_local(k)+Hzero(sitei,sitej)*Amat(sitei,sitej,k)
        enddo
     enddo
     do sitei=Nsite+1,2*Nsite,1
        do sitej=Nsite+1,2*Nsite,1
           Kinm_local(k)=Kinm_local(k)+Hzero(sitei,sitej)*Amat(sitei,sitej,k)
        enddo
     enddo
!     write(*,*)'Kinm_local(k) ',Kinm_local(k)
!     stop 'DEBUG'
   end do
   
end if
end subroutine get_klocal


subroutine bcs_get_vlocal(Vinm_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(OUT)::Vinm_local(Dtot)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite)!,Bmat(2*Nsite,2*Nsite,Dtot)
integer::k,sitei,sitej,ip,jp
complex(kind=8)::Vnormal !,Vanomalous
integer::i_deb,j_deb

k=1

if(Nbands.eq.1)then
  Vnormal=zero
  do sitei=1,Nsite,1
    Vnormal=Vnormal+Amat(sitei,sitei)*Amat(sitei+Nsite,sitei+Nsite)      &
     &             -Amat(sitei+Nsite,sitei)*Amat(sitei,sitei+Nsite)
  enddo
  Vinm_local(k)=dcmplx(onsitU)*Vnormal
elseif(Nbands.eq.3)then
   do sitei=1,Nbravais,1
          Vinm_local(k)=Vinm_local(k)+dcmplx(onsitUd)*                    &
    &                   (Amat(sitei,sitei)*Amat(sitei+Nsite,sitei+Nsite)  &
    &                   -Amat(sitei+Nsite,sitei)*Amat(sitei,sitei+Nsite))
   end do
   do sitei=Nbravais+1,3*Nbravais,1
          Vinm_local(k)=Vinm_local(k)+dcmplx(onsitUp)*                        &
    &                   (Amat(sitei,sitei)*Amat(sitei+Nsite,sitei+Nsite)      &
                        -Amat(sitei+Nsite,sitei)*Amat(sitei,sitei+Nsite))
   end do
endif


end subroutine bcs_get_vlocal


!--------------------------------------------
!Get the Vinm_local(k)=<phR|V|phL>/ <phR|phL>
!the Vinm_local is should from 1~Dtot
!--------------------------------------------
subroutine get_vlocal(Vinm_local,Amat)
use param
use phiT_param
use lattice_param
use model_param
implicit none
complex(kind=8),intent(OUT)::Vinm_local(Dtot)
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
integer::k,sitei,sitej
Vinm_local=zero


if(Nbands.eq.1)then
  if(dtype.EQ.'c') then
    do k=1,Dtot,1
       do sitei=1,Nsite,1
          Vinm_local(k)=Vinm_local(k)+Amat(sitei,sitei,k)*Amat(sitei+Nsite,sitei+Nsite,k) &
                     & -Amat(sitei,sitei+Nsite,k)*Amat(sitei+Nsite,sitei,k)
       end do
       Vinm_local(k)=dcmplx(onsitU)*Vinm_local(k)
    end do
  else if(dtype.EQ.'d') then
    do k=1,Dtot,1
       do sitei=1,Nsite,1
          Vinm_local(k)=Vinm_local(k)+Amat(sitei,sitei,k)*Amat(sitei+Nsite,sitei+Nsite,k)
       end do
       Vinm_local(k)=dcmplx(onsitU)*Vinm_local(k)
    end do
  end if
elseif(Nbands.eq.3)then
  if(dtype.EQ.'c') then
    do k=1,Dtot,1
       do sitei=1,Nbravais,1
          Vinm_local(k)=Vinm_local(k)+dcmplx(onsitUd)*                          &
                     &    (Amat(sitei,sitei,k)*Amat(sitei+Nsite,sitei+Nsite,k)  &
                     &    -Amat(sitei,sitei+Nsite,k)*Amat(sitei+Nsite,sitei,k))
       end do
       do sitei=Nbravais+1,3*Nbravais,1
         Vinm_local(k)=Vinm_local(k)+dcmplx(onsitUp)*                           &
                     &    (Amat(sitei,sitei,k)*Amat(sitei+Nsite,sitei+Nsite,k)  &
                     &    -Amat(sitei,sitei+Nsite,k)*Amat(sitei+Nsite,sitei,k))
       end do
    end do
  else if(dtype.EQ.'d') then
    do k=1,Dtot,1
       do sitei=1,Nbravais,1
          Vinm_local(k)=Vinm_local(k)+dcmplx(onsitUd)*                          &
                     &     (Amat(sitei,sitei,k)*Amat(sitei+Nsite,sitei+Nsite,k))
       end do
!       write(*,*)'Vinm_local(k) d-orbitals ',Vinm_local(k)
       do sitei=Nbravais+1,3*Nbravais,1
          Vinm_local(k)=Vinm_local(k)+dcmplx(onsitUp)*                          &
                     &     (Amat(sitei,sitei,k)*Amat(sitei+Nsite,sitei+Nsite,k))
       end do
!       write(*,*)'Vinm_local(k) total ',Vinm_local(k)
    end do
  end if
endif
end subroutine get_vlocal


!-------------------------------------------
!measure the nu_local and nd_local from Amat
!-------------------------------------------
subroutine bcs_get_nlocal(nu_local,nd_local,Amat)
use param
use lattice_param
use phiT_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite)
complex(kind=8),intent(OUT)::nu_local(Dtot),nd_local(Dtot)
integer::i,j,k,m,n

nu_local=zero;nd_local=zero
do k=1,Dtot,1
   do i=1,Nsite,1
      nu_local(k)=nu_local(k)+Amat(i,i)
      nd_local(k)=nd_local(k)+Amat(i+Nsite,i+Nsite)
   end do
end do
end subroutine bcs_get_nlocal

!-------------------------------------------
!measure the nu_local and nd_local from Amat
!-------------------------------------------
subroutine get_nlocal(nu_local,nd_local,Amat)
use param
use lattice_param
use phiT_param
implicit none
complex(kind=8),intent(IN)::Amat(2*Nsite,2*Nsite,Dtot)
complex(kind=8),intent(OUT)::nu_local(Dtot),nd_local(Dtot)
integer::i,j,k,m,n

nu_local=zero;nd_local=zero
do k=1,Dtot,1
   do i=1,Nsite,1
      nu_local(k)=nu_local(k)+Amat(i,i,k)
      nd_local(k)=nd_local(k)+Amat(i+Nsite,i+Nsite,k)
   end do
end do
end subroutine get_nlocal


!---------------------------------------------------------------
!We add the numberator to addm when measuring different quantity
!---------------------------------------------------------------
subroutine add_numerator(i,imp_local,tot_local,coe,w_meas,tot_meas,add_local,addm)
use param
use phiT_param
use phi_param
use method_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(IN)::imp_local(Dtot),tot_local,coe(Dtot)
complex(kind=8),intent(IN)::w_meas,tot_meas
complex(kind=8),intent(IN)::add_local(Dtot)
complex(kind=8),intent(INOUT)::addm
complex(kind=8)::rat_h
integer::k

if(I_wavefun.eq.1)then
!rat_h=<phL|hat{O}|phR>
  rat_h=zero
  do k=1,Dtot,1
     rat_h=rat_h+conjg(coe(k))*imp_local(k)*add_local(k)
  end do

!addm:
!rcpmc,fpmc: w_meas*rx_i*tot_meas*<phL|hat{O}|phR>/<phL|phR>
!cpmc: w_meas*rx_i*<phL|hat{O}|phR>/<phL|phR>
  if(crn.GT.0.d0) then
    addm=addm+w_meas*rx(i)*tot_meas*rat_h/tot_local
  else
    addm=addm+w_meas*rx(i)*rat_h/tot_local
  end if

elseif(I_wavefun.eq.2)then

  addm=addm+w_meas*rx(i)*add_local(1)

endif


!DEBUG
!write(*,*)'Measuring, walker  ',i
!write(*,*)
!write(*,*)'rat_h= ',rat_h
!write(*,*)'rx(i) = ',rx(i)
!write(*,*)'w_meas = ',w_meas
!write(*,*)'tot_imp(i) = ',tot_imp(i)
!write(*,*)'imp_local(1) ',imp_local(1)
!write(*,*)'tot_local = ',tot_local
!write(*,*)
!write(*,*)' = ',add_local(1)
!write(*,*)
!write(*,*)

end subroutine add_numerator



!-------------------------------------------
!We add the denominator to addm in measuring
!-------------------------------------------
subroutine add_denominator(i,tot_meas,w_meas,addm)
use phi_param
use method_param
implicit none
integer,intent(IN)::i
complex(kind=8),intent(IN)::tot_meas,w_meas
complex(kind=8),intent(INOUT)::addm
if(crn.GT.0.d0) then
   addm=addm+w_meas*rx(i)*tot_meas
else
   addm=addm+w_meas*rx(i)
end if
end subroutine add_denominator
