subroutine get_phi()
use param
use lattice_param
use model_param
use phiT_param
use phi_param
use method_param
use mc_loop_param
use mpi_serial_param
implicit none
integer::Nw
complex(kind=8),allocatable::phi_init(:,:,:)
complex(kind=8),allocatable::phi_read(:,:,:)
complex(kind=8),allocatable::phi_coe(:)
real(kind=8),allocatable::prob_init(:)
complex(kind=8),allocatable::r_init(:) !The complex phase
complex(kind=8)::dummy,dummy2

real(kind=8)::eps
integer::table(Nwalkers)
integer::i,j,k,j1,j2,npartu,npartd
real(kind=8)::norm
real(kind=8),external::dasum
complex(kind=8)::tot_tmp

call allocate_phi()
call initial_w_imp()


!---------------------------------------------------------------------
!Read wave function from file
!---------------------------------------------------------------------
if(PP.EQ.0) then   

   open(unit=10,file='phi_N.dat',status='old')
   read(10,*) Nw
   close(10)

#ifdef MPI
   if(Nw.EQ.MNwalkers) then
#else
  if(Nw.EQ.Nwalkers) then
#endif

     call read_phi_rank()
     if(rank.eq.0) write(*,*) "read phi from file directly."

#ifdef MPI
  else if(Nw.LT.MNwalkers) then
#else
  else if(Nw.LT.Nwalkers) then
#endif

     allocate(phi_init(2*Nsite,Ntot,Nw),phi_coe(Nw),prob_init(Nw),r_init(Nw))
     call read_phi(Nw,Nsite,phi_coe,phi_init)
     if(rank.eq.0) write(*,*) "read phi from file and pop contr."
     
     do k=1,Nw,1
        if(crn.GT.0.d0) then  !fp
           prob_init(k)=abs(phi_coe(k))
           r_init(k)=phi_coe(k)/abs(phi_coe(k)) !The first rx(i)
        else   !CPMC
           call get_tot_tmp(tot_tmp,phi_init(1,1,k))
           prob_init(k)=abs(phi_coe(k)*tot_tmp)
           r_init(k)=phi_coe(k)*tot_tmp/prob_init(k)
           !if(rank.eq.0) write(*,*) r_init(k)
        end if
     end do

     norm=1.d0/dasum(Nw,prob_init,1)
     call dscal(Nw,norm,prob_init,1)

     call initial_phi(Nw,prob_init,phi_init,table)

     do i=1,Nwalkers,1
        rx(i)=r_init(table(i))
     end do

     if(allocated(phi_init)) deallocate(phi_init)
     if(allocated(phi_coe)) deallocate(phi_coe)
     if(allocated(prob_init)) deallocate(prob_init)
     if(allocated(r_init)) deallocate(r_init)

  else

     if(rank.eq.0) then
        write(*,*) "Nw can not be larger than number of walkers",Nw
        call mystop
     end if

  end if

!----------------------------------------------------------------------
!Set phi from the trial wave function
!----------------------------------------------------------------------
else if (PP.EQ.1) then  

   if(I_wavefun.eq.1) then
     
      Nw=Dtot

#ifdef MPI
      if(Nw.EQ.MNwalkers) then
#else
      if(Nw.EQ.Nwalkers) then
#endif
         write(*,*) "Nw.EQ.(M)Nwalkers: we should read from the file."
         call mystop

#ifdef MPI
      else if(Nw.LT.MNwalkers) then
#else
      else if(Nw.LT.Nwalkers) then
#endif

      if(rank.eq.0) write(*,*) "read phi from phiT and pop contr."
      allocate(phi_init(2*Nsite,Ntot,Nw),phi_coe(Nw),prob_init(Nw),r_init(Nw))
      call zcopy(2*Nsite*Ntot*Nw,phiT(1,1,1),1,phi_init(1,1,1),1)
      call zcopy(Nw,coe_multi(1),1,phi_coe(1),1)      

      do k=1,Nw,1
         if(crn.GT.0.d0) then  !fp
            prob_init(k)=abs(phi_coe(k))
            r_init(k)=phi_coe(k)/abs(phi_coe(k)) !The first rx(i)
         else   !CPMC
            call get_tot_tmp(tot_tmp,phi_init(1,1,k))
            prob_init(k)=abs(phi_coe(k)*tot_tmp)
            r_init(k)=phi_coe(k)*tot_tmp/prob_init(k)
         end if
      end do

      norm=1.d0/dasum(Nw,prob_init,1)
      call dscal(Nw,norm,prob_init,1)
      
      call initial_phi(Nw,prob_init,phi_init,table)

      do i=1,Nwalkers,1
         rx(i)=r_init(table(i))
      end do

      if(allocated(phi_init)) deallocate(phi_init)
      if(allocated(phi_coe)) deallocate(phi_coe)
      if(allocated(prob_init)) deallocate(prob_init)
      if(allocated(r_init)) deallocate(r_init)
    
   else

      if(rank.eq.0) then
         write(*,*) "Nw can not be larger than number of walkers",Nw
         call mystop
      end if

   end if

    
 elseif(I_wavefun.eq.2)then

! BCS case, sampling the trial wave function

    if(rank.eq.0) then
       write(*,*) "Sampling BCS wave functions with ",Nwalkers, "   walkers"
    endif

    Nw=Nwalkers
    call sample_bcs(Nw)  
    !call sample_bcs_svd(Nw)    try to use SVD, still unclear how to use k --> -k

 endif

!----------------------------------------------------------------------
!Use non-interacting wavefunction as initial phi
!----------------------------------------------------------------------
else if(PP.eq.2) then

   Nw=1
   allocate(phi_init(2*Nsite,Ntot,Nw),phi_coe(Nw),prob_init(Nw),r_init(Nw))
   phi_coe(1)=cmplx(1.d0,0.d0)

   allocate(phi_read(2*Nsite,Ntot,Nw))
   call input_phi_d(phi_read)

   eps=1.d-3

   !normalize
   do i=1,Nw,1
      do k=1,Ntot,1
         dummy=zero
         do j=1,2*Nsite,1
            dummy=dummy+conjg(phi_read(j,k,i))*phi_read(j,k,i)
         end do
         do j=1,2*Nsite,1
            phi_read(j,k,i)=phi_read(j,k,i)/dsqrt(dble(dummy))
         end do
      enddo
   enddo

   if(dtype.eq.'c')then
      phi_init=phi_read
   elseif(dtype.eq.'d'.or.dtype.eq.'m')then
      do k=1,Nw,1
         npartu=0
         npartd=0
         do j2=1,Ntot
            dummy=zero
            do j1=1,Nsite,1
               dummy=dummy+conjg(phi_read(j1,j2,k))*phi_read(j1,j2,k)
            enddo
            if(abs(dummy-one).lt.eps)then
               npartu=npartu+1
               do j1=1,Nsite,1
                  phi_init(j1,npartu,k)=phi_read(j1,j2,k)
               enddo
               do j1=Nsite+1,2*Nsite,1
                  phi_init(j1,npartu,k)=zero
               enddo
            else
               dummy2=zero
               do j1=Nsite+1,2*Nsite,1
                  dummy2=dummy2+conjg(phi_read(j1,j2,k))*phi_read(j1,j2,k)
               enddo
               if(abs(dummy2-one).lt.eps)then
                  npartd=npartd+1
                  do j1=1,Nsite,1
                     phi_init(j1,Nspin(1)+npartd,k)=zero
                  enddo
                  do j1=Nsite+1,2*Nsite,1
                     phi_init(j1,Nspin(1)+npartd,k)=phi_read(j1,j2,k)
                  enddo
               endif
            endif
         enddo
      enddo
   endif

   deallocate(phi_read)

     
   if(crn.GT.0.d0) then  !fp
      prob_init(1)=abs(phi_coe(1))
      r_init(1)=phi_coe(1)/abs(phi_coe(1)) !The first rx(i)
   else   !CPMC
      call get_tot_tmp(tot_tmp,phi_init(1,1,1))
      prob_init(1)=abs(phi_coe(1)*tot_tmp)
      r_init(1)=phi_coe(1)*tot_tmp/prob_init(1)
      !if(rank.eq.0) write(*,*) r_init(k)
   end if

   if(rank.eq.0) write(*,*) 'imp function calculated'

   norm=1.d0/dasum(Nw,prob_init,1)
   call dscal(Nw,norm,prob_init,1)

   call initial_phi(Nw,prob_init,phi_init,table)

   rx(1)=r_init(table(1))

   if(allocated(phi_init)) deallocate(phi_init)
   if(allocated(phi_coe)) deallocate(phi_coe)
   if(allocated(prob_init)) deallocate(prob_init)
   if(allocated(r_init)) deallocate(r_init)
     
   !if(rank.eq.0) then 
   !   if(dtype.EQ.'c') then
   !      call input_phiT(2*Nsite,Ntot,Hzero,phi(1,1,1))
   !   else if(dtype.EQ.'d') then
   !      call input_phi_d()
   !   else
   !      if(rank.eq.0)write(*,*)'Please use c or d method'
   !      call mystop
   !   end if
   !end if
    
   !call get_tot_tmp(tot_tmp,phi(1,1,1))
   !rx(1)=tot_tmp/abs(tot_tmp)
!#ifdef MPI
!    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!    call MPI_BCAST(phi(1,1,1),1,phtype2,0,MPI_COMM_WORLD,IERR)
!#endif

 else

    write(*,*) "Something is wrong with PP input"
    call mystop
    
 end if

if(rank.eq.0) then
  write(*,*) "Nw of the phi:",Nw
end if


if(crn.LT.0)  call cal_imp_ovlap()


if(rank.eq.0) then
  write(*,*)
  write(*,*)
  write(*,*)'Initial importance function ',tot_imp
  write(*,*)
  write(*,*)
endif

end subroutine get_phi



!-----------------------------
!read the phi_coe and phi_init
!-----------------------------
subroutine read_phi(Nw,Nsite,phi_coe,phi_init)
use param
use io_module
use model_param
implicit none
integer,intent(IN)::Nw
integer,intent(IN)::Nsite
complex(kind=8),intent(OUT)::phi_coe(Nw)
complex(kind=8),intent(OUT)::phi_init(2*Nsite,Ntot,Nw)
integer::i,j,k,j1,j2,npartu,npartd
real(kind=8)::eps
complex(kind=8)::dummy,dummy2
complex(kind=8),dimension(:,:,:),allocatable::phi_read
character(len=300)::filen

eps=1.d-3

!Get phi_coe:
open(unit=10,file='phi_coe.dat',status='old')
  do i=1,Nw,1
     read(10,*) phi_coe(i)
  end do
close(10)
!Read phi_init
allocate(phi_read(2*Nsite,Ntot,Nw))
open(unit=10,file='phi.dat',status='old')
do i=1,Nw,1
  do k=1,Ntot,1
     do j=1,2*Nsite,1
        read(10,*) phi_read(j,k,i)
     end do
  end do
end do
!normalize
do i=1,Nw,1
 do k=1,Ntot,1
  dummy=zero
  do j=1,2*Nsite,1
    dummy=dummy+conjg(phi_read(j,k,i))*phi_read(j,k,i)
  end do
  do j=1,2*Nsite,1
    phi_read(j,k,i)=phi_read(j,k,i)/dsqrt(dble(dummy))
  end do
 enddo
enddo
close(10)


if(dtype.eq.'c')then
  phi_init=phi_read
elseif(dtype.eq.'d'.or.dtype.eq.'m')then
  do k=1,Nw,1
    npartu=0
    npartd=0
    do j2=1,Ntot
      dummy=zero
      do j1=1,Nsite,1
        dummy=dummy+conjg(phi_read(j1,j2,k))*phi_read(j1,j2,k)
      enddo
      if(abs(dummy-one).lt.eps)then
        npartu=npartu+1
        do j1=1,Nsite,1
           phi_init(j1,npartu,k)=phi_read(j1,j2,k)
        enddo
        do j1=Nsite+1,2*Nsite,1
           phi_init(j1,npartu,k)=zero
        enddo
      else
        dummy2=zero
        do j1=Nsite+1,2*Nsite,1
          dummy2=dummy2+conjg(phi_read(j1,j2,k))*phi_read(j1,j2,k)
        enddo
        if(abs(dummy2-one).lt.eps)then
          npartd=npartd+1
          do j1=1,Nsite,1
            phi_init(j1,Nspin(1)+npartd,k)=zero
          enddo
          do j1=Nsite+1,2*Nsite,1
            phi_init(j1,Nspin(1)+npartd,k)=phi_read(j1,j2,k)
          enddo
        endif
      endif
    enddo
  enddo
endif

deallocate(phi_read)

end subroutine read_phi



!---------------------------------------------------
!write the coe and phi to the file according to rank
!---------------------------------------------------
subroutine write_phi_rank()
use param
use mpi_serial_param
use io_module
use mc_loop_param
use phi_param
use lattice_param
use model_param
use method_param
use project_param
implicit none
character(len=300)::phi_name,coe_name
complex(kind=8)::phtmp(2*Nsite,Ntot)
integer::i,j,k
call createFileName(phi_name,'phi')
call appendBaseName(phi_name,'_',rank)
call appendBaseName(phi_name,'.dat')
call createFileName(coe_name,'coe')
call appendBaseName(coe_name,'_',rank)
call appendBaseName(coe_name,'.dat')
!call openUnit(phi_name,10,'R')
!call openUnit(coe_name,20,'R')
call openUnit(phi_name,10,'R')
call openUnit(coe_name,20,'R')
do i=1,Nwalkers,1

   !call zgemm('N','N',2*Nsite,Ntot,2*Nsite,one,exp_halfK,2*Nsite,phi(1,1,i),2*Nsite,zero,phtmp,2*Nsite)
   call k_to_ph_dc(exp_halfK,phi(1,1,i),phtmp(1,1))

   do k=1,Ntot,1
      do j=1,2*Nsite,1
         write(10,*) phtmp(j,k)
      end do
   end do
   !write(10) phtmp(1:2*Nsite,1:Ntot)

   if(crn.GT.0.d0) then !fp
     write(20,*) rx(i)*weight(i)
   else
     write(20,*) rx(i)*weight(i)/tot_imp(i)
   end if

end do
close(10)
close(20)
end subroutine write_phi_rank



!--------------------------------------------------
!read the coe and phi to the file according to rank
!--------------------------------------------------
subroutine read_phi_rank()
use mpi_serial_param
use io_module
use mc_loop_param
use phi_param
use lattice_param
use model_param
use method_param
implicit none
character(len=300)::phi_name,coe_name
complex(kind=8)::tmp,tot_tmp
integer::i,j,k

call createFileName(phi_name,'phi')
call appendBaseName(phi_name,'_',rank)
call appendBaseName(phi_name,'.dat')
call createFileName(coe_name,'coe')
call appendBaseName(coe_name,'_',rank)
call appendBaseName(coe_name,'.dat')
call openUnit(phi_name,10,'B')
call openUnit(coe_name,20,'B')
do i=1,Nwalkers,1

   read(10) phi(1:2*Nsite,1:Ntot,i)

   read(20) tmp
   if(crn.GT.0.d0) then !fp
     weight(i)=abs(tmp)
     rx(i)=tmp/weight(i)
   else  !CPMC
     call get_tot_tmp(tot_tmp,phi(1,1,i))
     weight(i)=abs(tmp*tot_tmp)
     rx(i)=tmp*tot_tmp/weight(i)
   end if

   !set dlogw
   if(weight(i).GT.0.d0) then
     dlogw(i)=dlog(weight(i))
   else if(weight(i).EQ.0.d0) then
     dlogw(i)=-1d100
   else
     write(*,*) "Something is wrong with weight init.",weight(i),i
     call mystop
   end if

end do
close(10)
close(20)
end subroutine read_phi_rank





!-----------------------
!We initial weight dlogw
!-----------------------
subroutine initial_w_imp
use param
use phi_param
implicit none
weight=1.d0   
dlogw=0.d0   !All the dlogw here is dlog(weight)
end subroutine initial_w_imp




!---------------------------------------------------------------
!We initial phi and inverse overlap function in this subroutine.
!---------------------------------------------------------------
subroutine initial_phi(Nw,prob_init,phi_init,table)
use param
use model_param
use lattice_param
use phi_param
use mc_loop_param
use rand_num
use method_param
use mpi_serial_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer,intent(IN)::Nw
real(kind=8),intent(IN)::prob_init(Nw)
complex(kind=8),intent(IN)::phi_init(2*Nsite,Ntot,Nw)
integer,intent(OUT)::table(Nwalkers)
integer,allocatable::Mtable(:)
real(kind=8)::eta!Get the probability of initial wave function
integer::i,j,k,npartu,npartd
integer::j1,j2
real(kind=8)::eps
complex(kind=8)::dummy,dummy2

eps=1.d-3

#ifdef MPI
 if(rank.eq.0) allocate(Mtable(MNwalkers))
#endif

#ifdef MPI
   if(rank.eq.0) call distri_p(Nw,prob_init,MNwalkers,Mtable)
#else
   call distri_p(Nw,prob_init,Nwalkers,table)
#endif


#ifdef MPI
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_SCATTER(Mtable(1),Nwalkers,MPI_INTEGER,table(1),Nwalkers,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
#endif

 do i=1,Nwalkers,1
   !---------------------------------------------------------
   !Set the initial wavefunction proportional to prob_init(k)
   !---------------------------------------------------------
   eta=1d-3
   !eta=100.d0 !test
   if(Nw.eq.1) eta=0.d0

   !Only use decoupled term here, for we do not want mixed term here.
   !Be careful here!
   if(dtype.EQ.'c') then
   
     do j1=1,2*Nsite,1
        do j2=1,Ntot,1
           phi(j1,j2,i)=phi_init(j1,j2,table(i))+eta*(2.d0*rndm()-1.d0)
        end do
     end do
   
   else if(dtype.EQ.'d'.or.dtype.EQ.'m') then

     !Just put in the up and dn term
     do j1=1,Nsite,1
        do j2=1,Nspin(1),1
           phi(j1,j2,i)=phi_init(j1,j2,table(i))+eta*(2.d0*rndm()-1.d0)
        end do
        !test
        do j2=Nspin(1)+1,Ntot,1
        !   write(*,*) rndm()
           phi(j1,j2,i)=zero
        end do
        !end test
     end do
     
     
     do j1=Nsite+1,2*Nsite,1
        !test
        do j2=1,Nspin(1),1
        !   write(*,*) rndm()
           phi(j1,j2,i)=zero
        end do
        !end test
        do j2=Nspin(1)+1,Ntot,1
           phi(j1,j2,i)=phi_init(j1,j2,table(i))+eta*(2.d0*rndm()-1.d0)
        end do
     end do

  end if
 end do


if(rank.eq.0)then
open(2,file='psi_init.info',status='unknown')
do i=1,1,1
  do k=1,Ntot,1
     do j=1,2*Nsite,1
        write(2,*) phi(j,k,i)
     end do
  end do
end do
close(2)
endif

!stop 'DEBUG'

#ifdef MPI
 if(rank.eq.0) deallocate(Mtable)
#endif

end subroutine initial_phi






!----------------------------------------------------------------
!This subroutine generate n walkers by the distribution of p(1:m)
!Similar to the subroutine reconfiguration, while here m/=n+++++
!----------------------------------------------------------------
subroutine distri_p(m,p,n,table)
use rand_num
implicit none
integer,intent(IN)::m,n
real(kind=8),intent(IN)::p(m)
integer,intent(OUT)::table(n)
real(kind=8)::wsum,z
integer::i,j
wsum=0.d0
do i=1,m,1
   wsum=wsum+p(i)
end do
if(abs(wsum-1.d0).GT.1d-10) then
  write(*,*) "p input is not normalized!"
  call mystop
end if

j=1
wsum=p(1)
do i=1,n,1
   z=dble((i-1)+rndm())/dble(n)
   do while(z>wsum)
      j=j+1
      wsum=wsum+p(j)
   enddo
   table(i)=j
enddo

end subroutine distri_p


!-------------------------------
!Get the tot_tmp=<phiT|phi_init>
!-------------------------------
subroutine get_tot_tmp(tot_tmp,phi_init)
use param
use lattice_param
use model_param
use phiT_param
implicit none
complex(kind=8),intent(OUT)::tot_tmp
complex(kind=8),intent(IN)::phi_init(2*Nsite,Ntot)
complex(kind=8)::imp
integer::m,i,k,j
tot_tmp=zero
do m=1,Dtot,1
   if(I_wavefun.eq.1)then
     call imp_fun_dc(phiT(1,1,m),phi_init(1,1),imp)
   elseif(I_wavefun.eq.2)then
     if(Nzeta.eq.0)then
       call bcs_imp_fun_dc(FPairing(1,1),phi_init(1,1),imp)
     else
       call unp_bcs_imp_fun_dc(FPairing(1,1),DUnpaired(1,1),phi_init(1,1),imp)
     endif
   endif
   if(abs(imp).LT.1d-18) then
      write(*,*)
      write(*,*) "Small overlap phiT with phi_init",m
      write(*,*)
      call mystop      
   end if
   tot_tmp=tot_tmp+conjg(coe_multi(m))*imp
end do

end subroutine get_tot_tmp


subroutine sample_bcs_svd(Nw)
use param
use lattice_param
use model_param
use phiT_param
use phi_param
use rand_num
use mpi_serial_param
implicit none

integer,intent(IN)::Nw

integer::iptrEk(Nsite),iorboccp(Nspin(2))
real(kind=8)   ,dimension(:)  ,allocatable::S
complex(kind=8),dimension(:,:),allocatable::A,B

integer::ioccnow,iorb_fr_indx,iorb_fr,iorb_to
integer::M,k,n,i,j,ip,iw,is,ir,imetrp,metrop
real(kind=8)::metrotest,NofS
complex(kind=8)::tot_tmp

real(kind=8), allocatable    ::  rwork(:)
complex(kind=8), allocatable ::  work(:)
INTEGER                      ::  INFO,LWORK
CHARACTER                    ::  JOBU, JOBVT

if(Nspin(1).ne.Nspin(2))then
  write(*,*)'ERROR: BCS sampling available only for the fully paired case'
  write(*,*)'Please read initial state from file'
  stop
endif

M=Nsite

allocate(S(M))
allocate(A(M,M))
allocate(B(M,M))

JOBU='A'
JOBVT='A'
LWORK=10*M
allocate(work(LWORK))
allocate(rwork(5*M))

call     ZGESVD( JOBU, JOBVT, M, M, FPairing, M, S, A, M, B, M,         &
     &           WORK, LWORK, RWORK, INFO )


if(rank.eq.0)then
  open(2,file='Singular_Values_Of_F.info',status='unknown')
    do i=1,M
      write(2,*)i,S(i)
    enddo
  close(2)
endif

!stop 'DEB'

metrop=Ntot*1000

do i=1,M
  iptrEk(i)=i
enddo


iw=1
do ip=1,Nspin(2) !number of paired orbitals
  iorboccp(ip)=iptrEk(ip)
enddo
phi(:,:,iw)=Zero
do ip=1,Nspin(1)
  is=iorboccp(ip) ! index of ip-th orbital
  do ir=1,Nsite
    phi(ir,      ip         ,iw)=A(ir,is)   !spin up
    phi(ir+Nsite,ip+Nspin(1),iw)=B(is,ir)   !spin dn
  enddo
enddo
call get_tot_tmp(tot_tmp,phi(1,1,iw))
rx(iw)=tot_tmp/abs(tot_tmp)


do iw=2,Nw
  do ip=1,Nspin(2) !number of paired orbitals
    iorboccp(ip)=iptrEk(ip)
  enddo
  do imetrp=1,metrop
    iorb_fr_indx=floor(Nspin(2)*rndm()+1.0)
    iorb_fr=iorboccp(iorb_fr_indx)
 7  iorb_to=floor(Nsite*rndm()+1.0)
    do ioccnow=1,Nspin(2)
      if(iorb_to.eq.iorboccp(ioccnow))then
        go to 7
      endif
    enddo
    metrotest=(dble(S(iorb_to))/dble(S(iorb_fr)))**2
    if(rndm().le.metrotest)then
      iorboccp(iorb_fr_indx)=iorb_to
    endif
  enddo

  phi(:,:,iw)=Zero
  do ip=1,Nspin(1)
    is=iptrEk(ip)     !iorboccp(ip) ! index of ip-th orbital
    do ir=1,Nsite
      phi(ir,      ip         ,iw)=A(ir,is)   !spin up
      phi(ir+Nsite,ip+Nspin(1),iw)=B(is,ir)   !spin dn
    enddo
  enddo

  call get_tot_tmp(tot_tmp,phi(1,1,iw))
  rx(iw)=tot_tmp/abs(tot_tmp)

enddo

deallocate(S)
deallocate(work)
deallocate(rwork)
deallocate(A)
deallocate(B)

end subroutine sample_bcs_svd


subroutine sample_bcs(Nw)
use param
use lattice_param
use model_param
use phiT_param
use phi_param
use rand_num
use mpi_serial_param
implicit none

integer,intent(IN)::Nw

complex(kind=8),dimension(:),allocatable::effek
real(kind=8),dimension(:),allocatable::eofk
real(kind=8)::kl(Dimen),kpth(Dimen)
integer::k,n,i,j,ip,iw,is,ir,imetrp,metrop
integer::ioccnow,iorb_fr_indx,iorb_fr,iorb_to
integer::iptrEk(Nsite),iorboccp(Nspin(2))
real(kind=8)::x,y,ph,rmrp,metrotest
complex(kind=8)::tot_tmp
complex(kind=8)::eband(Nbands,Nbands)
real(kind=8)::eigenvalues(Nbands)

if(Nspin(1).ne.Nspin(2))then
  write(*,*)'ERROR: BCS sampling available only for the fully paired case'
  write(*,*)'Please read initial state from file'
  stop
endif

allocate(effek(Nsite))
allocate(eofk(Nsite))

metrop=Ntot*1000

effek=zero
do k=1,Nbravais,1

  do n=1,Dimen,1
    kl(n)=2.d0*Pi*dble(coor(k,n)-1)/dble(Nl(n))
    if(kl(n).gt.pi)kl(n)=kl(n)-2.d0*pi
    kpth(n)=kl(n)+2.d0*Pi*kbound(n)/dble(Nl(n))
  enddo

!get dispersion relation

  if(Nbands.eq.1)then
    eofk(k)=0.d0
    do n=1,Dimen,1
      eofk(k)=eofk(k)-2.d0*dble(t1)*dcos(kpth(n))
    end do
  elseif(Nbands.eq.3)then
    eband(1,1)= epsd
    eband(1,2)=-2.d0*Xi*tpd*dsin(0.5d0*kpth(1))
    eband(1,3)= 2.d0*Xi*tpd*dsin(0.5d0*kpth(2))
    eband(2,1)= 2.d0*Xi*tpd*dsin(0.5d0*kpth(1))
    eband(2,2)= epsp 
    eband(2,3)= 4.d0*tpp*dsin(0.5d0*kpth(1))*dsin(0.5d0*kpth(2))
    eband(3,1)=-2.d0*Xi*tpd*dsin(0.5d0*kpth(2))
    eband(3,2)= 4.d0*tpp*dsin(0.5d0*kpth(1))*dsin(0.5d0*kpth(2))
    eband(3,3)= epsp

    call eigen(eband,Nbands,eigenvalues)

    eofk(k)           =eigenvalues(1)
    eofk(k+Nbravais)  =eigenvalues(2)
    eofk(k+2*Nbravais)=eigenvalues(3)
  endif

!get fourier transform of Fpairing

  do i=1,Nbravais
    do j=1,Nbravais
      ph=0.d0
      do n=1,Dimen,1
        rmrp=dble(coor(i,n)-coor(j,n))
        ph=ph+rmrp*kl(n)
      enddo
      effek(k)=effek(k)+exp(Xi*ph)*FPairing(i,j)/dble(Nbravais)
      if(Nbands.eq.3)then
        effek(k+Nbravais)  =effek(k+Nbravais)  +exp(Xi*ph)*FPairing(i+Nbravais,j+Nbravais)/dble(Nbravais)
        effek(k+2*Nbravais)=effek(k+2*Nbravais)+exp(Xi*ph)*FPairing(i+2*Nbravais,j+2*Nbravais)/dble(Nbravais)
      endif
    enddo
  enddo

enddo


if(rank.eq.0)then
  open(2,file='Fourier_Transform_Of_F.info',status='unknown')
    do k=1,Nbravais

      if(Nbands.eq.1)then 
        write(2,*)effek(k)
        flush(2)
      elseif(Nbands.eq.3)then
        write(2,*)effek(k),effek(k+Nbravais),effek(k+2*Nbravais)
        flush(2)
      endif 
    enddo
  close(2)
  open(2,file='Free_Particles_Dispersion.info',status='unknown')
    do k=1,Nbravais,1
      if(Nbands.eq.1)then
        write(2,*)eofk(k)
        flush(2)
      elseif(Nbands.eq.3)then
        write(2,*)eofk(k),eofk(k+Nbravais),eofk(k+2*Nbravais)
        flush(2)
      endif
    enddo
  close(2)
endif


call hsort(Nsite,eofk,iptrEk)

do iw=1,Nw
  do ip=1,Nspin(2) !number of paired orbitals
    iorboccp(ip)=iptrEk(ip)
  enddo


!  do imetrp=1,metrop
!    iorb_fr_indx=floor(Nspin(2)*rndm()+1.0)
!    iorb_fr=iorboccp(iorb_fr_indx)
! 1  iorb_to=floor(Nsite*rndm()+1.0)
!    do ioccnow=1,Nspin(2)
!      if(iorb_to.eq.iorboccp(ioccnow))then
!        go to 1
!      endif
!    enddo
!    metrotest=(dble(effek(iorb_to))/dble(effek(iorb_fr)))**2
!    if(rndm().le.metrotest)then
!      iorboccp(iorb_fr_indx)=iorb_to
!    endif
!  enddo

  phi(:,:,iw)=Zero
  do ip=1,Nspin(1)
    is=iorboccp(ip) ! index of ip-th orbital
    do n=1,Dimen,1
      kl(n)=2.d0*Pi*dble(coor(is,n)-1)/dble(Nl(n))
      if(kl(n).gt.pi)kl(n)=kl(n)-2.d0*pi
    enddo   
    do ir=1,Nbravais
      ph=0.d0
      do n=1,Dimen,1
        rmrp=dble(coor(ir,n)-1)
        ph=ph+rmrp*kl(n)
      enddo
      if(Nbands.eq.1)then
        phi(ir,      ip         ,iw)=exp( Xi*ph)/dsqrt(dble(Nbravais))   !spin up
        phi(ir+Nsite,ip+Nspin(1),iw)=conjg(phi(ir,ip,iw))                !spin dn
      elseif(Nbands.eq.3)then
        phi(ir           ,ip,iw)=eband(1,1)*exp( Xi*ph)/dsqrt(dble(Nbravais)) !d
        phi(ir+  Nbravais,ip,iw)=eband(2,1)*exp( Xi*ph)/dsqrt(dble(Nbravais)) !px
        phi(ir+2*Nbravais,ip,iw)=eband(3,1)*exp( Xi*ph)/dsqrt(dble(Nbravais)) !py
        phi(ir           +Nsite,ip+Nspin(1),iw)=phi(ir           ,ip,iw)
        phi(ir+  Nbravais+Nsite,ip+Nspin(1),iw)=phi(ir+  Nbravais,ip,iw)
        phi(ir+2*Nbravais+Nsite,ip+Nspin(1),iw)=phi(ir+2*Nbravais,ip,iw)
      endif
    enddo   
  enddo

  call get_tot_tmp(tot_tmp,phi(1,1,iw))
  rx(iw)=tot_tmp/abs(tot_tmp)

enddo


deallocate(effek)
deallocate(eofk)

end subroutine sample_bcs



      SUBROUTINE HSORT(N,A,IND)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(*),IND(*)
      DO 1 J=1,N
    1 IND(J)=J
      IF(N.LE.1)RETURN
      L=N/2+1
      IR=N
   10 CONTINUE
      IF(L.GT.1)THEN
       L=L-1
       INDXT=IND(L)
       Q=A(INDXT)
      ELSE
       INDXT=IND(IR)
       Q=A(INDXT)
       IND(IR)=IND(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
        IND(1)=INDXT
        RETURN
       ENDIF
      ENDIF
      I=L
      J=L+L
   20 IF(J.LE.IR)THEN
      IF(J.LT.IR)THEN
       IF(A(IND(J)).LT.A(IND(J+1)))J=J+1
      ENDIF
      IF(Q.LT.A(IND(J)))THEN
       IND(I)=IND(J)
       I=J
       J=J+J
      ELSE
       J=IR+1
      ENDIF
      GOTO20
      ENDIF
      IND(I)=INDXT
      GOTO10
      END

