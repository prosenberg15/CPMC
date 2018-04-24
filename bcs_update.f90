subroutine BCS_U_phi(i,ovlp)
use param
use rand_num
use mc_loop_param
use phi_param
use lattice_param
use model_param
use phiT_param
use method_param
use project_param
implicit none
#ifdef MPI
include "mpif.h"
#endif
integer,intent(IN)::i
complex(kind=8),intent(INOUT)::ovlp(Ntot,Ntot,Dtot)
complex(kind=8)::avn_up(Nsite),avn_dn(Nsite)
real(kind=8)::p_tilde(2)
integer::j,aux,k
real(kind=8)::norm,x,bias,cap,tmp
complex(kind=8)::tot_tmp,explr_up,explr_dn
real(kind=8)::dummy,sommau,sommad
real(kind=8)::t0,t2,t3

logical first
data first /.true./

dummy=1.d0

if(weight(i).le.0.d0) return


cap=1.d0/(dt**(1.d0/4.d0))

tot_tmp=tot_imp(i)


call compute_FB(ovlp(1:Nspin(1),1:Nspin(1),1),i,avn_up,avn_dn)

do j=1,Nsite,1
  if(dcp.eq.'S') then
! 1 --- x=-1
! 2 --- x=+1
    bias=dble(gama*(avn_up(j)-avn_dn(j)))
    p_tilde(1)=max(0.5d0*(1.d0-bias),0.d0)
    p_tilde(2)=max(0.5d0*(1.d0+bias),0.d0)
  else if(dcp.eq.'C') then
    bias=dble(avn_up(j)+avn_dn(j)-one)
    if(dabs(bias).ge.cap)then
      bias=cap
    endif
    p_tilde(1)=0.5d0*(1.d0-dble(gama)*bias)
    if(p_tilde(1).lt.0.d0)p_tilde(1)=0.d0
    p_tilde(2)=0.5d0*(1.d0+dble(gama)*bias)
    if(p_tilde(2).lt.0.d0)p_tilde(2)=0.d0
  endif

!Choose the auxiliary field for site j, sampling p_tilde
  if(p_tilde(1).ge.rndm()) then
    aux=1
    x=1.d0
  else
    aux=2
    x=2.d0
  end if

  if(back_pro) then
    back_store(j,i_back,i)=x
  end if

!sampled one particle propagator for site j
  if(Nbands.eq.1)then
    explr_up=expln_up(aux)
    explr_dn=expln_dn(aux)
  elseif(Nbands.eq.3)then
    if(j.le.Nbravais)then
      explr_up=explnd_up(aux)
      explr_dn=explnd_dn(aux)
    elseif(j.le.2*Nbravais)then
      explr_up=explnx_up(aux)
      explr_dn=explnx_dn(aux)
    else
      explr_up=explny_up(aux)
      explr_dn=explny_dn(aux)
    endif
  endif

!update of the determinant
  do k=1,Nspin(1),1
    phi(j,k,i)=phi(j,k,i)*(explr_up+one)
  end do
  do k=Nspin(1)+1,Ntot,1
    phi(j+Nsite,k,i)=phi(j+Nsite,k,i)*(explr_dn+one)
  end do

!contribution p(x_j)/ptilde(x_j) to the weight of the walker
  weight(i)=weight(i)*0.5d0/p_tilde(aux)  !factor 0.5 = p(x) from HS
  dlogw(i)=dlogw(i)+dlog(0.5d0/p_tilde(aux))

  if(dcp.eq.'C') then
    dummy=dummy*exp(dt*OnsitU/2.d0-dble(gama)*dble((-1)**aux))
    if(Nbands.eq.3)then
      write(*,*)'ERROR: please use spin decomposition for three band'
      stop
    endif
  endif
 
end do

call get_imp_inv(i,ovlp)
if(dble(tot_imp(i)/tot_tmp).LE.0.d0) then
   weight(i)=0.d0
   dlogw(i)=-1d100
   return
end if

!update the weight of the walker i
weight(i)=weight(i)*dble(dummy*tot_imp(i)/tot_tmp)  
dlogw(i)=dlogw(i)+dlog(dble(dummy*tot_imp(i)/tot_tmp))

if(weight(i).lt.0.d0)  then
  write(*,*) "something is wrong in BCS_U_phi."
  write(*,*)'weight(i) = ',weight(i)
  write(*,*)'tot_imp(i) = ',tot_imp(i)
  write(*,*)'tot_tmp = ',tot_tmp
  call mystop
end if


end subroutine BCS_U_phi


subroutine compute_FB(Am1,i,nu,nd)
use param
use phi_param
use lattice_param
use model_param
use phiT_param
use method_param
complex(kind=8),intent(IN)::Am1(Nspin(1),Nspin(1))
integer,intent(IN)::i
complex(kind=8),intent(OUT)::nu(Nsite),nd(Nsite)
integer::k,p,q,j
complex(kind=8)::phu(Nsite,Nspin(1)),phd(Nsite,Nspin(2))
complex(kind=8)::Am1phut(Nspin(1),Nsite),Fstarphd(Nsite,Nspin(2))
complex(kind=8)::FDphu(Nsite,Nspin(1)),FDphuAm1t(Nsite,Nspin(1))

phu=phi(1:Nsite,1:Nspin(1),i)
phd=phi(Nsite+1:2*Nsite,Nspin(1)+1:Ntot,i)


nu=zero
call zgemm('N','T',Nspin(1),Nsite,Nspin(1),one,Am1,Nspin(1),phu,Nsite,zero,Am1phut,Nspin(1))
call zgemm('N','N',Nsite,Nspin(2),Nsite,one,conjg(Fpairing(:,:)),Nsite,phd,Nsite,zero,Fstarphd,Nsite)
if(Nzeta.eq.0)then
  do j=1,Nsite
    do q=1,Nspin(2)
      nu(j)=nu(j)+Fstarphd(j,q)*Am1phut(q,j)
    enddo
  enddo
else 
  do j=1,Nsite
    do q=1,Nzeta
      nu(j)=nu(j)+conjg(Dunpaired(j,q))*Am1phut(q,j)
    enddo
    do q=1,Nspin(2)
      nu(j)=nu(j)+Fstarphd(j,q)*Am1phut(Nzeta+q,j)
    enddo
  enddo
endif

nd=zero
call zgemm('T','N',Nsite,Nspin(1),Nsite,one,conjg(Fpairing(:,:)),Nsite,phu,Nsite,zero,FDphu,Nsite)
call zgemm('N','T',Nsite,Nspin(1),Nspin(1),one,FDphu,Nsite,Am1,Nspin(1),zero,FDphuAm1t,Nsite)
if(Nzeta.eq.0)then
  do j=1,Nsite
    do q=1,Nspin(2)
      nd(j)=nd(j)+FDphuAm1t(j,q)*phd(j,q)
    enddo
  enddo
else
  do j=1,Nsite
    do q=1,Nspin(2)
      nd(j)=nd(j)+FDphuAm1t(j,Nzeta+q)*phd(j,q)
    enddo
  enddo
endif

end subroutine compute_FB 
