program read_data

implicit none
integer::Nbravais,Nbands,L,Lx,Ly,ibc,Dist,sample_size,Nsample,Nbeta,Npair,I_twob
real(kind=8)::elle,tpielx,tpiely,dt
real(kind=8),dimension(:,:),allocatable::r,q
complex(kind=8), dimension(2,2) :: sigmax,sigmay,sigmaz
complex(kind=8), dimension(:,:,:), allocatable::gofr,G
real(kind=8),dimension(:,:,:), allocatable::av_G,err_G
complex(kind=8), dimension(:,:), allocatable::Gt


real(kind=8)::pi
parameter(pi=3.14159265358979d0)
complex(kind=8)::zero,one,Xi,half
parameter(zero=dcmplx(0.d0,0.d0))
parameter(one=dcmplx(1.d0,0.d0))
parameter(Xi=dcmplx(0.d0,1.d0))
parameter(half=dcmplx(0.5d0,0.d0))

integer::counter
integer::i,j,ib,jb,sitei,i1,i2,i3,momi,alpha,beta,idist,isample,i_beta
character*3:: sfix
real(kind=8)::x,y,e,qr,tpiel,dr,rx,ry,dx,ex,ey,ez
real(kind=8)::edummyx,edummyy,edummyz,spin,espin
complex(kind=8)::dummy,dummyx,dummyy,dummyz,sx,sy,sz,Gq


sigmax(1,1)=zero
sigmax(2,1)=one
sigmax(1,2)=one
sigmax(2,2)=zero



sigmay(1,1)=zero
sigmay(2,1)=Xi
sigmay(1,2)=-Xi
sigmay(2,2)=zero


sigmaz(1,1)=one
sigmaz(2,1)=zero
sigmaz(1,2)=zero
sigmaz(2,2)=-one

write(*,*)'Nbands? '
read(*,*)Nbands
write(*,*)'Lx, Ly? '
read(*,*)Lx,Ly
write(*,*)'Open or closed boundary conditions [0/1]?'
read(*,*)ibc
if(ibc.eq.1)then
   Dist=Ly/2
else
   Dist=Ly
endif
write(*,*)'Npair?'
read(*,*)Npair
write(*,*)'Nbeta?'
read(*,*)Nbeta
write(*,*)'dt?'
read(*,*)dt


!Nbands=3
Nbravais=Lx*Ly
write(*,*)'Number of unit cells Lx x Ly  = ',Nbravais

allocate(r(2,Nbands*Nbravais))
allocate(q(2,Nbands*Nbravais))
allocate(gofr(Nbands*Nbravais,2,2))
allocate(G(Nbravais,Nbands,Nbands))
allocate(av_G(Nbravais,3,3))
allocate(err_G(Nbravais,3,3))
allocate(Gt(6*Nbravais,0:2*Nbeta))

open(2,file='lattice.info',status='old')
do i=1,Nbravais
  read(2,*)r(1,i),r(2,i)
enddo
if(Nbands.eq.3)then
do i=Nbravais+1,2*Nbravais
  r(1,i)=r(1,i-Nbravais)
  r(2,i)=r(2,i-Nbravais)
enddo
do i=2*Nbravais+1,3*Nbravais
  r(1,i)=r(1,i-2*Nbravais)
  r(2,i)=r(2,i-2*Nbravais)
enddo
endif
tpielx=2.d0*pi/dble(Lx)
tpiely=2.d0*pi/dble(Ly)
do i=1,Nbands*Nbravais
  q(1,i)=r(1,i)*tpielx
  q(2,i)=r(2,i)*tpiely
enddo

close(2)


open(10,file='spin_correlation',status='old')


sample_size=Nbravais

counter=0
av_G=0.d0
err_G=0.d0

do isample=1,1000
  do i=1,sample_size
    read(10,*,end=2)sitei,e,x,y
    alpha=1
    beta=1
    G(sitei,alpha,beta)=dcmplx(x,y)
    av_G(sitei,alpha,beta)=av_G(sitei,alpha,beta)+x
    err_G(sitei,alpha,beta)=err_G(sitei,alpha,beta)+x*x
  enddo
  counter=counter+1  
  write(*,*)'Reading spin correlation ',counter
  call xifs(sfix,isample)
!  open(11,file='read_spin_correlation.'//sfix,status='unknown')
!    do i=1,Nbravais
!     write(11,*)r(1,i),r(2,i),dble(G(i,1,1)),aimag(G(i,1,1))
!    enddo
!  close(11)
enddo
2 continue
close(10)

write(*,*)
write(*,*)'I read ',counter,' estimations'
write(*,*)
write(*,*)

do sitei=1,Nbravais
  ib=1
  jb=1
  av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
  if(counter.eq.1)then
    err_G(sitei,jb,ib)=0.d0
  else
    err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
  endif
enddo

open(11,file='average_spin_correlation.result',status='unknown')
  do sitei=1,Nbravais
    ib=1
    jb=1
    write(11,*)r(1,sitei),r(2,sitei),av_G(sitei,jb,ib),err_G(sitei,jb,ib)
  enddo
close(11)

open(9,file='density_correlation',status='old')

sample_size=Nbands*Nbands*Nbravais

counter=0
av_G=0.d0
err_G=0.d0

do isample=1,1000
  do i=1,sample_size
    read(9,*,end=3)sitei,jb,ib,x
    G(sitei,jb,ib)=dcmplx(x,0.d0)
  enddo
  counter=counter+1
  write(*,*)'Reading density correlations ',counter
  call xifs(sfix,isample)
!  open(21,file='read_density_correlation.'//sfix,status='unknown')
  do ib=1,Nbands
    do jb=ib,Nbands
      do i=1,Nbravais
        x=0.5d0*dble(G(i,jb,ib)+G(i,ib,jb))
        y=0.5d0*aimag(G(i,jb,ib)+G(i,ib,jb))
!        write(21,'(2f10.4,2f15.8,2i2)')r(1,i),r(2,i),x,y,jb,ib
        av_G(i,jb,ib)=av_G(i,jb,ib)+x
        err_G(i,jb,ib)=err_G(i,jb,ib)+x*x
      enddo       
    enddo
  enddo 
!  close(21)
enddo
3 continue
close(9)

write(*,*)
write(*,*)'I read ',counter,' estimations'
write(*,*)
write(*,*)


do ib=1,Nbands
    do jb=ib,Nbands
      do sitei=1,Nbravais  
        av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
        err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
        err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
        if(counter.eq.1)then
          err_G(sitei,jb,ib)=0.d0
        else
          err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
        endif
    enddo
  enddo
enddo

open(11,file='average_density_correlation.result',status='unknown')
do ib=1,Nbands
    do jb=ib,Nbands
       do sitei=1,Nbravais
         write(11,'(2f10.4,2f15.8,2i2)')r(1,sitei),r(2,sitei),av_G(sitei,jb,ib),err_G(sitei,jb,ib),jb,ib
       enddo
       write(11,*)
       write(11,*)
    enddo
enddo
close(11)


open(9,file='szsz_correlation',status='old')

sample_size=Nbands*Nbands*Nbravais

counter=0
av_G=0.d0
err_G=0.d0

do isample=1,1000
  do i=1,sample_size
    read(9,*,end=31)sitei,jb,ib,x
    G(sitei,jb,ib)=dcmplx(x,0.d0)
  enddo
  counter=counter+1
  write(*,*)'Reading szsz correlations ',counter
  call xifs(sfix,isample)
!  open(21,file='read_szsz_correlation.'//sfix,status='unknown')
  do ib=1,Nbands
    do jb=ib,Nbands
      do i=1,Nbravais
        x=0.5d0*dble(G(i,jb,ib)+G(i,ib,jb))
        y=0.5d0*aimag(G(i,jb,ib)+G(i,ib,jb))
!        write(21,'(2f10.4,2f15.8,2i2)')r(1,i),r(2,i),x,y,jb,ib
        av_G(i,jb,ib)=av_G(i,jb,ib)+x
        err_G(i,jb,ib)=err_G(i,jb,ib)+x*x
      enddo
    enddo
  enddo
!  close(21)
enddo
31 continue
close(9)

write(*,*)
write(*,*)'I read ',counter,' estimations'
write(*,*)
write(*,*)


do ib=1,Nbands
    do jb=ib,Nbands
      do sitei=1,Nbravais
        av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
        err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
        err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
        if(counter.eq.1)then
          err_G(sitei,jb,ib)=0.d0
        else
          err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
        endif
    enddo
  enddo
enddo

open(11,file='average_szsz_correlation.result',status='unknown')
do ib=1,Nbands
    do jb=ib,Nbands
       do sitei=1,Nbravais
         write(11,'(2f10.4,2f15.8,2i2)')r(1,sitei),r(2,sitei),av_G(sitei,jb,ib),err_G(sitei,jb,ib),jb,ib
       enddo
       write(11,*)
       write(11,*)
    enddo
enddo
close(11)



if(Npair.ne.0)then

open(9,file='pairing.1',status='old')
sample_size=Nbravais

counter=0
av_G=0.d0
err_G=0.d0

do isample=1,1000
  do i=1,sample_size
    read(9,*,end=33)sitei,x,y
    alpha=1
    beta=1
    gofr(sitei,alpha,beta)=dcmplx(x,y)
  enddo
  counter=counter+1
  call xifs(sfix,isample)
!  open(21,file='read_pairing_swave.'//sfix,status='unknown')
    do i=1,Nbravais
     jb=1
     ib=1
     x=dble(gofr(i,1,1))
     y=aimag(gofr(i,1,1))
!     write(21,*)r(1,i),r(2,i),x,y
     av_G(i,jb,ib)=av_G(i,jb,ib)+x
     err_G(i,jb,ib)=err_G(i,jb,ib)+x*x
    enddo
!  close(21)
enddo

33 continue
close(9)


do sitei=1,Nbravais
  ib=1
  jb=1
  av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
  if(counter.eq.1)then
    err_G(sitei,jb,ib)=0.d0
  else
    err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
  endif
enddo

open(11,file='average_on_site_swave_pairing.result',status='unknown')
  do sitei=1,Nbravais
    ib=1
    jb=1
    write(11,*)r(1,sitei),r(2,sitei),av_G(sitei,jb,ib),err_G(sitei,jb,ib)
  enddo
close(11)


if(Npair.ge.2)then
open(10,file='pairing.2',status='old')

counter=0
av_G=0.d0
err_G=0.d0


do isample=1,1000
  do i=1,sample_size
    read(10,*,end=34)sitei,x,y
    alpha=1
    beta=1
    gofr(sitei,alpha,beta)=dcmplx(x,y)
  enddo
  counter=counter+1
  call xifs(sfix,isample)
!  open(21,file='read_pairing_dwave.'//sfix,status='unknown')
    do i=1,Nbravais
      jb=1
      ib=1
      x=dble(gofr(i,1,1))
      y=aimag(gofr(i,1,1))
!      write(21,*)r(1,i),r(2,i),x,y
      av_G(i,jb,ib)=av_G(i,jb,ib)+x
      err_G(i,jb,ib)=err_G(i,jb,ib)+x*x
    enddo
!  close(21)
enddo

34 continue
close(10)


do sitei=1,Nbravais
  ib=1
  jb=1
  av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
  if(counter.eq.1)then
    err_G(sitei,jb,ib)=0.d0
  else
    err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
  endif
enddo

open(11,file='average_extended_swave_pairing.result',status='unknown')
  do sitei=1,Nbravais
    ib=1
    jb=1
    write(11,*)r(1,sitei),r(2,sitei),av_G(sitei,jb,ib),err_G(sitei,jb,ib)
  enddo
close(11)



if(Npair.ge.3)then


open(12,file='pairing.3',status='old')

counter=0
av_G=0.d0
err_G=0.d0


do isample=1,1000
  do i=1,sample_size
    read(12,*,end=341)sitei,x,y
    alpha=1
    beta=1
    gofr(sitei,alpha,beta)=dcmplx(x,y)
  enddo
  counter=counter+1
  call xifs(sfix,isample)
!  open(21,file='read_pairing_dwave.'//sfix,status='unknown')
    do i=1,Nbravais
      jb=1
      ib=1
      x=dble(gofr(i,1,1))
      y=aimag(gofr(i,1,1))
!      write(21,*)r(1,i),r(2,i),x,y
      av_G(i,jb,ib)=av_G(i,jb,ib)+x
      err_G(i,jb,ib)=err_G(i,jb,ib)+x*x
    enddo
!  close(21)
enddo

341 continue
close(12)

do sitei=1,Nbravais
  ib=1
  jb=1
  av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
  if(counter.eq.1)then
    err_G(sitei,jb,ib)=0.d0
  else
    err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
  endif
enddo

open(11,file='average_dwave_pairing.result',status='unknown')
  do sitei=1,Nbravais
    ib=1
    jb=1
    write(11,*)r(1,sitei),r(2,sitei),av_G(sitei,jb,ib),err_G(sitei,jb,ib)
  enddo
close(11)


endif


endif


endif

if(Nbeta.eq.0)return

go to 100

open(10,file='green_p_dynamical',status='unknown')
sample_size=Nbands*Nbravais*2*Nbeta
write(*,*)'Now care for dynamics',sample_size

do isample=1,1000
   
  do i=1,sample_size
    read(10,*,end=4)sitei,i_beta,x,y
    write(*,*)i
    Gt(sitei,i_beta)=dcmplx(x,y)
  enddo
  call xifs(sfix,isample)
  write(*,*)'Now open file'
  open(22,file='gqt.'//sfix,status='unknown')
  do i_beta=1,Nbeta
    do momi=1,Nbravais
        Gq=zero
        do sitei=1,Nbravais
          qr=q(1,momi)*r(1,sitei)+q(2,momi)*r(2,sitei)
          Gq=Gq+exp(Xi*qr)*Gt(sitei,i_beta)
        enddo
      write(22,*)i_beta*dt,momi,dble(Gq),aimag(Gq)
    enddo
  enddo
  close(22)
enddo
4 continue
close(10)


100 continue

deallocate(av_G,err_G)

allocate(av_G(0:Nbeta,3,3))
allocate(err_G(0:Nbeta,3,3))


open(15,file='cnuc+mu_dynamical',status='old')
open(16,file='c+nucmu_dynamical',status='old')

sample_size=Nbeta+1


counter=0
av_G=0.d0
err_G=0.d0

do isample=1,1000
  do i=1,sample_size
    read(15,*,end=32)sitei,i_beta,x,y
    Gt(sitei,i_beta)=dcmplx(x,y)
  enddo
  counter=counter+1
  call xifs(sfix,isample)
  write(*,*)'Now open file'
  open(22,file='read_cnuc+mu_beta.'//sfix,status='unknown')
  do i_beta=0,Nbeta
    do sitei=1,1
      jb=1
      ib=1
      Gq=Gt(sitei,i_beta)
      x=dble(Gq)
      write(22,*)i_beta*dt,sitei,dble(Gq),aimag(Gq)
      av_G(i_beta,jb,ib)=av_G(i_beta,jb,ib)+x
      err_G(i_beta,jb,ib)=err_G(i_beta,jb,ib)+x*x
    enddo
  enddo
  close(22)
enddo
32 continue

do sitei=0,Nbeta
  ib=1
  jb=1
  av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
  if(counter.eq.1)then
    err_G(sitei,jb,ib)=0.d0
  else
    err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
  endif
enddo

open(11,file='average_cnuc+mu_beta.result',status='unknown')
  do sitei=0,Nbeta
    ib=1
    jb=1
    write(11,*)sitei*dt,av_G(sitei,jb,ib),err_G(sitei,jb,ib)
  enddo
close(11)


counter=0
av_G=0.d0
err_G=0.d0

do isample=1,1000
  do i=1,sample_size
    read(16,*,end=93)sitei,i_beta,x,y
    Gt(sitei,i_beta)=dcmplx(x,y)
  enddo
  counter=counter+1
  call xifs(sfix,isample)
  write(*,*)'Now open file'
  open(22,file='read_c+nucmu_beta.'//sfix,status='unknown')
  do i_beta=0,Nbeta
    do sitei=1,1
      ib=1
      jb=1
      Gq=Gt(sitei,i_beta)
      x=dble(Gq)
      write(22,*)i_beta*dt,sitei,dble(Gq),aimag(Gq)
      av_G(i_beta,jb,ib)=av_G(i_beta,jb,ib)+x
      err_G(i_beta,jb,ib)=err_G(i_beta,jb,ib)+x*x
    enddo
  enddo
  close(22)
enddo
93 continue


do sitei=0,Nbeta
  ib=1
  jb=1
  av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
  if(counter.eq.1)then
    err_G(sitei,jb,ib)=0.d0
  else
    err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
  endif
enddo


open(11,file='average_c+nucmu_beta.result',status='unknown')
  do sitei=0,Nbeta
    ib=1
    jb=1
    write(11,*)sitei*dt,av_G(sitei,jb,ib),err_G(sitei,jb,ib)
  enddo
close(11)




deallocate(av_G,err_G)


write(*,*)'Did you compute two body dynamical responses? [1=yes, 0=no]'
read(*,*)I_twob

if(I_twob.eq.0)go to 76

open(13,file='nupnup_dynamical',status='old')
open(12,file='ndnnup_dynamical',status='old')


allocate(av_G(0:((Nbeta+1)*Nbravais),3,3))
allocate(err_G(0:((Nbeta+1)*Nbravais),3,3))


counter=0
av_G=0.d0
err_G=0.d0



sample_size=Nbravais*(Nbeta+1)
do isample=1,1000

  do i=1,sample_size
    read(13,*,end=5)sitei,i_beta,x,y
    write(*,*)i
    Gt(sitei,i_beta)=dcmplx(x,y)
  enddo
  counter=counter+1
  call xifs(sfix,isample)
  write(*,*)'Now open file'
  open(22,file='read_nupnup_beta.'//sfix,status='unknown')
  i3=0
  do i_beta=0,Nbeta
        do sitei=1,Nbravais
           ib=1
           jb=1
           Gq=Gt(sitei,i_beta)
           x=dble(Gq)
           write(22,*)i_beta*dt,sitei,dble(Gq),aimag(Gq)
           av_G(i3,jb,ib)=av_G(i3,jb,ib)+x
           err_G(i3,jb,ib)=err_G(i3,jb,ib)+x*x
           i3=i3+1
        enddo
  enddo
  close(22)
enddo
5 continue

do sitei=0,Nbravais*(Nbeta+1)
  ib=1
  jb=1
  av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
  if(counter.eq.1)then
    err_G(sitei,jb,ib)=0.d0
  else
    err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
  endif
enddo


open(11,file='average_nupnup_beta.result',status='unknown')

  sitei=0
  do i_beta=0,Nbeta
    do i=1,Nbravais
      ib=1
      jb=1
      write(11,*)i_beta*dt,i,av_G(sitei,jb,ib),err_G(sitei,jb,ib)
      sitei=sitei+1
    enddo
  enddo
close(11)

counter=0
av_G=0.d0
err_G=0.d0

do isample=1,1000

  do i=1,sample_size
    read(12,*,end=6)sitei,i_beta,x,y
    write(*,*)i
    Gt(sitei,i_beta)=dcmplx(x,y)
  enddo
  counter=counter+1
  call xifs(sfix,isample)
  write(*,*)'Now open file'
  open(22,file='read_ndnnup_beta.'//sfix,status='unknown')
  i3=0
  do i_beta=0,Nbeta
        do sitei=1,Nbravais
           ib=1
           jb=1
           Gq=Gt(sitei,i_beta)
           write(22,*)i_beta*dt,sitei,dble(Gq),aimag(Gq)
           x=dble(Gq)
           av_G(i3,jb,ib)=av_G(i3,jb,ib)+x
           err_G(i3,jb,ib)=err_G(i3,jb,ib)+x*x
           i3=i3+1
        enddo
  enddo
  close(22)
enddo
6 continue

do sitei=0,Nbravais*(Nbeta+1)
  ib=1
  jb=1
  av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
  if(counter.eq.1)then
    err_G(sitei,jb,ib)=0.d0
  else
    err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
  endif
enddo



open(11,file='average_ndnnup_beta.result',status='unknown')

  sitei=0
  do i_beta=0,Nbeta
    do i=1,Nbravais
      ib=1
      jb=1
      write(11,*)i_beta*dt,i,av_G(sitei,jb,ib),err_G(sitei,jb,ib)
      sitei=sitei+1
    enddo
  enddo
close(11)


76 continue

write(*,*)
write(*,*)'Superfluid density? [0 = no, 1 = yes]'
read(*,*)i3

if(i3.eq.0)go to 199

deallocate(av_G,err_G)

allocate(av_G(0:Nbeta,3,3))
allocate(err_G(0:Nbeta,3,3))

open(18,file='superfluid_dynamical',status='old')

sample_size=Nbeta+1


counter=0
av_G=0.d0
err_G=0.d0

do isample=1,1000
  do i=1,sample_size
    read(18,*,end=98)sitei,i_beta,x,y
    Gt(sitei,i_beta)=dcmplx(x,y)
  enddo
  counter=counter+1
  call xifs(sfix,isample)
  write(*,*)'Now open file'
  open(22,file='read_superfluid_beta.'//sfix,status='unknown')
  do i_beta=0,Nbeta
    do sitei=1,1
      ib=1
      jb=1
      Gq=Gt(sitei,i_beta)
      x=dble(Gq)
      write(22,*)i_beta*dt,sitei,dble(Gq),aimag(Gq)
      av_G(i_beta,jb,ib)=av_G(i_beta,jb,ib)+x
      err_G(i_beta,jb,ib)=err_G(i_beta,jb,ib)+x*x
    enddo
  enddo
  close(22)
enddo
98 continue


do sitei=0,Nbeta
  ib=1
  jb=1
  av_G(sitei,jb,ib)=av_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dble(counter)
  err_G(sitei,jb,ib)=dsqrt(abs(err_G(sitei,jb,ib)-av_G(sitei,jb,ib)**2))
  if(counter.eq.1)then
    err_G(sitei,jb,ib)=0.d0
  else
    err_G(sitei,jb,ib)=err_G(sitei,jb,ib)/dsqrt(dble(counter-1))
  endif
enddo


open(11,file='average_superfluid_beta.result',status='unknown')
  do sitei=0,Nbeta
    ib=1
    jb=1
    write(11,*)sitei*dt,av_G(sitei,jb,ib),err_G(sitei,jb,ib)
  enddo
close(11)





199 continue

stop
end


subroutine xifs(sfix,i)
implicit none
character*3 sfix
integer i
if(i.lt.10)then
  write(sfix,'(i1)')i
elseif(i.lt.100)then
  write(sfix,'(i2)')i
elseif(i.lt.1000)then
  write(sfix,'(i3)')i
endif
return
end
