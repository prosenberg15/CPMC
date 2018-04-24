!main program of code
program main
use mpi_serial_param
use timing_module
use rand_num
use phiT_param
use lattice_param
implicit none

!--------------------------------------------------
!Initial rank and Nsize, if mpi then inital the MPI
!--------------------------------------------------
 call start_code()

!-----------------------------------------------
!BeginTiming and EndTiming get the running time.
!-----------------------------------------------
 if(rank.eq.0) write(*,*) "Start to run the routine."
 call BeginTiming()
 call PrintTiming()


!----------------------------------------------
!initial parameter:random number and read param
!----------------------------------------------
 call init_genrand()
 call readparam()


!-----------------------------------------------------------
!Get the lattice information, Get Hzero matrix of the lattice
!-----------------------------------------------------------
 call set_lattice()


!--------------------------------------------------
!Set the mpi htype and phtype after set the lattice
!--------------------------------------------------
 call init_mpi_type()


!-------------------------------------
!Inital the projection part of K and V
!-------------------------------------
 call inital_k_v()


!-----------------------------------
!Get the trial wave function
!-----------------------------------
 call get_phiT()


!-----------------------------------------------
!Method we want to use:CPMC,RCPMC,FPMC
!-----------------------------------------------
 call chose_method()


!--------------------
!Prepare for the stop
!--------------------
 call clean

 
end program main




!---------------------------------------------------
!Print the method we want to use: CPMC, RCPMC, FPMC
!---------------------------------------------------
subroutine chose_method()
use mpi_serial_param
use method_param
implicit none


if(max_crn.LT.0) then        !Step by step run
  if(rank.eq.0) write(*,*)   "Running step by step QMC measurement."
  if(crn.GT.0.d0) then
    if(rank.eq.0) write(*,*) "This is a free projection run."
  else if(crn.LT.0.d0) then
    if(rank.eq.0) write(*,*) "This is a constraint run."
  else
    if(rank.eq.0) write(*,*) "Something is wrong with the crn input:",crn
    call mystop
  end if
  call step_qmc()
else if(max_crn.EQ.0) then   !CPMC run
  if(crn.LT.0.d0) then
    if(rank.eq.0) write(*,*) "Running CPQMC for the result."
    call cpmc()
  else
    if(rank.eq.0) write(*,*) "When max_crn EQ 0, we are using cpmc, crn must be LT 0:",crn
    call mystop
  end if  
else if(max_crn.GT.0) then   !RCPMC run
  if(crn.GT.0.d0) then
    if(rank.eq.0) write(*,*) "Running Released Constraint for the result."
    if(rank.eq.0) write(*,*) "Release step is:",max_crn
    crn=-1.d0 !Run CPQMC first
    call rcpmc()  
  else
    if(rank.eq.0) write(*,*) "When max_crn GT 0, we are use rcpmc, crn must be GT 0",crn
    call mystop
  end if
end if


if(rank.eq.0) write(*,*) ""
if(rank.eq.0) write(*,*) ""
if(rank.eq.0) write(*,*) ""

end subroutine chose_method




!-------------------------------------------------------
!The therm subroutine in the QMC for both FPMC CPMC RPMC
!-------------------------------------------------------
subroutine therm(i_pop,i_GS)
use param
use mc_loop_param
use mpi_serial_param
use adET_param
use meas_param
use method_param
!DEBUG
use phi_param
implicit none
integer,intent(INOUT)::i_pop,i_GS
logical::ad_ET
integer::i_therm,max_therm

max_therm=Thermblock*blockstep
ad_ET=.false.;E_sum=0.d0;E_sum_old=0.d0;i_ad=0

!to avoid problems if PopContrlstep=1
if(PopContrlstep.eq.1)then 
  if(rank.eq.0) write(*,*)"WARNING: it is better to choose PopContrlstep > 1"
  m_w=1.d0;m_w_old=m_w
endif

call K_phi_W(1)

call init_ET()

if(rank.eq.0) write(*,*) "Start therm:"
do i_therm=1,max_therm,1
   if(rank.eq.0) write(*,*) i_therm

   i_pop=i_pop+1;i_GS=i_GS+1

   !propagate one step in MC
   call K_V_update()

   !Do periodically Modified GS
   !When free-projection it is done
   !in the population control step
   if(crn.LT.0.d0) then
     if(i_GS.EQ.StepforGram) then
       call Modified_GS()
       i_GS=0
     end if
   end if

   !Do periodically population control and adjust ET
   !If crn.GT.0.d0  do modified gs
   if(i_pop.EQ.PopContrlstep) then
     if(crn.GT.0.d0) then
       call Modified_GS()
       i_GS=0
     end if
     call PopControl()
     i_pop=0
     ad_ET=.true.
   end if
   call adjustET(ad_ET)

end do

!write phi_rank
!if(rank.eq.0) then
!  write(*,*) "write the phi of different rank into file"
!  write(*,*) "We'd better write the phi just after population control"
!  write(*,*) "The next run with the phi might not do pop contrl at the beginning."
!  write(*,*) ""
!  write(*,*) ""
!  write(*,*) ""
!end if
!call write_phi_rank()
!call init_ET()
!call mystop

end subroutine therm
