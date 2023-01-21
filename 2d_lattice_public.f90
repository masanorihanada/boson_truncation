!###########################################################  
!######              2d lattice scalar QFT          ########
!######           written by Masanori Hanada        ########
!###########################################################
!Mersenne twister.
include 'mt19937.f90'
program Truncation_effect

  use mtmod !Mersenne twistor
  implicit none

  include 'size_2d.inc'
  !---------------------------------
  character(150) input_config,data_output,output_config
  !----------------------------------------
  double precision a_dig,a_lat 
  double precision LatticeSpacing_delta !lattice spacing \Delta for Trotterization. Sorry for using a confusing name.
  double precision temperature !temperature T
  double precision IR_cutoff !IR cutoff R
  double precision mass2 !m_{lat}^2
  double precision lambda !lambda=0 in the paper
  !----------------------------------------
  !     Iteration
  integer iteration
  integer NoSkip
  !    initial configuration
  !    IniConfig=0 -> new start
  !    IniConfig=1 -> use old configuration
  integer IniConfig
  !----------------------------------------
  !   nseed=0 -> use mersenne_seed in the input file (always so if iniconfig=0)
  integer nseed
  !----------------------------------------
  !"nx" is "n" in the paper.
  integer nx(1:nsite_t,1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  !----------------------------------------
  integer nx_check_start(1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  integer nx_check_end(1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  integer nxp1_check(1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  integer nxm1_check(1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  integer seedsize,mersenne_seed
  integer,allocatable:: seed(:)
  integer isite_t,isite_x,isite_y,iter,niter,isite_x_2,isite_y_2
  integer nx_old(1:nsite_t,1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  !------------------------------------------
  double precision metropolis,ransu!"ransu" is Japanese word for "random number"
  double precision action_init,action_fin
  double precision pot,x2_1,x2_2

  integer naccept,trial,NumCalObservable,rejection,CheckHam,ite,&
       &info,count_diff_p,count_diff_m,step,t_start,t_end,t,temp
  !***************************************
  !*** Read parameters from input file ***
  !***************************************
  open(unit=10,status='OLD',file='input_2d_public.dat',action='READ')
  read(10,*)input_config
  read(10,*)data_output
  read(10,*)output_config
  read(10,*)IniConfig
  read(10,*)nseed
  read(10,*)temperature
  read(10,*)mass2
  read(10,*)lambda
  read(10,*)a_lat
  read(10,*)IR_cutoff
  read(10,*)niter
  read(10,*)NoSkip
  read(10,*)mersenne_seed
  close(10)
  a_dig=2d0*IR_cutoff/dble(nlevel-1)
  LatticeSpacing_delta=1d0/temperature/dble(Nsite_T)
  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  if(IniConfig.EQ.0)then !new start  
     call sgrnd(mersenne_seed)!use the seed in the input file
     do isite_t=1,Nsite_t
        do isite_x=1,nsite_x
           do isite_y=1,nsite_y
              nx(isite_t,isite_x,isite_y)=nlevel/2
           end do
        end do
     end do
  else if(IniConfig.EQ.1)then!read from input configuration file
     open(unit=9,status='OLD',file=input_config,action='READ')
     read(9,*) nx
     if(nseed.eq.1)then
        call mtgetu(9)!use the seed in the config file
     else
        call sgrnd(mersenne_seed)!use the seed in the input file
     end if
     close(9)
  end if
  !*******************************
  !***    Make the output file ***
  !******************************* 
  open(unit=10,status='REPLACE',file=data_output,action='WRITE')
  write(10,*) "#numbers of sites along t,x,y directions=",nsite_t,nsite_x,nsite_y
  write(10,*) "#a_dig=",a_dig
  write(10,*) "#a_lat=",a_lat
  write(10,*) "#delta=",LatticeSpacing_delta
  write(10,*) "#delta/a^2=",LatticeSpacing_delta/a_dig/a_dig
  write(10,*) "#temperature=",temperature
  write(10,*) "#mass^2=",mass2
  write(10,*) "#lambda=",lambda
  write(10,*) "# iter,potential energy, acceptance rate"
  write(10,*)'#------------------------------------------------'
  !**********************
  !*** Reset counters ***
  !**********************
  naccept=0
  trial=0
  NumCalObservable=0
  !*****************
  !*** main loop ***
  !*****************
  do iter=1,niter
     do isite_t=1,nsite_t 
        do isite_x=1,nsite_x
           do isite_y=1,nsite_y
              step=int(grnd()*dble(nsite_t)*0.5d0)+1 !step is "B" in the paper.
              t_start=isite_t
              t_end=mod(isite_t+step,nsite_t)
              if(t_end.EQ.0)then
                 t_end=nsite_t
              end if
              !backup nx (in case the proposal is rejected)
              nx_old=nx

              call calc_action(action_init,nx_old,a_lat,a_dig,LatticeSpacing_delta,mass2,lambda,IR_cutoff)
              !************************************************
              !*** n->n' (note that nx is "n" in the paper) ***
              !************************************************
              ransu=grnd()
              info=0
              do t=t_start,t_start+step
                 if(t.LE.nsite_t)then
                    temp=t
                 else
                    temp=mod(t,nsite_t)
                 end if
                 if(nx(temp,isite_x,isite_y).EQ.0)then
                    if(ransu.LT.0.5d0)then
                       !automatic reject
                       info=info+1
                    else
                       nx(temp,isite_x,isite_y)=nx(temp,isite_x,isite_y)+1
                       !info=0
                    end if
                    
                 else if(nx(temp,isite_x,isite_y).EQ.nlevel-1)then
                    if(ransu.LT.0.5d0)then
                       nx(temp,isite_x,isite_y)=nx(temp,isite_x,isite_y)-1
                       !info=0
                    else
                       !automatic reject
                       info=info+1
                    end if
                 else
                    if(ransu.LT.0.5d0)then
                       nx(temp,isite_x,isite_y)=nx(temp,isite_x,isite_y)-1
                       !info=0
                    else
                       nx(temp,isite_x,isite_y)=nx(temp,isite_x,isite_y)+1
                       !info=0
                    end if
                 end if
              end do
              !****************************
              !*** test the constraint. ***
              !****************************
              do isite_x_2=1,nsite_x
                 do isite_y_2=1,nsite_y
                    nx_check_start(isite_x_2,isite_y_2)=nx(t_start,isite_x_2,isite_y_2)
                    nx_check_end(isite_x_2,isite_y_2)=nx(t_end,isite_x_2,isite_y_2)
                    if(t_end.NE.nsite_t)then
                       nxp1_check(isite_x_2,isite_y_2)=nx(t_end+1,isite_x_2,isite_y_2)
                    else if(t_end.EQ.nsite_t)then
                       nxp1_check(isite_x_2,isite_y_2)=nx(1,isite_x_2,isite_y_2)
                    end if
                    if(t_start.NE.1)then
                       nxm1_check(isite_x_2,isite_y_2)=nx(t_start-1,isite_x_2,isite_y_2)
                    else if(t_start.EQ.1)then
                       nxm1_check(isite_x_2,isite_y_2)=nx(nsite_t,isite_x_2,isite_y_2)
                    end if
                 end do
              end do
              count_diff_p=0
              count_diff_m=0
              do isite_x_2=1,nsite_x
                 do isite_y_2=1,nsite_y
                    count_diff_p=count_diff_p+&
                         &iabs(nx_check_end(isite_x_2,isite_y_2)-nxp1_check(isite_x_2,isite_y_2))
                    count_diff_m=count_diff_m+&
                            &iabs(nx_check_start(isite_x_2,isite_y_2)-nxm1_check(isite_x_2,isite_y_2))
                 end do
              end do
              
              if((count_diff_p.LE.1).AND.(count_diff_m.LE.1).AND.(info.EQ.0))then
                 call calc_action(action_fin,nx,a_lat,a_dig,LatticeSpacing_delta,mass2,lambda,IR_cutoff)
                 !***********************
                 !*** Metropolis test ***
                 !***********************
                 metropolis=grnd()
                 if(metropolis.LT.dexp(action_init-action_fin))then!accept
                    naccept=naccept+1
                 else
                    !reject.
                    nx=nx_old
                 end if
              else
                 !automatic reject due to constraint violation
                 nx=nx_old
              end if
              
           end do
        end do
     end do
     !*******************
     !*** Measurement ***
     !*******************
     if(MOD(iter,NoSkip).EQ.0)then
        NumCalObservable=NumCalObservable+1         
        call calc_potential(pot,nx,mass2,lambda,IR_cutoff,a_dig)
        call x2_4by4_lattice(x2_1,x2_2,nx,IR_cutoff,a_dig,a_lat)
        write(10,'(i10,4f15.8)')iter,pot,x2_1,x2_2,dble(naccept)/dble(nsite_t*nsite_x*nsite_y*iter)
        write(*,'(i10,4f15.8)')iter,pot,x2_1,x2_2,dble(naccept)/dble(nsite_t*nsite_x*nsite_y*iter)
     end if
     
  end do
  !************************
  !*** End of main loop ***
  !************************
  close(10)
  !**************************
  !*** Save configuration ***
  !**************************
  open(UNIT = 22, File = output_config, STATUS = "REPLACE", ACTION = "WRITE")
  write(22,*) nx
  call mtsaveu(22)
  close(22)

  
end program Truncation_effect

!***************************************************************
!*** Box-Muller method for generating Gaussian random number ***
!***************************************************************
SUBROUTINE BoxMuller(p,q)  
  
  use mtmod !Mersenne twistor
  implicit none 
  
  doubleprecision p,q,r,s,Pi
  
  Pi=2d0*DASIN(1d0)
  r=grnd()
  s=grnd()
  p=dsqrt(-2d0*dlog(r))*DSIN(2d0*Pi*s)
  q=dsqrt(-2d0*dlog(r))*DCOS(2d0*Pi*s)
  
  return
  
END SUBROUTINE BoxMuller
!****************************
!*** Calculate the action ***
!****************************
subroutine calc_action(action,nx,a_lat,a_dig,LatticeSpacing_delta,mass2,lambda,IR_cutoff)

  implicit none

  include 'size_2d.inc'

  integer nx(1:nsite_t,1:nsite_x,1:nsite_y),nxj,nxjp1,nxjm1,nbos
  INTEGER isite_t,isite_x,isite_y
  double precision action,a_dig,LatticeSpacing_delta,mass2,lambda,IR_cutoff,xj,xjm1,xjp1
  DOUBLE PRECISION a_lat
  integer count_diff_p,count_diff_m
  integer nx_check(1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  integer nxp1_check(1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  integer nxm1_check(1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1

  !**************
  !*** V(phi) ***
  !**************
  call calc_potential(action,nx,mass2,lambda,IR_cutoff,a_dig)
  action=action*dble(nsite_t*nsite_x*nsite_y)*LatticeSpacing_delta
  !***********
  !*** P^2 ***
  !***********
  do isite_t=1,nsite_t
     do isite_x=1,nsite_x
        do isite_y=1,nsite_y
           nx_check(isite_x,isite_y)=nx(isite_t,isite_x,isite_y)
           if(isite_t.NE.nsite_t)then
              nxp1_check(isite_x,isite_y)=nx(isite_t+1,isite_x,isite_y)
           else if(isite_t.EQ.nsite_t)then
              nxp1_check(isite_x,isite_y)=nx(1,isite_x,isite_y)
           end if
        end do
     end do
     count_diff_p=0
     do isite_x=1,nsite_x
        do isite_y=1,nsite_y
           
           count_diff_p=count_diff_p+&
                &iabs(nx_check(isite_x,isite_y)-nxp1_check(isite_x,isite_y))
  
           
        end do
     end do
     
     nbos=nsite_x*nsite_y

     if((count_diff_p.EQ.0))then
        action=action-dlog(1d0-LatticeSpacing_delta*dble(nbos)/(a_dig*a_dig))
     else if(count_diff_p.EQ.1)then
        action=action-dlog(0.5d0*LatticeSpacing_delta/(a_dig*a_dig))
     end if
  end do
  !********************************
  !*** Spatial derivative terms ***
  !********************************
  do isite_t=1,nsite_t
     do isite_x=1,nsite_x
        do isite_y=1,nsite_y
           !Add spatial derivative. x-direction
           nxj=nx(isite_t,isite_x,isite_y)
           if(isite_x.EQ.nsite_x)then
              nxjp1=nx(isite_t,1,isite_y)
           else
              nxjp1=nx(isite_t,isite_x+1,isite_y)
           end if
           xj=-IR_cutoff+a_dig*dble(nxj)
           xjp1=-IR_cutoff+a_dig*dble(nxjp1)
           action=action+0.5d0/a_lat/a_lat*(xj-xjp1)*(xj-xjp1)*LatticeSpacing_delta
           !Add spatial derivative. y-direction
           nxj=nx(isite_t,isite_x,isite_y)
           if(isite_y.EQ.nsite_y)then
              nxjp1=nx(isite_t,isite_x,1)
           else
              nxjp1=nx(isite_t,isite_x,isite_y+1)
           end if
           xj=-IR_cutoff+a_dig*dble(nxj)
           xjp1=-IR_cutoff+a_dig*dble(nxjp1)
           action=action+0.5d0/a_lat/a_lat*(xj-xjp1)*(xj-xjp1)*LatticeSpacing_delta
        end do
     end do
  end do
   
  return

END subroutine Calc_Action
!************************
!*** Calculate V(phi) ***
!************************
subroutine calc_potential(pot,nx,mass2,lambda,IR_cutoff,a_dig)

  implicit none

  include 'size_2d.inc'

  double precision pot,a_dig,IR_cutoff,mass2,lambda,xj
  integer nx(1:nsite_t,1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  
  integer isite_t,isite_x,isite_y

  pot=0d0
  do isite_t=1,nsite_t
     do isite_x=1,nsite_x
        do isite_y=1,nsite_y
           xj=-IR_cutoff+a_dig*dble(nx(isite_t,isite_x,isite_y))
           pot=pot+0.5d0*mass2*xj*xj+0.25d0*lambda*xj*xj*xj*xj
        end do
     end do
  end do
  pot=pot/dble(nsite_t*nsite_x*nsite_y)
  
  return
  
END subroutine Calc_potential
!****************************************************************************
!*** Calculate <|\tilde{phi}_{\vec{q}}|^2> for \vec{q}=(0,0) and (pi,pi). ***
!*** This subroutine works only for 4*4 lattice.           ******************
!****************************************************************************
subroutine x2_4by4_lattice(x2_1,x2_2,nx,IR_cutoff,a_dig,a_lat)

  implicit none

  include 'size_2d.inc'

  double precision a_dig,a_lat,IR_cutoff,xj,x_1,x_2,x2_1,x2_2
  integer nx(1:nsite_t,1:nsite_x,1:nsite_y)!runs from 0 to nlevel-1
  
  integer isite_t,isite_x,isite_y

  x2_1=0d0!\vec{q}=(0,0)
  x2_2=0d0!\vec{q}=(pi,pi)
  do isite_t=1,nsite_t
     x_1=0d0
     x_2=0d0
     do isite_x=1,nsite_x
        do isite_y=1,nsite_y
           xj=-IR_cutoff+a_dig*dble(nx(isite_t,isite_x,isite_y))
           x_1=x_1+xj
           if(mod(isite_x+isite_y,2).EQ.0)then
              x_2=x_2+xj
           else
              x_2=x_2-xj
           end if
        end do
     end do
     x2_1=x2_1+x_1*x_1/(dble(nsite_x*nsite_y))
     x2_2=x2_2+x_2*x_2/(dble(nsite_x*nsite_y))
  end do
  x2_1=x2_1/dble(nsite_t)
  x2_2=x2_2/dble(nsite_t)
  
  return
  
END subroutine x2_4by4_lattice
