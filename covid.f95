program Covid_model_Spain
  !SEIR model with Commuting
  !Commuting coded in force of infection following (1).
  !All-to-all interactions.
  !The time evolution can or cannot be divided in two periods: with-without lockdown.
  !Variables are commented when initilized for first time
  !-----------------------------------------------------------------------------
  !------------Code developed by Javier Aguilar---------------------------------
  !------------Help/comments: Sandro Meloni, Jose J. Ramasco and Pere Colet-----
  !------------https://ifisc.uib-csic.es/en/people/-----------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !------------For doubts/comments javieraguilar@ifisc.uib-csic.es--------------
  !-----------------------------------------------------------------------------
  !REFERENCES:
  !(1)  BALCAN, Duygu, et al. Modeling the spatial spread of infectious diseases:
  ! The GLobal Epidemic and Mobility computational model. Journal of computational
  ! science, 2010, vol. 1, no 3, p. 132-145.
  !-----------------------------------------------------------------------------
  !Main text, preprint version:
  !Aguilar, Javier, et al. "Impact of urban structure on COVID-19 spread."
  !arXiv preprint arXiv:2007.15367 (2020).
  !-----------------------------------------------------------------------------
  !---ACKNOWLEDGEMENTS----------------------------------------------------------
  !Richard Chandler and Paul Northrop, authors of RANDGEN.F (source of Random2.f)
  !https://www.ucl.ac.uk/~ucakarc/work/software/randgen.txt

  implicit none
  integer*8 :: j,k,M,L,dim,seed,it,Nseed,i,control
  integer*8 :: DR,totalN,dI,dRa,dS,clock,lockday
  integer*8 ::steps,realiz,ii,jj,kk,dum,itmax,dum1,dum2,jjcontrol
  integer*8 ::Iseed,totalI
  double precision :: dt1,ba,tE,tI,Returnt
  double precision ::pIR,pEI,pSE,dtSIR,lockS
  integer*8 , dimension(:), allocatable :: Ia,Ip,Ips,H,ICU,Ih,Po,R,S,E,N,D,inc,Rdoc,cumul1,cumul2
  integer*8 , dimension(:), allocatable :: dIaLS,dSLS,DRLS,DELS,dPoLS,dIpsLS
  double precision , dimension(:), allocatable :: Neff,Neff2
  integer*8 :: sumS,sumE,sumIa,sumR
  double precision , dimension(:), allocatable ::Iaseed,Rseed,Sseed,Eseed
  double precision , dimension(:), allocatable ::dj1,dj2,outbreak
  double precision , dimension(:,:), allocatable ::dji1,dji2
  double precision , dimension(:,:), allocatable ::In
  integer*8 , dimension(:,:), allocatable ::con,con2
  double precision , dimension(:,:,:), allocatable ::Out
  double precision , dimension(4) :: CPUt
  double precision ts,tf
  double precision fpSE2,fpSE
  integer*8 i_dran
  double precision, dimension(:, :), ALLOCATABLE :: dumnet1,dumnet2
  !WARNING:: if len=50 is changed, check for coherence in code with subroutines.
  character(len=50) str,outfile,histfile,site,Readpopfile,Formato 
  character(len=50) nametripcom2,nametripcom1,popfile
  character*50 order
  INTEGER*4 status, system
  logical check,logdoc,lockdown
  integer ZBQLBIN !Binomial distribution
  double precision inj,outji

  !Read Parameters from file!
  read(*,*) jjcontrol !Realization number: each realiztion is labeled with an integer that comes from keyboard/script.
  open(unit=1002,&
  file="Parameters.dat",&
  status="old")
  read(1002,*) dim !total number of populations (e.g. zip codes,counties,countries,...)
  read(1002,*) dt1 !Integration time: In one iteration, the time moves dt1 forward. In general, dt1= 1 day. So one iteration= one day.
  read(1002,*) itmax !Maximum iteration
  read(1002,*) realiz !realizations
  read(1002,*) seed !Initial position of desease
  print *,M,dt1,itmax,realiz
  clock=100
  lockday=3200
  dtSIR=dt1/dble(clock) !Integration time for epidemic dynamics

  close(1002)

  !--------------------------Define (allocate) arrays-------------------------------
  allocate(In(4,dim))
  allocate(Out(4,dim,dim))
  allocate(con(dim,dim+1))
  allocate(con2(dim,dim+1))
  allocate(Ia(dim))
  allocate(R(dim))
  allocate(S(dim))
  allocate(E(dim))
  allocate(inc(dim))
  allocate(outbreak(dim))
  allocate(cumul1(dim))
  allocate(cumul2(dim))

  allocate(N(dim))
  allocate(Neff(dim))
  allocate(Neff2(dim))
  allocate(dj1(dim))
  allocate(dj2(dim))
  allocate(dji1(dim,dim))
  allocate(dji2(dim,dim))

  !-------------------------------------------------------------------------------------
  dum=time()+jjcontrol
  call dran_ini(dum) !Seed for random number generator
  call ZBQLINI(dum) !Initialize seed, if 0 then seed=time()
  write(*,*) "Seed of random numbers=",dum

  nametripcom1="networks/Commuting.dat" !File with commuting number of travelers per day period 1 (no lockdown)
  nametripcom2="networks/Commuting.dat" !File with commuting number of travelers per day period 2 (lockdown)
  popfile="networks/Population.dat" !File with number of inhabitants

  check=.false. !if .true. is detected->program ends

  !Read population file
  N=0
  call read_populations(dim,popfile,N)
  write(*,*) "Total number of inhabitants=",sum(N)

  !first read integers, then operate with them and convert to real
  allocate(dumnet1(dim,dim))
  allocate(dumnet2(dim,dim))
  call Create_network_dble(dim,dumnet1,nametripcom1) !commuting first period
  call Create_network_dble(dim,dumnet2,nametripcom2) !commuting second period


  !Look for negative or NAN values in mobility data:
  call check_topology_dble(dumnet1,dim,nametripcom1)
  call check_topology_dble(dumnet2,dim,nametripcom2)

  !Dji are the "(sigmas ji)/t" in reference (1).This is, (commutibg rate j->i)/ (return rate)
  Returnt=1.0d0/3.0d0 !mean return time (for commuting) (inverse of return time)
  do j = 1, dim, 1
    do i = 1, dim, 1
      dji1(j,i)=dble(dumnet1(j,i))*Returnt/(dble(N(j))) !period 1: Mobility without lockdown
      dji2(j,i)=dble(dumnet2(j,i))*Returnt/(dble(N(j))) !period 2: Mobility with lockdown
    end do
  end do

  deallocate(dumnet1)
  deallocate(dumnet2)

  !Compute total commuting out (sigma j in ref. (1))
  dj1=0.0d0
  dj2=0.0d0
  con=0
  con2=0
  do i = 1, dim, 1 !compute total exit rate
    do j = 1, dim, 1
      if ( dji1(i,j).ne.0 ) then !Connections without lockdown
        con(i,1)=con(i,1)+1 !degree of node i
        con(i,con(i,1)+1)=j !connections of i
      end if
      if ( dji2(i,j).ne.0 ) then !Connections with lockdown
        con2(i,1)=con2(i,1)+1 !degree of node i
        con2(i,con2(i,1)+1)=j !connections of i
      end if
      dj1(i)=dj1(i)+dji1(i,j)
      dj2(i)=dj2(i)+dji2(i,j)
    end do
  end do
  write(*,*) "Total flux out before/after lock",sum(dj1)/sum(dj2)

  !-------------------------------------------------------------------------------
  !Infection parameters------------------
  tI=3.8 !<Time to become R from I> (Infectious period)
  tE=3.7!average time in exposed regime
  k=5 !mean number of contacts (IRRELEVANT IN THIS MODEL AS LONG AS K>0)

  ba=1.6/(k*tI) !!Transmission rate (S->E rate)


  pIR=1.0d0-exp(-dtSIR/tI) !prob. thet I-individual becomes R in time dtSIR
  pEI=1.0d0-exp(-dtSIR/tE) !prob. thet E-individual becomes Po in time dtSIR
!-------------------------------------------------------------------------------

  CPUt=0.0d0 !Measure mean CPU time spent in different segments of the code

  control=0 !Number of super-critical realization (global outbreaks)

  lockS=0.7 !% of S individuals that can interact during lockdown


   do jj = 1, 1, 1
     lockdown=.False.
     !initialize measures
     logdoc=.true. !if false there are detected-infection cases
     !Initial condition
     call read_populations(dim,popfile,N) !I need to do this each realization as the number of inhabitants can fluctuate
     totalN=sum(N) !Total population
     S=N !Susceptible
     E=0 !Exposed (undocumented)
     inc=0 !Incidence
     Ia=0 !Infected
     R=0 !Recovered

     !Initial condition
     Ia(seed)=10
     do i = 1, 50, 1
       seed=i_dran(dim)
       Ia(seed)=1

       if ( Ia(seed).gt.N(seed) ) then !Super small populations!
         Ia(seed)=1
       end if
     end do

     S=S-E-Ia-R

     totalI=SUM(E)+SUM(Ia)
     Nseed=totalI

     it=0 !total number of iterations

     do i = 1, dim, 1 !check if there are negative initial populations
       call isitneg(S(i),check,11)
       call isitneg(E(i),check,22)
       call isitneg(Ia(i),check,44)
       call isitneg(R(i),check,88)

       if ( S(i).gt.N(i) ) then
         write(*,*) "fatal error! S>N","i=",i
         write(*,*) "Si=",S(i),"N(i)=",N(i),i
         write(*,*) "realiz=",jj,"it=",it
         goto 101
       end if
     end do

     if ( check ) then
       write(*,*) "The fatal error was in IC!..."
       write(*,*) "realization=",jj
       write(*,*) "Exiting program"
       goto 101 !Exit program
     end if


    open(unit=2002,& !Save seed
    file="Results/T/t_n"// &
    trim(str(jjcontrol))//".dat",&
    status="unknown")
    write(2002,*) "#it  S I R E inc inc/N ID"

    open(unit=3002,& !Save seed
    file="Results/sum/sumonzips_n"// &
    trim(str(jjcontrol))//".dat",&
    status="unknown")
    write(3002,*) "#it  S/N I/N R/N E/N inc/N"

    totalI=SUM(E)+SUM(Ia)

       !Compute effective populations: people remaining+visitors
       call ComputeNeffSEIR(dim,dj1,dji1,Neff,N)
       call ComputeNeffSEIR(dim,dj2,dji2,Neff2,N)

       if ( sum(Neff-N).gt.1000 ) then !if Neff lower than 98% of N
         write(*,*) "fatal error, sum(Neff).ne.sum(N)"
         write(*,*) "t=",it,"realiz=",jj
         write(*,*) sum(Neff),sum(N)
         write(*,*) shape(dji1)
         dum=0
         do j = 1, dim, 1
            dum=dum+outji(N(1),dj1(1),dji1(1,j))
          end do

          write(*,*) N(1),inj(N(1),dj1(1)),dum,dum+inj(N(1),dj1(1))
         goto 101
       end if

       if ( sum(Neff2-N).gt.1000 ) then !if Neff2 lower than 98% of N
         write(*,*) "fatal error, sum(Neff2).ne.sum(N)"
         write(*,*) "t=",it,"realiz=",jj
         write(*,*) sum(Neff2),sum(N)
         write(*,*) shape(dji1)
         dum=0
         do j = 1, dim, 1
            dum=dum+outji(N(1),dj1(1),dji1(1,j))
          end do

          write(*,*) N(1),inj(N(1),dj1(1)),dum,dum+inj(N(1),dj1(1))
         goto 101
       end if

     do while (( totalI.gt.0 ).and.(it.lt.itmax)) !Epidemic dynamics + update averages/print history

       it=it+1

       call cpu_time(ts)
       inc=0
       do ii=1,clock,1 !Integrate SEIR
         !Compute effective groups
         if ( .not.lockdown ) then !Effective  populations having into account commuting, see reference (1)
           !Compute effective groups
           do j=1,dim,1
             In(1,j)=S(j)/(1.0d0+dj1(j)) !Effected Susceptible people form j in j
             In(2,j)=Ia(j)/(1.0d0+dj1(j)) !Effected infected people form j in j

             do dum = 1, con(j,1), 1
                 l=con(j,dum+1)
                 Out(1,j,l) =S(j)*dji1(j,l)/(1.0d0+dj1(j)) !Effective Susceptible people from j in l
                 Out(2,j,l) =Ia(j)*dji1(j,l)/(1.0d0+dj1(j)) !Effective infective people from j in l

             end do
           enddo
         else
           !Compute effective groups
           do j=1,dim,1
             In(1,j)=S(j)/(1.0d0+dj2(j)) !Effected Susceptible people form j in j
             In(2,j)=Ia(j)/(1.0d0+dj2(j)) !Effected infected people form j in j

             do dum = 1, con2(j,1), 1
                 l=con2(j,dum+1)
                 Out(1,j,l) =S(j)*dji2(j,l)/(1.0d0+dj2(j)) !Effective Susceptible people from j in l
                 Out(2,j,l) =Ia(j)*dji2(j,l)/(1.0d0+dj2(j)) !Effective infective people from j in l

             end do
           enddo
         end if

         do j = 1, dim, 1
           DRa=0
           DS=0
           DI=0

           !---------___EXTRACT USING BINOMIALS ------------------------
           DRa=ZBQLBIN(Ia(j),pIR)!Binomials to extract possitive contribution of R, and negative contributions of I
           DI=ZBQLBIN(E(j),pEI) !Binomial to extract negativw contribution of E, negative contribution of I

           !Extract possitive contriburion to E, negative (and only) contribution to S
           if ( S(j).gt.0 ) then !Though the process of computing pSE I have to divide by S(j)
             if ( .not.lockdown ) then
               pSE=fpSE2(dim,dji1,dj1,S(j),In,Out,j,Neff,ba,k,dtSIR,it,jj,con)

               DS=ZBQLBIN(S(j),pSE)
             else
               pSE=fpSE2(dim,dji2,dj2,S(j),In,Out,j,Neff2,ba,k,dtSIR,it,jj,con2)
               DS=ZBQLBIN(int(lockS*S(j)),pSE)
             endif

           end if
           !----------------------------------------------------------------
           !UPDATE GROUPS---------------------------------------------------
           !update susceptible
           S(j)=S(j)-DS
           !update Exposed
           E(j)=E(j)+DS-DI
           !update Infected
           Ia(j)=Ia(j)+DI-DRa !Prevalence
           inc(j)=inc(j)+DI !Incidence
           !Update Recovered
           R(j)=R(j)+DRa

         end do
       end do

       totalI=SUM(E)+SUM(Ia) !If totalI=0 this is the last update

       if ( SUM(Ia)/dble(sum(N)).gt.0.005) then !There is an outbreak if prevalence reaches certain value
         lockdown=.True.
       end if

       do i = 1, dim, 1 !check if there are negative initial populations
         call isitneg(S(i),check,511)
         call isitneg(E(i),check,522)
         call isitneg(Ia(i),check,544)
         call isitneg(R(i),check,588)
       end do

       call cpu_time(tf)

       CPUt(1)=CPUt(1)+tf-ts !CPU time for SIR step

       !------------------------------------------------------------------------
       !Update averages/print/check-------------------------------------------
       call cpu_time(ts)

      do i = 1, dim, 1 !Print history
        write(2002,'(6(I15),ES16.6E3,I15)') it,S(i),Ia(i),R(i),E(i),inc(i),inc(i)/dble(N(i)),i
      end do
      sumS=sum(S)
      sumE=sum(E)
      sumIa=sum(Ia)
      sumR=sum(R)
      dum=sum(N)
      !write(3002,*) it*dt1,sumS/dble(dum),sumIa/dble(dum),sumR/dble(dum),&
      Formato='(I15,5(ES16.6E3))'
      write(3002,Formato) it,sumS/dble(dum),sumIa/dble(dum),sumR/dble(dum),&
      sumE/dble(dum),sum(inc)/dble(dum)
      !Check quantities (conservation, sign), exit program if Fatal error--------
      !--------------------------------------------------------------------------
      call check_negative_total(sumS,sumE,sumIa,sumR,dum,check) !Check sign of sum over groups

      if ( sum(N).ne.totalN ) then !Check if total number of inhabuitants is conserved
        write(*,*) "FATAL ERROR!"
        write(*,*) "Total number of inhabitants not conserved"
        write(*,*) "Ninitial=",totalN,"Nfinal=",dum
        write(*,*) "realization=",jj,"iteration=",it
        write(*,*) "Exiting program"
        goto 101 !Exit program
      end if


      if ( S(seed).gt.N(seed)) then
        write(*,*) "Fatal error, more S than population 1"
        write(*,*) "realization=",jj,"iteration=",it
        write(*,*) "Exiting program"
        write(*,*) S(seed)/dble(N(seed)),S(seed),N(seed)
        goto 101 !Exit program
      end if

      if ( check ) then
        write(*,*) "The fatal error was in..."
        write(*,*) "realization=",jj,"iteration=",it
        write(*,*) "Exiting program"
        goto 101 !Exit program
      end if
      !--------------------------------------------------------------------------------
      !------------------------End checking..........................................
      call cpu_time(tf)
      CPUt(3)=CPUt(3)+tf-ts !CPU time for updating averages and check variables

      end do !loop for Gillespie

      if ( sum(R)+totalI.lt.Nseed ) then
        write(*,*)"Fatal error! sum(R)+sum(D)+totalI.lt.Nseed "
        write(*,*) "realiz=",jj
        write(*,*) "R+D=",sum(R)+sum(D),"Nseed=",Nseed
        goto 101
      end if

      close(2002)
      close(3002)

  end do !end loop over different realizations


  CPUt=CPUt/dble(realiz*itmax)
  write(*,*) "CPU time on SIR:",CPUt(1)
  write(*,*) "CPU time on LS step (multinomials):",CPUt(2)
  write(*,*) "CPU time in averaging and checks",CPUt(3)
  !Open files to save averages over realizations
  open(unit=1,& !Averaged history of seed population
   file="Measures.dat",&
   status="unknown",&
   access="append")
  write(1,*) "lockdown=",lockdown,"realiz=",jjcontrol,"site=",site

  close (1)


  !--------------------------Deallocate (allocatable) arrays-------------------------------
  deallocate(Ia)
  deallocate(inc)
  deallocate(R)
  deallocate(S)
  deallocate(E)
  deallocate(N)
  deallocate(Neff)
  deallocate(Neff2)
  deallocate(dj1)
  deallocate(dj2)
  deallocate(dji1)
  deallocate(dji2)

  deallocate(outbreak)

  deallocate(In)
  deallocate(Out)
  deallocate(con)
  deallocate(con2)
  !-------------------------------------------------------------------------------------
  !System commands:
  order="gzip -f Results/T/t_n"// &
  trim(str(jjcontrol))//".dat"
  status = system( order )
  if ( status .ne. 0 ) then
     stop 'system: error'
  end if

  order="gzip -f Results/sum/sumonzips_n"// &
  trim(str(jjcontrol))//".dat"
  status = system( order )
  if ( status .ne. 0 ) then
     stop 'system: error'
  end if

  write(*,*) "Finished!!!:))"
101 end program Covid_model_Spain
