program Covid_model_Spain
  !SEIIIIR model with two kind of mobility patterns: Commuting and long stays
  !Commuting coded in force of infection following (1).
  !Long stays are simulated as physical movement of agents.
  !All-to-all interactions.
  !The time evolution is divided in two periods of the same length. The 1st (2nd) period has the mobility computed from the week of the 4th (18th) of March
  !Variables are commented when initilized for first time
  !-------------------------------------------------------------------------------
  !------------Code developed by Javier Aguilar----------------------------------
  !------------Help/comments: Sandro Meloni, Jose J. Ramasco and Pere Colet-------
  !-------------------------------------------------------------------------------
  !REFERENCES:
  !(1)  BALCAN, Duygu, et al. Modeling the spatial spread of infectious diseases:
  ! The GLobal Epidemic and Mobility computational model. Journal of computational
  ! science, 2010, vol. 1, no 3, p. 132-145.

  implicit none
  integer*8 :: j,k,M,L,dim,seed,it,Nseed,i,control,a
  integer*8 :: DR,dIh,dIp,dIa,DPo,dLSji,DH,DD,DICU,DRps,DRH,DRICU,DDh,DDICU,dIps
  integer*8 :: totalN,dI,dRa,dRIh,dS,clock,lockday
  integer*8 ::steps,realiz,ii,jj,kk,dum,itmax,dum1,dum2,jjcontrol
  integer*8 ::Iseed,SLseed,RLseed,kmaxL,totalI
  double precision :: dt1,bh,ba,bp,bpo,PUnDoc,Pmild,tE,tI,tPo,Returnt,tdoc
  double precision ::pIR,pEPo,pSE,pPoI,dtSIR,pSElam,pSEBin,pSE2,lockS,Picu,ppauci
  double precision :: pHD,pHR,pICUR,pICUD,rHD,rHR,rICUR,rICUD,bps
  integer*8 , dimension(:), allocatable :: Ia,Ip,Ips,H,ICU,Ih,Po,R,S,E,N,D,inc,Rdoc,cumul1,cumul2
  integer*8 , dimension(:), allocatable :: dIaLS,dSLS,DRLS,DELS,dPoLS,dIpsLS
  integer*8 , dimension(:,:), allocatable :: tripLS1,tripLS2
  double precision , dimension(:), allocatable :: Neff,Neff2
  integer*8 sumS,sumE,sumIa,sumIp,sumIh,sumR,sumPo,sumH,sumICU,sumD,SumIps
  double precision , dimension(:), allocatable ::Iaseed,Ipseed,Ihseed,Rseed,Sseed,Eseed,Poseed,Hseed,ICUseed
  double precision , dimension(:), allocatable :: Ipsseed,Dseed
  double precision , dimension(:), allocatable ::dumIaseed,dumIpseed,dumIhseed,dumRseed,dumSseed,dumEseed
  double precision , dimension(:), allocatable :: dumPoseed,dumICUseed,dumHseed,dumDseed,dumIpsseed
  double precision , dimension(:), allocatable ::Iat,Ipt,Iht,Rt,St,Et,Pot,ICUt,Ht,Dt,Ipst
  double precision , dimension(:), allocatable ::dj1,dj2,outbreak
  double precision , dimension(:,:), allocatable ::Iax,Ipx,Ihx,Rx,Sx,Ex,Pox,Hx,ICUx,Dx,Ipsx
  double precision , dimension(:,:), allocatable :: dumIax,dumIpx,dumIhx,dumRx,dumSx,dumEx
  double precision , dimension(:,:), allocatable :: dumPox,dumHx,dumICUx,dumDx,dumIpsx !To store mean history all provinces
  double precision , dimension(:), allocatable ::dumIat,dumIpt,dumIht,dumRt,dumSt,dumEt,dumPot,dumDt,dumIpst !To store mean history of total population
  double precision , dimension(:), allocatable :: dumICUt,dumHt
  double precision , dimension(:,:), allocatable ::dji1,dji2
  double precision , dimension(:,:), allocatable ::In
  integer*8 , dimension(:,:), allocatable ::con,con2
  double precision , dimension(:,:,:), allocatable ::Out
  double precision , dimension(4) :: CPUt
  double precision ts,tf
  double precision fpSE2,fpSE
  integer*8 i_dran
  double precision, dimension(:, :), ALLOCATABLE :: dumnet1,dumnet2
  character(len=50) str,outfile,histfile,site,Readpopfile,Formato
  character(len=50) nametripLS1,nametripLS2,nametripcom2,nametripcom1,popfile
  character*50 order
  INTEGER*4 status, system
  logical check,logdoc,lockdown
  integer ZBQLBIN !Binomial distribution
  double precision inj,outji

  !read(*,"(I4.4,I2,A3,A11,A21,A11)") jjcontrol,a,site,nametripcom1,Readpopfile,nametripcom2
  read(*,*) jjcontrol,a,site,nametripcom1,Readpopfile,nametripcom2
  write(*,*) jjcontrol,a,site,nametripcom1,Readpopfile,nametripcom2
  site=trim(ADJUSTL(site))
  nametripcom1=trim(ADJUSTL(nametripcom1))
  nametripcom2=trim(ADJUSTL(nametripcom2))
  readpopfile=trim(ADJUSTL(readpopfile))
  !Read Parameters from file!
  open(unit=1002,&
  file="Parameters_"//trim(site)//".dat",&
  status="old")
  read(1002,*) dim !total number of populations
  read(1002,*) dt1 !Integration time: Each iteration time moves dt1 forward
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
  allocate(Ip(dim))
  allocate(Ih(dim))
  allocate(R(dim))
  allocate(S(dim))
  allocate(E(dim))
  allocate(inc(dim))
  allocate(Po(dim))
  allocate(H(dim))
  allocate(ICU(dim))
  allocate(Rdoc(dim))
  allocate(D(dim))
  allocate(Ips(dim))
  allocate(outbreak(dim))
  allocate(cumul1(dim))
  allocate(cumul2(dim))

  allocate(N(dim))
  allocate(dIaLS(dim))
  allocate(dSLS(dim))
  allocate(dELS(dim))
  allocate(dRLS(dim))
  allocate(dPoLS(dim))
  allocate(dIpsLS(dim))
  allocate(Neff(dim))
  allocate(Neff2(dim))
  allocate(dj1(dim))
  allocate(dj2(dim))
  allocate(tripLS1(dim,dim))
  allocate(tripLS2(dim,dim))
  allocate(dji1(dim,dim))
  allocate(dji2(dim,dim))

  !-------------------------------------------------------------------------------------
  dum=time()+jjcontrol
  call dran_ini(dum) !Seed for random number generator
  call ZBQLINI(dum) !Initialize seed, if 0 then seed=time()
  write(*,*) "Seed of random numbers=",dum

!  nametripLS1="LStrips_average_4_march.dat" !File with long stay number of travelers per day period 1
!  nametripLS2="LStrips_18_march.dat" !File with long stay number of travelers per day period 2

  nametripcom1="networks/"//trim(nametripcom1) !File with commuting number of travelers per day period 1

  nametripcom2="networks/"//trim(nametripcom2) !File with commuting number of travelers per day period 2
  popfile="networks/Population_PS_"//trim(site)//".dat" !File with number of inhabitants

  check=.false. !if .true. is detected->program ends

  !Read population file
  N=0
  call read_populations(dim,popfile,N)
  write(*,*) "Total number of inhabitants=",sum(N)

!dumnet is not the array that will contain the commuting information, I have to proceed like this since
!I have problems to read integer data as real (Segmentation fault (core dumped)). So I first read integers, then operate with them and convert to real
  allocate(dumnet1(dim,dim))
  allocate(dumnet2(dim,dim))
  call Create_network_dble(dim,dumnet1,nametripcom1) !commuting first period
  call Create_network_dble(dim,dumnet2,nametripcom2) !commuting second period
!  call Create_network(dim,tripLS1,nametripLS1) !Long stays first period
  !call Create_network(dim,tripLS2,nametripLS2) !long stays second period

  !Look for negative or NAN values in mobility data:
  !call check_topology(tripLS1,dim,nametripLS1)
  !call check_topology(tripLS2,dim,nametripLS2)
  call check_topology_dble(dumnet1,dim,nametripcom1)
  call check_topology_dble(dumnet2,dim,nametripcom2)
  tripLS1=0
  tripLS2=0
  !Dji are the "(sigmas ji)/t" in reference (1).This is, (commutibg rate j->i)/ (return time)
  Returnt=1.0d0/3.0d0 !mean return time (for commuting)
  do j = 1, dim, 1
    do i = 1, dim, 1
      dji1(j,i)=dble(dumnet1(j,i))*Returnt/(dble(N(j))) !period 1
      dji2(j,i)=dble(dumnet2(j,i))*Returnt/(dble(N(j))) !period 2
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
      if ( dji1(i,j).ne.0 ) then
        con(i,1)=con(i,1)+1 !degree of node i
        con(i,con(i,1)+1)=j !connections of i
      end if
      if ( dji2(i,j).ne.0 ) then
        con2(i,1)=con2(i,1)+1 !degree of node i
        con2(i,con2(i,1)+1)=j !connections of i
      end if
      dj1(i)=dj1(i)+dji1(i,j)
      dj2(i)=dj2(i)+dji2(i,j)
    end do
  end do
  write(*,*) "Total flux out before/after lock",sum(dj1)/sum(dj2)
  !goto 101
  !-------------------------------------------------------------------------------

  !Infection parameters------------------
  tI=3.8 !<Time to become R from I> (Infectious period)
  k=5 !mean number of contacts
  bp=1.6/(k*tI) !Transmission rate from S to E through Hospitalized (documented) infections
  bh=1.6/(k*tI) ! infection rate through Ih
  bps=1.6/(k*tI) ! infection rate through Ips
  ba=1.6/(k*tI)! infection rate through Ia
  bpo=1.6/(k*tI)! infection rate through Po

  PunDoc=1.0 !Prob of becoming  aymptotic UnDocumented infected from E
  ppauci=0.2 !Prob of becoming paucisintomatic undocumented from E
  Pmild=0.6 ! Prob of becoming mild documented infection from E
  Picu=0.28 !if severe, prob of going to ICU

  tE=3.7!average time in exposed regime
  tPo=0.00000001 !average rime in podromic phase


  !Nseed=50 !original number of aymptotic infected+Exposed

  rHR=0.047 !rate H->R
  rHD=0.009 !rate H->D
  rICUR=0.048 !rate ICU->R
  rICUD=0.018 !rate ICU->D

  !realiz=10  !number of realizations

  pIR=1.0d0-exp(-dtSIR/tI) !prob. thet I-individual becomes R in time dtSIR
  pPoI=1.0d0-exp(-dtSIR/tPo)!prob. thet Po-individual becomes I in time dtSIR
  pEPo=1.0d0-exp(-dtSIR/tE) !prob. thet E-individual becomes Po in time dtSIR
  pHR=1.0d0-exp(-dtSIR*rHR) !prob. H->R in dtSIR (provided that H is constant in dtSR)
  pHD=1.0d0-exp(-dtSIR*rHD) !prob. H->D in dtSIR
  pICUR=1.0d0-exp(-dtSIR*rICUR) !prob. ICU->R in dtSIR (provided that ICU is constant in dtSR)
  pICUD=1.0d0-exp(-dtSIR*rICUD) !prob. ICU->D in dtSIR (provided that ICU is constant in dtSR)

  tdoc=0.0 !mean time for first detection (creation of Ih or Ip)
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
     inc=0
     Po=0 !Infected individuals in Podromic phase (UnDocumented)
     Ia=0 !Asymptomatic infected (undocumented, can move)
     Ips=0 !paucisymptomatic infected (undocumented, they can move)
     Ih=0 !Household infected (documented, can't move)
     Ip=0 !Severe infected (documented, can't move)
     H=0 !Hospitalized not ICU
     ICU=0 !Hospitalized in ICU
     D=0 ! Dead
     R=0 !Recovered
     cumul1=0 !cumulative of active
     cumul2=0 !cumulative of active
     !seed=20 !Seed in madrid
     !seed=i_dran(dim) ! We assume that the disease starts in one random population (seed) on undocumented groups

     Ia(seed)=10
     do i = 1, 50, 1
       seed=i_dran(dim)
       Ia(seed)=1

       if ( Ia(seed).gt.N(seed) ) then !Super small populations!
         Ia(seed)=1
       end if
     end do
     write(*,*) "realiz=",jj,seed
     S=S-E-Po-Ips-Ia-Ih-Ip-H-ICU-D-R

     totalI=SUM(E)+SUM(Po)+SUM(Ip)+SUM(Ia)+SUM(Ih)+SUM(Ips)+SUM(H)+SUM(ICU)
     Nseed=totalI

     it=0 !total number of iterations

     do i = 1, dim, 1 !check if there are negative initial populations
       call isitneg(S(i),check,11)
       call isitneg(E(i),check,22)
       call isitneg(Po(i),check,33)
       call isitneg(Ia(i),check,44)
       call isitneg(Ih(i),check,55)
       call isitneg(Ip(i),check,66)
       call isitneg(Ips(i),check,77)
       call isitneg(R(i),check,88)
       call isitneg(D(i),check,99)
       call isitneg(ICU(i),check,1010)

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
    file=trim(site)//"/T_"//trim(str(a))//"/t_zip_n"// &
    trim(str(jjcontrol))//".dat",&
    status="unknown")
    write(2002,*) "#it  S I R E inc inc/N ID"

    open(unit=3002,& !Save seed
    file=trim(site)//"/sum/sumonzip_a"//trim(str(a))//"_n"// &
    trim(str(jjcontrol))//".dat",&
    status="unknown")
    write(3002,*) "#it  S/N I/N R/N E/N inc/N"

    totalI=SUM(E)+SUM(Po)+SUM(Ip)+SUM(Ia)+SUM(Ih)+SUM(Ips)+SUM(H)+SUM(ICU)


       call ComputeNeff(dim,dj1,dji1,Neff,N,Ip,Ih,ICU,D,H)

       call ComputeNeff(dim,dj2,dji2,Neff2,N,Ip,Ih,ICU,D,H)


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

     do while (( totalI.gt.0 ).and.(it.lt.itmax)) !Epidemic dynamics + diffusion + updating averages

       it=it+1
       !write(*,*) "it=",it,"Iseed=",Ia(seed),"Itotal=",sum(Ia)

       !write(*,*) "totalN-sum(N)=",totalN-sum(N)
       call cpu_time(ts)
       inc=0
       do ii=1,clock,1 !Integrate SEIIIIR with more integration steps than mobility
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
           DRIh=0
           DRps=0
           pSE=0
           DS=0
           DPo=0
           DH=0
           DDICU=0
           DDH=0
           DICU=0
           DRH=0
           DRICU=0
           !Binomials to extract possitive contribution of R, and negative contributions of I
           DRa=ZBQLBIN(Ia(j),pIR)
           DRps=ZBQLBIN(Ips(j),pIR)
           DRIh=ZBQLBIN(Ih(j),pIR)
           !Binomials to extract possitive contribution to H and ICU, and another negative contribution to I
           DH=ZBQLBIN(Ip(j),pIR)
           DICU=ZBQLBIN(DH,pICU)
           DH=DH-DICU
           !Binomials to extract positive contribution to documented R
           DRH=ZBQLBIN(H(j),pHR)
           DRICU=ZBQLBIN(ICU(j),pICUR)
           !Binomials to extract positive contribution to D
           DDICU=ZBQLBIN(ICU(j)-DRICU,pICUD)
           DDH=ZBQLBIN(H(j)-DRH,pHD)
           !--------------------------------------------------------------------
           DPo=ZBQLBIN(E(j),pEpo) !Binomial to extract positive contribution of Po, negative contribution of E
           DI=ZBQLBIN(Po(j),pPoI) !Binomial to extract positive contribution of I, negative contribution of Po

           !Extract possitive contriburion to E, negative (and only) contribution to S
           if ( S(j).gt.0 ) then !Though the process of computing pSE I have to divide by S(j)
             if ( .not.lockdown ) then
               pSE=fpSE2(dim,dji1,dj1,S(j),In,Out,j,Neff,ba,bp,bh,bpo,bps,k,dtSIR,it,jj,con)

               DS=ZBQLBIN(S(j),pSE)
             else
               pSE=fpSE2(dim,dji2,dj2,S(j),In,Out,j,Neff2,ba,bp,bh,bpo,bps,k,dtSIR,it,jj,con2)
               DS=ZBQLBIN(int(lockS*S(j)),pSE)
             endif
              !pSE2=1-exp(-bp*(Ia(j)+Po(j))*dtSIR/N(j))

           end if
           !----------------------------------------------------------------
           !Update groups
           !update susceptible
           S(j)=S(j)-DS
           !update Exposed
           E(j)=E(j)+DS-DPo
           inc(j)=inc(j)+DPo
           !Update podromic
           !Po(j)=Po(j)-DI+DPo
           DI=DPo
           !UPDATE Infected------------------------
           !infected asymptomatic
           DIa=ZBQLBIN(DI,PUnDoc)
           Ia(j)=Ia(j)+DIa-DRa
           !Infected paucisintomatic
           DIps=ZBQLBIN(DI-DIa,Ppauci)
           Ips(j)=Ips(j)+DIps-DRps

           !Infected mild (household)
           DIh=ZBQLBIN(DI-DIa-DIps,Pmild/(1.0-Ppauci)) !Possitive contribution to Documented household (mild infected I)

           Ih(j)=Ih(j)+DIh-DRIh
           !infected severe HOSPITAL AND ICU
           Ip(j)=Ip(j)+(DI-DIh-DIa-DIps)-DICU-DH
           !-----------------------------------------
           !update ICU
           ICU(j)=ICU(j)+DICU-DRICU-DDICU
           !Update H
           H(j)=H(j)+DH-DRH-DDH
           !Update D
           D(j)=D(j)+DDH+DDICU
           !Update recovers
           R(j)=R(j)+DRa+DRIh+DRps+DRICU+DRH
           Rdoc(j)=Rdoc(j)+DRIh+DRICU+DRH

         end do
       end do

       totalI=SUM(E)+SUM(Po)+SUM(Ip)+SUM(Ia)+SUM(Ih)+SUM(Ips)+SUM(H)+SUM(ICU) !If totalI=0 this is the last update
       !print *,(totalI/dble(sum(N)))
       if ( SUM(Ia)/dble(sum(N)).gt.0.005) then
 !        write(*,*) "YEEEEEEEEEEEEEEES"
         lockdown=.True.
       end if

       do i = 1, dim, 1 !check if there are negative initial populations
         call isitneg(S(i),check,511)
         call isitneg(E(i),check,522)
         call isitneg(Po(i),check,533)
         call isitneg(Ia(i),check,544)
         call isitneg(Ih(i),check,555)
         call isitneg(Ip(i),check,566)
         call isitneg(Ips(i),check,577)
         call isitneg(R(i),check,588)
         call isitneg(D(i),check,599)
         call isitneg(ICU(i),check,51010)

       end do

       call cpu_time(tf)

       CPUt(1)=CPUt(1)+tf-ts !CPU time for SIR step

       !------------------------------------------------------------------------
       !Update averages-------------------------------------------
       call cpu_time(ts)



       if ( (sumIh+sumIp).gt.0 ) then !measure first-detection time
         if ( logdoc ) then !One measure per realization
           logdoc=.false.
           tdoc=tdoc+it*dt1
         end if
       end if


      !Save trajectories-------------------------------------------
      cumul1=cumul1+Ip+Ih+ICU+H
      cumul2=cumul2+Ip+ICU+H

      do i = 1, dim, 1
        write(2002,'(6(I15),ES16.6E3,I15)') it,S(i),Ia(i),R(i),E(i),inc(i),inc(i)/dble(N(i)),i
      end do
      sumS=sum(S)
      sumE=sum(E)
      sumIa=sum(Ia)
      sumIh=sum(Ih)
      sumIp=sum(Ip)
      sumR=sum(R)
      sumPo=sum(Po)
      sumICU=sum(ICU)
      sumH=sum(h)
      sumIps=sum(Ips)
      sumD=sum(D)
      dum=sum(N)
      !write(3002,*) it*dt1,sumS/dble(dum),sumIa/dble(dum),sumR/dble(dum),&
      Formato='(I15,5(ES16.6E3))'
      write(3002,Formato) it,sumS/dble(dum),sumIa/dble(dum),sumR/dble(dum),&
      sumE/dble(dum),sum(inc)/dble(dum)
      !Check quantities (conservation, sign), exit program if Fatal error--------
      !--------------------------------------------------------------------------
      call check_negative_total(sumS,sumE,& !Check sign of sum over groups
          sumIa,sumIh,sumIp,sumR,sumPo,sumIps,sum(N),check)


      if ( sum(N).ne.totalN ) then !Check if total number of inhabuitants is conserved
        write(*,*) "FATAL ERROR!"
        write(*,*) "Total number of inhabitants not conserved"
        write(*,*) "Ninitial=",totalN,"Nfinal=",dum
        write(*,*) "realization=",jj,"iteration=",it
        write(*,*) "Exiting program"
        goto 101 !Exit program
      end if

      do i = 1, dim, 1 !check sign of groups that can't move (the rest of them arechecked in distribution step)
        call isitneg(Ip(i),check,555) !556 is the label of this error
        call isitneg(Ih(i),check,444) !445 is the label of this error
      end do

      if ( S(seed).gt.N(seed)) then
        write(*,*) "Fatal error, more S than population 1"
        write(*,*) "realization=",jj,"iteration=",it
        write(*,*) "Exiting program"
        write(*,*) S(seed)/dble(N(seed)),S(seed),N(seed)
        write(*,*) dSLS(seed),dELS(seed),dIaLS(seed),dRLS(seed),dPoLS(seed),DIpsLS(seed)
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

      if ( sum(R)+sum(D)+totalI.lt.Nseed ) then
        write(*,*)"Fatal error! sum(R)+sum(D)+totalI.lt.Nseed "
        write(*,*) "realiz=",jj
        write(*,*) "R+D=",sum(R)+sum(D),"Nseed=",Nseed
        goto 101
      end if

      close(2002)
      close(3002)

      if ( dble(R(seed)).ge.20 ) then !Average over super-critical samples. Super-critical if R>threshold.
        control=control+1
      end if
  end do !end loop over different realizations

  write(*,*) "Probability of outbreak=",dble(control)/realiz

  !compute averages on total population
  call read_populations(dim,popfile,N)


  tdoc=tdoc/realiz

  CPUt=CPUt/dble(realiz*itmax)
  write(*,*) "CPU time on SIR:",CPUt(1)
  write(*,*) "CPU time on LS step (multinomials):",CPUt(2)
  write(*,*) "CPU time in averaging and checks",CPUt(3)
  !Open files to save averages over realizations
  open(unit=1,& !Averaged history of seed population
   file="Measures.dat",&
   status="unknown",&
   access="append")
  write(1,*) "lockdown=",lockdown,"realiz=",jjcontrol,"a=",a,"site=",site

  close (1)


  !--------------------------Deallocate (allocate) arrays-------------------------------
  deallocate(Ia)
  deallocate(inc)
  deallocate(Ip)
  deallocate(Ih)
  deallocate(R)
  deallocate(S)
  deallocate(E)
  deallocate(Po)
  deallocate(Ips)
  deallocate(H)
  deallocate(D)
  deallocate(ICU)
  deallocate(N)
  deallocate(dIaLS)
  deallocate(dSLS)
  deallocate(dELS)
  deallocate(dRLS)
  deallocate(dPoLS)
  deallocate(dIpsLS)
  deallocate(Neff)
  deallocate(Neff2)
  deallocate(dj1)
  deallocate(dj2)
  deallocate(tripLS1)
  deallocate(tripLS2)
  deallocate(dji1)
  deallocate(dji2)

  deallocate(outbreak)
  deallocate(cumul1)
  deallocate(cumul2)
  deallocate(In)
  deallocate(Out)
  deallocate(con)
  deallocate(con2)
  !-------------------------------------------------------------------------------------
  !System commands:
  order="gzip -f "//trim(site)//"/T_"//trim(str(a))//"/t_zip_n"// &
  trim(str(jjcontrol))//".dat"
  status = system( order )
  if ( status .ne. 0 ) then
     stop 'system: error'
  end if

  order="gzip -f "//trim(site)//"/sum/sumonzip_a"//trim(str(a))//"_n"// &
  trim(str(jjcontrol))//".dat"
  status = system( order )
  if ( status .ne. 0 ) then
     stop 'system: error'
  end if

  write(*,*) "Finished!!!:))"
101 end program Covid_model_Spain
