!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!LAST UPDATE: wed-15-05-2020!!!!!!!!!!!!!!!

character(len=50) function str(k)
  implicit none
!   "Convert an integer to string."
  integer, intent(in) :: k
  write (str, *) k !write to a string
  str = adjustl(str)
end function str

subroutine initialize_meas(dim,S,E,Ia,Ih,Ip,R,Po,Ips,H,D,ICU)
  !Initialize averages of 1-d arrays
  implicit none
  integer*8, intent(in) :: dim
  double precision,dimension (0:dim), intent(inout) :: S,E,Ia,Ih,Ip,R,Po,Ips,H,D,ICU
  S=0.0d0
  E=0.0d0
  Ia=0.0d0
  Ih=0.0d0
  Ip=0.0d0
  R=0.0d0
  Po=0.0d0
  Ips=0.0d0
  H=0.0d0
  D=0.0d0
  ICU=0.0d0

  return
end subroutine initialize_meas

subroutine initialize_meas_2d(dim1,dim2,S,E,Ia,Ih,Ip,R,Po,Ips,H,D,ICU)
  !Initialize averages of 2d arrays
  implicit none
  integer*8, intent(in) :: dim1,dim2
  double precision,dimension (dim1,0:dim2), intent(inout) :: S,E,Ia,Ih,Ip,R,Po,Ips,H,D,ICU

  S=0.0d0
  E=0.0d0
  Ia=0.0d0
  Ih=0.0d0
  Ip=0.0d0
  R=0.0d0
  Po=0.0d0
  Ips=0.0d0
  H=0.0d0
  D=0.0d0
  ICU=0.0d0


  return
end subroutine initialize_meas_2d

subroutine Update_average_seed(S,E,Ia,Ih,Ip,R,Po,Ips,H,D,ICU,dS,dE,dIa,dIh,dIp,dR,dPo,dIps,dH,dD,dICU)
  !cummulative of DENSITY (double precision) of seed along realization in order to compute average
  implicit none
  double precision, intent(in) :: dS,dE,dIa,dIh,dIp,dR,dPo,dIps,dH,dD,dICU
  double precision, intent(inout) :: S,E,Ia,Ih,Ip,R,Po,Ips,H,D,ICU

  S=S+dS
  E=E+dE
  Ia=Ia+dIa
  Ih=Ih+dIh
  Ip=Ip+dIp
  R=R+dR
  Po=Po+dPo
  Ips=Ips+dIps
  H=H+Dh
  D=D+dD
  ICU=ICU+dICU

  return
end subroutine Update_average_seed

subroutine Update_average(S,E,Ia,Ih,Ip,R,Po,Ips,H,D,ICU,dS,dE,dIa,dIh,dIp,dR,dPo,dIps,dH,dD,dICU)
  !cummulative of sum of goups along realization in order to compute average
  implicit none
  integer*8, intent(in) :: dS,dE,dIa,dIh,dIp,dR,dPo,dIps,dH,dD,dICU
  double precision, intent(inout) :: S,E,Ia,Ih,Ip,R,Po,Ips,H,D,ICU

  S=S+dS
  E=E+dE
  Ia=Ia+dIa
  Ih=Ih+dIh
  Ip=Ip+dIp
  R=R+dR
  Po=Po+dPo
  Ips=Ips+dIps
  H=H+Dh
  D=D+dD
  ICU=ICU+dICU

  return
end subroutine Update_average

double precision function Inj(Xj,dj)
  implicit none
  integer*8, intent(in) :: Xj
  double precision, intent(in):: dj
  !Efective Compartiment population Xj that stays in j
  Inj = Xj/(1+dj)
end function Inj

double precision function Outji(Xj,dj,dji)
!Effective compartiment population Xj that go from j to i
  implicit none
  double precision, intent(in) :: dj,dji
  integer*8, intent(in) :: Xj
  Outji= Xj*dji/(1+dj)
end function Outji

subroutine ComputeNeff(dim,dj,dji,Neff,N,Ip,Ih,ICU,D,H)
  !Compute effective populations (people staying+incomming people)
  implicit none
  integer*8, intent(in) :: dim
  integer*8 :: ii,jj,dum
  integer*8, intent(in) :: Ip(dim),Ih(dim),N(dim),ICU(dim),D(dim),H(dim)
  double precision, intent(out) :: Neff(dim)
  double precision, dimension (:),intent(in) :: dj(dim),dji(dim,dim)
  double precision :: Inj,Outji

  Neff=0
  do ii = 1, dim, 1
    dum=Ip(ii)+Ih(ii)+H(ii)+D(ii)+ICU(ii) !people that will not go out of population ii
    Neff(ii)=dble(dum)+Inj(N(ii)-dum,dj(ii)) !People from ii staying at ii
    do jj = 1, dim, 1 !people from jj going to ii
      dum=N(jj)-Ip(jj)-Ih(jj)-H(jj)-D(jj)-ICU(jj)! people from jj that could travel
      Neff(ii)=Neff(ii)+Outji(dum,dj(jj),dji(jj,ii))
    end do
  end do
end subroutine ComputeNeff

subroutine ComputeNeffSEIR(dim,dj,dji,Neff,N)
  !Compute effective populations (people staying+incomming people)
  implicit none
  integer*8, intent(in) :: dim
  integer*8 :: ii,jj
  integer*8, intent(in) :: N(dim)
  double precision, intent(out) :: Neff(dim)
  double precision, dimension (:),intent(in) :: dj(dim),dji(dim,dim)
  double precision :: Inj,Outji

  Neff=0
  do ii = 1, dim, 1
    Neff(ii)=Inj(N(ii),dj(ii)) !People from ii staying at ii
    do jj = 1, dim, 1 !people from jj going to ii
      Neff(ii)=Neff(ii)+Outji(N(jj),dj(jj),dji(jj,ii))
    end do
  end do
end subroutine ComputeNeffSEIR

double precision function fpSE(dim,dji,dj,Sj,Ip,Ih,Ia,Po,Ips,j,Neff,ba,bp,bh,bpo,bps,k,dt,t,realiz)
!Compute DS on population j having into account commuting, this is:
!Having into account the Sjj that stays in j and get infected by Ijj(I-native of j) and Iij(I-visitors in j)
!and the Sji that leave j for a certain time and come back than can get infected from Iii (I-native of i) and Imi (I-visitors in i).
  implicit none
  integer*8 :: i,l
  integer*8, intent(in) :: dim,Sj,j,k,t,realiz
  double precision, intent(in) :: Neff(dim)
  integer*8,dimension(:), intent(in) :: Ia(dim),Ip(dim),Ih(dim),Po(dim),Ips(dim)
  double precision :: lji,ljj,expo,dum
  double precision, intent(in) :: ba,bp,bh,bps,dt,bpo
  double precision,dimension(:), intent(in) :: dj(dim)
  double precision,dimension(:,:), intent(in) :: dji(dim,dim)
  double precision :: Inj,Outji

  expo=0 !  prob=1-exp(-k*expo)

  do i = 1, dim, 1
    lji=0
    if ( dji(j,i).gt.0 ) then
      lji=lji+Inj(Ia(i),dj(i))*ba+Inj(Po(i),dj(i))*bpo+Inj(Ips(i),dj(i))*bps+Ip(i)*bp+Ih(i)*bh !Infecttions in neighbourhood of j
      do l = 1, dim, 1
        if ( dji(l,i).gt.0 ) then !Infecttions in neighbourhood of j due to visits to the neighborhood
          dum=Outji(Ia(l),dj(l),dji(l,i))*ba+Outji(Po(l),dj(l),dji(l,i))*bpo+Outji(Ips(l),dj(l),dji(l,i))*bps
          lji=lji+dum
        endif
      end do
      expo=expo+lji*Outji(Sj,dj(j),dji(j,i))/(Sj*Neff(i))
    end if
  end do

  !Infections in j
  ljj=Inj(Ia(j),dj(j))*ba+Inj(Po(j),dj(j))*bpo+Inj(Ips(j),dj(j))*bps+Ip(j)*bp+Ih(j)*bh
  do i = 1, dim, 1 ! I have to do second loop and if becouse dij could not be symmetric (i could have dij=0 and dji!=0)
    if ( dji(i,j).gt.0 ) then
      ljj=ljj+Outji(Ia(i),dj(i),dji(i,j))*ba+Outji(Po(i),dj(i),dji(i,j))*bpo+Outji(Ips(i),dj(i),dji(i,j))*bps
    endif
  end do

  expo=expo+ljj*Inj(Sj,dj(j))/(Sj*Neff(j))
  fpSE=1.0-exp(-k*expo*dt)

  if ( fpSE.ne.fpSE ) then
    write(*,*) "FATAL ERROR pSE=NaN"
    write(*,*) "iteration=",t
    write(*,*) "realization",realiz
  else if (fpSE.lt.0) then
    write(*,*) "FATAL ERROR pSE<0"
    write(*,*) "iteration=",t
    write(*,*) "realization",realiz
  else if (fpSE.gt.1) then
    write(*,*) "FATAL ERROR pSE>1"
    write(*,*) "iteration=",t
    write(*,*) "realization",realiz
  end if

end function fpSE

double precision function fpSE2(dim,dji,dj,Sj,In,Out,j,Neff,ba,k,dt,t,realiz,con)
!Compute DS on population j having into account commuting, this is:
!Having into account the Sjj that stays in j and get infected by Ijj(I-native of j) and Iij(I-visitors in j)
!and the Sji that leave j for a certain time and come back than can get infected from Iii (I-native of i) and Imi (I-visitors in i).
  implicit none
  integer*8 :: i,l,dum1,dum2
  integer*8, intent(in) :: dim,Sj,j,k,t,realiz
  double precision, intent(in) :: dt !*
  double precision, intent(in) :: Neff(dim)
  double precision :: lji,ljj,expo,dum
  double precision, intent(in) :: ba
  double precision,dimension(:), intent(in) :: dj(dim)
  double precision,dimension(:,:), intent(in) :: dji(dim,dim)
  integer*8,dimension(:,:),intent(in) :: con(dim,dim+1)
  double precision,dimension(:,:), intent(in) :: In(4,dim)
  double precision,dimension(:,:,:), intent(in) :: Out(4,dim,dim)
  double precision :: Inj,Outji

  expo=0 !  prob=1-exp(-k*expo)

  do dum1 = 1, con(j,1), 1
    i=con(j,dum1+1)
    lji=0
        lji=lji+In(2,i)*ba!Infecttions in neighbourhood of j
        do dum2 = 1, con(i,1), 1
            l=con(i,dum2+1)
           !Infecttions in neighbourhood of j due to visits to the neighborhood
            dum=Out(2,l,i)*ba
            lji=lji+dum

        end do
        expo=expo+lji*Out(1,j,i)/(Sj*Neff(i))

  end do

  !Infections in j
  ljj=In(2,j)*ba
  do dum1 = 1, con(j,1), 1
    i=con(j,dum1+1) !label of neighbour
    ljj=ljj+Out(2,i,j)*ba

  end do

  expo=expo+ljj*In(1,j)/(Sj*Neff(j))
  fpSE2=1.0-exp(-k*expo*dt)

  if ( fpSE2.ne.fpSE2 ) then
    write(*,*) "FATAL ERROR pSE=NaN"
    write(*,*) "iteration=",t
    write(*,*) "realization",realiz
  else if (fpSE2.lt.0) then
    write(*,*) "FATAL ERROR pSE<0"
    write(*,*) "iteration=",t
    write(*,*) "realization",realiz
  else if (fpSE2.gt.1) then
    write(*,*) "FATAL ERROR pSE>1"
    write(*,*) "iteration=",t
    write(*,*) "realization",realiz
  end if

end function fpSE2

subroutine check_topology(net,N,name)
  !Check if mobility data is negative or NAN
  !Entries are integer
  implicit none
  integer*4 :: i,j
  integer*8, intent(in):: N
  integer*8,dimension(:,:), intent(in) :: net(N,N)
  character(len=50), intent(in):: name

  write(*,*) "dim=",N,name
  do i = 1, N, 1
    do j = 1, N, 1
      if (net(i,j).ne.net(i,j)) then
        write(*,*) "FATAL ERROR (net(i,j)=NaN)",i,j,name
      else if (net(i,j).lt.0) then
        write(*,*) "FATAL ERROR (net(i,j)<0)",i,j,name,net(i,j)
      endif
    end do
  end do

end subroutine check_topology

subroutine check_topology_dble(net,N,name)
  !Check if mobility data is negative or NAN
  !Entries are double precision
  implicit none
  integer*4 :: i,j
  integer*8, intent(in):: N
  double precision,dimension(:,:), intent(in) :: net(N,N)
  character(len=50), intent(in):: name

  write(*,*) "dim=",N,name
  do i = 1, N, 1
    do j = 1, N, 1
      if (net(i,j).ne.net(i,j)) then
        write(*,*) "FATAL ERROR (net(i,j)=NaN)",i,j,name
      else if (net(i,j).lt.0) then
        write(*,*) "FATAL ERROR (net(i,j)<0)",i,j,name,net(i,j)
      endif
    end do
  end do

end subroutine check_topology_dble

integer*8 function map2D (x,y,N)
  integer*8, intent(in) :: x,y,N
  !x in [1,N]
  !y in [1,M] M>=N

  map2D=x+(y-1)*N !in [1,M*N]
  !x=resto(label-1,N)+1 in [1,N]
  !y=((label-1)/N)+1 in [1,M], M>=N
end function

subroutine read_populations(dim,name,N)
  !Read number of inhabitants of the dim populations
  !Format of file: dim rows, each row with one population number
  implicit none
  integer*8, intent(in) :: dim
  integer*8,dimension(:), intent(out) :: N(dim)
  character(len=50), intent(in) :: name
  !N=0
  open(unit=1001, file=name,status="old", action="read")
  read(1001,*) N
  close (1001)
end subroutine read_populations

subroutine Create_network(N,net,name)
!Read mobility data di-->j.
!Data is integer
!format: the quantity di-->j is in row i column j of the file
  implicit none

  integer*8 :: i,j
  integer*8,intent(in) :: N
  integer*8,intent(out) :: net(N,N)
  character(len=50),intent(in) :: name

  net=0

  open(unit=1001, file=name,status="old", action="read")

  read(1001,*) net
  net=transpose(net)

  close(1001)

end subroutine Create_network

subroutine Create_network_dble(N,net,name)
!Read mobility data di-->j.
!Data is real
!format: the quantity di-->j is in row i column j of the file
  implicit none

  integer*8 :: i,j
  integer*8,intent(in) :: N
  double precision,intent(out) :: net(N,N)
  character(len=50),intent(in) :: name

  net=0

  open(unit=1001, file=name,status="old", action="read")

  read(1001,*) net
  net=transpose(net)

  close(1001)

end subroutine Create_network_dble

subroutine check_negative_mov(S,E,Ia,R,Po,Ips,N,out)
  !Checks if groups that can move have negative value after distributing travelers
  !if so, print message error
  implicit none
  integer*8 i
  integer*8,intent(in) :: S,E,Ia,R,N,Po,Ips
  logical, intent(out) :: out

  call isitneg(S,out,771)
  call isitneg(E,out,772)
  call isitneg(Ia,out,773)
  call isitneg(R,out,774)
  call isitneg(N,out,775)
  call isitneg(Po,out,776)
  call isitneg(Ips,out,777)
  if ( out ) then
    write(*,*) "------------------------------"
    write(*,*) "Subroutine check_negative_mov:"
    write(*,*) "Fatal Error in Multinomials!!"
    write(*,*) "------------------------------"
  end if

  return
end subroutine check_negative_mov

subroutine check_negative_total(S,E,Ia,R,N,out)
  !Checks if sum over groups is negative (e.g. is S=sum(S)<0?)
  !if so, print message error
  implicit none
  integer*8 i
  integer*8,intent(in) :: S,E,Ia,R,N
  logical, intent(out) :: out

  call isitneg(S,out,881) !LABEL DOES NOT WORK!!
  call isitneg(E,out,882)
  call isitneg(Ia,out,883)
  call isitneg(R,out,886)
  call isitneg(N,out,887)
  if ( out ) then
    write(*,*) "Fatal Error in sum!! Negative people!!",i
  end if

  return
end subroutine check_negative_total

subroutine isitneg(a,out,label)
  !checks if a is negative, if so, write error message and error tracker is set to true
  !isitneg is called several times, there is a label to identificate the origin of the error
  implicit none
  logical, intent(out) :: out
  integer*8, intent(in) :: a
  integer*4, intent(in) ::label

  if ( a.lt.0 ) then
    write(*,*) "-------------------"
    write(*,*) "Subroutine isitneg:"
    write(*,*) "Fatal error definite possitive quantity is < 0",a
    write(*,*) "Label of the error:",label
    write(*,*) "-------------------"
    out=.true.
  end if
  return
end subroutine isitneg

subroutine fixmax(N,NMax)
  !Checks if N>NMax if so, N=NMax
  implicit none
  integer*8,intent(in) :: NMax
  integer*8,intent(inout) :: N

  if ( N>NMax ) then
    N=NMax
  end if

  return
end subroutine fixmax

subroutine  multinomial(D,S,E,Ia,R,Po,Ips,dim,N,i,j,dS,dE,dIa,dR,dPo,dIps)
  !Distribute D individuals in 4 boxes according to
  !multinomial distribution
  !different probabilities (X/N) for each box
  implicit none
  integer*8 dum,Neff
  integer*8,intent(in) :: dim,i,j,S,E,Ia,R,Po,Ips,N
  integer*8,intent(inout) :: D
  double precision :: p,ps,pE,pIa,pBin,pR,pPo,pIps,norm
  integer*8,dimension(:) ,intent(inout) :: dS(dim),dE(dim),dIa(dim),dR(dim),dPo(dim),dIps(dim)
  integer ZBQLBIN !Binomial distribution

  Neff=S+E+Ia+R+Po+Ips !People that can move
  if ( D.gt.Neff ) then
    write(*,*) "FATAL ERROR, FLUX>POPULATION"
    D=Neff !maximum value
  end if

  !Probabilities of choosing a type of individual
  ps=dble(S)/dble(Neff)
  pIa=dble(Ia)/dble(Neff)
  pE=dble(E)/dble(Neff)
  pR=dble(R)/dble(Neff)
  pPo=dble(Po)/dble(Neff)
  pIps=dble(Ips)/dble(Neff)
  norm=ps+pIa+pE+pR

  if ( pR.lt.0 ) then
    write(*,*) "fatal error pmultin<0",pR,S,E,Ia,R
  end if

  if ( norm.eq.pR ) then
    call fixmax(D,R+dR(j))
    dR(i)=dR(i)+D
    dR(j)=dR(j)-D
    D=0

  else
    dum=ZBQLBIN(D,pR)
    call fixmax(dum,R+dR(j))
    dR(i)=dR(i)+dum
    dR(j)=dR(j)-dum
    D=D-dum


    if ( D.gt.0 ) then !If D=0 I generate NaN in pBin
      pBin=pIa/(1-pR)
      if ( pBin.gt.1.0d0 ) then
        pBin=1.0d0
        call fixmax(D,Ia+dIa(j))
        dIa(i)=dIa(i)+D
        dIa(j)=dIa(j)-D
        D=0
      else
        dum=ZBQLBIN(D,pBin)
        call fixmax(dum,Ia+dIa(j))
        dIa(i)=dIa(i)+dum
        dIa(j)=dIa(j)-dum
        D=D-dum
      end if

      if (D.gt.0) then!If D=0 I generate NaN in pBin
        pBin=pE/(1.0d0-pR-pIa)
        if ( pBin.gt.1 ) then

          pBin=1.0d0
          call fixmax(D,E+dE(j))
          dE(i)=dE(i)+D
          dE(j)=dE(j)-D
          D=0
        else
          dum=ZBQLBIN(D,pBin)
          call fixmax(dum,E+dE(j))
          dE(i)=dE(i)+dum
          dE(j)=dE(j)-dum
          D=D-dum
        end if

        if (D.gt.0) then!If D=0 I generate NaN in pBin
          pBin=pPo/(1.0d0-pR-pIa-pE)
          if ( pBin.gt.1 ) then

            pBin=1.0d0
            call fixmax(D,Po+dPo(j))
            dPo(i)=dPo(i)+D
            dPo(j)=dPo(j)-D
            D=0
          else
            dum=ZBQLBIN(D,pBin)
            call fixmax(dum,Po+dPo(j))
            dPo(i)=dPo(i)+dum
            dPo(j)=dPo(j)-dum
            D=D-dum
          end if
          if ( D.gt.0 ) then
            pBin=pIps/(1.0d0-pR-pIa-pE-pPo)
            if ( pBin.gt.1 ) then
              pBin=1.0d0
              call fixmax(D,Ips+dIps(j))
              dIps(i)=dIps(i)+D
              dIps(j)=dIps(j)-D
              D=0
            else
              dum=ZBQLBIN(D,pBin)
              call fixmax(dum,Ips+dIps(j))
              dIps(i)=dIps(i)+dum
              dIps(j)=dIps(j)-dum
              D=D-dum
            end if
          end if
        endif
      end if
    end if
  end if

  call fixmax(D,S+dS(j))
  if ( D.lt.0 ) then
     write(*,*)"FAAAAATAL ERROR. D<0!!!!!!"
  end if

  if ( S.gt.0 ) then
    dS(i)=dS(i)+D !positive increment for distributiion step
    dS(j)=dS(j)-D !people extracted
  end if

  return
end subroutine multinomial
