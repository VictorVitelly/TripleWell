program main
    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none

    integer(i4), parameter :: N=64,thermalization=5000,eachsweep=500,Nmsrs=200,Nmsrs2=120
    integer(i4), parameter :: Mbin(5)=(/4,5,10,15,20/),bins=201
    real(dp), parameter :: lambda0=1._dp, maxx=3._dp,minn=-3._dp,dphi=0.5_dp
    real(dp), parameter :: binwidth=(maxx-minn)/real(bins,dp)
    !call vary_mu(0.0_dp,-3.0_dp,21)
    call make_histogram(-5._dp)

contains

  function iv(i)
    integer(i4), intent(in) :: i
    integer(i4) :: iv
    if(i==N+1) then
      iv=1
    else if(i==0) then
      iv=N
    else
      iv=i
    end if
  end function

  function lagrangian(m02,phi,i1)
    real(dp), intent(in) :: m02
    real(dp), dimension(:), intent(in) :: phi
    integer(i4), intent(in) :: i1
    real(dp) :: lagrangian
    lagrangian=( (phi(iv(i1+1))-phi(i1) )**2  &
              &-m02*(phi(i1)**2+0.1_dp)**2 *((phi(i1)-1.4_dp)*(phi(i1)+1.4_dp) )**2)/2._dp
  end function lagrangian

  function S(m02,phi)
    real(dp), intent(in) :: m02
    real(dp), dimension(:), intent(in) :: phi
    real(dp) :: S
    integer(i4) :: i1,Narr
    Narr=size(phi,dim=1)
    S=0._dp
    do i1=1,Narr
        S=S+lagrangian(m02,phi,i1)
    end do
  end function S

  function DeltaS(m02,phi,i1,phi2)
    real(dp), intent(in) :: m02
    real(dp), dimension(:), intent(in) :: phi
    integer(i4), intent(in) :: i1
    real(dp), intent(in) :: phi2
    real(dp), dimension(size(phi)) :: phiy
    real(dp) :: DeltaS
    real(dp) :: DSa,DSb
    phiy(:)=phi(:)
    phiy(i1)=phi2
    DSa=lagrangian(m02,phiy,i1)+0.5_dp*(phi2-phi(iv(i1-1)))**2 
    DSb=lagrangian(m02,phi,i1)+0.5_dp*(phi(i1)-phi(iv(i1-1)))**2
    DeltaS=DSa -DSb
  end function DeltaS

  function mean(phi)
    real(dp), dimension(:), intent(in) :: phi
    integer(i4):: i1,Narr
    real(dp) :: mean
    Narr=size(phi,dim=1)
    mean=0._dp
    do i1=1,Narr
        mean=mean+phi(i1)
    end do
  end function mean

  subroutine random_phi(x,bound)
    real(dp),intent(out) :: x
    real(dp), intent(in) :: bound
    real(dp) :: y
    call random_number(y)
    x = 2._dp*bound*y -bound
  end subroutine random_phi
  
  subroutine hot_start(x,dphi)
  real(dp), dimension(:), intent(out) :: x
  real(dp), intent(in) :: dphi
  integer(i4) :: i
  do i=1,size(x)
    call random_phi(x(i),2._dp*dphi)
  end do
  end subroutine

  subroutine montecarlo(m0,dphi,phi,AR)
    real(dp), intent(in) :: m0,dphi
    real(dp), dimension(N), intent(inout) :: phi
    real(dp), intent(out) :: AR
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i4) :: i1,i2
    AR=0._dp
    do i1=1,N
        call random_phi(deltaphi,dphi)
        phi2=phi(i1)+deltaphi
        DS=DeltaS(m0,phi,i1,phi2)
        if(DS .le. 0._dp) then
          phi(i1)=phi2
          AR=AR+1._dp
        else
          call random_number(r)
          p=Exp(-DS)
          AR=AR+p
          if(r < p ) then
            phi(i1)=phi2
          end if
        end if
    end do
    AR=AR/real(N,dp)
  end subroutine montecarlo

  subroutine standard_error(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: variance
    integer(i4) :: k,Narr
    Narr=size(x)
    deltay=0._dp
    variance=0._dp
    do k=1,Narr
      variance=variance+(x(k) -y)**2
    end do
    variance=variance/real(Narr-1,dp)
    deltay=Sqrt(variance/real(Narr,dp))
  end subroutine standard_error

  subroutine jackknife(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(in) :: y
    real(dp), intent(out) :: deltay
    real(dp) :: jackk
    real(dp), allocatable :: xmean(:), delta_y(:)
    integer(i4) :: k,Narr,i,j
      Narr=size(x)
      allocate(delta_y(size(Mbin)))
      do j=1,size(Mbin)
        allocate(xmean(Mbin(j)))
        jackk=0._dp
        xmean=0._dp
        do i=1,Mbin(j)
          do k=1,Narr
            if(k .le. (i-1)*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            else if(k > i*Narr/Mbin(j)) then
              xmean(i)=xmean(i)+x(k)
            end if
          end do
          xmean(i)=xmean(i)/(real(Narr,dp) -real(Narr/Mbin(j),dp))
        end do
        do k=1,Mbin(j)
          jackk=jackk+(xmean(k)-y )**2
        end do
        delta_y(j)=Sqrt(real(Mbin(j)-1,dp)*jackk/real(Mbin(j),dp))
        deallocate(xmean)
      end do
      deltay=maxval(delta_y)
  end subroutine jackknife

  subroutine mean_0(x,y)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y
    integer(i4) :: k,Narr
    Narr=size(x)
    y=0._dp
    do k=1,Narr
      y=y+x(k)
    end do
    y=y/real(Narr,dp)
  end subroutine mean_0

  subroutine mean_scalar(x,y,deltay)
    real(dp), dimension(:), intent(in) :: x
    real(dp), intent(out) :: y,deltay
    integer(i4) :: k
    call mean_0(x,y)
    !call standard_error(x,y,deltay)
    call jackknife(x,y,deltay)
  end subroutine mean_scalar

  subroutine correlation(phi,corr1,corr2)
    real(dp), dimension(N), intent(in) :: phi
    real(dp), dimension(N), intent(inout) :: corr1
    real(dp), dimension(N,N), intent(inout) :: corr2
    real(dp) :: xx
    integer(i4) :: i1,i2
    do i1=1,N
      !corr1(i1)=corr1(i1)+abs(phi(i1))
      do i2=1,N
        corr2(i1,i2)=corr2(i1,i2)+(phi(i1)*phi(i2))
      end do
    end do
  end subroutine correlation

  subroutine vary_mu(mi,mf,Nps)
    integer(i4), intent(in) :: Nps
    real(dp), intent(in) :: mi,mf
    real(dp) :: phi(N),AR,m0
    integer(i4) :: i,i1,j,k
    real(dp) :: magnet(Nmsrs2),magnetps(Nmsrs2),action(Nmsrs2),arate(Nmsrs2)
    real(dp) :: magnet_ave,magnet_err,action_ave,action_err,arate_ave,arate_err
    real(dp) :: magnetps_ave,magnetps_err
    real(dp), allocatable :: corr1(:),corr2(:,:),CF(:,:),CF_ave(:,:),CF_delta(:,:)
    open(10, file = 'data/history.dat', status = 'replace')
    open(20, file = 'data/action.dat', status = 'replace')
    open(30, file = 'data/magnet.dat', status = 'replace')
    open(40, file = 'data/magnetps.dat', status = 'replace')
    open(60, file = 'data/corrfunc.dat', status = 'replace')
    allocate(corr1(N))
    allocate(corr2(N,N))
    allocate(CF(N,Nmsrs2))
    allocate(CF_ave(N,Nps))
    allocate(CF_delta(N,Nps))

    do k=1,Nps
      phi(:)=0._dp
      arate(:)=0._dp
      action(:)=0._dp
      magnet(:)=0._dp
      magnetps(:)=0._dp
      CF(:,:)=0._dp
      m0=mi+(mf-mi)*real(k-1,dp)/real(Nps-1,dp)
      do i=1,thermalization
          call montecarlo(m0,dphi,phi,AR)
          !write(10,*) i, S(m0,phi)
      end do
      do i=1,Nmsrs2
        corr1(:)=0._dp
        corr2(:,:)=0._dp
        phi(:)=-phi(:)
        do i1=1,Nmsrs
          do j=1,eachsweep
            call montecarlo(m0,dphi,phi,AR)
          end do
          arate(i)=arate(i)+AR
          action(i)=action(i)+S(m0,phi)
          magnet(i)=magnet(i)+abs(mean(phi))
          magnetps(i)=magnetps(i)+abs(phi(1))
          call correlation(phi,corr1,corr2)
        end do
        corr1(:)=corr1(:)/real(Nmsrs,dp)
        corr2(:,:)=corr2(:,:)/real(Nmsrs,dp)
        do i1=1,N
          CF(i1,i)=corr2(i1,1)!-(corr1(1)*corr1(i1))
        end do
      end do
      arate(:)=arate(:)/real(Nmsrs,dp)
      action(:)=action(:)/real(Nmsrs,dp)
      magnet(:)=magnet(:)/real(Nmsrs,dp)
      magnetps(:)=magnetps(:)/real(Nmsrs,dp)
      do i=1,N
          call mean_scalar(CF(i,:),CF_ave(i,k),CF_delta(i,k))
      end do

      call mean_scalar(arate,arate_ave,arate_err)
      call mean_scalar(action,action_ave,action_err)
      call mean_scalar(magnet,magnet_ave,magnet_err)
      call mean_scalar(magnetps,magnetps_ave,magnetps_err)
      write(*,*) m0,arate_ave,arate_err
      write(20,*) m0,action_ave/real(N,dp), action_err/real(N,dp)
      write(30,*) m0,magnet_ave/real(N,dp), magnet_err/real(N,dp)
      write(40,*) m0,magnetps_ave, magnetps_err
    end do

    do k=1,N+1
      write(60,*) abs(k-1), CF_ave(iv(k),:), CF_delta(iv(k),:)
    end do
    close(10)
    close(20)
    close(30)
    close(40)
    close(60)
    deallocate(corr1,corr2,CF,CF_ave,CF_delta)
  end subroutine vary_mu


  subroutine histogram(x,A1,A2)
  real(dp), dimension(:), intent(in) :: x
  integer(i4), dimension(bins), intent(inout) :: A1
  real(dp), dimension(bins), intent(in) :: A2
  integer(i4) :: i,j
  do i=1,bins
    do j=1,size(x,dim=1)
        if(x(j) .le. real(A2(i),dp)+binwidth/2._dp .and. x(j)>real(A2(i),dp)-binwidth/2._dp ) then
          A1(i)=A1(i)+1
          cycle
        end if
    end do
  end do
  end subroutine histogram

  subroutine histogram2(x,A1,A2)
  real(dp), intent(in) :: x
  integer(i4), dimension(bins), intent(inout) :: A1
  real(dp), dimension(bins), intent(in) :: A2
  integer(i4) :: i
  do i=1,bins
    if(x .le. real(A2(i),dp)+binwidth/2._dp .and. x>real(A2(i),dp)-binwidth/2._dp ) then
      A1(i)=A1(i)+1
    cycle
    end if
  end do
  end subroutine histogram2

  subroutine make_histogram(m0)
    real(dp), intent(in) :: m0
    real(dp) :: phi(N),norm,AR
    integer(i4) :: i,j,k
    real(dp), allocatable :: A2(:)
    integer(i4), allocatable :: A1(:)
    real(dp) :: arate(Nmsrs2), ar_ave, ar_err
    open(50, file = 'data/histogram.dat', status = 'replace')
    allocate(A1(bins))
    allocate(A2(bins))
    !phi(:)=0._dp
    call hot_start(phi,dphi)
    arate(:)=0._dp
    do i=1,bins
      A2(i)=minn+binwidth/2._dp+real(i-1,dp)*binwidth
    end do
    A1=0

    do i=1,thermalization
      call montecarlo(m0,dphi,phi,AR)
    end do

    do i=1,Nmsrs2
      do j=1,Nmsrs
        do k=1,eachsweep
          call montecarlo(m0,dphi,phi,AR)
        end do
        arate(i)=arate(i)+AR
        !call histogram(phi,A1,A2)
        !call histogram2(mean(phi)/real(N,dp),A1,A2)
        call histogram2(phi(1),A1,A2)
      end do
    end do

    norm=0._dp
    do i=1,bins
      norm=norm+A1(i)
    end do
    norm=norm*(real(maxx-minn,dp) )/real(bins,dp)
    do i=1,bins
      write(50,*) A2(i), A1(i)/norm, sqrt( real(A1(i),dp) )/norm
    end do
    
    arate(:)=arate(:)/real(Nmsrs,dp)
    call mean_scalar(arate,ar_ave,ar_err)
    write(*,*) m0,ar_ave,ar_err
    
    deallocate(A1,A2)
    close(50)
  end subroutine make_histogram

end program main
