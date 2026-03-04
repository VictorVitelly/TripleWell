module statistics
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  implicit none

contains
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

  subroutine metropolis(m0,dphi,phi)
    real(dp), intent(in) :: m0,dphi
    real(dp), dimension(N), intent(inout) :: phi
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i4) :: i1
    do i1=1,N
        call random_phi(deltaphi,dphi)
        phi2=phi(i1)+deltaphi
        DS=DeltaS(m0,phi,i1,phi2)
        if(DS .le. 0._dp) then
          phi(i1)=phi2
        else
          call random_number(r)
          p=Exp(-DS)
          if(r < p ) then
            phi(i1)=phi2
          end if
        end if
    end do
  end subroutine metropolis

  subroutine montecarlo(m0,dphi,phi,AR)
    real(dp), intent(in) :: m0,dphi
    real(dp), dimension(N), intent(inout) :: phi
    real(dp), intent(out) :: AR
    real(dp) :: deltaphi,phi2,DS,r,p
    integer(i4) :: i1
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

  subroutine cluster(phi)
    real(dp), dimension(N),intent(inout) :: phi
    integer(i4), dimension(N) :: spin
    logical, dimension(N) :: bond_x
    integer(i4) :: i,label(N),parent(N),next_label,left_label
    logical, allocatable :: flip_cluster(:)
    real(dp) :: beta,r,p
    spin(:)=nint(sign(1._dp,phi(:)),i4)
    do i=1,N
        if(spin(i)==spin(mod(i,N)+1) ) then
          beta=abs(phi(i))*abs(phi(mod(i,N)+1))
          p=1._dp-exp(-2._dp*beta )
          call random_number(r)
          bond_x(i)=(r<p)
        else
          bond_x(i)=.false.
        end if
    end do

    label(:)=0
    do i=1,N
      parent(i)=i
    end do
    next_label=1
    left_label=0
    do i=1,N
        left_label=0
        if(i>1 .and. bond_x(i-1) ) then
          left_label=label(i-1)
        end if
        if(left_label==0) then
          label(i)=next_label
          next_label=next_label+1
        else if(left_label /= 0) then
          label(i)=left_label
        end if
    end do
    if(bond_x(N) ) then
      call union(label(1),label(N),parent )
    end if

    allocate(flip_cluster(next_label) )
    flip_cluster(:)=.false.
    do i=1,next_label-1
      call random_number(r)
      flip_cluster(i)=(r<0.5_dp)
    end do
    do i=1,N
      if(flip_cluster(label(i))) then
        phi(i)=-phi(i)
      end if
    end do
    deallocate(flip_cluster)

  end subroutine cluster

  subroutine cycles(m0,dphi,phi,AR)
    integer(i4) :: i
    real(dp), intent(in) :: m0,dphi
    real(dp), dimension(N), intent(inout) :: phi
    real(dp), intent(out) :: AR
    do i=1,3
      call metropolis(m0,dphi,phi)
    end do
    call montecarlo(m0,dphi,phi,AR)
    !call cluster(phi)
  end subroutine cycles
!ERROR STATISTICS

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
    call mean_0(x,y)
    !call standard_error(x,y,deltay)
    call jackknife(x,y,deltay)
  end subroutine mean_scalar

end module statistics
