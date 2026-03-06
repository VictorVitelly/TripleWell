module functions
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  implicit none

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

  function potential(m02,phi,i1)
    real(dp), intent(in) :: m02
    real(dp), dimension(:), intent(in) :: phi
    integer(i4), intent(in) :: i1
    real(dp) :: potential
    potential=(m02*phi(i1)**2+0.5_dp*c*phi(i1)**4)/2._dp
    !potential=m02*(phi(i1)**2) *((phi(i1)-c)*(phi(i1)+c) )**2 /2._dp
    !potential=(phi(i1)**6+2._dp*phi(i1)**4-2._dp*(2._dp*m02+1)*phi(i1)**2) /2._dp
    !potential=m02*(phi(i1)**2+1._dp) *((phi(i1)-c)*(phi(i1)+c) )**2 /2._dp
  end function potential

  function lagrangian(m02,phi,i1)
    real(dp), intent(in) :: m02
    real(dp), dimension(:), intent(in) :: phi
    integer(i4), intent(in) :: i1
    real(dp) :: lagrangian
    lagrangian= ( (phi(iv(i1+1))-phi(i1))/dt )**2 /2._dp &
              &+potential(m02,phi,i1)
  end function lagrangian

  function S(m02,phi)
    real(dp), intent(in) :: m02
    real(dp), dimension(:), intent(in) :: phi
    real(dp) :: S
    integer(i4) :: i1,Narr
    Narr=size(phi,dim=1)
    S=0._dp
    do i1=1,Narr
        S=S+dt*lagrangian(m02,phi,i1)
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
    DeltaS=dt*(DSa -DSb)
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

  recursive function find(x,parent) result(out)
    integer(i4), intent(in) :: x
    integer(i4), intent(inout) :: parent(:)
    integer(i4) :: out
    if(parent(x) /= x) then
      parent(x)=find(parent(x),parent )
    end if
    out=parent(x)
  end function find

  subroutine union(x,y,parent)
    integer(i4),intent(in) :: x,y
    integer(i4),intent(inout) :: parent(:)
    integer :: root_x, root_y
    root_x=find(x,parent)
    root_y=find(y,parent)
    if (root_x /= root_y) then
      parent(root_y)=root_x
    end if
  end subroutine union

end module functions
