program test
    use iso_fortran_env, only : dp => real64, i4 => int32
    implicit none
    integer(i4),parameter :: N=6
    real(dp), dimension(N) :: phi

    write(*,*) 'Configuración inicial'
    call hot_start(phi,1._dp)
    write(*,*) phi(:)
    call cluster(phi)
    write(*,*) 'Configuración final'
    write(*,*) phi(:)

contains

  subroutine random_phi(x,bound)
    implicit none
    real(dp),intent(out) :: x
    real(dp), intent(in) :: bound
    real(dp) :: y
    call random_number(y)
    x = 2._dp*bound*y -bound
  end subroutine random_phi

  subroutine hot_start(phi,hotphi)
    real(dp), dimension(N), intent(out) :: phi
    real(dp), intent(in) :: hotphi
    integer(i4) :: i1
    do i1=1,N
        call random_phi(phi(i1),hotphi)
    end do
  end subroutine hot_start

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
          p=1._dp!-exp(-2._dp*beta )
          call random_number(r)
          bond_x(i)=(r<p)
        else
          bond_x(i)=.false.
        end if
    end do
    
    write(*,*) 'Enlaces'
    write(*,*) bond_x(:)

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
    
    write(*,*) 'Labels before union'
    write(*,*) label(:)

    if(bond_x(N) ) then
      call union(label(1),label(N),parent )
    end if
    do i=1,N
      label(i)=find(label(i),parent)
    end do
    
    write(*,*) 'Labels after union'
    write(*,*) label(:)

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

end program test
