module parameters
  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none
    integer(i4), parameter :: N=64,thermalization=5000,eachsweep=100,Nmsrs=100,Nmsrs2=120
    integer(i4), parameter :: Mbin(5)=(/4,5,10,15,20/),bins=201
    real(dp), parameter :: dt=1._dp,a=3.5_dp, c=1._dp, maxx=3._dp,minn=-3._dp,dphi=0.5_dp
    real(dp), parameter :: binwidth=(maxx-minn)/real(bins,dp)

end module parameters
