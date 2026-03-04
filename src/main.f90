program main
  use iso_fortran_env, only : dp => real64, i4 => int32
  use parameters
  use functions
  use statistics
  use measurements
  implicit none

    call vary_mu(0.0_dp,-3.0_dp,21)
    !call make_histogram(a)

end program main
