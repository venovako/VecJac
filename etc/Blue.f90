! see la_constants.f90 from LAPACK
program Blue
  implicit none

  integer, parameter :: dp = kind(0d0)
  real(dp), parameter :: dtsml = real(radix(0._dp), dp)**ceiling((minexponent(0._dp) - 1) * 0.5_dp)
  real(dp), parameter :: dtbig = real(radix(0._dp), dp)**floor((maxexponent(0._dp) - digits(0._dp) + 1) * 0.5_dp)

  integer :: c, f

  c = ceiling((minexponent(0._dp) - 1) * 0.5_dp)
  print *, 'ceiling((', minexponent(0._dp), ' - 1) / 2) =', c
  f = floor((maxexponent(0._dp) - digits(0._dp) + 1) * 0.5_dp)
  print *, 'floor((', maxexponent(0._dp), ' -', digits(0._dp), ' + 1) / 2) =', f

  print *, 'dtsml =', dtsml, ' = 2**', c
  print *, 'dtbig =', dtbig, ' = 2**', f
end program Blue
