module parameters

! initialize to default values
integer :: nx = 100
double precision :: t = 0.0d+0
double precision :: dt = 0.5d-3
double precision :: tend = 100d+0
double precision :: dumpinterval = 1d+0
double precision :: dx,dy
double precision :: Omega, R, gamma, lambda, L1, L2
double precision :: pi
character(20) :: filename
character(75) :: outfile
character(20) :: spatial_disc
character(20) :: pistvel
character(20) :: prefix
character(20) :: bc
logical :: use_upwind

! Also include data here
double precision, dimension(:), allocatable :: x,y,v_hat
double precision, dimension(:,:), allocatable :: u, v, rho, P, Temp
double precision, dimension(:,:), allocatable :: tau_xx,tau_yy,tau_xy,tau_yx,qx,qy
double precision, dimension(:,:), allocatable :: rho_u, rho_v, Et, rho_new, Et_new

end module parameters


! =========================================== !
! subroutine to get  inputs from file 'input' !
! =========================================== !

subroutine read_params

use parameters

implicit none

integer :: iunit

iunit = 4

! read variables 
open(UNIT=iunit, FILE=trim(filename))
read(iunit, *)
read(iunit, *)
read(iunit, *) Omega
read(iunit, *) R
read(iunit, *) lambda
read(iunit, *) gamma
read(iunit, *) L1
read(iunit, *) L2
read(iunit, *)
read(iunit, *) 
read(iunit, *) nx
read(iunit, *) dt
read(iunit, *) tend
read(iunit, *) dumpinterval
read(iunit, *) spatial_disc
read(iunit, *) pistvel
read(iunit, *) prefix
read(iunit, *) bc
close(iunit)
use_upwind = (trim(spatial_disc).eq.'upwind')
   

! print variables
print *, 'nx', nx
print *, 'dt', dt
print *, 'tend', tend
print *, 'dump interval', dumpinterval
print *, 'Omega', Omega
print *, 'R', R
print *, 'lambda', lambda
print *, 'gamma', gamma
print *, 'L1', L1
print *, 'L2', L2

! name output file
outfile = trim(prefix) // '_' // trim(spatial_disc) 
write(outfile,"(A,I0)"), trim(outfile) // '_nx', nx
write(outfile,"(A,ES8.2E2)") trim(outfile) // '_dt', dt
write(outfile,"(A,ES8.2E2)") trim(outfile) // '_tend', tend

!outfile = 'testing.out'
print *, 'outfile  ', outfile


! calculate some other variables
pi = 4*atan(1d+0)
dx = 1d+0 / dble(nx)
dy = L1 / dble(nx)

print *, 'where'

end subroutine read_params

