module parameters

! initialize to default values
integer :: nx = 100
double precision :: t = 0.0d+0
double precision :: dt = 0.5d-3
double precision :: tend = 100d+0
double precision :: dumpinterval = 1d+0
double precision :: dx
double precision :: Omega,Re,gamma,Ma,Pr
double precision :: pi
character(20) :: filename
character(75) :: outfile
character(20) :: spatial_disc
character(20) :: prefix
character(20) :: pistvel
character(20) :: bc
logical :: use_upwind

! Also include data here
double precision, dimension(:,:), allocatable :: u, v, rho, P, Temp
double precision, dimension(:,:), allocatable :: tau_xx,tau_yy,tau_xy,tau_yx,qx,qy
double precision, dimension(:,:), allocatable :: rho_u, rho_v, Et, Et_new, rho_new

! For moving immersed boundary
logical, dimension(:,:), allocatable :: coveredcells, freshlycleared, ghostcells
double precision, dimension(:), allocatable :: boundaryloc, xx, yy, boundary_u, boundary_v, boundary_rho
double precision :: F
integer :: ny

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
read(iunit, *) Re 
read(iunit, *) gamma
read(iunit, *) Ma
read(iunit, *) Pr
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
print *, 'Re', Re
print *, 'gamma', gamma
print *, 'Ma', Ma
print *, 'Pr', Pr

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
F = 1d+0 / Omega
ny = nx*(1+2*F)+1
print *, nx, F, ny

end subroutine read_params

