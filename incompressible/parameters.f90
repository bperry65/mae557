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
character(20) :: solver
character(20) :: prefix

! Also include data here
double precision, dimension(:,:), allocatable :: u,v,H1,H2,DP
double precision, dimension(:), allocatable :: p,dH
integer :: N, KL, KU, NRHS, LDAB, LD, INFO, LDB
integer, dimension(:), allocatable :: IPIV
doubleprecision, dimension(:,:), allocatable :: DP1

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
read(iunit, *) solver
read(iunit, *) prefix
close(iunit)

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
outfile = trim(prefix) // '_' // trim(solver) 
write(outfile,"(A,I0)"), trim(outfile) // '_nx', nx
write(outfile,"(A,ES8.2E2)") trim(outfile) // '_dt', dt
write(outfile,"(A,ES8.2E2)") trim(outfile) // '_tend', tend

!outfile = 'testing.out'
print *, 'outfile  ', outfile


! calculate some other variables
pi = 4*atan(1d+0)
dx = 1d+0 / dble(nx)

end subroutine read_params

