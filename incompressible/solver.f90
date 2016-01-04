program solver

  use parameters
  
  implicit none
  doubleprecision :: tdump, tprint
  integer :: i,j
  
  ! read inputs and get ready
  call getarg(1,filename)
  call read_params
  call initialize

  tdump = dumpinterval
  tprint = tend/1000
  
  ! run simulation
  do while (t.lt.(tend - 1d-10))
     call timestep

     if (t.ge.tprint) then
        print *, 't', t
        print *, 'u', u(10,19)
        tprint = tprint + tend/1000
     end if
     if (t.ge.tdump - 1d-10) then
        call dumpdata
        tdump = tdump + dumpinterval
     end if
  end do
  
 ! do i=1,nx
 !    do j =1,nx
 !       print *, 'divg u at i,j', u(i+1,j+1) - u(i,j+1) + v(i+1,j+1) - v(i+1,j), i,j
 !    end do
 ! end do

!  print *, 'H1', H1(nx,nx+1), H1(nx+1,nx+1)
!  print *, 'u', u(nx,nx+1), u(nx+1,nx+1)
!  print *, 'H2', H2(nx+1,nx), H2(nx+1,nx+1)
!  print *, 'v', v(nx+1,nx), v(nx+1,nx+1)
!  print *, 'p', p(nx*(nx-1)), p(nx*nx), p(nx*nx-1)
  
  
endprogram solver


! ============================================================!
!                        Initialization                       !
! ============================================================!

subroutine initialize

  use parameters
  implicit none

  integer :: i, iunit=8

  
  allocate(u(nx+1,nx+2))
  allocate(v(nx+2,nx+1))
  allocate(H1(nx+1,nx+2))
  allocate(H2(nx+2,nx+1))
  allocate(dH(nx*nx))
  allocate(p(nx*nx))
  
  ! initialize and apply initial conditions
  t = 0d+0
  u = 0d+0
  v = 0d+0
  H1 = 0d+0
  H2 = 0d+0
  dH = 0d+0
  p = 0d+0
  DP = 0d+0

  ! generate DP matrix
  NRHS = 1
  KL = nx
  KU = nx
  N = nx*nx
  LDAB = 2*KL+KU+1
  LDB = N
  INFO = 0
  allocate(DP(LDAB,N))
  allocate(DP1(LDAB,N))
  allocate(IPIV(N))
  DP(KL+1,:) = 1d+0
  DP(KL+KL,:) = 1d+0
  DP(KL*2+1,:) = -4d+0
  DP(KL*2+2,:) = 1d+0
  DP(KL*2+1+KU,:) = 1d+0
  DP(KL+KL+1,1) = 0d+0 ! Check and maybe delete later
  do i=1,nx
     DP(KL+KL,(i-1)*nx+1) = 0d+0
     DP(KL+KL+2,i*nx) = 0d+0
     DP(KL+KL+1,(i-1)*nx+1) = DP(KL+KL+1,(i-1)*nx+1) + 1d+0
     DP(KL+KL+1,i*nx) = DP(KL+KL+1,i*nx) + 1d+0
     DP(KL+KL+1,i) = DP(KL+KL+1,i) + 1d+0
     DP(KL+KL+1,nx*nx+1-i) = DP(KL+KL+1,nx*nx+1-i) + 1d+0
  end do
  DP = DP/dx/dx
  DP1 = DP
     
  ! apply boundary condition on top wall
  call apply_vel_bc
 ! print *, 'u', u
 ! print *, 'v', v
  
  ! intialize output file
  open(iunit,FILE=trim(outfile))
  write(iunit, *) nx
  close(iunit)
  
end subroutine initialize



! ============================================================!
!                        Run Timestep                         !
! ============================================================!

subroutine timestep

  use parameters
  implicit none

  integer :: i,j
  
  t = t + dt
  
  call calculateH

  call DGBSV(n, KL, KU, NRHS, DP, LDAB, IPIV, dH, LDB, INFO)

  DP = DP1
  p = dH
  
  ! Advance U in time
  do j=2,nx+1
     do i=2,nx
        u(i,j) = u(i,j) + dt/Omega*(H1(i,j) - (p(nx*(j-2)+i) - p(nx*(j-2)+i-1))/dx )
     end do
  end do

  do i=2,nx+1
     do j=2,nx
        v(i,j) = v(i,j) + dt/Omega*(H2(i,j) - (p(nx*(j-1)+i-1) - p(nx*(j-2)+i-1))/dx )
     end do
  end do
  
  call apply_vel_bc

end subroutine timestep



! ============================================================= !
!                           velocity BCs                        !
! ============================================================= !

subroutine apply_vel_bc

  use parameters
  implicit none
  integer :: i
  
  do i = 2,nx
     u(i,1) = -u(i,2)
     u(i,nx+2) = -u(i,nx+1) + 2d+0 *sin(t)
     v(1,i) = -v(2,i)
     v(nx+2,i) = -v(nx+1,i)
  end do
  
end subroutine apply_vel_bc


! ============================================================= !
!                    subroutine to calculate H
! ============================================================= !

subroutine calculateH

  use parameters
  implicit none

  integer :: i,j
  doubleprecision :: d1,d2,d3,d4

  ! Calculate H1
  do j = 2,nx+1
     do i = 2,nx
        if (i.eq.1) then
           d1 = ( u(i+1,j)**2d+0 - u(i,j)**2d+0 ) /dx
	   d3 = ( u(i+2,j) - 2d+0*u(i+1,j) + u(i,j) ) / dx**2d+0
        else if (i.eq.nx+1) then
           d1 = ( u(i,j)**2d+0 - u(i-1,j)**2d+0 ) /dx
	   d3 = ( u(i-2,j) - 2d+0*u(i-1,j) + u(i,j) ) / dx**2d+0
        else
           d1 = ( u(i+1,j)**2d+0 - u(i-1,j)**2d+0 ) /(2d+0*dx)
	   d3 = ( u(i+1,j) - 2d+0*u(i,j) + u(i-1,j) ) / dx**2d+0
        end if
        d2 = ( (u(i,j+1)+u(i,j))*(v(i,j)+v(i+1,j)) - (u(i,j-1)+u(i,j))*(v(i,j-1)+v(i+1,j-1)) ) / (4d+0*dx)
	d4 = ( u(i,j+1) - 2d+0*u(i,j) + u(i,j-1) ) / dx**2d+0
        H1(i,j) = -d1 - d2 + 1d+0/Re*(d3 + d4)
     end do
  end do

  !Calculate H2
  do i = 2,nx+1
     do j = 2,nx ! double check this
        if (j.eq.1) then
           d1 = ( v(i,j+1)**2d+0 - v(i,j)**2d+0 ) /dx
           d3 = ( v(i,j+2) - 2d+0*v(i,j+1) + v(i,j) ) / dx**2d+0
        else if (j.eq.nx+1) then
           d1 = ( v(i,j)**2d+0 - v(i,j-1)**2d+0 ) /dx
           d3 = ( v(i,j-2) - 2d+0*v(i,j-1) + v(i,j) ) / dx**2d+0
        else
           d1 = ( v(i,j+1)**2d+0 - v(i,j-1)**2d+0 ) /(2d+0*dx)
           d3 = ( v(i,j+1) - 2d+0*v(i,j) + v(i,j-1) ) / dx**2d+0
        end if
        d2 = ( (v(i+1,j)+v(i,j))*(u(i,j)+u(i,j+1)) - (v(i-1,j)+v(i,j))*(u(i-1,j)+u(i-1,j+1)) ) / (4d+0*dx)
        d4 = ( v(i+1,j) - 2d+0*v(i,j) + v(i-1,j) ) / dx**2d+0
        H2(i,j) = -d1 - d2 + 1d+0/Re*(d3 + d4)
     end do
  end do

  !Calculate dH
  do j = 1,nx
     do i = 1,nx
        dH((j-1)*nx + i) = (H1(i+1,j+1) - H1(i,j+1) + H2(i+1,j+1) - H2(i+1,j) ) / dx
     end do
  end do

  
end subroutine calculateH




! Subroutine to dump data

subroutine dumpdata

  use parameters
  implicit none

  integer :: i,j,iunit=8
  doubleprecision :: pp, uu, vv
  open(iunit, file=trim(outfile),position='APPEND')
  write(iunit,*) 'time', t, 'nondimensional'
  write(iunit,*) 'P ','U ', 'V '
  do i = 1,nx+1
     do j = 1,nx+1
        uu = 0.5d+0 * (u(i,j) + u(i,j+1))
        vv = 0.5d+0 * (v(i,j) + v(i+1,j))
        if ((i.ne.1) .and. (i.ne.nx+1) .and. (j.ne.1) .and. (j.ne.nx+1)) then
           pp = 0.25d+0 * ( p(nx*(j-2)+i-1) + p(nx*(j-2)+i) + p(nx*(j-1)+i-1) + p(nx*(j-1)+i) )
        else if ( (i.eq.1) .and. (j.ne.1) .and. (j.ne.nx+1)) then
           pp = 0.5d+0 * (p(nx*(j-2)+i) + p(nx*(j-1)+i) )
        else if ( (i.eq.nx+1) .and. (j.ne.1) .and. (j.ne.nx+1)) then
           pp = 0.5d+0 * (p(nx*(j-2)+i-1) + p(nx*(j-1)+i-1) )
        else if ( (j.eq.1) .and. (i.ne.1) .and. (i.ne.nx+1)) then
           pp = 0.25d+0 * (p(nx*(j-1)+i-1) + p(nx*(j-1)+i) )
        else if ( (j.eq.nx+1) .and. (i.ne.1) .and. (i.ne.nx+1)) then
           pp = 0.25d+0 * (p(nx*(j-2)+i-1) + p(nx*(j-2)+i) )
        else if ( (j.eq.1) .and. (i.eq.1) ) then
           pp = p(1)
        else if ( (j.eq.nx+1) .and. (i.eq.1) ) then
           pp = p(nx*(nx-1)+1)
        else if ( (j.eq.1) .and. (i.eq.nx+1) ) then
           pp = p(nx)
        else if ( (j.eq.nx+1) .and. (i.eq.nx+1) ) then
           pp = p(nx*nx)
        end if
        write(iunit,*) pp, uu, vv
     end do
  end do
  close(iunit)
  
end subroutine dumpdata

