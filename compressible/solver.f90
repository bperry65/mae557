program solver

  use parameters
  
  implicit none
  integer :: j
  doubleprecision :: tdump, tprint
  
  ! read inputs and get ready

  call getarg(1,filename)
  call read_params
  call initialize

  tdump = dumpinterval
  tprint = tend / 1000d+0
  
  ! run simulation
  do while (t.lt.(tend - 1d-10))
     
     call timestep

     call mass_conservation
     
     if (t.ge.tprint) then
        !print *, 't', t
        tprint = tprint + tend / 1000d+0
        !print *, 'v', v(10,15)
        !print *, 'vbound', boundary_v(3), 'ybound', boundaryloc(3)
        print *, t, total_mass
        !print *, 'P  ', P(5,:)
        !print *, 'rho', rho(5,:)
        !print *, 'v  ', v(5,:)
        !print *, 'T  ', Temp(5,:)
     end if
     
     if (t.ge.tdump - 1d-10) then
        call dumpdata
        tdump = tdump + dumpinterval
     end if

  end do

  call dumpdata
  
endprogram solver


! ============================================================!
!                        Initialization                       !
! ============================================================!

subroutine initialize

  use parameters
  implicit none

  integer :: i, iunit=8

  allocate(u(nx+1,ny+1))
  allocate(v(nx+1,ny+1))
  allocate(Temp(nx+1,ny+1))
  allocate(P(nx+1,ny+1))
  allocate(rho(nx+1,ny+1))
  allocate(qy(nx+1,ny+1))
  allocate(qx(nx+1,ny+1))
  allocate(tau_xy(nx+1,ny+1))
  allocate(tau_yx(nx+1,ny+1))
  allocate(tau_yy(nx+1,ny+1))
  allocate(tau_xx(nx+1,ny+1))
  allocate(Et(nx+1,ny+1))
  allocate(Et_new(nx+1,ny+1))
  allocate(rho_v(nx+1,ny+1))
  allocate(rho_u(nx+1,ny+1))
  allocate(rho_new(nx+1,ny+1))

  allocate(coveredcells(nx+1,ny+1))
  allocate(freshlycleared(nx+1,ny+1))
  allocate(ghostcells(nx+1,ny+1))
  allocate(boundaryloc(nx+1))
  allocate(boundary_u(nx+1))
  allocate(boundary_v(nx+1))
  allocate(xx(nx+1))
  allocate(yy(ny+1))

  if (trim(bc).eq.'density') then
     allocate(boundary_rho(nx+1))
     boundary_rho(i) = 1d+0
  end if
  
  ! initialize and apply initial conditions
  t = 0d+0
  u = 0d+0
  v = 0d+0
  Temp = 1d+0
  P = 1d+0
  rho = 1d+0
  qx = 0d+0
  qy = 0d+0
  tau_xx = 0d+0
  tau_yy = 0d+0
  tau_xy = 0d+0
  tau_yx = 0d+0
  rho_u = 0d+0
  rho_v = 0d+0
  Et = 1d+0

  do i = 1,nx+1
     xx(i) = dble(i-1)/dble(nx)
  end do
  do i = 1,ny+1
     yy(i) = dble(i-1)/dble(nx)
  end do
  
  ! apply boundary condition on top wall
  call update_moving_boundary
  call update_mask
  freshlycleared = .false.

  call update_ghost

  ! intialize output file
  open(iunit,FILE=trim(outfile))
  write(iunit, *) nx
  write(iunit, *) ny
  close(iunit)
  
end subroutine initialize



! ============================================================!
!                        Run Timestep                         !
! ============================================================!

subroutine timestep

  use parameters
  implicit none

  integer :: i,j,v_up,u_up
  
  t = t + dt

  call calc_fluxes

  do i = 2,nx
     do j = 2,ny
        if (.not. coveredcells(i,j)) then
           rho_new(i,j) = rho(i,j) - dt/Omega * ( &
                0.5d+0 * (rho(i+1,j)*u(i+1,j) - rho(i-1,j)*u(i-1,j))/dx + &
                0.5d+0 * (rho(i,j+1)*v(i,j+1) - rho(i,j-1)*v(i,j-1))/dx )
           rho_u(i,j) = rho_u(i,j) - dt/Omega * ( &
                0.5d+0 / dx * &
                ( rho(i+1,j)*u(i+1,j)*u(i+1,j) - rho(i-1,j)*u(i-1,j)*u(i-1,j) &
                + 1d+0/gamma/Ma**2d+0 * (P(i+1,j) - P(i-1,j)) &
                - 2d+0 * (tau_xx(i,j) - tau_xx(i-1,j))  &
                + rho(i,j+1)*u(i,j+1)*v(i,j+1) - rho(i,j-1)*u(i,j-1)*v(i,j-1) &
                - 2d+0 * (tau_yx(i,j) - tau_yx(i,j-1)) ))
           rho_v(i,j) = rho_v(i,j) - dt/Omega * ( &
                0.5d+0 / dx * &
                ( rho(i+1,j)*u(i+1,j)*v(i+1,j) - rho(i-1,j)*u(i-1,j)*v(i-1,j) &
                - 2d+0 * (tau_xy(i,j) - tau_xy(i-1,j)) &
                + rho(i,j+1)*v(i,j+1)*v(i,j+1) - rho(i,j-1)*v(i,j-1)*v(i,j-1) &
                + 1d+0/gamma/Ma**2d+0 * (P(i,j+1) - P(i,j-1)) &
                - 2d+0 * (tau_yy(i,j) - tau_yy(i,j-1)) ))
           Et_new(i,j) = Et(i,j) - dt/Omega * ( &
                0.5d+0 / dx * &
                ( Et(i+1,j)*u(i+1,j) - Et(i-1,j)*u(i-1,j) &
                + (gamma-1d+0) * (P(i+1,j)*u(i+1,j) - P(i-1,j)*u(i-1,j)) &
                - Ma**2 * (gamma-1d+0) * gamma * ( (u(i+1,j)+u(i,j))*tau_xx(i,j) - (u(i-1,j)+u(i,j))*tau_xx(i-1,j) ) &
                - Ma**2 * (gamma-1d+0) * gamma * ( (v(i+1,j)+v(i,j))*tau_xy(i,j) - (v(i-1,j)+v(i,j))*tau_xy(i-1,j) ) &
                + 2d+0 * (qx(i,j) - qx(i-1,j)) &
                + Et(i,j+1)*v(i,j+1) - Et(i,j-1)*v(i,j-1) &
                + (gamma-1d+0) * (P(i,j+1)*v(i,j+1) - P(i,j-1)*v(i,j-1)) &
                - Ma**2 * (gamma-1d+0) * gamma * ( (u(i,j+1)+u(i,j))*tau_yx(i,j) - (u(i,j-1)+u(i,j))*tau_yx(i,j-1) ) &
                - Ma**2 * (gamma-1d+0) * gamma * ( (v(i,j+1)+v(i,j))*tau_yy(i,j) - (v(i,j-1)+v(i,j))*tau_yy(i,j-1) ) &
                + 2d+0 * (qy(i,j) - qy(i,j-1)) ))
        else
           rho_new(i,j) = 1d+0
           rho_u(i,j) = 0d+0
           rho_v(i,j) = 0d+0
           Et_new(i,j) = 0d+0
        end if
     end do
  end do

  call calc_primatives

  call update_moving_boundary
  call update_mask
  
  ! Apply boundary conditions
  select case (bc)
  case('density')
     call apply_P_bc_rho
  case('zpg')
     call apply_P_bc_zpg
  case('cpg')
     call apply_P_bc_cpg
  end select

  call update_ghost ! update ghost updates the velocity and pressure BCs at the moving boundary

  
end subroutine timestep



! ============================================================!
!                        Calculate Fluxes                         !
! ============================================================!

subroutine calc_fluxes

  use parameters
  implicit none

  integer :: i,j
  doubleprecision :: dudx, dudy, dvdx, dvdy

  ! Fluxes located midway between the points where other info is stored (staggered grid)
  ! tau_xx, tau_xy, qx located at x-direction midpoints
  ! tau_yy, tau_yx, qy located at y-direction midpoints
  ! All calculated using 2nd order central difference
  
  do i =1,nx
     do j=2,ny
        !if (.not. coveredcells(i,j)) then
           dudx = (u(i+1,j) - u(i,j))/dx
           dvdx = (v(i+1,j) - v(i,j))/dx
           if (i.eq.1) then
              dudy = 0.25d+0 *(u(i+1,j+1) - u(i+1,j-1))/dx
              dvdy = 0.25d+0 *(v(i+1,j+1) - v(i+1,j-1))/dx
           else if (i.eq.nx) then
              dudy = 0.25d+0 *(u(i,j+1) - u(i,j-1))/dx
              dvdy = 0.25d+0 *(v(i,j+1) - v(i,j-1))/dx
           else
              dudy = 0.25d+0 *(u(i+1,j+1) + u(i,j+1) - u(i+1,j-1) - u(i,j-1))/dx
              dvdy = 0.25d+0 *(v(i+1,j+1) + v(i,j+1) - v(i+1,j-1) - v(i,j-1))/dx
           end if
           qx(i,j) = -gamma/Re/Pr * (Temp(i+1,j) - Temp(i,j))/dx
           tau_xx(i,j) = 1d+0 / Re * 2d+0 / 3d+0 * (2d+0 * dudx - dvdy)
           tau_xy(i,j) = 1d+0 / Re *               (dudy + dvdx)
        !end if
     end do
  end do
  
  do i = 2,nx
     do j= 1,ny
        !if (.not. coveredcells(i,j)) then
           dudy = (u(i,j+1) - u(i,j))/dx
           dvdy = (v(i,j+1) - v(i,j))/dx
           if (j.eq.1) then
              dudx = 0.25d+0 * (u(i+1,j+1) - u(i-1,j+1))/dx 
              dvdx = 0.25d+0 * (v(i+1,j+1) - v(i-1,j+1))/dx 
           else
              dudx = 0.25d+0 * (u(i+1,j+1) + u(i+1,j) - u(i-1,j+1) - u(i-1,j))/dx 
              dvdx = 0.25d+0 * (v(i+1,j+1) + v(i+1,j) - v(i-1,j+1) - v(i-1,j))/dx 
           end if
           qy(i,j) = -gamma/Re/Pr * (Temp(i,j+1) - Temp(i,j))/dx
           tau_yy(i,j) = 1d+0 / Re * 2d+0 / 3d+0 * (2d+0 * dvdy - dudx)
           tau_yx(i,j) = 1d+0 / Re *               (dudy + dvdx)
        !end if
     end do
  end do
  
end subroutine calc_fluxes

! ============================================================= !
!                   calculate primitives                        !
! ============================================================= !

subroutine calc_primatives

  use parameters
  implicit none

  integer :: i,j

  do i=2,nx
     do j=2,ny
        rho(i,j) = rho_new(i,j)
        Et(i,j) = Et_new(i,j)
        u(i,j) = rho_u(i,j) / rho(i,j)
        v(i,j) = rho_v(i,j) / rho(i,j)
        P(i,j) = Et(i,j) - gamma*(gamma-1d+0)* Ma**2d+0 * rho(i,j)* 0.5d+0 * (u(i,j)**2d+0 + v(i,j)**2d+0)
        Temp(i,j) = P(i,j)/rho(i,j)
     end do
  end do
  
end subroutine calc_primatives


! ============================================================= !
!                           BCs                                 !
! ============================================================= !

! Pressure
! density determined from continuity, pressure from ideal gas law
subroutine apply_P_bc_rho

  use parameters
  implicit none
  integer :: i

  do i = 2,nx
     rho(i,1)    = rho(i,1)    - dt/Omega/dx*rho(i,2)*v(i,2)
     P(i,1) = rho(i,1)*Temp(i,1)
  end do

  do i = 2,ny
     P(1,i) = rho(1,i)*Temp(1,i)
     P(nx+1,i) = rho(nx+1,i)*Temp(nx+1,i)
     rho(1,i)    = rho(1,i)    - dt/Omega/dx*rho(2,i)*u(2,i)
     rho(nx+1,i) = rho(nx+1,i) + dt/Omega/dx*rho(nx,i)*u(nx,i)
  end do
  
end subroutine apply_P_bc_rho

! Zero Pressure Gradient
! apply to only moving boundaries
subroutine apply_P_bc_zpg

  use parameters
  implicit none
  integer :: i

  do i = 2,nx
     P(i,1)    = P(i,2)    
     rho(i,1) = P(i,1)/Temp(i,1)
  end do

  do i = 2,ny
     P(1,i)    = P(2,i)    
     P(nx+1,i) = P(nx,i)
     rho(1,i) = P(1,i)/Temp(1,i)
     rho(nx+1,i) = P(nx+1,i)/Temp(nx+1,i)
  end do
     
end subroutine apply_P_bc_zpg


! Constant Pressure Gradient
! (pressure determined by value at adjacent cell and gradient at adjacent cell * dx)
! apply to only moving boundaries
subroutine apply_P_bc_cpg

  use parameters
  implicit none
  integer :: i

  do i = 2,nx
     P(i,1)    = P(i,2)  * 2d+0 - P(i,3)  
     rho(i,1) = P(i,1)/Temp(i,1)
  end do

  do i = 2,ny
     rho(1,i) = P(1,i)/Temp(1,i)
     rho(nx+1,i) = P(nx+1,i)/Temp(nx+1,i)
     P(1,i)    = P(2,i)  * 2d+0 - P(3,i)  
     P(nx+1,i) = P(nx,i) * 2d+0 - P(nx-1,i)
  end do
  
end subroutine apply_P_bc_cpg

  
! ============================================================= !
!                           Data Dumping                        !
! ============================================================= !

subroutine dumpdata

  use parameters
  implicit none

  integer :: i,j,iunit=8

  open(iunit, file=trim(outfile),position='APPEND')
  write(iunit,*) 'time', t, 'nondimensional ', 'a ', 'b '
  write(iunit,*) 'P ','U ', 'V ', 'T ', 'rho'
  do i=1,nx+1
     do j = 1,ny+1
        write(iunit,"(E20.10$)") P(i,j), U(i,j), V(i,j), Temp(i,j), rho(i,j), xx(i), yy(j)
        write(iunit,*) ''
     end do
  end do
  close(iunit)
  
end subroutine dumpdata

! ============================================================= !
!                 Immersed Boundary - masked cells              !
! ============================================================= !

subroutine update_mask

  use parameters
  implicit none

  integer :: i,j
  logical :: notcovered

  ghostcells = .false.
  
  do i = 1,nx+1
     do j = 1,ny+1
        if (boundaryloc(i) .le. yy(j)) then
           coveredcells(i,j) = .true.
           freshlycleared(i,j) = .false.
        else
           if (coveredcells(i,j) .eq. .true.) then
              coveredcells(i,j) = .false.
              freshlycleared(i,j) = .true.
           else
              freshlycleared(i,j) = .false.
              coveredcells(i,j) = .false.
           end if
        end if
     end do
  end do

  do i = 1,nx+1
     j = 0
     notcovered = .true.
     do while (notcovered)
        j = j + 1
        if (coveredcells(i,j) .eq. .true.) then
           notcovered = .false.
           ghostcells(i,j) = .true.
        end if
     end do
  end do
  
  
end subroutine update_mask

! ============================================================= !
!                 Immersed Boundary - boundary location          !
! ============================================================= !

subroutine update_moving_boundary

  use parameters
  implicit none

  integer :: i

  select case (pistvel)
  case('constantplus')
     do i = 1,nx+1
        boundaryloc(i) = 1d+0 + F*t
        boundary_u(i) = 0d+0
        boundary_v(i) = 1d+0
     end do
  case('constantminus')
     do i = 1,nx+1
        boundaryloc(i) = 1d+0 - F*t
        boundary_u(i) = 0d+0
        boundary_v(i) = -1d+0
     end do     
  case('constant')
     do i = 1, nx+1
        if (mod(t+pi/2d+0,2d+0*pi) .lt. pi) then
           boundaryloc(i) = 1d+0 - pi/2d+0*F + F * (mod(t+pi/2d+0,2d+0*pi))
           boundary_u(i) = 0d+0
           boundary_v(i) = 1d+0
        else
           boundaryloc(i) = 1d+0 + pi/2d+0*F - F * (mod(t+pi/2d+0,2d+0*pi) - pi)
           boundary_u(i) = 0d+0
           boundary_v(i) = -1d+0
        end if
     end do
  case('sinusoid')
     do i = 1, nx+1
        boundaryloc(i) = 1d+0 + F*(1d+0 - cos(t))
        boundary_u(i) = 0d+0
        boundary_v(i) = +sin(t)
     end do
  end select

end subroutine update_moving_boundary


! ============================================================= !
!                 Immersed Boundary - update ghost cells         !
! ============================================================= !

subroutine update_ghost

  ! this subroutine assumes a flat piston head, moving in the direction normal to its surface
  
  use parameters
  implicit none

  integer :: i,j, nfreshcleared
  double precision :: a_int ! coefficient for interpolation

  ! Where are the ghost cells?
  ghostrow = 0
  do j = 1, ny+1
     if (ghostrow.eq.0 .and. ghostcells(2,j).eq.(.true.)) ghostrow = j
  end do
  
  ! Are there freshly cleared cells?
  nfreshcleared = 0
  j = 0
  do while (freshlycleared(2,ghostrow-j-1) .eq. .true.)
     j = j+1
     nfreshcleared = nfreshcleared + 1
  end do
  a_int = (boundaryloc(2) - yy(ghostrow-1-nfreshcleared)) / (yy(ghostrow) - yy(ghostrow-1-nfreshcleared))
  
  ! determine value at ghost cell that results in interpolation at boundary location to meet the boundary condition
  do i = 2,nx
     ! Interpolation
     u(i,ghostrow) = ( boundary_u(i) - u(i,ghostrow-1-nfreshcleared)*(1d+0 - a_int) )/a_int
     v(i,ghostrow) = ( boundary_v(i) - v(i,ghostrow-1-nfreshcleared)*(1d+0 - a_int) )/a_int
     Temp(i,ghostrow) = ( 1d+0          - Temp(i,ghostrow-1-nfreshcleared)*(1d+0 - a_int) )/a_int

   !  if (a_int .le. 1d-4) then
   !     u(i,ghostrow) = boundary_u(i)
   !     v(i,ghostrow) = boundary_v(i)
   !     Temp(i,ghostrow) = 1d+0
   !  end if

     ! pressure boundary condition
     select case(bc)
     case('density') ! calculate change in density following boundary, extrapolate to ghost cells
        boundary_rho(i) = boundary_rho(i) - dt*boundary_rho(i)*(v(i,ghostrow)-v(i,ghostrow-1-nfreshcleared))/(dx+nfreshcleared)
        rho(i,ghostrow) = ( boundary_rho(i) - rho(i,ghostrow-1-nfreshcleared)*(1d+0 - a_int) )/a_int
        P(i,ghostrow) = rho(i,ghostrow)*Temp(i,ghostrow)
     case('cpg')
        P(i,ghostrow) = P(i,ghostrow-1-nfreshcleared) &
             + dble(nfreshcleared+1) * (P(i,ghostrow-1-nfreshcleared) - P(i,ghostrow-2-nfreshcleared))
        rho(i,ghostrow) = P(i,ghostrow)/Temp(i,ghostrow)
     case('zpg')
        P(i,ghostrow) = P(i,ghostrow-1-nfreshcleared)
        rho(i,ghostrow) = P(i,ghostrow)/Temp(i,ghostrow)
     end select

     ! calculation of other variables

     rho_u(i,ghostrow) = rho(i,ghostrow)*u(i,ghostrow) 
     rho_v(i,ghostrow) = rho(i,ghostrow)*v(i,ghostrow)
     Et(i,ghostrow) = P(i,ghostrow) + rho(i,ghostrow) * gamma*(gamma-1d+0) * Ma**2d+0 * &
          0.5d+0 * (u(i,ghostrow)**2d+0 + v(i,ghostrow)**2d+0)
  end do

  ! Deal with the freshly cleared cells: interpolate from boundary and adjacent cells
  if (nfreshcleared .ne. 0) then
     do i = 2,nx
        do j = ghostrow-1, ghostrow-nfreshcleared
           ! Interpolation
           u(i,j) = u(i,ghostrow-1-nfreshcleared) &
                + dble(j-ghostrow+nfreshcleared+1)/(dble(nfreshcleared+1)*a_int) &
                * (boundary_u(i) - u(i,ghostrow-1-nfreshcleared))
           v(i,j) = v(i,ghostrow-1-nfreshcleared) &
                + dble(j-ghostrow+nfreshcleared+1)/(dble(nfreshcleared+1)*a_int) &
                * (boundary_v(i) - v(i,ghostrow-1-nfreshcleared))
           Temp(i,j) = Temp(i,ghostrow-1-nfreshcleared) &
                + dble(j-ghostrow+nfreshcleared+1)/(dble(nfreshcleared+1)*a_int) &
                * (1d+0          - Temp(i,ghostrow-1-nfreshcleared))

           ! Pressure boundary condition
           select case(bc)
           case('density')
              P(i,j) = P(i,ghostrow-1-nfreshcleared)
           case('cpg')
              P(i,j) = P(i,ghostrow-1-nfreshcleared) &
                   + dble(j-ghostrow+nfreshcleared+1)*(P(i,ghostrow-1-nfreshcleared) - P(i,ghostrow-2-nfreshcleared))
           case('zpg')
              P(i,j) = P(i,ghostrow-1-nfreshcleared)
           end select
        
           !calculation of other variables
           rho(i,j) = P(i,j)/Temp(i,j)
           rho_u(i,j) = rho(i,j)*u(i,j) 
           rho_v(i,j) = rho(i,j)*v(i,j)
           Et(i,j) = P(i,j) + rho(i,j) * gamma*(gamma-1d+0) * Ma**2d+0 * &
                0.5d+0 * (u(i,j)**2d+0 + v(i,j)**2d+0)
        end do
     end do
  end if
  
end subroutine update_ghost



! ============================================================= !
!                     Check mass conservation                   !
! ============================================================= !

subroutine mass_conservation

  use parameters
  implicit none

  integer :: i,j

  total_mass = 0d+0
  
  do i = 2,nx
     do j=2,ny
        if (.not. coveredcells(i,j)) then
           total_mass = total_mass + rho(i,j)*dx*dx
        end if
     end do
     total_mass = total_mass + &
          (rho(i,ghostrow-1) + (boundaryloc(i) - yy(ghostrow-1))/dx * (rho(i,ghostrow)-rho(i,ghostrow))) &
          * dx *(boundaryloc(i) - yy(ghostrow-1))  
     total_mass = total_mass + rho(i,1)*dx*dx/2
  end do

  do j = 2,ny
     total_mass = total_mass + (rho(1,j) + rho(nx+1,j))*dx*dx/2
  end do

  total_mass = total_mass + (rho(1,1) + rho(nx+1,1))*dx*dx/4
  
end subroutine mass_conservation
