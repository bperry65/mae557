program solver

  use parameters
  
  implicit none
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
     
     if (t.ge.tprint) then
        print *, 't', t
        tprint = tprint + tend / 1000d+0
        print *, 'u', u(10,nx-1)
       ! print *, 'u-top',rho(:,nx+1)
       ! print *, 'u-bottom',rho(:,1)
       ! print *, 'u-left',rho(1,:)
       ! print *, 'u-right',rho(nx+1,:)
     end if
     
     if (t.ge.tdump - 1d-10) then
        call dumpdata
        tdump = tdump + dumpinterval
     end if
     
  end do

  call dumpdata
  !print *, 'u', u
  !print *
  ! print *, 'v', v
  ! print *
  !print *, 'P', P
  !print *
  ! print *, 'T', Temp
  ! print *
  ! print *, 'rho', rho
  
endprogram solver


! ============================================================!
!                        Initialization                       !
! ============================================================!

subroutine initialize

  use parameters
  implicit none

  integer :: i, iunit=8

  allocate(u(nx+1,nx+1))
  allocate(v(nx+1,nx+1))
  allocate(Temp(nx+1,nx+1))
  allocate(P(nx+1,nx+1))
  allocate(rho(nx+1,nx+1))
  allocate(qy(nx+1,nx+1))
  allocate(qx(nx+1,nx+1))
  allocate(tau_xy(nx+1,nx+1))
  allocate(tau_yx(nx+1,nx+1))
  allocate(tau_yy(nx+1,nx+1))
  allocate(tau_xx(nx+1,nx+1))
  allocate(Et(nx+1,nx+1))
  allocate(rho_v(nx+1,nx+1))
  allocate(rho_u(nx+1,nx+1))
  allocate(rho_new(nx+1,nx+1))
  
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
  rho_u = 0d+0
  rho_v = 0d+0
  Et = 1d+0
  
  ! apply boundary condition on top wall
  call apply_vel_bc
  
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

  call calc_fluxes
  
  ! Calculate rho, u, v, T at next time step using Forward Euler
  do i = 2,nx
     do j = 2,nx
        rho_new(i,j) = rho(i,j) - dt/Omega * ( &
             0.5d+0 * (rho_u(i+1,j) - rho_u(i-1,j))/dx + &
             0.5d+0 * (rho_v(i,j+1) - rho_v(i,j-1))/dx )
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
        Et(i,j) = Et(i,j) - dt/Omega * ( &
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
     end do
  end do

  call calc_primatives
  
  ! Apply boundary conditions
  select case (bc)
  case('nscbc')
     call apply_P_bc_nscbc
  case('rho')
     call apply_P_bc
  case('zpg')
     call apply_P_bc_zpg
  case('constP')
     if (t.le.1e-12) print *, 'using const P bc'
  end select
  
  call apply_vel_bc
  
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
     do j=2,nx
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
        tau_xy(i,j) = 1d+0 / Re * 2d+0 / 3d+0 * (dudy + dvdx)
     end do
  end do
  
  do i = 2,nx
     do j= 1,nx
        dudy = (u(i,j+1) - u(i,j))/dx
        dvdy = (v(i,j+1) - v(i,j))/dx
        if (j.eq.1) then
           dudx = 0.25d+0 * (u(i+1,j+1) - u(i-1,j+1))/dx 
           dvdx = 0.25d+0 * (v(i+1,j+1) - v(i-1,j+1))/dx 
        else if (j.eq.nx) then
           dudx = 0.25d+0 * (u(i+1,j) - u(i-1,j))/dx 
           dvdx = 0.25d+0 * (v(i+1,j) - v(i-1,j))/dx 
        else
           dudx = 0.25d+0 * (u(i+1,j+1) + u(i+1,j) - u(i-1,j+1) - u(i-1,j))/dx 
           dvdx = 0.25d+0 * (v(i+1,j+1) + v(i+1,j) - v(i-1,j+1) - v(i-1,j))/dx 
        end if
        
        qy(i,j) = -gamma/Re/Pr * (Temp(i,j+1) - Temp(i,j))/dx
        tau_yy(i,j) = 1d+0 / Re * 2d+0 / 3d+0 * (2d+0 * dvdy - dudx)
        tau_yx(i,j) = 1d+0 / Re * 2d+0 / 3d+0 * (dudy + dvdx)        
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
     do j=2,nx
        rho(i,j) = rho_new(i,j)
        u(i,j) = rho_u(i,j) / rho(i,j)
        v(i,j) = rho_v(i,j) / rho(i,j)
        P(i,j) = Et(i,j) - gamma*(gamma-1d+0)*rho(i,j)* 0.5d+0 * (u(i,j)**2d+0 + v(i,j)**2d+0)
        Temp(i,j) = P(i,j)/rho(i,j)
     end do
  end do
  
end subroutine calc_primatives


! ============================================================= !
!                           velocity BCs                        !
! ============================================================= !

subroutine apply_vel_bc

  use parameters
  implicit none
  integer :: i
  
  do i = 2,nx
     u(i,nx+1) = 1.0d+0 *sin(t)
     rho_u(i,nx+1) = u(i,nx+1)*rho(i,nx+1)
     Et(i,nx+1) = P(i,nx+1) + rho(i,nx+1) * gamma*(gamma-1d+0) * Ma**2d+0 * 0.5d+0 * u(i,nx+1)**2d+0
  end do
  
end subroutine apply_vel_bc

! Pressure
subroutine apply_P_bc

  use parameters
  implicit none
  integer :: i

  ! density determined from continuity, pressure from ideal gas law
  do i = 2,nx
     rho(i,nx+1) = rho(i,nx+1) + dt/Omega/dx*rho(i,nx)*v(i,nx)
     rho(i,1)    = rho(i,1)    - dt/Omega/dx*rho(i,2)*v(i,2)
     rho(1,i)    = rho(1,i)    - dt/Omega/dx*rho(2,i)*u(2,i)
     rho(nx+1,i) = rho(nx+1,i) + dt/Omega/dx*rho(nx,i)*u(nx,i)
     P(i,nx+1) = rho(i,nx+1)*Temp(i,nx+1)
     P(i,1) = rho(i,1)*Temp(i,1)
     P(1,i) = rho(1,i)*Temp(1,i)
     P(nx+1,i) = rho(nx+1,i)*Temp(nx+1,i)
  end do
  
end subroutine apply_P_bc


! Zero Pressure Gradient
subroutine apply_P_bc_zpg

  use parameters
  implicit none
  integer :: i

  ! density determined from continuity, pressure from ideal gas law
  do i = 2,nx
     P(i,nx+1) = rho(i,nx) 
     P(i,1)    = rho(i,2)    
     P(1,i)    = rho(2,i)    
     P(nx+1,i) = rho(nx,i) 
     rho(i,nx+1) = P(i,nx+1)/Temp(i,nx+1)
     rho(i,1) = P(i,1)/Temp(i,1)
     rho(1,i) = P(1,i)/Temp(1,i)
     rho(nx+1,i) = P(nx+1,i)/Temp(nx+1,i)
  end do
  
end subroutine apply_P_bc_zpg


! Pressure NSCBC
subroutine apply_P_bc_nscbc

  use parameters
  implicit none
  integer :: i
  double precision :: L5,L1
  
  do i = 2,nx
     ! Top Boundary
     L5 = (v(i,nx) + Temp(i,nx)**0.5d+0 / Ma) * &
          ( (P(i,nx) - P(i,nx-1))/dx &
          + Ma * gamma * rho(i,nx) * Temp(i,nx)**0.5d+0 * (v(i,nx)-v(i,nx-1))/dx)
     P(i,nx+1) = P(i,nx) &
          + L5 * dx / (v(i,nx)-Temp(i,nx)**0.5d+0 / Ma) &
          + Ma * gamma * rho(i,nx) * Temp(i,nx)**0.5d+0 * v(i,nx)
     ! Right Boundary
     L5 = (u(nx,i) + Temp(nx,i)**0.5d+0 / Ma) * &
          ( (P(nx,i) - P(nx-1,i))/dx &
          + Ma * gamma * rho(nx,i) * Temp(nx,i)**0.5d+0 * (u(nx,i)-u(nx-1,i))/dx)
     P(nx+1,i) = P(nx,i) &
          + L5 * dx / (u(nx,i)-Temp(nx,i)**0.5d+0 / Ma) &
          + Ma * gamma * rho(nx,i) * Temp(nx,i)**0.5d+0 * u(nx,i)
     ! Bottom Boundary
     L1 = (v(i,2) - Temp(i,2)**0.5d+0 / Ma) * &
          ( (P(i,3) - P(i,2))/dx &
          - Ma * gamma * rho(i,2) * Temp(i,2)**0.5d+0 * (v(i,3)-v(i,2))/dx)
     P(i,1) = P(i,2) &
          - L1 * dx / (v(i,2) + Temp(i,2)**0.5d+0 / Ma) &
          + Ma * gamma * rho(i,2) * Temp(i,2)**0.5d+0 * v(i,2)
     ! Left Boundary
     L1 = (u(2,i) - Temp(2,i)**0.5d+0 / Ma) * &
          ( (P(3,i) - P(2,i))/dx &
          - Ma * gamma * rho(2,i) * Temp(2,i)**0.5d+0 * (u(3,i)-u(2,i))/dx)
     P(1,i) = P(2,i) &
          - L1 * dx / (u(2,i) + Temp(2,i)**0.5d+0 / Ma) &
          + Ma * gamma * rho(2,i) * Temp(2,i)**0.5d+0 * u(2,i)
     ! Update Densities
     rho(i,nx+1) = P(i,nx+1) / Temp(i,nx+1)
     rho(i,1) = P(i,1) / Temp(i,1)
     rho(nx+1,i) = P(nx+1,i) / Temp(nx+1,i)
     rho(1,i) = P(1,i) / Temp(1,i)
  end do

end subroutine apply_P_bc_nscbc
  
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
     do j = 1,nx+1
        write(iunit,"(E20.10$)") P(i,j), U(i,j), V(i,j), Temp(i,j), rho(i,j)
        write(iunit,*) ''
     end do
  end do
  close(iunit)
  
end subroutine dumpdata

