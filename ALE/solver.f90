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

  allocate(x(nx+1))
  allocate(y(nx+1))
  allocate(v_hat(nx+1))
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
  allocate(Et_new(nx+1,nx+1))
  
  ! initialize and apply initial conditions
  t = 0d+0
  u = 0d+0
  v = 0d+0
  Temp = 300d+0
  P = 100000d+0
  rho = 100000d+0/(R*300d+0)
  qx = 0d+0
  qy = 0d+0
  tau_xx = 0d+0
  tau_yy = 0d+0
  tau_xy = 0d+0
  tau_yx = 0d+0
  rho_u = 0d+0
  rho_v = 0d+0
  Et = 1/(gamma-1d+0)*100000d+0
  do i = 1,nx+1 
     x(i) = 1.0d+0*(i-1)/dble(nx)
     y(i) = 1.0d+0*(i-1)/dble(nx)*L1
  end do
  
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

  integer :: i,j,v_up,u_up
  
  t = t + dt

  call calc_fluxes

  call grid_adv

  ! First Order Upwinding
  if( use_upwind ) then
     if (t.lt.1d-12+dt) print *, 'Using upwind'
     if (u(i,j).lt.0d+0) then
        u_up = 1
     else
        u_up = 0
     end if
     if (v(i,j).lt.0d+0) then
        v_up = 1
     else
        v_up = 0
     end if
     
     do i = 2,nx
        do j = 2,nx     
           rho_new(i,j) = rho(i,j) - dt * ( &
                (u(i,j)) * (rho(i+u_up,j) - rho(i-1+u_up,j))/dx + &
                rho(i,j) * (u(i+u_up,j) - u(i-1+u_up,j))/dx + &
                (v(i,j)-v_hat(j)) * (rho(i,j+v_up) - rho(i,j-1+v_up))/dy + &
                rho(i,j) * (v(i,j+v_up) - v(i,j-1+v_up))/dy )
           rho_u(i,j) = rho_u(i,j) - dt * &
                (1d+0/dx *( u(i,j)*(rho(i+u_up,j)*u(i+u_up,j) - rho(i-1+u_up,j)*u(i-1+u_up,j)) &
                + rho(i,j)*u(i,j)*(u(i+u_up,j) - u(i-1+u_up,j)) &
                + 1d+0 * (P(i+u_up,j) - P(i-1+u_up,j))  &
                - 1d+0 * (tau_xx(i,j) - tau_xx(i-1,j)))  &
                + 1d+0/dy *( (v(i,j)-v_hat(j))*(rho(i,j+v_up)*u(i,j+v_up) - rho(i,j-1+v_up)*u(i,j-1+v_up)) &
                + rho(i,j)*u(i,j)*(v(i,j+v_up) - v(i,j-1+v_up)) &
                - 1d+0 * (tau_yx(i,j) - tau_yx(i,j-1)) ))
           rho_v(i,j) = rho_v(i,j) - dt * &
                (1d+0/dx *( u(i,j)*(rho(i+u_up,j)*v(i+u_up,j) - rho(i-1+u_up,j)*v(i-1+u_up,j)) &
                + rho(i,j)*v(i,j) * (u(i+u_up,j) - u(i-1+u_up,j)) &
                - 1d+0 * (tau_xy(i,j) - tau_xy(i-1,j))) &
                + 1d+0/dy *( (v(i,j)-v_hat(j))*(rho(i,j+v_up)*v(i,j+v_up) - rho(i,j-1+v_up)*v(i,j-1+v_up)) &
                + rho(i,j)*v(i,j) * (v(i,j+v_up) - v(i,j-1+v_up)) &
                + 1d+0 * (P(i,j+v_up) - P(i,j-1+v_up)) &
                - 1d+0 * (tau_yy(i,j) - tau_yy(i,j-1)) ))
           Et_new(i,j) = Et(i,j) - dt * &
                (1d+0/dx * (u(i,j) * ( Et(i+u_up,j) - Et(i-1+u_up,j)) &
                +  Et(i,j) * (u(i+u_up,j) - u(i-1+u_up,j)) &
                + 1d+0 * (P(i+u_up,j)*u(i+u_up,j) - P(i-1+u_up,j)*u(i-1+u_up,j)) &
                - 0.5d+0 * ( (u(i+1,j)+u(i,j))*tau_xx(i,j) - (u(i-1,j)+u(i,j))*tau_xx(i-1,j) ) &
                - 0.5d+0 * ( (v(i+1,j)+v(i,j))*tau_xy(i,j) - (v(i-1,j)+v(i,j))*tau_xy(i-1,j) ) &
                + 1d+0 * (qx(i,j) - qx(i-1,j))) &
                + 1d+0/dy * ((v(i,j)-v_hat(j)) * (Et(i,j+v_up) - Et(i,j-1+v_up)) &
                + Et(i,j) * (v(i,j+v_up) - v(i,j-1+v_up)) &
                + 1d+0 * (P(i,j+v_up)*v(i,j+v_up) - P(i,j-1+v_up)*v(i,j-1+v_up)) &
                - 0.5d+0 * ( (u(i,j+1)+u(i,j))*tau_yx(i,j) - (u(i,j-1)+u(i,j))*tau_yx(i,j-1) ) &
                - 0.5d+0 * ( (v(i,j+1)+v(i,j))*tau_yy(i,j) - (v(i,j-1)+v(i,j))*tau_yy(i,j-1) ) &
                + 1d+0 * (qy(i,j) - qy(i,j-1))))
        end do
     end do

  ! Central Difference
  else
     do i = 2,nx
        do j = 2,nx     
           rho_new(i,j) = rho(i,j) - dt * ( &
                0.5d+0 * (u(i,j) * (rho(i+1,j) - rho(i-1,j)) &
                + rho(i,j) * (u(i+1,j) - u(i-1,j)))/dx &
                + 0.5d+0 * ((v(i,j)-v_hat(j))*(rho(i,j+1) - rho(i,j-1)) &
                + rho(i,j)*(v(i,j+1) - v(i,j-1)))/dy )
           rho_u(i,j) = rho_u(i,j) - dt * ( &
                0.5d+0 / dx * &
                ( rho(i,j)*u(i,j)*(u(i+1,j) - u(i-1,j)) &
                + u(i,j)*(rho(i+1,j)*u(i+1,j) - rho(i-1,j)*u(i-1,j)) &
                + 1d+0 * (P(i+1,j) - P(i-1,j)) &
                - 2d+0 * (tau_xx(i,j) - tau_xx(i-1,j)))  &
                + 0.5d+0 / dy *(rho(i,j)*u(i,j)*(v(i,j+1) - v(i,j-1)) &
                + (v(i,j+1)-v_hat(j))*(rho(i,j+1)*u(i,j+1) - rho(i,j-1)*u(i,j-1)) &
                - 2d+0 * (tau_yx(i,j) - tau_yx(i,j-1)) ))
           rho_v(i,j) = rho_v(i,j) - dt * ( &
                0.5d+0 / dx * &
                (u(i,j)*( rho(i+1,j)*v(i+1,j) - rho(i-1,j)*v(i-1,j)) &
                + rho(i,j)*v(i,j)*(u(i+1,j) - u(i-1,j)) &
                - 2d+0 * (tau_xy(i,j) - tau_xy(i-1,j))) &
                + 0.5d+0 / dy *(rho(i,j)*v(i,j)*(v(i,j+1) - v(i,j-1)) &
                + (v(i,j)-v_hat(j))*(rho(i,j+1)*v(i,j+1) - rho(i,j-1)*v(i,j-1)) &
                + 1d+0 * (P(i,j+1) - P(i,j-1)) &
                - 2d+0 * (tau_yy(i,j) - tau_yy(i,j-1)) ))
           Et_new(i,j) = Et(i,j) - dt * ( &
                0.5d+0 / dx * &
                (u(i,j) * ( Et(i+1,j) - Et(i-1,j)) & !
                + Et(i,j) * (u(i+1,j) - u(i-1,j)) & !
                + (P(i+1,j)*u(i+1,j) - P(i-1,j)*u(i-1,j)) & !
                - ( (u(i+1,j)+u(i,j))*tau_xx(i,j) - (u(i-1,j)+u(i,j))*tau_xx(i-1,j) ) &
                - ( (v(i+1,j)+v(i,j))*tau_xy(i,j) - (v(i-1,j)+v(i,j))*tau_xy(i-1,j) ) &
                + 2d+0 * (qx(i,j) - qx(i-1,j))) &
                + 0.5d+0 / dy *((v(i,j)-v_hat(j)) * (Et(i,j+1) - Et(i,j-1)) & !
                + Et(i,j) * (v(i,j+1) - v(i,j-1)) & !
                + (P(i,j+1)*v(i,j+1) - P(i,j-1)*v(i,j-1)) & !
                - ( (u(i,j+1)+u(i,j))*tau_yx(i,j) - (u(i,j-1)+u(i,j))*tau_yx(i,j-1) ) &
                - ( (v(i,j+1)+v(i,j))*tau_yy(i,j) - (v(i,j-1)+v(i,j))*tau_yy(i,j-1) ) &
                + 2d+0 * (qy(i,j) - qy(i,j-1)) ))
        end do
     end do
  end if
  
  call calc_primatives
  
  ! Apply boundary conditions
  select case (bc)
!!$  case('nscbc')
!!$     call apply_P_bc_nscbc
  case('rho')
     call apply_P_bc
  case('zpg')
     call apply_P_bc_zpg
  case('cpg')
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
           dudy = 0.25d+0 *(u(i+1,j+1) - u(i+1,j-1))/dy
           dvdy = 0.25d+0 *(v(i+1,j+1) - v(i+1,j-1))/dy
        else if (i.eq.nx) then
           dudy = 0.25d+0 *(u(i,j+1) - u(i,j-1))/dy
           dvdy = 0.25d+0 *(v(i,j+1) - v(i,j-1))/dy
        else
           dudy = 0.25d+0 *(u(i+1,j+1) + u(i,j+1) - u(i+1,j-1) - u(i,j-1))/dy
           dvdy = 0.25d+0 *(v(i+1,j+1) + v(i,j+1) - v(i+1,j-1) - v(i,j-1))/dy
        end if
        qx(i,j) = lambda * (Temp(i+1,j) - Temp(i,j))/dx
        tau_xx(i,j) = 2d+0 / 3d+0 * visc * (2d+0 * dudx - dvdy)
        tau_xy(i,j) =               visc * (dudy + dvdx)
     end do
  end do
  
  do i = 2,nx
     do j= 1,nx
        dudy = (u(i,j+1) - u(i,j))/dy
        dvdy = (v(i,j+1) - v(i,j))/dy
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
        
        qy(i,j) = lambda * (Temp(i,j+1) - Temp(i,j))/dy
        tau_yy(i,j) = 2d+0 / 3d+0 * visc * (2d+0 * dvdy - dudx)
        tau_yx(i,j) =               visc * (dudy + dvdx)        
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
        Et(i,j) = Et_new(i,j)
        u(i,j) = rho_u(i,j) / rho(i,j)
        v(i,j) = rho_v(i,j) / rho(i,j)
        P(i,j) = (gamma-1d+0)* Et(i,j) + rho(i,j)* 0.5d+0 * (u(i,j)**2d+0 + v(i,j)**2d+0)
        Temp(i,j) = P(i,j)/(rho(i,j)*R)
     end do
  end do

  ! Fold in corner calculations
  rho(1,1) = rho(2,1)
  rho(1,nx+1) = rho(2,nx+1)
  rho(nx+1,1) = rho(nx,1)
  rho(nx+1,nx+1) = rho(nx,nx+1)

  Et(1,1) = Et(2,1)
  Et(1,nx+1) = Et(2,nx+1)
  Et(nx+1,1) = Et(nx,1)
  Et(nx+1,nx+1) = Et(nx,nx+1)

  u(1,1) = u(2,1)
  u(1,nx+1) = u(2,nx+1)
  u(nx+1,1) = u(nx,1)
  u(nx+1,nx+1) = u(nx,nx+1)

  v(1,1) = v(2,1)
  v(1,nx+1) = v(2,nx+1)
  v(nx+1,1) = v(nx,1)
  v(nx+1,nx+1) = v(nx,nx+1)

  P(1,1) = P(2,1)
  P(1,nx+1) = P(2,nx+1)
  P(nx+1,1) = P(nx,1)
  P(nx+1,nx+1) = P(nx,nx+1)

  Temp(1,1) = P(2,1)/(rho(2,1)*R)
  Temp(1,nx+1) = P(2,nx+1)/(rho(2,nx+1)*R)
  Temp(nx+1,1) = P(nx,1)/(rho(nx,1)*R)
  Temp(nx+1,nx+1) = P(nx,nx+1)/(rho(nx,nx+1)*R)

  
end subroutine calc_primatives


! ============================================================= !
!                           velocity BCs                        !
! ============================================================= !

subroutine apply_vel_bc

  use parameters
  implicit none
  integer :: i

  double precision :: v_w
  select case (pistvel)
  ! Constant advancement and withdrawal wall BCs
  case('constant')
     v_w = 2d+0*(L2-L1)*Omega !2d+0*Omega*(L1-L2)*cos(2d+0*pi*Omega*(t+1d+0/4d+0))/sqrt(1d+0-sin(2d+0*pi*Omega*(t+1/4d+0))**2d+0)
     do i = 2,nx
        v(i,nx+1) = v_w
        rho_v(i,nx+1) = v(i,nx+1)*rho(i,nx+1)
        Et(i,nx+1) = P(i,nx+1)/(gamma-1d+0) + rho(i,nx+1) * 0.5d+0 * v(i,nx+1)**2d+0
     end do
     ! Grid BC
     do i = 2,nx
        v_hat(i) = v_w*(i-1)/dble(nx)
     end do
     

  case('sinusoid')
     ! Full Cycle BC
     v_w = pi*Omega*(L2-L1)*sin(2d+0*pi*Omega*t)
     do i = 2,nx
        v(i,nx+1) = v_w
        rho_v(i,nx+1) = v(i,nx+1)*rho(i,nx+1)
        Et(i,nx+1) = P(i,nx+1)/(gamma-1d+0) + rho(i,nx+1) * 0.5d+0 * v(i,nx+1)**2d+0
     end do
     do i = 2,nx
        v_hat(i) = v_w*(i-1d+0)/dble(nx)
     end do

  end select

end subroutine apply_vel_bc

! Pressure
subroutine apply_P_bc

  use parameters
  implicit none
  integer :: i

  ! density determined from continuity, pressure from ideal gas law
  do i = 2,nx
     rho(i,nx+1) = rho(i,nx+1) + dt/dy*rho(i,nx)*v(i,nx) !/Omega
     rho(i,1)    = rho(i,1)    - dt/dy*rho(i,2)*v(i,2) !/Omega
     rho(1,i)    = rho(1,i)    - dt/dx*rho(2,i)*u(2,i) !/Omega
     rho(nx+1,i) = rho(nx+1,i) + dt/dx*rho(nx,i)*u(nx,i) !/Omega
     P(i,nx+1) = rho(i,nx+1)*R*Temp(i,nx+1)
     P(i,1) = rho(i,1)*R*Temp(i,1)
     P(1,i) = rho(1,i)*R*Temp(1,i)
     P(nx+1,i) = rho(nx+1,i)*R*Temp(nx+1,i)
  end do
  
end subroutine apply_P_bc


! Zero Pressure Gradient
subroutine apply_P_bc_zpg

  use parameters
  implicit none
  integer :: i

  ! density determined from continuity, pressure from ideal gas law
  do i = 2,nx
     P(i,nx+1) = P(i,nx)
     P(i,1)    = P(i,2)    
     P(1,i)    = P(2,i)    
     P(nx+1,i) = P(nx,i) 
     rho(i,nx+1) = P(i,nx+1)/(R*Temp(i,nx+1))
     rho(i,1) = P(i,1)/(R*Temp(i,1))
     rho(1,i) = P(1,i)/(R*Temp(1,i))
     rho(nx+1,i) = P(nx+1,i)/(R*Temp(nx+1,i))
  end do
  
end subroutine apply_P_bc_zpg

! Constant Pressure Gradient
! (pressure determined by value at adjacent cell and gradient at adjacent cell * dx)
subroutine apply_P_bc_cpg

  use parameters
  implicit none
  integer :: i

  ! density determined from continuity, pressure from ideal gas law
  do i = 2,nx
     P(i,nx+1) = P(i,nx) * 2d+0 - P(i,nx-1)
     P(i,1)    = P(i,2)  * 2d+0 - P(i,3)  
     P(1,i)    = P(2,i)  * 2d+0 - P(3,i)  
     P(nx+1,i) = P(nx,i) * 2d+0 - P(nx-1,i)
     rho(i,nx+1) = P(i,nx+1)/(R*Temp(i,nx+1))
     rho(i,1) = P(i,1)/(R*Temp(i,1))
     rho(1,i) = P(1,i)/(R*Temp(1,i))
     rho(nx+1,i) = P(nx+1,i)/(R*Temp(nx+1,i))
  end do
  
end subroutine apply_P_bc_cpg


!!$! Pressure NSCBC
!!$subroutine apply_P_bc_nscbc
!!$
!!$  use parameters
!!$  implicit none
!!$  integer :: i
!!$  double precision :: L5,L4
!!$  
!!$  do i = 2,nx
!!$     ! Top Boundary
!!$     L5 = (v(i,nx) + Temp(i,nx)**0.5d+0) * &
!!$          ( (P(i,nx) - P(i,nx-1))/dy &
!!$          +  gamma * rho(i,nx) * Temp(i,nx)**0.5d+0 * (v(i,nx)-v(i,nx-1))/dy)
!!$     P(i,nx+1) = P(i,nx) &
!!$          + L5 * dy / (v(i,nx)-Temp(i,nx)**0.5d+0 ) &
!!$          +  gamma * rho(i,nx) * Temp(i,nx)**0.5d+0 * v(i,nx)
!!$     ! Right Boundary
!!$     L5 = (u(nx,i) + Temp(nx,i)**0.5d+0) * &
!!$          ( (P(nx,i) - P(nx-1,i))/dx &
!!$          +  gamma * rho(nx,i) * Temp(nx,i)**0.5d+0 * (u(nx,i)-u(nx-1,i))/dx)
!!$     P(nx+1,i) = P(nx,i) &
!!$          + L5 * dx / (u(nx,i)-Temp(nx,i)**0.5d+0 ) &
!!$          +  gamma * rho(nx,i) * Temp(nx,i)**0.5d+0 * u(nx,i)
!!$     ! Bottom Boundary
!!$     L4 = (v(i,2) - Temp(i,2)**0.5d+0 ) * &
!!$          ( (P(i,3) - P(i,2))/dy &
!!$          -  gamma * rho(i,2) * Temp(i,2)**0.5d+0 * (v(i,3)-v(i,2))/dy)
!!$     P(i,1) = P(i,2) &
!!$          - L4 * dy / (v(i,2) + Temp(i,2)**0.5d+0 ) &
!!$          +  gamma * rho(i,2) * Temp(i,2)**0.5d+0 * v(i,2)
!!$     ! Left Boundary
!!$     L4 = (u(2,i) - Temp(2,i)**0.5d+0 / Ma) * &
!!$          ( (P(3,i) - P(2,i))/dx &
!!$          - Ma * gamma * rho(2,i) * Temp(2,i)**0.5d+0 * (u(3,i)-u(2,i))/dx)
!!$     P(1,i) = P(2,i) &
!!$          - L4 * dx / (u(2,i) + Temp(2,i)**0.5d+0 / Ma) &
!!$          + Ma * gamma * rho(2,i) * Temp(2,i)**0.5d+0 * u(2,i)
!!$     ! Update Densities
!!$     rho(i,nx+1) = P(i,nx+1) /(R* Temp(i,nx+1))
!!$     rho(i,1) = P(i,1) /(R* Temp(i,1))
!!$     rho(nx+1,i) = P(nx+1,i) /(R* Temp(nx+1,i))
!!$     rho(1,i) = P(1,i) /(R* Temp(1,i))
!!$  end do
!!$
!!$end subroutine apply_P_bc_nscbc

subroutine grid_adv
  
  use parameters
  implicit none

  integer :: i
  doubleprecision :: y_w
  
  ! Tracks grid for given case to create images later
  select case (pistvel)
  case('constant')
     do i = 1,nx+1
        y_w =  1d+0 + 2d+0*(L2-L1)*t !(L2+1d+0)/2d+0 !+ (L1-L2)/pi*asin(sin(2d+0*pi*Omega*(t+0.25d+0)))
        y(i) = y_w*1.0d+0*(i-1)/ dble(nx)
     end do
  case('sinusoid')
     do i = 1,nx+1
        y_w = (L1-L2)/2*cos(2d+0*pi*Omega*t) + (L2+1d+0)/2d+0
        y(i) = y_w*1.0d+0*(i-1)/ dble(nx)
     end do
  end select

  dy = y_w / dble(nx)


end subroutine grid_adv
  
! ============================================================= !
!                           Data Dumping                        !
! ============================================================= !

subroutine dumpdata

  use parameters
  implicit none

  integer :: i,j,iunit=8,junit=9
  double precision :: mass

  mass = 0d+0

  open(iunit, file=trim(outfile),position='APPEND')
  write(iunit,*) 'time', t, 'nondimensional ', 'a ', 'b ', 'c'
  write(iunit,*) 'P ','U ', 'V ', 'T ', 'rho ', 'y_grid'
  do i=1,nx+1
     do j = 1,nx+1
        write(iunit,"(E20.10$)") P(i,j), U(i,j), V(i,j), Temp(i,j), rho(i,j), y(j)
        write(iunit,*) ''
     end do
  end do
  close(iunit)

  ! Mass conservation
  do i=2,nx
     do j=2,nx
        mass = mass + rho(i,j)*dx*dy*nx**2
     end do
  end do
  do i=2,nx
     mass = mass + 0.5d+0*rho(1,i)*dx*dy*nx**2
  end do
  do i=2,nx
     mass = mass + 0.5d+0*rho(nx+1,i)*dx*dy*nx**2
  end do
  do i=2,nx
     mass = mass + 0.5d+0*rho(i,1)*dx*dy*nx**2
  end do
  do i=2,nx
     mass = mass + 0.5d+0*rho(i,nx+1)*dx*dy*nx**2
  end do
  mass = mass + 0.25d+0*rho(1,nx+1)*dx*dy*nx**2
  mass = mass + 0.25d+0*rho(nx+1,nx+1)*dx*dy*nx**2
  mass = mass + 0.25d+0*rho(1,1)*dx*dy*nx**2
  mass = mass + 0.25d+0*rho(nx+1,1)*dx*dy*nx**2
  

  open(junit, file=trim('mass'),position='APPEND')
  !write(junit,*) 'time', 'mass'
  write(junit,"(E20.10$)") t, mass
  write(junit,*) ''
  close(junit)

  
end subroutine dumpdata

