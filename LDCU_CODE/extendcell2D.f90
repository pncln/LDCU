    subroutine extendcell2D(rho,m1,m2,E)
    implicit none
    integer, parameter :: Nmax=1200
    real*8, dimension(0:Nmax,0:Nmax) ::rho,m1,m2,E
    integer i,j,Nx,Ny,Nx5,Ny5
    real*8,dimension (1:Nmax) :: x,y
    real*8 xb,yb,r2,c1,u,v,pi,TT,dTT,SS
    character(20) BoundC
    real*8 gamma,dx,dy,x0,y0,xlength,ylength,dt,dt0,Tfinal,t,theta
    common/parameter1/dx,dy,x0,y0,xlength,ylength,dt,dt0,Tfinal,t,theta,gamma
    common/parameter2/Nx,Ny,Nx5,Ny5
    common/parameter3/x,y
    common/parameter4/BoundC
    
    pi=4.d0*datan(1.d0)
    select case (BoundC)
    case('fixed')
         do i=1,Nx
            xb=x0+(i-0.5d0)*dx-5.d0-t
            yb=y0-0.5d0*dy-5.d0-t
            r2=xb**2+yb**2
            rho(i,0)=(1.d0-(gamma-1.d0)*25.d0/(8.d0*gamma*pi**2)*dexp(1.d0-r2))**(1.d0/(gamma-1.d0))
            c1=2.5d0/pi*dexp(0.5d0*(1.d0-r2))
            u=1.d0-c1*yb
            v=1.d0+c1*xb
            m1(i,0)=rho(i,0)*u
            m2(i,0)=rho(i,0)*v
            E(i,0)=(rho(i,0))**gamma/(gamma-1.d0)+0.5d0*(m1(i,0)**2+m2(i,0)**2)/rho(i,0)
            xb=x0+(i-0.5d0)*dx-5.d0-t
            yb=y0+ylength+0.5d0*dy-5.d0-t
            r2=xb**2+yb**2
            rho(i,Ny+1)=(1.d0-(gamma-1.d0)*25.d0/(8.d0*gamma*pi**2)*dexp(1.d0-r2))**(1.d0/(gamma-1.d0))
            c1=2.5d0/pi*dexp(0.5d0*(1.d0-r2))
            u=1.d0-c1*yb
            v=1.d0+c1*xb
            m1(i,Ny+1)=rho(i,Ny+1)*u
            m2(i,Ny+1)=rho(i,Ny+1)*v
            E(i,Ny+1)=(rho(i,Ny+1))**gamma/(gamma-1.d0)+0.5d0*(m1(i,Ny+1)**2+m2(i,Ny+1)**2)/rho(i,Ny+1)
         enddo
          do j=1,Ny
            xb=x0-0.5d0*dx-5.d0-t
            yb=y0+(j-0.5d0)*dy-5.d0-t
            r2=xb**2+yb**2
            rho(0,j)=(1.d0-(gamma-1.d0)*25.d0/(8.d0*gamma*pi**2)*dexp(1.d0-r2))**(1.d0/(gamma-1.d0))
            c1=2.5d0/pi*dexp(0.5d0*(1.d0-r2))
            u=1.d0-c1*yb
            v=1.d0+c1*xb
            m1(0,j)=rho(0,j)*u
            m2(0,j)=rho(0,j)*v
            E(0,j)=(rho(0,j))**gamma/(gamma-1.d0)+0.5d0*(m1(0,j)**2+m2(0,j)**2)/rho(0,j)
            xb=x0+xlength+0.5d0*dx-5.d0-t
            yb=y0+(j-0.5d0)*dy-5.d0-t
            r2=xb**2+yb**2
            rho(Nx+1,j)=(1.d0-(gamma-1.d0)*25.d0/(8.d0*gamma*pi**2)*dexp(1.d0-r2))**(1.d0/(gamma-1.d0))
            c1=2.5d0/pi*dexp(0.5d0*(1.d0-r2))
            u=1.d0-c1*yb
            v=1.d0+c1*xb
            m1(Nx+1,j)=rho(Nx+1,j)*u
            m2(Nx+1,j)=rho(Nx+1,j)*v
            E(Nx+1,j)=(rho(Nx+1,j))**gamma/(gamma-1.d0)+0.5d0*(m1(Nx+1,j)**2+m2(Nx+1,j)**2)/rho(Nx+1,j)
        enddo
    case('periodic')
       ! !==periodic boundary====
        do i=1,Nx
            rho(i,0)=rho(i,Ny)
            rho(i,Ny+1)=rho(i,1)
            m1(i,0)=m1(i,Ny)
            m1(i,Ny+1)=m1(i,1)
            m2(i,0)=m2(i,Ny)
            m2(i,Ny+1)=m2(i,1)
            E(i,0)=E(i,Ny)
            E(i,Ny+1)=E(i,1)
        enddo
        do j=1,Ny
            rho(0,j)=rho(Nx,j)
            rho(Nx+1,j)=rho(1,j)
            m1(0,j)=m1(Nx,j)
            m1(Nx+1,j)=m1(1,j)
            m2(0,j)=m2(Nx,j)
            m2(Nx+1,j)=m2(1,j)
            E(0,j)=E(Nx,j)
            E(Nx+1,j)=E(1,j)
        enddo
    case('reflectfree')
       !! ==reflecting boundary(solid wall)+free====
        do i=1,Nx
            rho(i,0)=rho(i,1)
            rho(i,Ny+1)=rho(i,Ny)
            m1(i,0)=m1(i,1)
            m1(i,Ny+1)=m1(i,Ny)
            m2(i,0)=-m2(i,1)
            m2(i,Ny+1)=m2(i,Ny)
            E(i,0)=E(i,1)
            E(i,Ny+1)=E(i,Ny)
        enddo
        do j=1,Ny
            rho(0,j)=rho(1,j)
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=-m1(1,j)
            m1(Nx+1,j)=m1(Nx,j)
            m2(0,j)=m2(1,j)
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=E(1,j)
            E(Nx+1,j)=E(Nx,j)
        enddo
    case('BCFFSP')
        do j=1,Ny
            rho(0,j)=1.4d0
            m1(0,j)=4.2d0
            m2(0,j)=0.d0
            E(0,j)=8.8d0
            if(j>Ny5)then
             rho(Nx+1,j)=rho(Nx,j)
             m1(Nx+1,j)=m1(Nx,j)
             m2(Nx+1,j)=m2(Nx,j)
             E(Nx+1,j)=E(Nx,j)
            else 
             rho(Nx5,j)=rho(Nx5-1,j)
             m1(Nx5,j)=-m1(Nx5-1,j)
             m2(Nx5,j)=m2(Nx5-1,j)
             E(Nx5,j)=E(Nx5-1,j)
            end if
        end do
        do i=1,Nx
            rho(i,Ny+1)=rho(i,Ny)
            m1(i,Ny+1)=m1(i,Ny)
            m2(i,Ny+1)=-m2(i,Ny)
            E(i,Ny+1)=E(i,Ny)
            if(i<Nx5)then
              rho(i,0)=rho(i,1)
              m1(i,0)=m1(i,1)
              m2(i,0)=-m2(i,1)
              E(i,0)=E(i,1) 
            else
              rho(i,Ny5+1)=rho(i,Ny5+2)
              m1(i,Ny5+1)=m1(i,Ny5+2)
              m2(i,Ny5+1)=-m2(i,Ny5+2)
              E(i,Ny5+1)=E(i,Ny5+2)
            end if    
        end do
     case('free')
        !!==free====
        do i=1,Nx
            rho(i,0)=rho(i,1)
            rho(i,Ny+1)=rho(i,Ny)
            m1(i,0)=m1(i,1)
            m1(i,Ny+1)=m1(i,Ny)
            m2(i,0)=m2(i,1)
            m2(i,Ny+1)=m2(i,Ny)
            E(i,0)=E(i,1)
            E(i,Ny+1)=E(i,Ny)
        enddo
        do j=1,Ny
            rho(0,j)=rho(1,j)
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=m1(1,j)
            m1(Nx+1,j)=m1(Nx,j)
            m2(0,j)=m2(1,j)
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=E(1,j)
            E(Nx+1,j)=E(Nx,j)
        enddo
      case('fixedleft1')
        !!==free====
        do i=1,Nx
            rho(i,0)=rho(i,1)
            rho(i,Ny+1)=rho(i,Ny)
            m1(i,0)=m1(i,1)
            m1(i,Ny+1)=m1(i,Ny)
            m2(i,0)=m2(i,1)
            m2(i,Ny+1)=m2(i,Ny)
            E(i,0)=E(i,1)
            E(i,Ny+1)=E(i,Ny)
        enddo
        do j=1,Ny
          if(dabs(y(j))<0.05d0)then  
            rho(0,j)=5.d0
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=150.d0
            m1(Nx+1,j)=m1(Nx,j)
            m2(0,j)=0.d0
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=2250.d0+0.4127d0/(gamma-1.d0)
            E(Nx+1,j)=E(Nx,j)
          else 
            rho(0,j)=rho(1,j)
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=m1(1,j)
            m1(Nx+1,j)=m1(Nx,j)
            m2(0,j)=m2(1,j)
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=E(1,j)
            E(Nx+1,j)=E(Nx,j)
          end if    
        enddo
      case('fixedleft2')
        !!==free====
        do i=1,Nx
            rho(i,0)=rho(i,1)
            rho(i,Ny+1)=rho(i,Ny)
            m1(i,0)=m1(i,1)
            m1(i,Ny+1)=m1(i,Ny)
            m2(i,0)=m2(i,1)
            m2(i,Ny+1)=m2(i,Ny)
            E(i,0)=E(i,1)
            E(i,Ny+1)=E(i,Ny)
        enddo
        do j=1,Ny
          if(dabs(y(j))<0.05d0)then  
            rho(0,j)=5.d0
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=4000.d0
            m1(Nx+1,j)=m1(Nx,j)
            m2(0,j)=0.d0
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=1.6d6+0.4127d0/(gamma-1.d0)
            E(Nx+1,j)=E(Nx,j)
          else 
            rho(0,j)=rho(1,j)
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=m1(1,j)
            m1(Nx+1,j)=m1(Nx,j)
            m2(0,j)=m2(1,j)
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=E(1,j)
            E(Nx+1,j)=E(Nx,j)
          end if    
        enddo
       case('fixedleft3')
        !!==free====
        do i=1,Nx
            rho(i,0)=rho(i,1)
            rho(i,Ny+1)=rho(i,Ny)
            m1(i,0)=m1(i,1)
            m1(i,Ny+1)=m1(i,Ny)
            m2(i,0)=m2(i,1)
            m2(i,Ny+1)=m2(i,Ny)
            E(i,0)=E(i,1)
            E(i,Ny+1)=E(i,Ny)
        enddo
       do j=1,Ny
            rho(0,j)=5.26829268d0
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=rho(0,j)*4.86111111d0
            m1(Nx+1,j)=m1(Nx,j)
            m2(0,j)=0.d0
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=29.88095238d0/(gamma-1.d0)+0.5d0*m1(0,j)**2/rho(0,j)
            E(Nx+1,j)=E(Nx,j)
        enddo
     case('solidwall')
       !! ==reflecting boundary(solid wall)====
        do i=1,Nx
            rho(i,0)=rho(i,1)
            rho(i,Ny+1)=rho(i,Ny)
            m1(i,0)=m1(i,1)
            m1(i,Ny+1)=m1(i,Ny)
            m2(i,0)=-m2(i,1)
            m2(i,Ny+1)=-m2(i,Ny)
            E(i,0)=E(i,1)
            E(i,Ny+1)=E(i,Ny)
        enddo
        do j=1,Ny
            rho(0,j)=rho(1,j)
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=-m1(1,j)
            m1(Nx+1,j)=-m1(Nx,j)
            m2(0,j)=m2(1,j)
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=E(1,j)
            E(Nx+1,j)=E(Nx,j)
        enddo
      case('solidwall1')
       !! ==reflecting boundary(solid wall)====
        do i=1,Nx
            rho(i,0)=2.d0
            rho(i,Ny+1)=1.d0
            m1(i,0)=0.d0
            m1(i,Ny+1)=0.d0
            m2(i,0)=0.d0
            m2(i,Ny+1)=0.d0
            E(i,0)=1.d0/(gamma-1.d0)
            E(i,Ny+1)=2.5d0/(gamma-1.d0)
        enddo
        do j=1,Ny
            rho(0,j)=rho(1,j)
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=-m1(1,j)
            m1(Nx+1,j)=-m1(Nx,j)
            m2(0,j)=m2(1,j)
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=E(1,j)
            E(Nx+1,j)=E(Nx,j)
        enddo
    end select
    
    return
    end subroutine