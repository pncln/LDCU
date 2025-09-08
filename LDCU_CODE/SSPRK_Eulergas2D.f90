    subroutine SSPRK_Eulergas2D(rho0,mx0,my0,E0,F1,F2,F3,F4,G1,G2,G3,G4,rhonew,mxnew,mynew,Enew)
    implicit none 
    integer,parameter :: Nmax=1200
    real*8,dimension (0:Nmax,0:Nmax) ::rho0,mx0,my0,E0,rho1,mx1,my1,E1,rho2,mx2,my2,E2,rhonew,mxnew,mynew,Enew
    real*8,dimension (0:Nmax,1:Nmax):: F1,F2,F3,F4
    real*8,dimension (1:Nmax,0:Nmax):: G1,G2,G3,G4
    real*8,dimension (1:Nmax) :: x,y
    integer i,j,Nx,Ny,newdtneed,Nx5,Ny5
    real*8 dx,dy,x0,y0,xlength,ylength,dt,dt0,dt1,Tfinal,t,theta,gamma,mu,lambda,grav
    common/parameter1/dx,dy,x0,y0,xlength,ylength,dt,dt0,Tfinal,t,theta,gamma
    common/parameter2/Nx,Ny,Nx5,Ny5
    common/parameter3/x,y
!105 dt=dt0
    lambda=dt/dx
    mu=dt/dy
    do j=1,Ny
        do i=1,Nx
            !if((i<Nx5) .or. (j>Ny5))then
            rho1(i,j)=rho0(i,j)-(lambda*(F1(i,j)-F1(i-1,j))+mu*(G1(i,j)-G1(i,j-1)))
            mx1(i,j)=mx0(i,j)-(lambda*(F2(i,j)-F2(i-1,j))+mu*(G2(i,j)-G2(i,j-1)))
            my1(i,j)=my0(i,j)-(lambda*(F3(i,j)-F3(i-1,j))+mu*(G3(i,j)-G3(i,j-1)))
            E1(i,j)=E0(i,j)-(lambda*(F4(i,j)-F4(i-1,j))+mu*(G4(i,j)-G4(i,j-1)))
            !end if
       enddo
    enddo
    
    dt1=dt
    call extendcell2D(rho1,mx1,my1,E1)
    newdtneed=0
    call numerical_flux_Euler2D1(rho1,mx1,my1,E1,newdtneed,F1,F2,F3,F4,G1,G2,G3,G4)
    !if(dt0<dt1)then
    !  goto 105
    !end if 
    do i=1,Nx
        do j=1,Ny
            !if((i<Nx5) .or. (j>Ny5))then
        rho2(i,j)=((3.d0*rho0(i,j)+rho1(i,j))-(lambda*(F1(i,j)-F1(i-1,j))+mu*(G1(i,j)-G1(i,j-1))))/4.d0
        mx2(i,j)=((3.d0*mx0(i,j)+mx1(i,j))-(lambda*(F2(i,j)-F2(i-1,j))+mu*(G2(i,j)-G2(i,j-1))))/4.d0
        my2(i,j)=((3.d0*my0(i,j)+my1(i,j))-(lambda*(F3(i,j)-F3(i-1,j))+mu*(G3(i,j)-G3(i,j-1))))/4.d0
        E2(i,j)=((3.d0*E0(i,j)+E1(i,j))-(lambda*(F4(i,j)-F4(i-1,j))+mu*(G4(i,j)-G4(i,j-1))))/4.d0 
        !end if
        enddo
    enddo
    
    call extendcell2D(rho2,mx2,my2,E2)
    newdtneed=0
    call numerical_flux_Euler2D1(rho2,mx2,my2,E2,newdtneed,F1,F2,F3,F4,G1,G2,G3,G4)
    !if(dt0<dt1)then
    !  goto 105
    !end if 
    do i=1,Nx
        do j=1,Ny
            !if((i<Nx5) .or. (j>Ny5))then
        rhonew(i,j)=((rho0(i,j)+2.d0*rho2(i,j))-2.d0*(lambda*(F1(i,j)-F1(i-1,j))+mu*(G1(i,j)-G1(i,j-1))))/3.d0
        mxnew(i,j)=((mx0(i,j)+2.d0*mx2(i,j))-2.d0*(lambda*(F2(i,j)-F2(i-1,j))+mu*(G2(i,j)-G2(i,j-1))))/3.d0
        mynew(i,j)=((my0(i,j)+2.d0*my2(i,j))-2.d0*(lambda*(F3(i,j)-F3(i-1,j))+mu*(G3(i,j)-G3(i,j-1))))/3.d0
        Enew(i,j)=((E0(i,j)+2.d0*E2(i,j))-2.d0*(lambda*(F4(i,j)-F4(i-1,j))+mu*(G4(i,j)-G4(i,j-1))))/3.d0
        !end if
        enddo
    enddo

    return
    end subroutine
    