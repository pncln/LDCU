    program Euler_gas_main_2D
    implicit none 
    integer,parameter :: Nmax=1200
    real*8,dimension (0:Nmax,0:Nmax) ::rho0,mx0,my0,E0,rho,mx,my,E 
    real*8,dimension (1:Nmax) :: x,y
    real*8,dimension (0:Nmax,1:Nmax):: F1,F2,F3,F4,alphaI
    real*8,dimension (1:Nmax,0:Nmax):: G1,G2,G3,G4,alphaJ
    integer i,j,Nx,Ny,newdtneed,nstep,istep,Nx5,Ny5
    character(40) filename,rhon,un,vn,En,alphan,betan
    character(20) BoundC,Extype,scaleSw,shrinkSw
    real*8 pstar,lambda_max,p1,p2,t1
    real*8 gamma,dx,dy,x0,y0,xlength,ylength,dt,dt0,Tfinal,t,theta,dtstep,dxinv,dyinv
    common/parameter1/dx,dy,x0,y0,xlength,ylength,dt,dt0,Tfinal,t,theta,gamma,dxinv,dyinv
    common/parameter2/Nx,Ny,Nx5,Ny5
    common/parameter3/x,y
    common/parameter4/BoundC,Extype,scaleSw,shrinkSw
    common/parameter5/alphaI,alphaJ
    !=======input the data========================
    
    Nx=600
    Ny=600
    tfinal=2.5d0
    x0=0.d0
    y0=0.d0
    xlength=0.3d0
    ylength=0.3d0
    !Nx5=Nx/5+1
    !Ny5=Ny/5
    !RT 
    !Nx=200
    !Ny=800
    !Tfinal=1.95d0
    !x0=-0.2d0
    !y0=0.d0
    !xlength=0.4d0
    !ylength=0.8d0

    dx=xlength/Nx
    dy=ylength/Ny
    dxinv=1.d0/dx
    dyinv=1.d0/dy
    theta=1.3d0
    gamma=1.4d0
    !gamma=5.d0/3.d0 !!RT example
    nstep=1
    dtstep=Tfinal/nstep
    
    !Extype='2DRP'
    !Extype='MCW'
    !Extype='Exp'
    !Extype='DMR'
    !Extype='Quirk'
    Extype='Imp'
    !Extype='accuracy'
    !Extype='KH'
    !Extype='RT'
    !Extype='FFSP'
    !Extype='HighMach2'
    !Extype='Cfg3'
   ! Extype='Vortex'
    !shrinkSw='off'
    !scaleSw='off'  
    !MLLF='off'
    
    !!===========================initial data=====================================================
    call initialset_2D(rho0,mx0,my0,E0)
    !call initialset_2DRT(rho0,mx0,my0,E0)
    !call ini_Riemann(rho0,mx0,my0,E0)
    
    t=0.0d0 
    do istep=0,nstep
        Tfinal=istep*dtstep
        !if(istep==1)then
        !Tfinal=0.001d0
        !end if
        !if(istep==2)then
        !Tfinal=0.0015d0
        !end if
        !if(istep==3)then
        !Tfinal=4.d0
        !end if
        do while (t<Tfinal-1.0d-12)
            newdtneed=1
            call  numerical_flux_Euler2D1(rho0,mx0,my0,E0,newdtneed,F1,F2,F3,F4,G1,G2,G3,G4)
            call  SSPRK_Eulergas2D(rho0,mx0,my0,E0,F1,F2,F3,F4,G1,G2,G3,G4,rho,mx,my,E)
            !call  numerical_flux_RT(rho0,mx0,my0,E0,newdtneed,F1,F2,F3,F4,G1,G2,G3,G4)
            !call  SSPRK_EulergasRT(rho0,mx0,my0,E0,F1,F2,F3,F4,G1,G2,G3,G4,rho,mx,my,E)
            call extendcell2d(rho,mx,my,E)
            t=t+dt
            write(*,*) 't=', t
            rho0=rho
            mx0=mx
            my0=my
            E0=E
        enddo
    
        write(rhon,101) istep
        open(11,file=rhon)
        !! .....Writing the solution.............................................
        do j=1,Ny
            do i=1,Nx
                write(11,100) x(i),y(j),rho0(i,j)
            enddo
        enddo
        close(11)
        
        !write(alphan,102) istep
        !open(12,file=alphan)
        !! .....Writing the solution.............................................
        !do j=1,Ny
        !    do i=1,Nx
        !        write(12,100) x(i),y(j),(E0(i,j)-0.5d0*(mx0(i,j)**2+my0(i,j)**2)/rho0(i,j))*(gamma-1.d0)
        !    enddo
        !enddo
        !close(12)
    !    !
        !write(betan,103) istep
        !open(13,file=betan)
        !!! .....Writing the solution.............................................
        !do j=1,Ny
        !    do i=1,Nx
        !        !write(13,*) x(i),y(j),alphaJ(i,j)
        !        write(13,*) x(i),y(j),E0(i,j)
        !    enddo
        !enddo
    !    close(13)
    !    write(betan,104) istep
    !    open(14,file=betan)
    !    !! .....Writing the solution.............................................
    !    do j=1,Ny
    !        do i=1,Nx
    !            !write(13,*) x(i),y(j),alphaJ(i,j)
    !            write(14,*) x(i),y(j),E0(i,j)
    !        enddo
    !    enddo
    !    close(14)
     call CPU_TIME(t1)
     write(*,*)t1
    enddo
!102 format("out1accmx400",i3.3)
!103 format("out1accmy400",i3.3)  
!104 format("out1accE400",i3.3)
!101 format("out1DMR",i3.3)
!101 format("out111MCWrhonewV2",i3.3)
!101 format("out12DRP",i3.3)
!101 format("out1CPU",i3.3)
!101 format("out1KH1rhooldV",i3.3)
101 format("outImp",i3.3)
!101 format("out1FFSPrhonewVVVVV",i3.3)
!101 format("out00Cfg3rho",i3.3)
!102 format("out1p",i3.3)   
!103 format("out12E",i3.3)
!102 format("out1KH1alpha",i3.3)   
!103 format("out1KH1beta",i3.3)
    
100 format (' ', f24.16,' ', f24.16,' ', f24.16,' ', f24.16,' ', f24.16,' ', f24.16)
    end program

   real*8 function dminmod(v1,v2,v3)
    implicit none
    real*8 v1,v2,v3,u
    dminmod=0.5d0*(dsign(1.d0,v1)+dsign(1.d0,v3))*dmin1(dabs(v1),dabs(v2),dabs(v3))
    return
    end function
    
    real*8 function dminmod4(a,b,c,d)
    implicit none
    real*8 a,b,c,d,e,f
    e = dmin1(a,b,c,d)
    f = dmax1(a,b,c,d)
    dminmod4 = 0.d0
    if(e>0.d0) then
        dminmod4 = e
    elseif(f<0.d0) then
        dminmod4 = f
    endif
    return
    end function
    
    real*8 function minmod(x,y)
    implicit none
    real*8 x,y
    minmod=0.5d0*(dsign(1.d0,x)+dsign(1.d0,y))*dmin1(dabs(x),dabs(y))
    return
    end function

