    subroutine initialset_2D(rho,m1,m2,E)
    implicit none
    integer,parameter :: Nmax=1200
    real*8,dimension (0:Nmax,0:Nmax) :: rho,m1,m2,E,p
    real*8,dimension (1:Nmax) :: x,y
    character(20) BoundC,Extype
    real*8 alpha1,alpha2,beta1,beta2,pi
    real*8 eps,J1,J2,I1,I2,Y1,Y2,u,v,ei
    real*8 grav,yint
    real*8 xW,xE,yN,yS
    real*8 ::Ysjn(3),Xwie(3)
    real*8,dimension (1:3,1:3,4) ::qk
    integer i,j,k1,k2,Nx,Ny,Nx5,Ny5
    real*8 dx,dy,x0,y0,xlength,ylength,dt,dt0,Tfinal,t,theta,gamma,r
    common/parameter1/dx,dy,x0,y0,xlength,ylength,dt,dt0,Tfinal,t,theta,gamma
    common/parameter2/Nx,Ny,Nx5,Ny5
    common/parameter3/x,y
    common/parameter4/BoundC,Extype

    do i=1,Nx
        x(i)=x0+(i-0.5d0)*dx
    enddo
    do j=1,Ny
        y(j)=y0+(j-0.5d0)*dy
    enddo
    pi=4.d0*atan(1.d0)

    select case(Extype)
    case ('MCW')
        !!=========Moving Contact Wave=========
        BoundC='free'
        !do j=1,Ny
        !    yN=y(j)+0.5d0*dy
        !    yS=y(j)-0.5d0*dy
        !    Ysjn=(/yS,y(j),yN/)
        !    do i=1,Nx
        !        xE=x(i)+0.5d0*dx
        !        xW=x(i)-0.5d0*dx
        !        Xwie=(/xW,x(i),xE/)
        !        do k2=1,3
        !            do k1=1,3
        !                if ((Xwie(k1)<=0.1d0 .and. Xwie(k1)>=-0.1d0 .and. Ysjn(k2)<=0.02d0) .or. &
        !                (Xwie(k1)<=0.02d0 .and. Xwie(k1)>=-0.02d0 .and. Ysjn(k2)>=0.02d0 .and. Ysjn(k2)<=0.1d0) .or. &
        !                (Xwie(k1)+0.02d0)**2+(Ysjn(k2)-0.02d0)**2<=0.08d0**2 .or. &
        !                (Xwie(k1)-0.02d0)**2+(Ysjn(k2)-0.02d0)**2<=0.08d0**2) then
        !                    qk(k1,k2,1)=1.4d0
        !                    !qk(k1,k2,2)=0.28d0
        !                    qk(k1,k2,2)=0.d0
        !                    qk(k1,k2,3)=0.28d0
        !                    qk(k1,k2,4)=1.d0/(gamma-1.d0)+0.5d0*(qk(k1,k2,2)**2+qk(k1,k2,3)**2)/qk(k1,k2,1)
        !                else
        !                    qk(k1,k2,1)=1.0d0
        !                    !qk(k1,k2,2)=0.2d0
        !                    qk(k1,k2,2)=0.d0
        !                    qk(k1,k2,3)=0.2d0
        !                    qk(k1,k2,4)=1.d0/(gamma-1.d0)+0.5d0*(qk(k1,k2,2)**2+qk(k1,k2,3)**2)/qk(k1,k2,1)
        !                endif
        !            enddo
        !        enddo
        !        rho(i,j)=(qk(1,1,1)+qk(3,1,1)+qk(1,3,1)+qk(3,3,1)+4.d0*(qk(2,1,1)+qk(1,2,1)+qk(3,2,1)+qk(2,3,1))+16.d0*qk(2,2,1))/36.d0
        !        m1(i,j)=(qk(1,1,2)+qk(3,1,2)+qk(1,3,2)+qk(3,3,2)+4.d0*(qk(2,1,2)+qk(1,2,2)+qk(3,2,2)+qk(2,3,2))+16.d0*qk(2,2,2))/36.d0
        !        m2(i,j)=(qk(1,1,3)+qk(3,1,3)+qk(1,3,3)+qk(3,3,3)+4.d0*(qk(2,1,3)+qk(1,2,3)+qk(3,2,3)+qk(2,3,3))+16.d0*qk(2,2,3))/36.d0
        !        E(i,j)=(qk(1,1,4)+qk(3,1,4)+qk(1,3,4)+qk(3,3,4)+4.d0*(qk(2,1,4)+qk(1,2,4)+qk(3,2,4)+qk(2,3,4))+16.d0*qk(2,2,4))/36.d0
        !    enddo
        !enddo

        do j=1,Ny
            do i=1,Nx
                if ((x(i)<=0.1d0 .and. x(i)>=-0.1d0 .and. y(j)<=0.02d0) .or. &
                (x(i)<=0.02d0 .and. x(i)>=-0.02d0 .and. y(j)>=0.02d0 .and. y(j)<=0.1d0) .or. &
                (x(i)+0.02d0)**2+(y(j)-0.02d0)**2<=0.08d0**2 .or. &
                (x(i)-0.02d0)**2+(y(j)-0.02d0)**2<=0.08d0**2) then
                    rho(i,j)=1.4d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.28d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=1.0d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.2d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo

        !do j=1,Ny
        !    do i=1,Nx
        !        if ((x(i)<=0.1d0 .and. y(j)<=0.02d0) .or. (x(i)<=0.02d0 .and. y(j)<=0.1d0) &
        !            .or. (x(i)-0.02d0)**2+(y(j)-0.02d0)**2<=0.08d0**2) then
        !            rho(i,j)=1.4d0
        !            m1(i,j)=0.28d0
        !            m2(i,j)=0.28d0
        !            E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
        !        else
        !            rho(i,j)=1.0d0
        !            m1(i,j)=0.2d0
        !            m2(i,j)=0.2d0
        !            E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
        !        endif
        !    enddo
        !enddo
    case ('Exp')
        !!=========Explosion=========
        BoundC='reflectfree'
        do j=1,Ny
            do i=1,Nx
                if (x(i)**2+y(j)**2<0.16d0-1.d-8) then
                    rho(i,j)=1.d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=1.d0/(gamma-1.d0)
                else
                    rho(i,j)=0.125d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=0.1d0/(gamma-1.d0)
                endif
            enddo
        enddo
    case('FFSP')    
        BoundC='BCFFSP'
        do j=1,Ny
         do i=1,Nx
           if((i<Nx5) .or. (j>Ny5))then
            rho(i,j)=1.4d0
            m1(i,j)=4.2d0
            m2(i,j)=0.d0
            E(i,j)=8.8d0 
           end if 
         end do       
        end do      
      case('accuracy')    
        BoundC='periodic'
        do j=1,Ny
         do i=1,Nx
            rho(i,j)=1.d0+0.5d0*sin(pi*(x(i)+y(j)))
            m1(i,j)=rho(i,j)
            m2(i,j)=-0.7d0*rho(i,j)
            E(i,j)=2.5d0+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
         end do       
        end do
    case('2DRP')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx
                if(x(i)>1.d0.and.y(j)>1.d0) then
                    rho(i,j)=1.5d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=1.5d0/(gamma-1.d0)
                else if(x(i)<1.d0.and.y(j)>1.d0)then
                    rho(i,j)=0.5323d0
                    m1(i,j)=0.5323d0*1.206d0
                    m2(i,j)=0.d0
                    E(i,j)=0.3d0/(gamma-1.d0)+0.5d0*0.5323d0*1.206d0**2
                else if(x(i)<1.d0.and.y(j)<1.d0)then
                    rho(i,j)=0.138d0
                    m1(i,j)=0.138d0*1.206d0
                    m2(i,j)=0.138d0*1.206d0
                    E(i,j)=0.029d0/(gamma-1.d0)+0.5d0*0.138d0*(1.206d0**2+1.206d0**2)
                else
                    rho(i,j)=0.5323d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.5323d0*1.206d0
                    E(i,j)=0.3d0/(gamma-1.d0)+0.5d0*0.5323d0*1.206d0**2
                end if
            end do
        end do    
    case ('Imp')
        !!=========implosion=========
        BoundC='solidwall'
        do j=1,Ny
            do i=1,Nx
                if (dabs(x(i))+dabs(y(j))<0.15d0-1.0d-8) then
                    rho(i,j)=0.125d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=0.14d0/(gamma-1.d0)
                else
                    rho(i,j)=1.d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=1.d0/(gamma-1.d0)
                endif
            enddo
        enddo
    case ('KH')
        !!=========Kelvin¨CHelmholtz=========
        BoundC='periodic'
        pi=4.d0*datan(1.d0)
        I1=0.00625d0
        do j=1,Ny
            do i=1,Nx
                if (y(j)<-0.25d0) then
                    rho(i,j)=1.d0
                    u=-0.5d0+0.5d0*dexp((y(j)+0.25d0)/I1)
                    v=0.01d0*dsin(4.d0*pi*x(i))
                    ei=3.75d0
                    m1(i,j)=rho(i,j)*u
                    m2(i,j)=rho(i,j)*v
                    E(i,j)=rho(i,j)*(ei+0.5d0*(u**2+v**2))
                elseif(y(j)<0.d0) then
                    rho(i,j)=2.d0
                    u=0.5d0-0.5d0*dexp((-y(j)-0.25d0)/I1)
                    v=0.01d0*dsin(4.d0*pi*x(i))
                    ei=1.875d0
                    m1(i,j)=rho(i,j)*u
                    m2(i,j)=rho(i,j)*v
                    E(i,j)=rho(i,j)*(ei+0.5d0*(u**2+v**2))
                elseif (y(j)<0.25d0) then
                    rho(i,j)=2.d0
                    u=0.5d0-0.5d0*dexp((y(j)-0.25d0)/I1)
                    v=0.01d0*dsin(4.d0*pi*x(i))
                    ei=1.875d0
                    m1(i,j)=rho(i,j)*u
                    m2(i,j)=rho(i,j)*v
                    E(i,j)=rho(i,j)*(ei+0.5d0*(u**2+v**2))
                else
                    rho(i,j)=1.d0
                    u=-0.5d0+0.5d0*dexp((-y(j)+0.25d0)/I1)
                    v=0.01d0*dsin(4.d0*pi*x(i))
                    ei=3.75d0
                    m1(i,j)=rho(i,j)*u
                    m2(i,j)=rho(i,j)*v
                    E(i,j)=rho(i,j)*(ei+0.5d0*(u**2+v**2))
                endif
            enddo
        enddo
    case ('HighMach1')    
        BoundC='fixedleft1'
        do j=1,Ny
            do i=1,Nx               
                    rho(i,j)=0.5d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=0.4127d0/(gamma-1.d0)
            enddo
        enddo
        
    case ('HighMach2')    
        BoundC='fixedleft2'
        do j=1,Ny
            do i=1,Nx               
                    rho(i,j)=0.5d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=0.4127d0/(gamma-1.d0)
            enddo
        enddo    

    case ('Vortex')
        !!=========accuracy test Vortex in Shu's paper=========
        BoundC='periodic'
        pi=4.d0*datan(1.d0)
        !do j=1,Ny
        !    do i=1,Nx
        !        r=(y(j)-5.d0)**2+(5.d0-x(i))**2
        !        rho(i,j)=(1.d0-(gamma-1.d0)*5.d0**2/(8.d0*gamma*pi**2)*dexp(1.d0-r**2))**(1.d0/(gamma-1.d0))
        !        u=1.d0+2.5d0/pi*dexp(0.5d0*(1.d0-r**2))*(x(i)-5.d0)
        !        v=1.d0+2.5d0/pi*dexp(0.5d0*(1.d0-r**2))*(y(j)-5.d0)
        !        m1(i,j)=rho(i,j)*u
        !        m2(i,j)=rho(i,j)*v
        !        E(i,j)=(rho(i,j))**gamma/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
        !    enddo
        !enddo
        do j=1,Ny
            yN=y(j)+0.5d0*dy
            yS=y(j)-0.5d0*dy
            Ysjn=(/yS,y(j),yN/)
            do i=1,Nx
                xE=x(i)+0.5d0*dx
                xW=x(i)-0.5d0*dx
                Xwie=(/xW,x(i),xE/)
                do k2=1,3
                    do k1=1,3
                        r=(Ysjn(k2)-5.d0)**2+(5.d0-Xwie(k1))**2
                        qk(k1,k2,1)=(1.d0-(gamma-1.d0)*25.d0/(8.d0*gamma*pi**2)*dexp(1.d0-r**2))**(1.d0/(gamma-1.d0))
                        u=1.d0-2.5d0/pi*dexp(0.5d0*(1.d0-r**2))*(Ysjn(k2)-5.d0)
                        v=1.d0+2.5d0/pi*dexp(0.5d0*(1.d0-r**2))*(Xwie(k1)-5.d0)
                        qk(k1,k2,2)=qk(k1,k2,1)*u
                        qk(k1,k2,3)=qk(k1,k2,1)*v
                        qk(k1,k2,4)=qk(k1,k2,1)**gamma/(gamma-1.d0)+0.5d0*(qk(k1,k2,2)**2+qk(k1,k2,3)**2)/qk(k1,k2,1)
                    enddo
                enddo
                rho(i,j)=(qk(1,1,1)+qk(3,1,1)+qk(1,3,1)+qk(3,3,1)+4.d0*(qk(2,1,1)&
                    +qk(1,2,1)+qk(3,2,1)+qk(2,3,1))+16.d0*qk(2,2,1))/36.d0
                m1(i,j)=(qk(1,1,2)+qk(3,1,2)+qk(1,3,2)+qk(3,3,2)+4.d0*(qk(2,1,2)&
                    +qk(1,2,2)+qk(3,2,2)+qk(2,3,2))+16.d0*qk(2,2,2))/36.d0
                m2(i,j)=(qk(1,1,3)+qk(3,1,3)+qk(1,3,3)+qk(3,3,3)+4.d0*(qk(2,1,3)&
                    +qk(1,2,3)+qk(3,2,3)+qk(2,3,3))+16.d0*qk(2,2,3))/36.d0
                E(i,j)=(qk(1,1,4)+qk(3,1,4)+qk(1,3,4)+qk(3,3,4)+4.d0*(qk(2,1,4)&
                    +qk(1,2,4)+qk(3,2,4)+qk(2,3,4))+16.d0*qk(2,2,4))/36.d0
            enddo
        enddo
    case ('DMR')
        BoundC='parsolidwall'
        do j=1,Ny
            do i=1,Nx
                if (x(i)<1.d0/6.d0+y(j)/dsqrt(3.d0)) then
                    rho(i,j)=8.d0
                    m1(i,j)=33.d0*dsqrt(3.d0)
                    m2(i,j)=-33.d0
                    E(i,j)=116.5d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=1.4d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=1.d0/(gamma-1.d0)
                endif
            enddo
        enddo
    case ('Quirk')
        BoundC='fixedleft3'
        do j=1,Ny
            do i=1,Nx
                    rho(i,j)=1.d0+1.d-3
                    m1(i,j)=0.d0+1.d-3
                    m2(i,j)=0.d0+1.d-3
                    E(i,j)=1.d0/(gamma*(gamma-1.d0))+1.d-3
            enddo
        enddo
    end select

    !pi=4.d0*datan(1.d0)
    !alpha1=1.d0
    !alpha2=1.d0
    !!call random_seed()
    !!call random_number(beta1)
    !!call random_number(beta2)
    !!beta1=pi*(2.d0*beta1-1.d0)
    !!beta2=pi*(2.d0*beta2-1.d0)
    !beta1=0.5d0*pi
    !beta2=0.5d0*pi
    !eps=5.0d-3
    !J1=0.25d0
    !J2=0.75d0
    !
    !do j=1,Ny
    !    do i=1,Nx
    !        Y1=alpha1*dcos(beta1+2.d0*pi*x(i))
    !        I1=J1+eps*Y1
    !        Y2=alpha2*dcos(beta2+2.d0*pi*x(i))
    !        I2=J2+eps*Y2
    !        if (y(j)<I1 .or. y(j)>I2) then
    !            rho(i,j)=1.d0
    !            !m1(i,j)=0.5d0*rho(i,j)
    !            m1(i,j)=rho(i,j)
    !            m2(i,j)=0.d0
    !            E(i,j)=2.5d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
    !        else
    !            rho(i,j)=2.d0
    !            !m1(i,j)=-0.5d0*rho(i,j)
    !            m1(i,j)=-rho(i,j)
    !            m2(i,j)=0.d0
    !            E(i,j)=2.5d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
    !        endif
    !    enddo
    !enddo




    !!!==RSC=====
    !do j=1,Ny
    !    do i=1,Nx
    !        if (x(i)<0.5d0 .and. y(j)<0.5d0) then
    !            rho(i,j)=1.d0
    !            m1(i,j)=-0.75d0*rho(i,j)
    !            m2(i,j)=0.5d0*rho(i,j)
    !            E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
    !        else if (x(i)<0.5d0 .and. y(j)>0.5d0) then
    !            rho(i,j)=2.d0
    !            m1(i,j)=0.75d0*rho(i,j)
    !            m2(i,j)=0.5d0*rho(i,j)
    !            E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
    !        else if (x(i)>0.5d0 .and. y(j)>0.5d0) then
    !            rho(i,j)=1.d0
    !            m1(i,j)=0.75d0*rho(i,j)
    !            m2(i,j)=-0.5d0*rho(i,j)
    !            E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
    !        else
    !            rho(i,j)=3.d0
    !            m1(i,j)=-0.75d0*rho(i,j)
    !            m2(i,j)=-0.5d0*rho(i,j)
    !            E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
    !        endif
    !    enddo
    !enddo



    !!!=========Moving Contact Wave=========
    !do j=1,Ny
    !    yN=y(j)+0.5d0*dy
    !    yS=y(j)-0.5d0*dy
    !    Ysjn=(/yS,y(j),yN/)
    !    do i=1,Nx
    !        xE=x(i)+0.5d0*dx
    !        xW=x(i)-0.5d0*dx
    !        Xwie=(/xW,x(i),xE/)
    !        do k2=1,3
    !            do k1=1,3
    !            if ((Xwie(k1)<=0.1d0 .and. Ysjn(k2)<=0.02d0) .or. (Xwie(k1)<=0.02d0 .and. Ysjn(k2)<=0.1d0) .or. (Xwie(k1)-0.02d0)**2+(Ysjn(k2)-0.02d0)**2<=0.08d0**2) then
    !                qk(k1,k2,1)=1.4d0
    !                !qk(k1,k2,2)=0.28d0
    !                qk(k1,k2,2)=0.d0
    !                qk(k1,k2,3)=0.28d0
    !                qk(k1,k2,4)=1.d0/(gamma-1.d0)+0.5d0*(qk(k1,k2,2)**2+qk(k1,k2,3)**2)/qk(k1,k2,1)
    !            else
    !                qk(k1,k2,1)=1.0d0
    !                !qk(k1,k2,2)=0.2d0
    !                qk(k1,k2,2)=0.d0
    !                qk(k1,k2,3)=0.2d0
    !                qk(k1,k2,4)=1.d0/(gamma-1.d0)+0.5d0*(qk(k1,k2,2)**2+qk(k1,k2,3)**2)/qk(k1,k2,1)
    !            endif
    !            enddo
    !        enddo
    !        rho(i,j)=(qk(1,1,1)+qk(3,1,1)+qk(1,3,1)+qk(3,3,1)+4.d0*(qk(2,1,1)+qk(1,2,1)+qk(3,2,1)+qk(2,3,1))+16.d0*qk(2,2,1))/36.d0
    !        m1(i,j)=(qk(1,1,2)+qk(3,1,2)+qk(1,3,2)+qk(3,3,2)+4.d0*(qk(2,1,2)+qk(1,2,2)+qk(3,2,2)+qk(2,3,2))+16.d0*qk(2,2,2))/36.d0
    !        m2(i,j)=(qk(1,1,3)+qk(3,1,3)+qk(1,3,3)+qk(3,3,3)+4.d0*(qk(2,1,3)+qk(1,2,3)+qk(3,2,3)+qk(2,3,3))+16.d0*qk(2,2,3))/36.d0
    !        E(i,j)=(qk(1,1,4)+qk(3,1,4)+qk(1,3,4)+qk(3,3,4)+4.d0*(qk(2,1,4)+qk(1,2,4)+qk(3,2,4)+qk(2,3,4))+16.d0*qk(2,2,4))/36.d0
    !    enddo
    !enddo


    call extendcell2D(rho,m1,m2,E)

    return
    end subroutine

