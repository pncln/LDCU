   subroutine numerical_flux_Euler2D1(rho,mx,my,E,newdtneed,F1,F2,F3,F4,G1,G2,G3,G4)
    implicit none 
    integer,parameter :: Nmax=1200
    real*8, dimension (0:Nmax,0:Nmax) :: rho,mx,my,E,uu,vv,pp
    real*8,dimension (0:Nmax,1:Nmax) :: ap,am,ap2,am2
    real*8,dimension (1:Nmax,0:Nmax) :: bp,bm,bp2,bm2
    real*8,dimension (0:Nmax,1:Nmax) :: rhoE,mxE,myE,EE,ppE,uuE,vvE
    real*8,dimension (1:Nmax,0:Nmax) :: rhoN,mxN,myN,EN,ppN,uuN,vvN
    real*8,dimension (1:Nmax,1:Nmax) :: rhoW,mxW,myW,EW,rhoS,mxS,myS,ES,ppW,uuW,vvW,ppS,uuS,vvS
    real*8,dimension (0:Nmax,1:Nmax):: F1,F2,F3,F4,alphaI
    real*8,dimension (1:Nmax,0:Nmax):: G1,G2,G3,G4,alphaJ
    real*8,dimension (1:Nmax) :: x,y
    real*8, external :: minmod, dminmod,dminmod4
    real*8 pE,uE,vE,pW,uW,vW,pS,uS,vS,pN,uN,vN,f1E,f2E,f3E,f4E,f1W,f2W,f3W,f4W,g1S,g2S,g3S,g4S,g1N,g2N,g3N,g4N
    character(20) BoundC,Extype,scaleSw,shrinkSw
    real*8 differ1,differ2,differ,prod_a,prod_a1,prod_b,prod_b1,alpha,ratio1,ratio2,ratio3,ratio4,ratiomax,ratiomin
    real*8 df1,df2,df3,df4,dg1,dg2,dg3,dg4,du1,du2,du3,du4,du1eps,du2eps,du3eps,du4eps
    real*8 d1,d2,d3,slx,sly,amax,bmax,dist_a,dist_b,drho,dmx,dmy,dE,rho_star,mx_star,my_star,E_star,wl,wr,u_star,v_star,delta,dn,p_star,dmmax,dp
    real*8 eps_tol,slx1,sly1,slx2,sly2,slx3,sly3,slx4,sly4,bpx,bpy,xslmx,xslmy,yslmx,yslmy,rhoL,rhoR
    real*8 xconserx,xconsery,yconserx,yconsery,xconstx,xconsty,yconstx,yconsty
    integer i,j,k,Nx,Ny,newdtneed,Nx5,Ny5
    real*8 xb,yb,r2,c1,u,v,pi,k0,ppp,eps,eps1,eps2,pmin,ap1,am1,bp1,bm1!,u_star1,v_star1
    real*8 dx,dy,x0,y0,xlength,ylength,dt,dt0,Tfinal,t,theta,gamma,dxinv,dyinv,urhox,urhoy,vrhox,vrhoy,sigmax,sigmay,coEW,coNS
    common/parameter1/dx,dy,x0,y0,xlength,ylength,dt,dt0,Tfinal,t,theta,gamma,dxinv,dyinv
    common/parameter2/Nx,Ny,Nx5,Ny5
    common/parameter3/x,y
    common/parameter4/BoundC,Extype,scaleSw,shrinkSw
    common/parameter5/alphaI,alphaJ
    !===using the linear appriximation or 2nd order approximation=================
    !eps=(gamma-1.d0)*(-0.5d0*(mx(1,1)**2+my(1,1)**2)/rho(1,1)+E(1,1))
    
    do j=0,Ny+1
        do i=0,Nx+1
          if(rho(i,j)>1.d-12)then  
           pp(i,j)=(gamma-1.d0)*(-0.5d0*(mx(i,j)**2+my(i,j)**2)/rho(i,j)+E(i,j))
           uu(i,j)=mx(i,j)/rho(i,j)
           vv(i,j)=my(i,j)/rho(i,j)
           if(pp(i,j)<0.d0)then
            write(*,*)pp(i,j)
           pause 
           end if
          end if 
        end do
    end do
    do j=1,Ny
        do i=1,Nx
            !if((i<Nx5) .or. (j>Ny5))then
            d1=theta*(rho(i+1,j)-rho(i,j))
            d2=0.5d0*(rho(i+1,j)-rho(i-1,j))
            d3=theta*(rho(i,j)-rho(i-1,j))   
            slx=0.5d0*dminmod(d1,d2,d3)
            rhoE(i,j)=rho(i,j)+slx
            rhoW(i,j)=rho(i,j)-slx
            d1=theta*(rho(i,j+1)-rho(i,j))
            d2=0.5d0*(rho(i,j+1)-rho(i,j-1))
            d3=theta*(rho(i,j)-rho(i,j-1))
            sly=0.5d0*dminmod(d1,d2,d3)  
            rhoN(i,j)=rho(i,j)+sly
            rhoS(i,j)=rho(i,j)-sly
            d1=theta*(mx(i+1,j)-mx(i,j))
            d2=0.5d0*(mx(i+1,j)-mx(i-1,j))
            d3=theta*(mx(i,j)-mx(i-1,j))   
            slx=0.5d0*dminmod(d1,d2,d3)
            mxE(i,j)=mx(i,j)+slx
            mxW(i,j)=mx(i,j)-slx
            d1=theta*(mx(i,j+1)-mx(i,j))
            d2=0.5d0*(mx(i,j+1)-mx(i,j-1))
            d3=theta*(mx(i,j)-mx(i,j-1))
            sly=0.5d0*dminmod(d1,d2,d3)  
            mxN(i,j)=mx(i,j)+sly
            mxS(i,j)=mx(i,j)-sly
            d1=theta*(my(i+1,j)-my(i,j))
            d2=0.5d0*(my(i+1,j)-my(i-1,j))
            d3=theta*(my(i,j)-my(i-1,j))   
            slx=0.5d0*dminmod(d1,d2,d3)
            myE(i,j)=my(i,j)+slx
            myW(i,j)=my(i,j)-slx
            d1=theta*(my(i,j+1)-my(i,j))
            d2=0.5d0*(my(i,j+1)-my(i,j-1))
            d3=theta*(my(i,j)-my(i,j-1))
            sly=0.5d0*dminmod(d1,d2,d3)  
            myN(i,j)=my(i,j)+sly
            myS(i,j)=my(i,j)-sly
            d1=theta*(E(i+1,j)-E(i,j))
            d2=0.5d0*(E(i+1,j)-E(i-1,j))
            d3=theta*(E(i,j)-E(i-1,j))   
            slx=0.5d0*dminmod(d1,d2,d3)
            EE(i,j)=E(i,j)+slx
            EW(i,j)=E(i,j)-slx
            d1=theta*(E(i,j+1)-E(i,j))
            d2=0.5d0*(E(i,j+1)-E(i,j-1))
            d3=theta*(E(i,j)-E(i,j-1))
            sly=0.5d0*dminmod(d1,d2,d3)  
            EN(i,j)=E(i,j)+sly
            ES(i,j)=E(i,j)-sly
            !end if
            ppp=(gamma-1.d0)*(E(i,j)-0.5d0*(mx(i,j)**2+my(i,j)**2)/rho(i,j))
            if(ppp<0.d0)then
             write(*,*)ppp
             pause
            end if
            eps=dmin1(1.d-12,ppp)
            pE=(gamma-1.d0)*(EE(i,j)-0.5d0*(mxE(i,j)**2+myE(i,j)**2)/rhoE(i,j))
            pW=(gamma-1.d0)*(EW(i,j)-0.5d0*(mxW(i,j)**2+myW(i,j)**2)/rhoW(i,j))
            pN=(gamma-1.d0)*(EN(i,j)-0.5d0*(mxN(i,j)**2+myN(i,j)**2)/rhoN(i,j))
            pS=(gamma-1.d0)*(ES(i,j)-0.5d0*(mxS(i,j)**2+myS(i,j)**2)/rhoS(i,j))
            pmin=dmin1(pE,pW,pN,pS)
            if (pmin.ge.eps) then
             delta=1.d0
            else   
             delta=(ppp-eps)/(ppp-pmin)   
            endif
            rhoE(i,j)=delta*rhoE(i,j)+(1.d0-delta)*rho(i,j)
            rhoW(i,j)=delta*rhoW(i,j)+(1.d0-delta)*rho(i,j)
            rhoN(i,j)=delta*rhoN(i,j)+(1.d0-delta)*rho(i,j)
            rhoS(i,j)=delta*rhoS(i,j)+(1.d0-delta)*rho(i,j)
            mxE(i,j)=delta*mxE(i,j)+(1.d0-delta)*mx(i,j)
            mxW(i,j)=delta*mxW(i,j)+(1.d0-delta)*mx(i,j)
            mxN(i,j)=delta*mxN(i,j)+(1.d0-delta)*mx(i,j)
            mxS(i,j)=delta*mxS(i,j)+(1.d0-delta)*mx(i,j)
            myE(i,j)=delta*myE(i,j)+(1.d0-delta)*my(i,j)
            myW(i,j)=delta*myW(i,j)+(1.d0-delta)*my(i,j)
            myN(i,j)=delta*myN(i,j)+(1.d0-delta)*my(i,j)
            myS(i,j)=delta*myS(i,j)+(1.d0-delta)*my(i,j)
            EE(i,j)=delta*EE(i,j)+(1.d0-delta)*E(i,j)
            EW(i,j)=delta*EW(i,j)+(1.d0-delta)*E(i,j)
            EN(i,j)=delta*EN(i,j)+(1.d0-delta)*E(i,j)
            ES(i,j)=delta*ES(i,j)+(1.d0-delta)*E(i,j)
        enddo
    enddo
    
   !do j=1,Ny
   !     do i=1,Nx
   !         d1=theta*(rho(i+1,j)-rho(i,j))
   !         d2=0.5d0*(rho(i+1,j)-rho(i-1,j))
   !         d3=theta*(rho(i,j)-rho(i-1,j))   
   !         slx1=0.5d0*dminmod(d1,d2,d3)
   !         rhoE(i,j)=rho(i,j)+slx1
   !         rhoW(i,j)=rho(i,j)-slx1
   !         d1=theta*(rho(i,j+1)-rho(i,j))
   !         d2=0.5d0*(rho(i,j+1)-rho(i,j-1))
   !         d3=theta*(rho(i,j)-rho(i,j-1))
   !         sly1=0.5d0*dminmod(d1,d2,d3)  
   !         rhoN(i,j)=rho(i,j)+sly1
   !         rhoS(i,j)=rho(i,j)-sly1
   !         d1=theta*(uu(i+1,j)-uu(i,j))
   !         d2=0.5d0*(uu(i+1,j)-uu(i-1,j))
   !         d3=theta*(uu(i,j)-uu(i-1,j))   
   !         slx2=0.5d0*dminmod(d1,d2,d3)
   !         uuE(i,j)=uu(i,j)+slx2
   !         uuW(i,j)=uu(i,j)-slx2
   !         d1=theta*(uu(i,j+1)-uu(i,j))
   !         d2=0.5d0*(uu(i,j+1)-uu(i,j-1))
   !         d3=theta*(uu(i,j)-uu(i,j-1))
   !         sly2=0.5d0*dminmod(d1,d2,d3)  
   !         uuN(i,j)=uu(i,j)+sly2
   !         uuS(i,j)=uu(i,j)-sly2
   !         d1=theta*(vv(i+1,j)-vv(i,j))
   !         d2=0.5d0*(vv(i+1,j)-vv(i-1,j))
   !         d3=theta*(vv(i,j)-vv(i-1,j))   
   !         slx3=0.5d0*dminmod(d1,d2,d3)
   !         vvE(i,j)=vv(i,j)+slx3
   !         vvW(i,j)=vv(i,j)-slx3
   !         d1=theta*(vv(i,j+1)-vv(i,j))
   !         d2=0.5d0*(vv(i,j+1)-vv(i,j-1))
   !         d3=theta*(vv(i,j)-vv(i,j-1))
   !         sly3=0.5d0*dminmod(d1,d2,d3)  
   !         vvN(i,j)=vv(i,j)+sly3
   !         vvS(i,j)=vv(i,j)-sly3
   !         d1=theta*(pp(i+1,j)-pp(i,j))
   !         d2=0.5d0*(pp(i+1,j)-pp(i-1,j))
   !         d3=theta*(pp(i,j)-pp(i-1,j))   
   !         slx4=0.5d0*dminmod(d1,d2,d3)
   !         ppE(i,j)=pp(i,j)+slx4
   !         ppW(i,j)=pp(i,j)-slx4
   !         d1=theta*(pp(i,j+1)-pp(i,j))
   !         d2=0.5d0*(pp(i,j+1)-pp(i,j-1))
   !         d3=theta*(pp(i,j)-pp(i,j-1))
   !         sly4=0.5d0*dminmod(d1,d2,d3)  
   !         ppN(i,j)=pp(i,j)+sly4
   !         ppS(i,j)=pp(i,j)-sly4
   !         slx1=2.d0*dxinv*slx1
   !         sly1=2.d0*dyinv*sly1
   !         mxE(i,j)=rhoE(i,j)*uuE(i,j)
   !         mxW(i,j)=2.d0*mx(i,j)-mxE(i,j)
   !         mxN(i,j)=rhoN(i,j)*uuN(i,j)
   !         mxS(i,j)=2.d0*mx(i,j)-mxN(i,j)
   !         myE(i,j)=rhoE(i,j)*vvE(i,j)
   !         myW(i,j)=2.d0*my(i,j)-myE(i,j)
   !         myN(i,j)=rhoN(i,j)*vvN(i,j)
   !         myS(i,j)=2.d0*my(i,j)-myN(i,j)
   !         if((mxW(i,j)/rhoW(i,j)-uuW(i,j))*slx2<0.d0)then
   !          mxW(i,j)=rhoW(i,j)*uuW(i,j)
   !          mxE(i,j)=2.d0*mx(i,j)-mxW(i,j)
   !         end if
   !         if((mxS(i,j)/rhoS(i,j)-uuS(i,j))*sly2<0.d0)then
   !          mxS(i,j)=rhoS(i,j)*uuS(i,j)
   !          mxN(i,j)=2.d0*mx(i,j)-mxS(i,j)
   !         end if
   !         if((myW(i,j)/rhoW(i,j)-vvW(i,j))*slx3<0.d0)then
   !          myW(i,j)=rhoW(i,j)*vvW(i,j)
   !          myE(i,j)=2.d0*my(i,j)-myW(i,j)
   !         end if
   !         if((myS(i,j)/rhoS(i,j)-vvS(i,j))*sly3<0.d0)then
   !          myS(i,j)=rhoS(i,j)*vvS(i,j)
   !          myN(i,j)=2.d0*my(i,j)-myS(i,j)
   !         end if
   !         coEW=dx**2*(gamma-1.d0)*rho(i,j)/(8.d0*rhoE(i,j)*rhoW(i,j))
   !         coNS=dy**2*(gamma-1.d0)*rho(i,j)/(8.d0*rhoN(i,j)*rhoS(i,j))
   !         eps=dmin1(1.d-13,pp(i,j))
   !         bpx=(pp(i,j)-eps)/coEW
   !         bpy=(pp(i,j)-eps)/coNS
   !         urhox=uu(i,j)*slx1
   !         urhoy=uu(i,j)*sly1
   !         vrhox=vv(i,j)*slx1
   !         vrhoy=vv(i,j)*sly1
   !         xslmx=2.d0*dxinv*(mxE(i,j)-mx(i,j))
   !         xslmy=2.d0*dxinv*(myE(i,j)-my(i,j))
   !         yslmx=2.d0*dyinv*(mxN(i,j)-mx(i,j))
   !         yslmy=2.d0*dyinv*(myN(i,j)-my(i,j))
   !         if(bpx>0.d0)then
   !          sigmax=((urhox-xslmx)**2+(vrhox-xslmy)**2)/bpx
   !         else
   !          sigmax=0.d0
   !          mxE(i,j)=mx(i,j)+0.5d0*dx*urhox
   !          myE(i,j)=my(i,j)+0.5d0*dx*vrhox
   !          mxW(i,j)=2.d0*mx(i,j)-mxE(i,j)
   !          myW(i,j)=2.d0*my(i,j)-myE(i,j)
   !          xslmx=2.d0*dxinv*(mxE(i,j)-mx(i,j))
   !          xslmy=2.d0*dxinv*(myE(i,j)-my(i,j))
   !         end if
   !         if(bpy>0.d0)then
   !          sigmay=((urhoy-yslmx)**2+(vrhoy-yslmy)**2)/bpy
   !         else
   !          sigmay=0.d0
   !          mxN(i,j)=mx(i,j)+0.5d0*dy*urhoy
   !          myN(i,j)=my(i,j)+0.5d0*dy*vrhoy
   !          mxS(i,j)=2.d0*mx(i,j)-mxN(i,j)
   !          myS(i,j)=2.d0*my(i,j)-myN(i,j)
   !          yslmx=2.d0*dyinv*(mxN(i,j)-mx(i,j))
   !          yslmy=2.d0*dyinv*(myN(i,j)-my(i,j))
   !         end if    
   !         if(sigmax>1.d0)then
   !          sigmax=dsqrt(1.d0/sigmax)
   !          mxE(i,j)=mx(i,j)+0.5d0*dx*((1.d0-sigmax)*urhox+sigmax*xslmx)
   !          myE(i,j)=my(i,j)+0.5d0*dx*((1.d0-sigmax)*vrhox+sigmax*xslmy)
   !          mxW(i,j)=2.d0*mx(i,j)-mxE(i,j)
   !          myW(i,j)=2.d0*my(i,j)-myE(i,j)
   !          xslmx=2.d0*dxinv*(mxE(i,j)-mx(i,j))
   !          xslmy=2.d0*dxinv*(myE(i,j)-my(i,j))
   !         end if
   !         if(sigmay>1.d0)then
   !          sigmay=dsqrt(1.d0/sigmay)
   !          mxN(i,j)=mx(i,j)+0.5d0*dy*((1.d0-sigmay)*urhoy+sigmay*yslmx)
   !          myN(i,j)=my(i,j)+0.5d0*dy*((1.d0-sigmay)*vrhoy+sigmay*yslmy)
   !          mxS(i,j)=2.d0*mx(i,j)-mxN(i,j)
   !          myS(i,j)=2.d0*my(i,j)-myN(i,j)
   !          yslmx=2.d0*dyinv*(mxN(i,j)-mx(i,j))
   !          yslmy=2.d0*dyinv*(myN(i,j)-my(i,j))
   !         end if
   !         coEW=coEW*((urhox-xslmx)**2+(vrhox-xslmy)**2)
   !         coNS=coNS*((urhoy-yslmx)**2+(vrhoy-yslmy)**2)
   !         if(dabs(slx4)>coEW)then
   !          if(slx4<0.d0)then
   !           EE(i,j)=ppE(i,j)/(gamma-1.d0)+0.5d0*(mxE(i,j)**2+myE(i,j)**2)/rhoE(i,j)
   !           EW(i,j)=2.d0*E(i,j)-EE(i,j)
   !          else
   !           EW(i,j)=ppW(i,j)/(gamma-1.d0)+0.5d0*(mxW(i,j)**2+myW(i,j)**2)/rhoW(i,j)  
   !           EE(i,j)=2.d0*E(i,j)-EW(i,j)
   !          end if 
   !         else 
   !          EE(i,j)=(pp(i,j)-coEW)/(gamma-1.d0)+0.5d0*(mxE(i,j)**2+myE(i,j)**2)/rhoE(i,j)
   !          EW(i,j)=2.d0*E(i,j)-EE(i,j)
   !         end if 
   !         if(dabs(sly4)>coNS)then
   !          if(sly4<0.d0)then
   !           EN(i,j)=ppN(i,j)/(gamma-1.d0)+0.5d0*(mxN(i,j)**2+myN(i,j)**2)/rhoN(i,j)
   !           ES(i,j)=2.d0*E(i,j)-EN(i,j)
   !          else
   !           ES(i,j)=ppS(i,j)/(gamma-1.d0)+0.5d0*(mxS(i,j)**2+myS(i,j)**2)/rhoS(i,j)
   !           EN(i,j)=2.d0*E(i,j)-ES(i,j)
   !          end if 
   !         else
   !          EN(i,j)=(pp(i,j)-coNS)/(gamma-1.d0)+0.5d0*(mxN(i,j)**2+myN(i,j)**2)/rhoN(i,j)
   !          ES(i,j)=2.d0*E(i,j)-EN(i,j)   
   !         end if  
   !     enddo
   ! enddo
    !  do j=1,Ny
    !    do i=1,Nx
    !        d1=theta*(rho(i+1,j)-rho(i,j))
    !        d2=0.5d0*(rho(i+1,j)-rho(i-1,j))
    !        d3=theta*(rho(i,j)-rho(i-1,j))   
    !        slx=0.5d0*dminmod(d1,d2,d3)
    !        rhoE(i,j)=rho(i,j)+slx
    !        rhoW(i,j)=rho(i,j)-slx
    !        d1=theta*(rho(i,j+1)-rho(i,j))
    !        d2=0.5d0*(rho(i,j+1)-rho(i,j-1))
    !        d3=theta*(rho(i,j)-rho(i,j-1))
    !        sly=0.5d0*dminmod(d1,d2,d3)  
    !        rhoN(i,j)=rho(i,j)+sly
    !        rhoS(i,j)=rho(i,j)-sly
    !        d1=theta*(uu(i+1,j)-uu(i,j))
    !        d2=0.5d0*(uu(i+1,j)-uu(i-1,j))
    !        d3=theta*(uu(i,j)-uu(i-1,j))   
    !        slx=0.5d0*dminmod(d1,d2,d3)
    !        uuE(i,j)=uu(i,j)+slx
    !        uuW(i,j)=uu(i,j)-slx
    !        d1=theta*(uu(i,j+1)-uu(i,j))
    !        d2=0.5d0*(uu(i,j+1)-uu(i,j-1))
    !        d3=theta*(uu(i,j)-uu(i,j-1))
    !        sly=0.5d0*dminmod(d1,d2,d3)  
    !        uuN(i,j)=uu(i,j)+sly
    !        uuS(i,j)=uu(i,j)-sly
    !        d1=theta*(vv(i+1,j)-vv(i,j))
    !        d2=0.5d0*(vv(i+1,j)-vv(i-1,j))
    !        d3=theta*(vv(i,j)-vv(i-1,j))   
    !        slx=0.5d0*dminmod(d1,d2,d3)
    !        vvE(i,j)=vv(i,j)+slx
    !        vvW(i,j)=vv(i,j)-slx
    !        d1=theta*(vv(i,j+1)-vv(i,j))
    !        d2=0.5d0*(vv(i,j+1)-vv(i,j-1))
    !        d3=theta*(vv(i,j)-vv(i,j-1))
    !        sly=0.5d0*dminmod(d1,d2,d3)  
    !        vvN(i,j)=vv(i,j)+sly
    !        vvS(i,j)=vv(i,j)-sly
    !        d1=theta*(pp(i+1,j)-pp(i,j))
    !        d2=0.5d0*(pp(i+1,j)-pp(i-1,j))
    !        d3=theta*(pp(i,j)-pp(i-1,j))   
    !        slx=0.5d0*dminmod(d1,d2,d3)
    !        ppE(i,j)=pp(i,j)+slx
    !        ppW(i,j)=pp(i,j)-slx
    !        d1=theta*(pp(i,j+1)-pp(i,j))
    !        d2=0.5d0*(pp(i,j+1)-pp(i,j-1))
    !        d3=theta*(pp(i,j)-pp(i,j-1))
    !        sly=0.5d0*dminmod(d1,d2,d3)  
    !        ppN(i,j)=pp(i,j)+sly
    !        ppS(i,j)=pp(i,j)-sly
    !        mxE(i,j)=rhoE(i,j)*uuE(i,j)
    !        mxW(i,j)=rhoW(i,j)*uuW(i,j)
    !        mxN(i,j)=rhoN(i,j)*uuN(i,j)
    !        mxS(i,j)=rhoS(i,j)*uuS(i,j)
    !        myE(i,j)=rhoE(i,j)*vvE(i,j)
    !        myW(i,j)=rhoW(i,j)*vvW(i,j)
    !        myN(i,j)=rhoN(i,j)*vvN(i,j)
    !        myS(i,j)=rhoS(i,j)*vvS(i,j)
    !        EE(i,j)=ppE(i,j)/(gamma-1.d0)+0.5d0*(mxE(i,j)**2+myE(i,j)**2)/rhoE(i,j)
    !        EW(i,j)=ppW(i,j)/(gamma-1.d0)+0.5d0*(mxW(i,j)**2+myW(i,j)**2)/rhoW(i,j)
    !        EN(i,j)=ppN(i,j)/(gamma-1.d0)+0.5d0*(mxN(i,j)**2+myN(i,j)**2)/rhoN(i,j)
    !        ES(i,j)=ppS(i,j)/(gamma-1.d0)+0.5d0*(mxS(i,j)**2+myS(i,j)**2)/rhoS(i,j)
    !    enddo
    !enddo
    select case (BoundC)
    case('periodic')
        !!===periodic boundary======
        do j=1,Ny
            rhoE(0,j)=rhoE(Nx,j)
            mxE(0,j)=mxE(Nx,j)
            myE(0,j)=myE(Nx,j)
            EE(0,j)=EE(Nx,j)
            rhoW(Nx+1,j)=rhoW(1,j)
            mxW(Nx+1,j)=mxW(1,j)
            myW(Nx+1,j)=myW(1,j)
            EW(Nx+1,j)=EW(1,j)
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoN(i,Ny)
            mxN(i,0)=mxN(i,Ny)
            myN(i,0)=myN(i,Ny)
            EN(i,0)=EN(i,Ny)
            rhoS(i,Ny+1)=rhoS(i,1)
            mxS(i,Ny+1)=mxS(i,1)
            myS(i,Ny+1)=myS(i,1)
            ES(i,Ny+1)=ES(i,1)
        enddo
    case('BCFFSP')
        do j=1,Ny
            rhoE(0,j)=1.4d0
            mxE(0,j)=4.2d0
            myE(0,j)=0.d0
            EE(0,j)=8.8d0
            if(j>Ny5)then
             rhoW(Nx+1,j)=rhoE(Nx,j)
             mxW(Nx+1,j)=mxE(Nx,j)
             myW(Nx+1,j)=myE(Nx,j)
             EW(Nx+1,j)=EE(Nx,j)
            else 
             rhoW(Nx5,j)=rhoE(Nx5-1,j)
             mxW(Nx5,j)=-mxE(Nx5-1,j)
             myW(Nx5,j)=myE(Nx5-1,j)
             EW(Nx5,j)=EE(Nx5-1,j)
            end if
        end do
        do i=1,Nx
            rhoS(i,Ny+1)=rhoN(i,Ny)
            mxS(i,Ny+1)=mxN(i,Ny)
            myS(i,Ny+1)=-myN(i,Ny)
            ES(i,Ny+1)=EN(i,Ny)
            if(i<Nx5)then
              rhoN(i,0)=rhoS(i,1)
              mxN(i,0)=mxS(i,1)
              myN(i,0)=-myS(i,1)
              EN(i,0)=ES(i,1) 
            else
              rhoN(i,Ny5+1)=rhoS(i,Ny5+2)
              mxN(i,Ny5+1)=mxS(i,Ny5+2)
              myN(i,Ny5+1)=-myS(i,Ny5+2)
              EN(i,Ny5+1)=ES(i,Ny5+2)
            end if    
        end do
        !write(*,*)rhoS(Nx5,Ny5+2)
        !stop
    case('reflectfree')
        !! ==reflecting boundary(solid wall)+free====
        do j=1,Ny
            rhoE(0,j)=rhoW(1,j)
            mxE(0,j)=-mxW(1,j)
            myE(0,j)=myW(1,j)
            EE(0,j)=EW(1,j)
            rhoW(Nx+1,j)=rhoE(Nx,j)
            mxW(Nx+1,j)=mxE(Nx,j)
            myW(Nx+1,j)=myE(Nx,j)
            EW(Nx+1,j)=EE(Nx,j) 
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoS(i,1)
            mxN(i,0)=mxS(i,1)
            myN(i,0)=-myS(i,1)
            EN(i,0)=ES(i,1)
            rhoS(i,Ny+1)=rhoN(i,Ny)
            mxS(i,Ny+1)=mxN(i,Ny)
            myS(i,Ny+1)=myN(i,Ny)
            ES(i,Ny+1)=EN(i,Ny)
        enddo
    case('free')
        ! !==free====
        do j=1,Ny
            rhoE(0,j)=rhoW(1,j)
            mxE(0,j)=mxW(1,j)
            myE(0,j)=myW(1,j)
            EE(0,j)=EW(1,j)
            rhoW(Nx+1,j)=rhoE(Nx,j)
            mxW(Nx+1,j)=mxE(Nx,j)
            myW(Nx+1,j)=myE(Nx,j)
            EW(Nx+1,j)=EE(Nx,j)
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoS(i,1)
            mxN(i,0)=mxS(i,1)
            myN(i,0)=myS(i,1)
            EN(i,0)=ES(i,1)
            rhoS(i,Ny+1)=rhoN(i,Ny)
            mxS(i,Ny+1)=mxN(i,Ny)
            myS(i,Ny+1)=myN(i,Ny)
            ES(i,Ny+1)=EN(i,Ny)
        enddo  
    case('solidwall')
        !! ==reflecting boundary(solid wall)====
        do j=1,Ny
            rhoE(0,j)=rhoW(1,j)
            mxE(0,j)=-mxW(1,j)
            myE(0,j)=myW(1,j)
            EE(0,j)=EW(1,j)
            rhoW(Nx+1,j)=rhoE(Nx,j)
            mxW(Nx+1,j)=-mxE(Nx,j)
            myW(Nx+1,j)=myE(Nx,j)
            EW(Nx+1,j)=EE(Nx,j)
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoS(i,1)
            mxN(i,0)=mxS(i,1)
            myN(i,0)=-myS(i,1)
            EN(i,0)=ES(i,1)
            rhoS(i,Ny+1)=rhoN(i,Ny)
            mxS(i,Ny+1)=mxN(i,Ny)
            myS(i,Ny+1)=-myN(i,Ny)
            ES(i,Ny+1)=EN(i,Ny)
        enddo
     case('fixedleft1')
        do j=1,Ny
          if(dabs(y(j))<0.05d0)then  
            rhoE(0,j)=5.d0
            mxE(0,j)=rhoE(0,j)*30.d0
            myE(0,j)=0.d0
            EE(0,j)=0.4127d0/(gamma-1.d0)+2250.d0
          else  
            rhoE(0,j)=rhoW(1,j)
            mxE(0,j)=mxW(1,j)
            myE(0,j)=myW(1,j)
            EE(0,j)=EW(1,j)
          end if  
            rhoW(Nx+1,j)=rhoE(Nx,j)
            mxW(Nx+1,j)=mxE(Nx,j)
            myW(Nx+1,j)=myE(Nx,j)
            EW(Nx+1,j)=EE(Nx,j)
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoS(i,1)
            mxN(i,0)=mxS(i,1)
            myN(i,0)=myS(i,1)
            EN(i,0)=ES(i,1)
            rhoS(i,Ny+1)=rhoN(i,Ny)
            mxS(i,Ny+1)=mxN(i,Ny)
            myS(i,Ny+1)=myN(i,Ny)
            ES(i,Ny+1)=EN(i,Ny)
        enddo  
     case('fixedleft2')
        do j=1,Ny
          if(dabs(y(j))<0.05d0)then  
            rhoE(0,j)=5.d0
            mxE(0,j)=rhoE(0,j)*800.d0
            myE(0,j)=0.d0
            EE(0,j)=0.4127d0/(gamma-1.d0)+1.6d6
          else  
            rhoE(0,j)=rhoW(1,j)
            mxE(0,j)=mxW(1,j)
            myE(0,j)=myW(1,j)
            EE(0,j)=EW(1,j)
          end if  
            rhoW(Nx+1,j)=rhoE(Nx,j)
            mxW(Nx+1,j)=mxE(Nx,j)
            myW(Nx+1,j)=myE(Nx,j)
            EW(Nx+1,j)=EE(Nx,j)
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoS(i,1)
            mxN(i,0)=mxS(i,1)
            myN(i,0)=myS(i,1)
            EN(i,0)=ES(i,1)
            rhoS(i,Ny+1)=rhoN(i,Ny)
            mxS(i,Ny+1)=mxN(i,Ny)
            myS(i,Ny+1)=myN(i,Ny)
            ES(i,Ny+1)=EN(i,Ny)
        enddo  
     case('fixedleft3')
        do j=1,Ny 
            rhoE(0,j)=5.26829268d0
            mxE(0,j)=rhoE(0,j)*4.86111111d0
            myE(0,j)=0.d0
            EE(0,j)= 29.88095238d0/(gamma-1.d0)+0.5d0*mxE(0,j)**2/rhoE(0,j)
            rhoW(Nx+1,j)=rhoE(Nx,j)
            mxW(Nx+1,j)=mxE(Nx,j)
            myW(Nx+1,j)=myE(Nx,j)
            EW(Nx+1,j)=EE(Nx,j)
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoS(i,1)
            mxN(i,0)=mxS(i,1)
            myN(i,0)=myS(i,1)
            EN(i,0)=ES(i,1)
            rhoS(i,Ny+1)=rhoN(i,Ny)
            mxS(i,Ny+1)=mxN(i,Ny)
            myS(i,Ny+1)=myN(i,Ny)
            ES(i,Ny+1)=EN(i,Ny)
        enddo  
     case('parsolidwall')
        !! ==reflecting boundary(solid wall)====
        k0=(4.125**dsqrt(3.d0)+dsqrt(gamma*116.5d0/8.d0))*0.5d0*dsqrt(3.d0)
        !k0=4.125**dsqrt(3.d0)*10.d0/6.6d0
        !k0=10.d0
        do j=1,Ny
            rhoE(0,j)=8.d0
            mxE(0,j)=mxW(1,j)
            myE(0,j)=myW(1,j)
            EE(0,j)=EW(1,j)
            rhoW(Nx+1,j)=rhoE(Nx,j)
            mxW(Nx+1,j)=mxE(Nx,j)
            myW(Nx+1,j)=myE(Nx,j)
            EW(Nx+1,j)=EE(Nx,j)  
        enddo
        do i=1,Nx
            if(x(i)>1.d0/6.d0)then
             rhoN(i,0)=rhoS(i,1)
             mxN(i,0)=mxS(i,1)
             myN(i,0)=-myS(i,1)
             EN(i,0)=ES(i,1)
            else
             rhoN(i,0)=8.d0
             mxN(i,0)=mxS(i,1)
             myN(i,0)=myS(i,1)
             EN(i,0)=ES(i,1)
            end if 
            if(x(i)<(1.d0+2.d0*dsqrt(3.d0))/6.d0+k0*t)then
             rhoS(i,Ny+1)=8.d0
             mxS(i,Ny+1)=33.d0*dsqrt(3.d0)
             myS(i,Ny+1)=-33.d0
             ES(i,Ny+1)=116.5d0/(gamma-1.d0)+0.5d0*(mxN(i,Ny)**2+myN(i,Ny)**2)/rhoN(i,Ny)
            else
             rhoS(i,Ny+1)=1.4d0
             mxS(i,Ny+1)=0.d0
             myS(i,Ny+1)=0.d0
             ES(i,Ny+1)=1.d0/(gamma-1.d0)         
            end if
        enddo
    case('fixed')
        pi=4.d0*datan(1.d0)
         do i=1,Nx
            xb=x0+(i-0.5d0)*dx-5.d0-t
            yb=y0-5.d0-t
            r2=xb**2+yb**2
            rhoN(i,0)=(1.d0-(gamma-1.d0)*25.d0/(8.d0*gamma*pi**2)*dexp(1.d0-r2))**(1.d0/(gamma-1.d0))
            c1=2.5d0/pi*dexp(0.5d0*(1.d0-r2))
            u=1.d0-c1*yb
            v=1.d0+c1*xb
            mxN(i,0)=rhoN(i,0)*u
            myN(i,0)=rhoN(i,0)*v
            EN(i,0)=(rhoN(i,0))**gamma/(gamma-1.d0)+0.5d0*(mxN(i,0)**2+myN(i,0)**2)/rhoN(i,0)
            xb=x0+(i-0.5d0)*dx-5.d0-t
            yb=y0+ylength-5.d0-t
            r2=xb**2+yb**2
            rhoS(i,Ny+1)=(1.d0-(gamma-1.d0)*25.d0/(8.d0*gamma*pi**2)*dexp(1.d0-r2))**(1.d0/(gamma-1.d0))
            c1=2.5d0/pi*dexp(0.5d0*(1.d0-r2))
            u=1.d0-c1*yb
            v=1.d0+c1*xb
            mxS(i,Ny+1)=rhoS(i,Ny+1)*u
            myS(i,Ny+1)=rhoS(i,Ny+1)*v
            ES(i,Ny+1)=(rhoS(i,Ny+1))**gamma/(gamma-1.d0)+0.5d0*(mxS(i,Ny+1)**2+myS(i,Ny+1)**2)/rhoS(i,Ny+1)
         enddo  
         do j=1,Ny
            xb=x0-5.d0-t
            yb=y0+(j-0.5d0)*dy-5.d0-t
            r2=xb**2+yb**2
            rhoE(0,j)=(1.d0-(gamma-1.d0)*25.d0/(8.d0*gamma*pi**2)*dexp(1.d0-r2))**(1.d0/(gamma-1.d0))
            c1=2.5d0/pi*dexp(0.5d0*(1.d0-r2))
            u=1.d0-c1*yb
            v=1.d0+c1*xb
            mxE(0,j)=rhoE(0,j)*u
            myE(0,j)=rhoE(0,j)*v
            EE(0,j)=(rhoE(0,j))**gamma/(gamma-1.d0)+0.5d0*(mxE(0,j)**2+myE(0,j)**2)/rhoE(0,j)
            xb=x0+xlength-5.d0-t
            yb=y0+(j-0.5d0)*dy-5.d0-t
            r2=xb**2+yb**2
            rhoW(Nx+1,j)=(1.d0-(gamma-1.d0)*25.d0/(8.d0*gamma*pi**2)*dexp(1.d0-r2))**(1.d0/(gamma-1.d0))
            c1=2.5d0/pi*dexp(0.5d0*(1.d0-r2))
            u=1.d0-c1*yb
            v=1.d0+c1*xb
            mxW(Nx+1,j)=rhoW(Nx+1,j)*u
            myW(Nx+1,j)=rhoW(Nx+1,j)*v
            EW(Nx+1,j)=(rhoW(Nx+1,j))**gamma/(gamma-1.d0)+0.5d0*(mxW(Nx+1,j)**2+myW(Nx+1,j)**2)/rhoW(Nx+1,j)
         enddo  
    end select   

    eps_tol=1.0d-12
    do j=1,Ny 
        do i=0,Nx
            uE=mxE(i,j)/rhoE(i,j)
            vE=myE(i,j)/rhoE(i,j)
            pE=(gamma-1.d0)*(-0.5d0*rhoE(i,j)*(uE**2+vE**2)+EE(i,j))
            if(pE<0.d0)then
             pE=0.d0
            end if
            uW=mxW(i+1,j)/rhoW(i+1,j)
            vW=myW(i+1,j)/rhoW(i+1,j)
            pW=(gamma-1.d0)*(-0.5d0*rhoW(i+1,j)*(uW**2+vW**2)+EW(i+1,j))
            if(pW<0.d0)then
             pW=0.d0
            end if
            ap(i,j)=dmax1(uW+dsqrt(gamma*pW/rhoW(i+1,j)),uE+dsqrt(gamma*pE/rhoE(i,j)),0.d0)            
            am(i,j)=dmin1(uW-dsqrt(gamma*pW/rhoW(i+1,j)),uE-dsqrt(gamma*pE/rhoE(i,j)),0.d0)
            f1E=mxE(i,j)
            f2E=mxE(i,j)*uE+pE
            f3E=mxE(i,j)*vE
            f4E=uE*(EE(i,j)+pE)
            f1W=mxW(i+1,j)
            f2W=mxW(i+1,j)*uW+pW
            f3W=vW*mxW(i+1,j)
            f4W=uW*(EW(i+1,j)+pW)
            dist_a=ap(i,j)-am(i,j)
            prod_a=ap(i,j)*am(i,j)
            if (dist_a>eps_tol) then
                rho_star=((ap(i,j)*rhoW(i+1,j)-am(i,j)*rhoE(i,j))-(f1W-f1E))/dist_a
                mx_star=((ap(i,j)*mxW(i+1,j)-am(i,j)*mxE(i,j))-(f2W-f2E))/dist_a
                my_star=((ap(i,j)*myW(i+1,j)-am(i,j)*myE(i,j))-(f3W-f3E))/dist_a
                E_star=((ap(i,j)*EW(i+1,j)-am(i,j)*EE(i,j))-(f4W-f4E))/dist_a
                p_star=(gamma-1.d0)*(E_star-0.5d0*(mx_star**2+my_star**2)/rho_star)
                dp=dmax1(p_star-1.d-13,0.d0)
                u_star=mx_star/rho_star
                v_star=my_star/rho_star
                ap1=ap(i,j)-u_star
                am1=am(i,j)-u_star
                wl=-am1*(rhoE(i,j)-rho_star)
                wr=ap1*(rho_star-rhoW(i+1,j))
                drho=minmod(wl,wr)
                wl=-am1*(myE(i,j)-my_star)
                wr=ap1*(my_star-myW(i+1,j))
                dmy=minmod(wl,wr)
                rhoL=am1*rho_star-drho
                rhoR=ap1*rho_star-drho 
                dmmax=dsqrt(-2.d0*rhoL*rhoR*dp/(rho_star*(gamma-1.d0)))
                if(dabs(dmy-v_star*drho)>dmmax)then
                 dmy=v_star*drho+sign(dmmax,dmy-v_star*drho)
                end if 
                dE=0.5d0*((my_star-dmy/ap1)**2/(rho_star-drho/ap1)-(my_star-dmy/am1)**2/(rho_star-drho/am1))/dist_a
                if(u_star>0.d0)then        
                 drho=drho*am(i,j)/am1
                 dmy=dmy*am(i,j)/am1
                 dE=dE*am(i,j)*ap1
                else 
                 drho=drho*ap(i,j)/ap1
                 dmy=dmy*ap(i,j)/ap1
                 dE=dE*ap(i,j)*am1
                end if
                if(dabs(u_star)>=dmax1(ap(i,j),-am(i,j)))then
                 drho=0.d0
                 dmy=0.d0
                 dE=0.d0
                 pause
                end if
                dmx=u_star*drho
                dE=dE+0.5d0*u_star**2*drho 
                !drho=0.d0
                ! dmx=0.d0
                ! dmy=0.d0
                ! dE=0.d0
                F1(i,j)=((ap(i,j)*f1E-am(i,j)*f1W)+prod_a*((rhoW(i+1,j)-rhoE(i,j))))/dist_a-drho
                F2(i,j)=((ap(i,j)*f2E-am(i,j)*f2W)+prod_a*((mxW(i+1,j)-mxE(i,j))))/dist_a-dmx
                F3(i,j)=((ap(i,j)*f3E-am(i,j)*f3W)+prod_a*((myW(i+1,j)-myE(i,j))))/dist_a-dmy
                F4(i,j)=((ap(i,j)*f4E-am(i,j)*f4W)+prod_a*((EW(i+1,j)-EE(i,j))))/dist_a-dE
                !if(dabs(prod_a)>1.d-12)then 
                ! F4(i,j)=((ap(i,j)*f4E-am(i,j)*f4W)+(prod_a*(EW(i+1,j)-EE(i,j))-0.5d0*prod_a1*((my_star-dmy/ap(i,j))**2/(rho_star-drho/ap(i,j))-(my_star-dmy/am(i,j))**2/(rho_star-drho/am(i,j)))))/dist_a-dE
                !else
                ! F4(i,j)=(ap(i,j)*f4E-am(i,j)*f4W)/dist_a                     
                !end if
            else
                F1(i,j)=0.5d0*(f1E+f1W)
                F2(i,j)=0.5d0*(f2E+f2W) 
                F3(i,j)=0.5d0*(f3E+f3W)
                F4(i,j)=0.5d0*(f4E+f4W)
            endif 
        enddo
    enddo    
    do i=1,Nx 
        do j=0,Ny
            uN=mxN(i,j)/rhoN(i,j)
            vN=myN(i,j)/rhoN(i,j)  
            pN=(gamma-1.d0)*(-0.5d0*rhoN(i,j)*(uN**2+vN**2)+EN(i,j))
            if(pN<0.d0)then
             pN=0.d0
            end if
            uS=mxS(i,j+1)/rhoS(i,j+1)
            vS=myS(i,j+1)/rhoS(i,j+1)   
            pS=(gamma-1.d0)*(-0.5d0*rhoS(i,j+1)*(uS**2+vS**2)+ES(i,j+1))
            if(pS<0.d0)then
             pS=0.d0
            end if
            bp(i,j)=dmax1(vS+dsqrt(gamma*pS/rhoS(i,j+1)),vN+dsqrt(gamma*pN/rhoN(i,j)),0.d0)
            bm(i,j)=dmin1(vS-dsqrt(gamma*pS/rhoS(i,j+1)),vN-dsqrt(gamma*pN/rhoN(i,j)),0.d0)
            g1N=myN(i,j)
            g2N=myN(i,j)*uN
            g3N=myN(i,j)*vN+pN
            g4N=vN*(EN(i,j)+pN)
            g1S=myS(i,j+1)
            g2S=myS(i,j+1)*uS
            g3S=myS(i,j+1)*vS+pS
            g4S=vS*(ES(i,j+1)+pS)
            dist_b=bp(i,j)-bm(i,j)
            prod_b=bp(i,j)*bm(i,j)
            if (dist_b>eps_tol) then
                rho_star=((bp(i,j)*rhoS(i,j+1)-bm(i,j)*rhoN(i,j))-(g1S-g1N))/dist_b
                mx_star=((bp(i,j)*mxS(i,j+1)-bm(i,j)*mxN(i,j))-(g2S-g2N))/dist_b
                my_star=((bp(i,j)*myS(i,j+1)-bm(i,j)*myN(i,j))-(g3S-g3N))/dist_b
                E_star=((bp(i,j)*ES(i,j+1)-bm(i,j)*EN(i,j))-(g4S-g4N))/dist_b
                p_star=(gamma-1.d0)*(E_star-0.5d0*(mx_star**2+my_star**2)/rho_star)
                dp=dmax1(p_star-1.d-13,0.d0)
                u_star=mx_star/rho_star
                v_star=my_star/rho_star
                bp1=bp(i,j)-v_star
                bm1=bm(i,j)-v_star
                wl=-bm1*(rhoN(i,j)-rho_star)
                wr=bp1*(rho_star-rhoS(i,j+1))
                drho=minmod(wl,wr)
                wl=-bm1*(mxN(i,j)-mx_star)
                wr=bp1*(mx_star-mxS(i,j+1))
                dmx=minmod(wl,wr)
                rhoL=bm1*rho_star-drho
                rhoR=bp1*rho_star-drho 
                dmmax=dsqrt(-2.d0*rhoL*rhoR*dp/(rho_star*(gamma-1.d0)))
                if(dabs(dmx-u_star*drho)>dmmax)then
                 dmx=u_star*drho+sign(dmmax,dmx-u_star*drho)
                end if   
                dE=0.5d0*((mx_star-dmx/bp1)**2/(rho_star-drho/bp1)-(mx_star-dmx/bm1)**2/(rho_star-drho/bm1))/dist_b
                if(v_star>0.d0)then        
                 drho=drho*bm(i,j)/bm1
                 dmx=dmx*bm(i,j)/bm1
                 dE=dE*bm(i,j)*bp1
                else 
                 drho=drho*bp(i,j)/bp1
                 dmx=dmx*bp(i,j)/bp1
                 dE=dE*bp(i,j)*bm1
                end if
                if(dabs(v_star)>=dmax1(bp(i,j),-bm(i,j)))then
                 drho=0.d0
                 dmx=0.d0
                 dE=0.d0
                 pause
                end if
                dmy=v_star*drho 
                dE=dE+0.5d0*v_star**2*drho
                G1(i,j)=((bp(i,j)*g1N-bm(i,j)*g1S)+prod_b*((rhoS(i,j+1)-rhoN(i,j))))/dist_b-drho
                G2(i,j)=((bp(i,j)*g2N-bm(i,j)*g2S)+prod_b*((mxS(i,j+1)-mxN(i,j))))/dist_b-dmx
                G3(i,j)=((bp(i,j)*g3N-bm(i,j)*g3S)+prod_b*((myS(i,j+1)-myN(i,j))))/dist_b-dmy
                G4(i,j)=((bp(i,j)*g4N-bm(i,j)*g4S)+prod_b*((ES(i,j+1)-EN(i,j))))/dist_b-dE
                !if(dabs(prod_b)>1.d-12)then 
                ! G4(i,j)=((bp(i,j)*g4N-bm(i,j)*g4S)+(prod_b*(ES(i,j+1)-EN(i,j))-0.5d0*prod_b1*((mx_star-dmx/bp(i,j))**2/(rho_star-drho/bp(i,j))-(mx_star-dmx/bm(i,j))**2/(rho_star-drho/bm(i,j)))))/dist_b-dE
                !else
                ! G4(i,j)=(bp(i,j)*g4N-bm(i,j)*g4S)/dist_b                       
                !end if
            else
                G1(i,j)=0.5d0*(g1N+g1S)
                G2(i,j)=0.5d0*(g2N+g2S) 
                G3(i,j)=0.5d0*(g3N+g3S)
                G4(i,j)=0.5d0*(g4N+g4S)
            endif 
            !end if
201        enddo
    enddo
    
    !CFL
    amax=0.d0
    do j=1,Ny 
        do i=0,Nx
            amax=dmax1(amax,ap(i,j),-am(i,j))
        enddo
    enddo
    bmax=0.d0
    do j=0,Ny
        do i=1,Nx 
            bmax=dmax1(bmax,bp(i,j),-bm(i,j)) 
        enddo
    enddo
        !amax=0.d0
        !do j=1,Ny 
        !    do i=1,Nx
        !        amax=dmax1(amax,ap(i-1,j)-am(i,j))
        !    enddo
        !enddo
        !bmax=0.d0
        !do j=1,Ny
        !    do i=1,Nx 
        !        bmax=dmax1(bmax,bp(i,j-1)-bm(i,j)) 
        !    enddo
        !enddo
        dt=0.475d0*dmin1(dx/amax,dy/bmax)
         !dt0=0.5d0/(amax/dx+bmax/dy)    
        if (t+dt>Tfinal)  then
            dt=Tfinal-t
        endif
    !if (newdtneed==1) then     
    !   dt=dt0 
    !endif
    
    return
    end subroutine