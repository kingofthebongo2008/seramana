      program panel
      common t_lift(10000)
c
c	smith-hess panel method for single
c	element lifting airfoil in 2-d
c	incompressible flow
c
      dimension z(10000)
      pi=3.1415926585
c
      call indata
      call setup
c

      alpha=5
      twist=5
      c_tip=3
      c_root=1
      s_span=5
c
	q_dyn=1.225*0.5*5**2
      alpha=alpha+twist
	d_s=s_span/10000.
      d_twist=twist/10000.
	d_chord=(c_root-c_tip)/10000.  
	
      tlift=0.
      


	  DO i=1,10000,1

c
      cosalf=cos(alpha*pi/180.)
      sinalf=sin(alpha*pi/180.)
      call cofish(sinalf,cosalf)
      call gauss(1)
      call veldis(sinalf,cosalf)
      call fandm(sinalf,cosalf,cl)

	  alpha=alpha-d_twist
	  chord=c_root-d_chord*i	  
	  area=d_s*chord

	  t_lift(i)=cl*q_dyn*area

	
	  ENDDO

      DO i=1,10000,1
        tlift=tlift+t_lift(i)
        ENDDO

      
c	print*, "Lift coef:",2.*t_lift/(q_dyn*(c_root+c_tip)*0.5*s_span)
	  print*, "Total Lift :",2.*tlift, "  [N]"
c      go to 100
c
  200 stop
      end
c
c**********************************************************
c
      subroutine setup
c
      common /bod/ nlower,nupper,nodtot,x(100),y(100),
     2		   costhe(100),sinthe(100)
      common /num/ pi,pi2inv
c
      pi=3.1415926585
      pi2inv=.5/pi
c
c  set coordinates of nodes on body surface
c
      write(6,1000)
 1000 format(///' body shape'//4x,'x',9x,'y'/)
      npoints=nlower
      sign=-1.0
      nstart=0
      do 110 nsurf=1,2
      do 100 n=1,npoints
      fract=float(n-1)/float(npoints)
      z=.5*(1.-cos(pi*fract))
      i=nstart+n
      call body(z,sign,x(i),y(i))
      write(6,1010)x(i),y(i)
 1010 format(f8.4,f10.4)
  100 continue
      npoints=nupper
      sign=1.0
      nstart=nlower
  110 continue
      nodtot=nupper+nlower
      x(nodtot+1)=x(1)
      y(nodtot+1)=y(1)
c
c  set slopes of panels
c
      do 200 i=1,nodtot
      dx=x(i+1)-x(i)
      dy=y(i+1)-y(i)
      dist=sqrt(dx*dx+dy*dy)
      sinthe(i)=dy/dist
      costhe(i)=dx/dist
  200 continue
c
      return
      end
c
c**********************************************************
c
      subroutine body(z,sign,x,y)
c
c  return coordinates of point on body surface
c
c     z = node spacing parameter
c     x,y = cartesian coordinates
c     sign = +1. for upper surface, -1. for lower surface
c
      common /par/ naca,tau,epsmax,ptmax
c
      if (sign.lt.0.) z=1.-z
      call naca45(z,thick,camber,beta)
      x=z-sign*thick*sin(beta)
      y=camber+sign*thick*cos(beta)
c
      return
      end
c
c**********************************************************
c
      subroutine cofish(sinalf,cosalf)
c
c  set coefficients of linear system
c
      common /bod/ nlower,nupper,nodtot,x(100),y(100),
     2		   costhe(100),sinthe(100)
      common /cof/ a(101,111),kutta
      common /num/ pi,pi2inv
c
      kutta=nodtot+1
c
c  initialize coefficients
c
      do 90 j=1,kutta
   90 a(kutta,j)=0.0
c
c  set vn=0. at midpoint of ith panel
c
      do 120 i=1,nodtot
      xmid=.5*(x(i)+x(i+1))
      ymid=.5*(y(i)+y(i+1))
      a(i,kutta)=0.0
c
c  find contribution of jth panel
c
      do 110 j=1,nodtot
      flog=0.0
      ftan=pi
      if (j.eq.i) go to 100
      dxj=xmid-x(j)
      dxjp=xmid-x(j+1)
      dyj=ymid-y(j)
      dyjp=ymid-y(j+1)
      flog=.5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
      ftan=atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
  100 ctimtj=costhe(i)*costhe(j)+sinthe(i)*sinthe(j)
      stimtj=sinthe(i)*costhe(j)-sinthe(j)*costhe(i)
      a(i,j)=pi2inv*(ftan*ctimtj+flog*stimtj)
      b=pi2inv*(flog*ctimtj-ftan*stimtj)
      a(i,kutta)=a(i,kutta)+b
      if ((i.gt.1).and.(i.lt.nodtot)) go to 110
c
c  if ith panel touches trailing edge, add contribution
c    to kutta condition
c
      a(kutta,j)=a(kutta,j)-b
      a(kutta,kutta)=a(kutta,kutta)+a(i,j)
  110 continue
c
c  fill in known sides
c
      a(i,kutta+1)=sinthe(i)*cosalf-costhe(i)*sinalf
  120 continue
      a(kutta,kutta+1)=-(costhe(1)+costhe(nodtot))*cosalf
     2		       -(sinthe(1)+sinthe(nodtot))*sinalf
c
      return
      end
c
c**********************************************************
c
      subroutine veldis(sinalf,cosalf)
c
c  compute and print out pressure distribution
c
      common /bod/ nlower,nupper,nodtot,x(100),y(100),
     2		   costhe(100),sinthe(100)
      common /cof/ a(101,111),kutta
      common /cpd/ cp(100)
      common /num/ pi,pi2inv
      dimension q(150)
c
c      write(6,1000)
 1000 format(///' pressure distribution'//4x,'x',8x,'cp'/)
c
c  retrieve solution from a-matrix
c
      do 50 i=1,nodtot
   50 q(i)=a(i,kutta+1)
      gamma=a(kutta,kutta+1)
c
c  find vt and cp at midpoint of ith panel
c
      do 130 i=1,nodtot
      xmid=.5*(x(i)+x(i+1))
      ymid=.5*(y(i)+y(i+1))
      vtang=cosalf*costhe(i)+sinalf*sinthe(i)
c
c  add contributions of jth panel
c
      do 120 j=1,nodtot
      flog=0.0
      ftan=pi
      if (j.eq.i) go to 100
      dxj=xmid-x(j)
      dxjp=xmid-x(j+1)
      dyj=ymid-y(j)
      dyjp=ymid-y(j+1)
      flog=.5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
      ftan=atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
  100 ctimtj=costhe(i)*costhe(j)+sinthe(i)*sinthe(j)
      stimtj=sinthe(i)*costhe(j)-sinthe(j)*costhe(i)
      aa=pi2inv*(ftan*ctimtj+flog*stimtj)
      b=pi2inv*(flog*ctimtj-ftan*stimtj)
      vtang=vtang-b*q(j)+gamma*aa
  120 continue
      cp(i)=1.-vtang*vtang
c      write(6,1010)xmid,cp(i)
 1010 format(f8.4,f10.4)
  130 continue
c
      return
      end
c
c**********************************************************
c
      subroutine fandm(sinalf,cosalf,cl)
c
c  compute and print out cd,cl,cmle
c
      common /bod/ nlower,nupper,nodtot,x(100),y(100),
     2		   costhe(100),sinthe(100)
      common /cpd/ cp(100) 
c
      cfx=0.0
      cfy=0.0
      cm=0.0
c
      do 100 i=1,nodtot
      xmid=.5*(x(i)+x(i+1))
      ymid=.5*(y(i)+y(i+1))
      dx=x(i+1)-x(i)
      dy=y(i+1)-y(i)
      cfx=cfx+cp(i)*dy
      cfy=cfy-cp(i)*dx
      cm=cm+cp(i)*(dx*xmid+dy*ymid)
  100 continue
      cd=cfx*cosalf+cfy*sinalf
      cl=cfy*cosalf-cfx*sinalf
      write(6,1000)cd,cl,cm
 1000 format(////'    cd =',f8.5,'    cl =',f8.5,'    cm =',f8.5)
c
      return
      end
c
c**********************************************************
c
      subroutine gauss(nrhs)
c
c  solution of linear algebraic system by
c  gaussian elimination with partial pivoting
c
c	  [a] = coefficient matrix
c	  neqns = number of equations
c	  nrhs = number of right-hand sides
c
c	  right-hand sides and solutions stored in
c	  columns neqns+1 thru neqns+nrhs of a
c
      common /cof/ a(101,111),neqns
c
      np=neqns+1
      ntot=neqns+nrhs
c
c  gauss reduction
c
      do 150 i=2,neqns
c
c  search for largest entry in (i-1)th column
c  on or below main diagonal
c
      im=i-1
      imax=im
      amax=abs(a(im,im))
      do 110 j=i,neqns
      if (amax.ge.abs(a(j,im))) go to 110
      imax=j
      amax=abs(a(j,im))
  110 continue
c
c  switch (i-1)th and imaxth equations
c
      if (imax.ne.im) go to 140
      do 130 j=im,ntot
      temp=a(im,j)
      a(im,j)=a(imax,j)
      a(imax,j)=temp
  130 continue
c
c  eliminate (i-1)th unknown from
c  ith thru neqnsth equations
c
  140 do 150 j=i,neqns
      r=a(j,im)/a(im,im)
      do 150 k=i,ntot
  150 a(j,k)=a(j,k)-r*a(im,k)
c
c  back substitution
c
      do 220 k=np,ntot
      a(neqns,k)=a(neqns,k)/a(neqns,neqns)
      do 210 l=2,neqns
      i=neqns+1-l
      ip=i+1
      do 200 j=ip,neqns
  200 a(i,k)=a(i,k)-a(i,j)*a(j,k)
  210 a(i,k)=a(i,k)/a(i,i)
  220 continue
c
      return
      end
c
c**********************************************************
c
      subroutine indata
c
c  set parameters of body shape, flow
c  situation, and node distribution
c
c  user must input:
c  nlower = number of nodes on lower surface
c  nupper = number of nodes on upper surface
c  plus data on body
c
      common /bod/ nlower,nupper,nodtot,x(100),y(100),
     2		   costhe(100),sinthe(100)
      common /par/ naca,tau,epsmax,ptmax
c
c      write(6,*)'input nlower,nupper'
c      read(5,*)nlower,nupper
c      write(6,*)'input naca number'
c      read(5,*)naca
c
      nlower=50
      nupper=50
      naca=2412
      ieps=naca/1000
      iptmax=naca/100-10*ieps
      itau=naca-1000*ieps-100*iptmax
      epsmax=ieps*0.01
      ptmax=iptmax*0.1
      tau=itau*0.01
      if (ieps.lt.10) return
      ptmax=0.2025
      epsmax=2.6595*ptmax**3
c
      return
      end
c
c**********************************************************
c
      subroutine naca45(z,thick,camber,beta)
c
      common /par/ naca,tau,epsmax,ptmax
c
c  evaluate thickness and camber
c  for naca 4- or 5-digit airfoil
c
      thick=0.0
      if (z.lt.1.e-10) go to 100
      thick=5.*tau*(.2969*sqrt(z)-z*(.126+z*(.3537
     2	    -z*(.2843-z*.1015))))
  100 if (epsmax.eq.0.) go to 130
      if (naca.gt.9999) go to 140
      if (z.gt.ptmax) go to 110
      camber=epsmax/ptmax/ptmax*(2.*ptmax-z)*z
      dcamdx=2.*epsmax/ptmax/ptmax*(ptmax-z)
      go to 120
  110 camber=epsmax/(1.-ptmax)**2*(1.+z-2.*ptmax)*(1.-z)
      dcamdx=2.*epsmax/(1.-ptmax)**2*(ptmax-z)
  120 beta=atan(dcamdx)
c
      return
c
  130 camber=0.0
      beta=0.0
c
      return
c
  140 if (z.gt.ptmax) go to 150
      w=z/ptmax
      camber=epsmax*w*((w-3.)*w+3.-ptmax)
      dcamdx=epsmax*3.*w*(1.-w)/ptmax
      go to 120
  150 camber=epsmax*(1.-z)
      dcamdx=-epsmax
      go to 120
c
      end
