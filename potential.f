c-----------------------------------------------------------------------
      double precision function potential(r,theta1,theta2,phi)
c-----------------------------------------------------------------------
c
c     Potential for C3-H2 [cm-1]
c     r: intermolecular distance [angstrom]
c     theta1: angle defined by intermolecular distance and the H2 molecule [degrees]
c     theta2: angle defined by intermolecular distance and the C3 molecule [degrees]
c     phi: relative orientation between H2 and C3 [degrees]
c
      implicit none 
      logical first
      integer:: ifirst1,ifb,k,nfb1
       double precision :: angles(3)
        double precision :: r,theta1,theta2,phi,x,t2,q_kernel
      double precision, allocatable :: rval(:),coeff(:,:),fonbas(:),
     >            icut(:),alpha(:),beta(:),tau(:),gamma(:),coekr(:,:)
      integer :: i,j, nr,nfb,l1max,l2max,l1sym,l2sym,mmax , nang
        
      integer, parameter :: iprint=0
       
        open(unit=10,file='COEFF_dummy',status='old')
        open(unit=87,file='w1',status='old')
        
        
         read(10,*) nr,nfb,l1max,l2max,l1sym,l2sym,mmax
        allocate(rval(nr),coekr(nr,nfb),fonbas(nfb))
         do k=1,nr
         read(10,*) rval(k)
         end do

        
        do ifb=1,nfb
          do j=1,nr
           read(87,*) coekr(j,ifb)
          end do 
        end do

        close(10)
        close(87)
        ifirst1 =1

      angles(1) = theta1
      angles(2) = theta2
      angles(3) = phi
        
       nang = nfb
      call primbasis(iprint,l1max,l2max,l1sym,l2sym,mmax,angles,
     >               nang,fonbas,nfb1)
           

      if (nfb1.ne.nfb) then
      print*,'NFB different de NFB1'
      stop 
      end if

      
        potential = 0.d0
       do ifb=1,nfb
         t2=0.d0
          do j=1,nr
           x=rval(j)
           t2 = t2 + coekr(j,ifb)*q_kernel(x,r)
          end do
        potential = potential + t2*fonbas(ifb)
        end do
        
      deallocate(fonbas,rval,coekr)  
      return
  101 format(/'Nr =',i3,', Nfb =',i3,', l1max =',i2,', l2max =',i2,
     >       ', l1sym =',i2,', l2sym =',i2,', mmax =',i2)
  102 format(i3,4e15.6)
  111 format(f10.4,50e12.4)
      end
c-----------------------------------------------------------------------
      subroutine primbasis(iprint,l1max,l2max,l1sym,l2sym,mmax,angles,
     >                     nang,fonbas,nfb)
c-----------------------------------------------------------------------
c     calcul de la valeur des fonctions:
c               P(l1,m,theta1)*P(l2,mtheta2)*cos(m*phi)
c     pour toutes les valeurs possibles de la triade (l1,l2,m) telles que:
c              l1 <= l1max, l2 <= l2max,  m <= min(l1,l2,mmax)
c     au point defini par:
c         angles(1)=theta1, angles(2)=theta2, angles(3)=phi
c
      implicit double precision (a-h,o-z)
      dimension angles(3), fonbas(nang)
      allocatable plmt1(:,:), plmt2(:,:)
      logical norm,first
      data norm/.true./,first/.true./
c
      allocate(plmt1(0:l1max,0:l1max),plmt2(0:l2max,0:l2max))
      xt1 = cosd(angles(1))
      call flegan(l1max,l1max,norm,xt1,plmt1)
      xt2 = cosd(angles(2))
      call flegan(l2max,l2max,norm,xt2,plmt2)
      ifb = 0
      do 10 l1=0,l1max,l1sym
        do 10 l2=0,l2max,l2sym
	  do 15 m=0,min(l1,l2,mmax)
	    ifb = ifb + 1
	    if (ifb.gt.nang) stop 'PRIMBASIS: IFB > NANG'
            if (first)
     >          write(6,103) ifb,l1,l2,m,plmt1(l1,m),plmt2(l2,m)
	    fonbas(ifb) = plmt1(l1,m) * plmt2(l2,m) * cosd(m*angles(3))
   15     continue
   10 continue
      nfb = ifb
      first = .false.
c
      return
  103 format(i3,', l1 =',i3,', l2 =',i3,', m =',i3,
     >       ', P1 =',f10.6,', P2 =',f10.6)
      end
c-----------------------------------------------------------------------
       double precision function q_kernel(x,r)
c-----------------------------------------------------------------------
        implicit none
        real*8 :: x,r,x_l,x_s
        x_l=r
           x_s=x
        if(r < x) then
            x_s=r
            x_l=x
        endif
   
		
c	    !for n=2 and m=4
            q_kernel = (2d0/(15d0*x_l**5)) * (1d0-5*x_s/(7*x_l))	
      end function q_kernel