c-----------------------------------------------------------------------
      subroutine flegan(lmax,mmax,norm,x,plm)
c-----------------------------------------------------------------------
c     Calcul par recurrence de la valeur des fonctions de Legendre
c     associees Plm(x) pour l=0,...lmax, et m=0,...min(l,mmax) au point 
c     x. Si NORM=.True., alors les Plm sont normalises.
c     Les valeurs sont placees dans le tableau plm(l,m).     
c     Attention, ces Plm ne contiennent pas le facteur (-1)**m, ce
c     qui correspond aux definitions de la plupart des auteurs:
c     Joachain, Arkfen, Ayant&Borg, mais pas du Numerical Recipes.
c     Ref: Numerical Recipes 2ed, page 246
c          Ayant & Borg, page 143
c          Arfken, page 669
c          Joachain, page 667
c     Note: en raison de la factorielle (2*m-1) qui apparait dans le
c           calcul, la valeur des Plm grandit tres vite en fonction
c           de m. On peut aller jusqu'a m=150, apres on depasse 
c           10**308, et on obtient +INF.
c     Voir aussi la routine folegas.f pour les fonctions de Legendre
c     associees non normalisees.
c
      implicit double precision (a-h,o-z)
      dimension plm(0:lmax,0:mmax)
      logical norm
      double precision, allocatable, dimension(:) :: ftab
      data argmax/708.d0/
c
c     print *,'entree dans flegan: lmax =',lmax,' mmax =',mmax
      if (mmax.gt.lmax) mmax=lmax
      if (abs(x).gt.1.d0) then
	write(6,100) x
	stop
      end if
c
c     Calcul des Plm non normalises
c     (existe-t-il une formule de recurrence sur des Plm normalises?)
c
      sinteta = dsqrt(1.d0-x*x)
      plm(0,0) = 1.d0
      if (lmax.gt.0) then
        do 10 m=1,mmax
   10   plm(m,m) = plm(m-1,m-1) * (2*m-1)*sinteta
c
        do 20 m=0,min(mmax,lmax-1)
	  plm(m+1,m) = x*(2*m+1)*plm(m,m)  
          do 20 l=m+2,lmax
	  plm(l,m) = (x*(2*l-1)*plm(l-1,m)-(l+m-1)*plm(l-2,m)) / (l-m)
   20   continue
      end if
c
c     Coefficients de normalisation. Attention a l'argument de 
c     l'exponentielle qui peut etre grand et negatif. La valeur limite
c     est -708. Au dela, dexp(argu) vaudra zero. 
c
      if (norm) then
        ifmax = lmax + mmax
        allocate(ftab(0:ifmax))
        call lnfactrl(ifmax,ftab)
c
        do 30 l=0,lmax
	  fl = dsqrt(dble(2*l+1)/2.d0)
	  do 30 m=0,min(l,mmax)
            argu = 0.5d0*(ftab(l-m)-ftab(l+m))
	    if (dabs(argu).gt.argmax) then
              write(6,102)
              stop
            else
              plm(l,m) = plm(l,m) * fl * dexp(argu)
	    end if
   30   continue
      end if
c
      return
  100 format(/'Erreur dans FLEGAN: X =',f22.16,' est plus grand',
     >       ' que 1.0 en valeur absolue')
  102 format(/'Erreur dans FLEGAN: Underflow ou Overflow dans le ',
     >      ' calcul de l''exponentielle du facteur de normalisation')
      end 
c-----------------------------------------------------------------------

      double precision function cosd(x)
      implicit none
      double precision :: x, pi
      
      pi=acos(-1.d0)
       cosd = cos(pi*x/180.d0)
      
      
      
      return
      end function

c-----------------------------------------------------------------------
      subroutine lnfactrl(n,ftab)
c-----------------------------------------------------------------------
c     FACTRL retourne dans le vecteur ftab(i), i=0,n les valeurs du 
c     log neperien de i!, c.a.d. ln(i!).  
c
      implicit double precision (a-h,o-z)
      dimension ftab(0:n)
c
      if (n.lt.0) then
	write(6,101) n
	stop
      end if
c
c     ln(i!) = ln(1) + ln(2) + ln(3) + .... + ln(i)
c
      ftab(0) = 0.d0  
      if (n.eq.0) return
      ftab(1) = 0.d0
      do 10 i=2,n
        ftab(i) = dlog(dble(i)) +  ftab(i-1)
c       write(6,100) i,ftab(i),dexp(ftab(i))
   10 continue
c
      return
  100 format(i3,2e17.6e3)
  101 format(/'Erreur dans FACTRL: N =',i5)
      end
