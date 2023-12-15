!
!  numerical constants
!
!
!  last revised: 15 January 2011
!
      module numconstants
      use mpidefs
      implicit none
      logical, target :: light_up
      integer :: print_intermediate_results,global_rank
      integer, allocatable :: monen(:)
      integer, private :: nmax=0
      real(8) :: pi
      real(8), allocatable :: bcof(:,:),fnr(:),vwh_coef(:,:,:,:)
      real(8), allocatable :: vcc_const(:,:,:),fnm1_const(:,:),fn_const(:,:),fnp1_const(:,:)
      real(8), allocatable :: tran_coef(:,:,:)
      data pi/3.1415926535897932385/
      data light_up/.false./

      contains

         subroutine init(notd)
         implicit none
         integer :: notd,l,n,ierr,nbc,m,mm1,mp1,np1,nm1,nn1
         real(8) :: fnorm1,fnorm2
!
!  bcof(n,l)=((n+l)!/(n!l!))^(1/2)
!
         if(notd.le.nmax) return
         nmax=max(nmax,notd)
         nbc=6*notd+6
         if(allocated(fnr)) deallocate(monen,fnr,bcof)
         allocate (monen(0:2*notd),bcof(0:nbc,0:nbc),fnr(0:2*nbc),stat=ierr)
!         write(*,'('' nmax, bcof status:'',2i5)') nmax,ierr
         do n=0,2*notd
            monen(n)=(-1)**n
         enddo
         fnr(0)=0.d0
         do n=1,2*nbc
            fnr(n)=dsqrt(dble(n))
         enddo
         bcof(0,0)=1.d0
         do n=0,nbc-1
            do l=n+1,nbc
               bcof(n,l)=fnr(n+l)*bcof(n,l-1)/fnr(l)
               bcof(l,n)=bcof(n,l)
            enddo
            bcof(n+1,n+1)=fnr(n+n+2)*fnr(n+n+1)*bcof(n,n)/fnr(n+1)/fnr(n+1)
         enddo
         if(allocated(vwh_coef)) deallocate(vwh_coef)
         allocate(vwh_coef(-notd:notd,1:notd,-1:1,-1:1))
!
!  constants used for calculation of svwf functions.
!
         do n=1,notd
            nn1=n*(n+1)
            np1=n+1
            nm1=n-1
            fnorm1=-.5d0/fnr(n+n+1)/fnr(n)/fnr(n+1)
            fnorm2=-.5d0*fnr(n+n+1)/fnr(n)/fnr(n+1)
            m=-n
            mp1=m+1
            mm1=m-1
            vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
            vwh_coef(m,n, 1,-1)=fnorm1*np1*fnr(n-m)*fnr(nm1-m)
            vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
            vwh_coef(m,n,-1,-1)=0.d0
            vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
            vwh_coef(m,n, 0,-1)=0.d0
            vwh_coef(m,n, 1, 0)=-fnorm2*fnr(n-m)*fnr(np1+m)
            vwh_coef(m,n,-1, 0)=-0.d0
            vwh_coef(m,n, 0, 0)=-fnorm2*m
            do m=-n+1,-1
               mp1=m+1
               mm1=m-1
               vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
               vwh_coef(m,n, 1,-1)=fnorm1*np1*fnr(n-m)*fnr(nm1-m)
               vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
               vwh_coef(m,n,-1,-1)=-fnorm1*np1*fnr(n+m)*fnr(nm1+m)
               vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
               vwh_coef(m,n, 0,-1)=fnorm1*np1*fnr(n+m)*fnr(n-m)
               vwh_coef(m,n, 1, 0)=-fnorm2*fnr(n-m)*fnr(np1+m)
               vwh_coef(m,n,-1, 0)=-fnorm2*fnr(n+m)*fnr(np1-m)
               vwh_coef(m,n, 0, 0)=-fnorm2*m
            enddo
            do m=0,n-1
               mp1=m+1
               mm1=m-1
               vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
               vwh_coef(m,n, 1,-1)=fnorm1*np1*fnr(n-m)*fnr(nm1-m)
               vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
               vwh_coef(m,n,-1,-1)=-fnorm1*np1*fnr(n+m)*fnr(nm1+m)
               vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
               vwh_coef(m,n, 0,-1)=fnorm1*np1*fnr(n+m)*fnr(n-m)
               vwh_coef(m,n, 1, 0)=-fnorm2*fnr(n-m)*fnr(np1+m)
               vwh_coef(m,n,-1, 0)=-fnorm2*fnr(n+m)*fnr(np1-m)
               vwh_coef(m,n, 0, 0)=-fnorm2*m
            enddo
            m=n
            mp1=m+1
            mm1=m-1
            vwh_coef(m,n, 1, 1)=-fnorm1*n*fnr(np1+m)*fnr(np1+mp1)
            vwh_coef(m,n, 1,-1)=0.d0
            vwh_coef(m,n,-1, 1)=fnorm1*n*fnr(np1-m)*fnr(np1-mm1)
            vwh_coef(m,n,-1,-1)=-fnorm1*np1*fnr(n+m)*fnr(nm1+m)
            vwh_coef(m,n, 0, 1)=fnorm1*n*fnr(np1+m)*fnr(np1-m)
            vwh_coef(m,n, 0,-1)=0.d0
            vwh_coef(m,n, 1, 0)=-0.d0
            vwh_coef(m,n,-1, 0)=-fnorm2*fnr(n+m)*fnr(np1-m)
            vwh_coef(m,n, 0, 0)=-fnorm2*m
         enddo
         end subroutine init

      end module numconstants
!
!  special function for the multiple sphere problem
!
      module specialfuncs
      implicit none
      contains

         subroutine timewrite(iunit,char1,time,line_break)
         use intrinsics
         implicit none
         integer :: iunit
         real(8) :: time,time2
         logical :: linebreak
         logical, optional :: line_break
         character(*) :: char1
         if(present(line_break)) then
            linebreak=line_break
         else
            linebreak=.true.
         endif
         if(time.gt.3600.d0) then
            time2=time/3600.d0
            if(linebreak) then
               write(iunit,'(a,f9.3,'' hours'')') char1,time2
            else
               write(iunit,'(a,f9.3,'' hours'',$)') char1,time2
            endif
         elseif(time.gt.60.d0) then
            time2=time/60.d0
            if(linebreak) then
               write(iunit,'(a,f9.2,'' min'')') char1,time2
            else
               write(iunit,'(a,f9.2,'' min'',$)') char1,time2
            endif
         else
            if(linebreak) then
               write(iunit,'(a,f9.2,'' sec'')') char1,time
            else
               write(iunit,'(a,f9.2,'' sec'',$)') char1,time
            endif
         endif
         if(linebreak) call flush(iunit)
         end subroutine timewrite
!
!  ricatti-bessel function psi(n), real argument
!
         subroutine ricbessel(n,ds,eps,nmax,psi)
         implicit none
         integer :: n,nmax,ns,i
         real(8) :: ds,dns,sn,psi(0:n),psit,ds2,sum,eps,err
         if(int(ds).lt.n) then
            ns=nint(ds+4.*(ds**.3333d0)+17)
            ns=max(n+10,ns)
            dns=0.d0
            do i=ns-1,n,-1
               sn=dble(i+1)/ds
               dns=sn-1.d0/(dns+sn)
            enddo
            psi(n)=dns
            psi(n-1)=dble(n)/ds-1.d0/(dns+dble(n)/ds)
            do i=n-2,1,-1
               sn=dble(i+1)/ds
               psi(i)=sn-1.d0/(psi(i+1)+sn)
            enddo
            psit=dsin(ds)
            psi(0)=psit
            ds2=ds*ds
            sum=psit*psit/ds2
            do i=1,n
               psit=psit/(dble(i)/ds+psi(i))
               sum=sum+dble(i+i+1)*psit*psit/ds2
               err=dabs(1.d0-sum)
               psi(i)=psit
               if(err.lt.eps) then
                  nmax=i
                  return
               endif
            enddo
            nmax=n
         else
            psi(0)=dsin(ds)
            psi(1)=psi(0)/ds-dcos(ds)
            do i=1,n-1
               sn=dble(i+i+1)/ds
               psi(i+1)=sn*psi(i)-psi(i-1)
            enddo
            nmax=n
         endif
         end subroutine ricbessel
!
!  ricatti-hankel function xi(n), real argument
!
!
!  last revised: 15 January 2011
!
         subroutine richankel(n,ds,xi)
         implicit none
         integer :: n,i,ns
         real(8) :: ds,dns,sn,chi0,chi1,chi2,psi,psi0,psi1
         complex(8) :: xi(0:n)
         if(int(ds).lt.n) then
            ns=nint(ds+4.*(ds**.3333)+17)
            ns=max(n+10,ns)
            dns=0.d0
            do i=ns-1,n,-1
               sn=dble(i+1)/ds
               dns=sn-1.d0/(dns+sn)
            enddo
            xi(n)=dns
            xi(n-1)=dble(n)/ds-1.d0/(dns+dble(n)/ds)
            do i=n-2,1,-1
               sn=dble(i+1)/ds
               xi(i)=sn-1.d0/(xi(i+1)+sn)
            enddo
            chi0=-dcos(ds)
            psi=dsin(ds)
            chi1=chi0/ds-psi
            xi(0)=dcmplx(psi,chi0)
            do i=1,n
               chi2=dble(i+i+1)/ds*chi1-chi0
               psi=psi/(dble(i)/ds+xi(i))
               xi(i)=dcmplx(psi,chi1)
               chi0=chi1
               chi1=chi2
            enddo
            return
         else
            chi0=-dcos(ds)
            psi0=dsin(ds)
            chi1=chi0/ds-psi0
            psi1=psi0/ds+chi0
            xi(0)=dcmplx(psi0,chi0)
            xi(1)=dcmplx(psi1,chi1)
            do i=1,n-1
               sn=dble(i+i+1)/ds
               xi(i+1)=sn*xi(i)-xi(i-1)
            enddo
            return
         endif
         end subroutine richankel
!
!  ricatti-bessel function psi(n), complex argument
!
!
!  last revised: 15 January 2011
!
         subroutine cricbessel(n,ds,psi)
         implicit none
         integer :: n,i
         complex(8) :: ds,psi(0:n),chi(0:n)
         call cspherebessel(n,ds,psi,chi)
         do i=0,n
            psi(i)=psi(i)*ds
         enddo
         return
         end subroutine cricbessel
!
!  ricatti-hankel function psi(n), complex argument
!
!
!  last revised: 15 January 2011
!  March 2013
!  The condition abs(xi(i))/abs(psi(0)) << 1.d-6
!  implies an argument with large imag part, and use of xi = psi + i chi will have
!  round off problems.   Upwards recurrence is used in this case.
!
         subroutine crichankel(n,ds,xi)
         implicit none
         integer :: n,i
         complex(8) :: ds,psi(0:n),chi(0:n),xi(0:n),ci,&
                       psi0
         data ci/(0.d0,1.d0)/
         xi(0)=-ci*cdexp(ci*ds)
         psi0=cdsin(ds)
         if(cdabs(xi(0))/cdabs(psi0).lt.1.d-6) then
            xi(1)=-cdexp(ci*ds)*(ci+ds)/ds
            do i=2,n
               xi(i)=dble(i+i+1)/ds*xi(i-1)-xi(i-2)
            enddo
         else
            call cspherebessel(n,ds,psi,chi)
            do i=1,n
               xi(i)=(psi(i)+ci*chi(i))*ds
            enddo
         endif
         end subroutine crichankel
!
!     ==========================================================
!     Purpose: Compute spherical Bessel functions jn(z) & yn(z)
!              for a complex argument
!     Input :  z --- Complex argument
!              n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
!     Output:  CSJ(n) --- jn(z)
!              CSY(n) --- yn(z)
!              NM --- Highest order computed
!     Routines called:
!              MSTA1 and MSTA2 for computing the starting
!              point for backward recurrence
!     ==========================================================
!
!    obtained from, and copywrited by, Jian-Ming Jin
!    http://jin.ece.uiuc.edu/
!
!
!  last revised: 15 January 2011
!
         subroutine cspherebessel(n,z,csj,csy)
         implicit none
         integer :: n,nm,k,m
         real(8) :: a0
         complex(8) :: z,csj(0:n),csy(0:n),csa,csb,cs,cf0,cf1,cf
         a0=cdabs(z)
         nm=n
         if (a0.lt.1.0d-60) then
            csj=(0.d0,0.d0)
            csy=(-1.d300,0.d0)
            csy(0)=(1.d0,0.d0)
            return
         endif
         csj=(0.d0,0.d0)
         csj(0)=cdsin(z)/z
         csj(1)=(csj(0)-cdcos(z))/z
         if (n.ge.2) then
            csa=csj(0)
            csb=csj(1)
            m=msta1(a0,200)
            if (m.lt.n) then
               nm=m
            else
               m=msta2(a0,n,15)
            endif
            cf0=0.0d0
            cf1=1.0d0-100
            do k=m,0,-1
               cf=(2.0d0*k+3.0d0)*cf1/z-cf0
               if (k.le.nm) csj(k)=cf
               cf0=cf1
               cf1=cf
            enddo
            if (cdabs(csa).gt.cdabs(csb)) cs=csa/cf
            if (cdabs(csa).le.cdabs(csb)) cs=csb/cf0
            do k=0,min(nm,n)
               csj(k)=cs*csj(k)
            enddo
         endif
         csy=(1.d200,0.d0)
         csy(0)=-cdcos(z)/z
         csy(1)=(csy(0)-cdsin(z))/z
         do k=2,min(nm,n)
            if (cdabs(csj(k-1)).gt.cdabs(csj(k-2))) then
               csy(k)=(csj(k)*csy(k-1)-1.0d0/(z*z))/csj(k-1)
            else
               csy(k)=(csj(k)*csy(k-2)-(2.0d0*k-1.0d0)/z**3)/csj(k-2)
            endif
         enddo
         end subroutine cspherebessel

         subroutine ch12n ( n, z, nm, chf1)
         use numconstants

         !*****************************************************************************80
         !
         !! CH12N computes Hankel functions of first and second kinds, complex argument.
         !
         !  Discussion:
         !
         !    Both the Hankel functions and their derivatives are computed.
         !
         !  Licensing:
         !
         !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
         !    they give permission to incorporate this routine into a user program
         !    provided that the copyright is acknowledged.
         !
         !  Modified:
         !
         !    26 July 2012
         !
         !  Author:
         !
         !    Shanjie Zhang, Jianming Jin
         !
         !  Reference:
         !
         !    Shanjie Zhang, Jianming Jin,
         !    Computation of Special Functions,
         !    Wiley, 1996,
         !    ISBN: 0-471-11963-6,
         !    LC: QA351.C45.
         !
         !  Parameters:
         !
         !    Input, integer ( kind = 4 ) N, the order of the functions.
         !
         !    Input, complex ( kind = 8 ) Z, the argument.
         !
         !    Output, integer ( kind = 4 ) NM, the highest order computed.
         !
         !    Output, complex ( kind = 8 ) CHF1(0:n), CHD1(0:n), CHF2(0:n), CHD2(0:n),
         !    the values of Hn(1)(z), Hn(1)'(z), Hn(2)(z), Hn(2)'(z).
         !
           implicit none

           integer ( kind = 4 ) n

           complex ( kind = 8 ) cbi(0:n+1)
           complex ( kind = 8 ) cbj(0:n+1)
           complex ( kind = 8 ) cbk(0:n+1)
           complex ( kind = 8 ) cby(0:n+1)
           complex ( kind = 8 ) cdi(0:n+1)
           complex ( kind = 8 ) cdj(0:n+1)
           complex ( kind = 8 ) cdk(0:n+1)
           complex ( kind = 8 ) cdy(0:n+1)
           complex ( kind = 8 ) cf1
           complex ( kind = 8 ) cfac
           complex ( kind = 8 ) chf1(0:n)
           complex ( kind = 8 ) ci
           integer ( kind = 4 ) k
           integer ( kind = 4 ) nm
           complex ( kind = 8 ) z
           complex ( kind = 8 ) zi

           ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
           if ( imag ( z ).le. 0.0D+00 ) then
             call cjynb ( n, z, nm, cbj, cdj, cby, cdy )
             nm=min(n,nm)
             do k = 0, nm
               chf1(k) = cbj(k) + ci * cby(k)
             end do
           else
             zi = - ci * z
             call ciknb ( n, zi, nm, cbi, cdi, cbk, cdk )
             cf1 = -ci
             cfac = 2.0D+00 / ( pi * ci )
             nm=min(n,nm)
             do k = 0, nm
               chf1(k) = cfac * cbk(k)
               cfac = cfac * cf1
             end do
           end if
           return
         end subroutine ch12n

         subroutine ciknb ( n, z, nm, cbi, cdi, cbk, cdk)
         use numconstants

         !*****************************************************************************80
         !
         !! CIKNB computes complex modified Bessel functions In(z) and Kn(z).
         !
         !  Discussion:
         !
         !    This procedure also evaluates the derivatives.
         !
         !  Licensing:
         !
         !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
         !    they give permission to incorporate this routine into a user program
         !    provided that the copyright is acknowledged.
         !
         !  Modified:
         !
         !    30 July 2012
         !
         !  Author:
         !
         !    Shanjie Zhang, Jianming Jin
         !
         !  Reference:
         !
         !    Shanjie Zhang, Jianming Jin,
         !    Computation of Special Functions,
         !    Wiley, 1996,
         !    ISBN: 0-471-11963-6,
         !    LC: QA351.C45.
         !
         !  Parameters:
         !
         !    Input, integer ( kind = 4 ) N, the order of In(z) and Kn(z).
         !
         !    Input, complex ( kind = 8 ) Z, the argument.
         !
         !    Output, integer ( kind = 4 ) NM, the highest order computed.
         !
         !    Output, complex ( kind = 8 ) CB((0:N), CDI(0:N), CBK(0:N), CDK(0:N),
         !    the values of In(z), In'(z), Kn(z), Kn'(z).
         !
           implicit none

           integer ( kind = 4 ) n

           real ( kind = 8 ) a0
           complex ( kind = 8 ) ca0
           complex ( kind = 8 ) cbi(0:n+1)
           complex ( kind = 8 ) cdi(0:n+1)
           complex ( kind = 8 ) cbkl
           complex ( kind = 8 ) cbs
           complex ( kind = 8 ) cbk(0:n+1)
           complex ( kind = 8 ) cdk(0:n+1)
           complex ( kind = 8 ) cf
           complex ( kind = 8 ) cf0
           complex ( kind = 8 ) cf1
           complex ( kind = 8 ) cg
           complex ( kind = 8 ) cg0
           complex ( kind = 8 ) cg1
           complex ( kind = 8 ) ci
           complex ( kind = 8 ) cr
           complex ( kind = 8 ) cs0
           complex ( kind = 8 ) csk0
           real ( kind = 8 ) el
           real ( kind = 8 ) fac
           integer ( kind = 4 ) k
           integer ( kind = 4 ) k0
           integer ( kind = 4 ) l
           integer ( kind = 4 ) m
!           integer ( kind = 4 ) msta1
!           integer ( kind = 4 ) msta2
           integer ( kind = 4 ) nm
           real ( kind = 8 ) vt
           complex ( kind = 8 ) z
           complex ( kind = 8 ) z1
           el = 0.57721566490153D+00
           a0 = abs ( z )
           nm = n
           if ( a0 < 1.0D-100 ) then
             do k = 0, n
               cbi(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
               cbk(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
               cdi(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
               cdk(k) = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
             end do
             cbi(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
             cdi(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
             return
           end if
           ci = cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
           if ( real ( z, kind = 8 ) < 0.0D+00 ) then
             z1 = -z
           else
             z1 = z
           end if
           if ( n == 0 ) then
             nm = 1
           end if
           m = msta1 ( a0, 200 )
           if ( m < nm ) then
             nm = m
           else
             m = msta2 ( a0, nm, 15 )
           end if
           cbs = 0.0D+00
           csk0 = 0.0D+00
           cf0 = 0.0D+00
           cf1 = 1.0D-100
           do k = m, 0, -1
             cf = 2.0D+00 * ( k + 1.0D+00 ) * cf1 / z1 + cf0
             if ( k <= nm ) then
               cbi(k) = cf
             end if
             if ( k /= 0 .and. k == 2 * int ( k / 2 ) ) then
               csk0 = csk0 + 4.0D+00 * cf / k
             end if
             cbs = cbs + 2.0D+00 * cf
             cf0 = cf1
             cf1 = cf
           end do
           cs0 = exp ( z1 ) / ( cbs - cf )
           do k = 0, nm
             cbi(k) = cs0 * cbi(k)
           end do
           if ( a0 <= 9.0D+00 ) then
             cbk(0) = - ( log ( 0.5D+00 * z1 ) + el ) * cbi(0) + cs0 * csk0
             cbk(1) = ( 1.0D+00 / z1 - cbi(1) * cbk(0) ) / cbi(0)
           else
             ca0 = sqrt ( pi / ( 2.0D+00 * z1 ) ) * exp ( -z1 )
             if ( a0 < 25.0D+00 ) then
               k0 = 16
             else if ( a0 < 80.0D+00 ) then
               k0 = 10
             else if ( a0 < 200.0D+00 ) then
               k0 = 8
             else
               k0 = 6
             end if
             do l = 0, 1
               cbkl = 1.0D+00
               vt = 4.0D+00 * l
               cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
               do k = 1, k0
                 cr = 0.125D+00 * cr &
                   * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
                 cbkl = cbkl + cr
               end do
               cbk(l) = ca0 * cbkl
             end do
           end if
           cg0 = cbk(0)
           cg1 = cbk(1)
           do k = 2, nm
             cg = 2.0D+00 * ( k - 1.0D+00 ) / z1 * cg1 + cg0
             cbk(k) = cg
             cg0 = cg1
             cg1 = cg
           end do
           if ( real ( z, kind = 8 ) < 0.0D+00 ) then
             fac = 1.0D+00
             do k = 0, nm
               if ( imag ( z ) < 0.0D+00 ) then
                 cbk(k) = fac * cbk(k) + ci * pi * cbi(k)
               else
                 cbk(k) = fac * cbk(k) - ci * pi * cbi(k)
               end if
               cbi(k) = fac * cbi(k)
               fac = - fac
             end do
           end if
           cdi(0) = cbi(1)
           cdk(0) = -cbk(1)
           do k = 1, nm
             cdi(k) = cbi(k-1) - k / z * cbi(k)
             cdk(k) = - cbk(k-1) - k / z * cbk(k)
           end do
           return
         end subroutine ciknb

         subroutine bessel_integer_complex(n,z,nmax,b)
         use numconstants
         implicit none
         integer :: n,nmax
         complex(8) :: z,b(0:n),cbj(0:n+1),cdj(0:n+1),cby(0:n+1),cdy(0:n+1)

         call cjynb ( n, z, nmax, cbj, cdj, cby, cdy )
         nmax=min(n,nmax)
         b(0:nmax)=cbj(0:nmax)
         end subroutine bessel_integer_complex

         subroutine cjynb ( n, z, nm, cbj, cdj, cby, cdy )
         use numconstants

         !*****************************************************************************80
         !
         !! CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
         !
         !  Licensing:
         !
         !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
         !    they give permission to incorporate this routine into a user program
         !    provided that the copyright is acknowledged.
         !
         !  Modified:
         !
         !    03 August 2012
         !
         !  Author:
         !
         !    Shanjie Zhang, Jianming Jin
         !
         !  Reference:
         !
         !    Shanjie Zhang, Jianming Jin,
         !    Computation of Special Functions,
         !    Wiley, 1996,
         !    ISBN: 0-471-11963-6,
         !    LC: QA351.C45.
         !
         !  Parameters:
         !
         !    Input, integer ( kind = 4 ) N, the order of Jn(z) and Yn(z).
         !
         !    Input, complex ( kind = 8 ) Z, the argument of Jn(z) and Yn(z).
         !
         !    Output, integer ( kind = 4 ) NM, the highest order computed.
         !
         !    Output, complex ( kind = 8 ) CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N),
         !    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
         !
           implicit none
           integer ( kind = 4 ) n
           real ( kind = 8 ), save, dimension ( 4 ) :: a = (/ &
             -0.7031250000000000D-01, 0.1121520996093750D+00, &
             -0.5725014209747314D+00, 0.6074042001273483D+01 /)
           real ( kind = 8 ) a0
           real ( kind = 8 ), save, dimension ( 4 ) :: a1 = (/ &
             0.1171875000000000D+00,-0.1441955566406250D+00, &
             0.6765925884246826D+00,-0.6883914268109947D+01 /)
           real ( kind = 8 ), save, dimension ( 4 ) :: b = (/  &
             0.7324218750000000D-01,-0.2271080017089844D+00, &
             0.1727727502584457D+01,-0.2438052969955606D+02 /)
           real ( kind = 8 ), save, dimension ( 4 ) :: b1 = (/ &
            -0.1025390625000000D+00,0.2775764465332031D+00, &
            -0.1993531733751297D+01,0.2724882731126854D+02 /)
           complex ( kind = 8 ) cbj(0:n+1)
           complex ( kind = 8 ) cbj0
           complex ( kind = 8 ) cbj1
           complex ( kind = 8 ) cbjk
           complex ( kind = 8 ) cbs
           complex ( kind = 8 ) cby(0:n+1)
           complex ( kind = 8 ) cby0
           complex ( kind = 8 ) cby1
           complex ( kind = 8 ) cdj(0:n+1)
           complex ( kind = 8 ) cdy(0:n+1)
           complex ( kind = 8 ) ce
           complex ( kind = 8 ) cf
           complex ( kind = 8 ) cf1
           complex ( kind = 8 ) cf2
           complex ( kind = 8 ) cp0
           complex ( kind = 8 ) cp1
           complex ( kind = 8 ) cq0
           complex ( kind = 8 ) cq1
           complex ( kind = 8 ) cs0
           complex ( kind = 8 ) csu
           complex ( kind = 8 ) csv
           complex ( kind = 8 ) ct1
           complex ( kind = 8 ) ct2
           complex ( kind = 8 ) cu
           complex ( kind = 8 ) cyy
           real ( kind = 8 ) el
           integer ( kind = 4 ) k
           integer ( kind = 4 ) m
!           integer ( kind = 4 ) msta1
!           integer ( kind = 4 ) msta2
           integer ( kind = 4 ) nm
           real ( kind = 8 ) r2p
           real ( kind = 8 ) y0
           complex ( kind = 8 ) z

           el = 0.5772156649015329D+00
           r2p = 0.63661977236758D+00
           y0 = abs ( imag ( z ) )
           a0 = abs ( z )
           nm = n
           if ( a0 < 1.0D-100 ) then
             do k = 0, n
               cbj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
               cdj(k) = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
               cby(k) = - cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
               cdy(k) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
             end do
             cbj(0) = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
             cdj(1) = cmplx ( 0.5D+00, 0.0D+00, kind = 8 )
             return
           end if
           if ( a0 <= 300.0D+00 .or. 80 < n ) then
             if ( n == 0 ) then
               nm = 1
             end if
             m = msta1 ( a0, 200 )
             if ( m < nm ) then
               nm = m
             else
               m = msta2 ( a0, nm, 15 )
             end if
             cbs = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
             csu = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
             csv = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
             cf2 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
             cf1 = cmplx ( 1.0D-30, 0.0D+00, kind = 8 )
             do k = m, 0, -1
               cf = 2.0D+00 * ( k + 1.0D+00 ) / z * cf1 - cf2
               if ( k <= nm ) then
                 cbj(k) = cf
               end if
               if ( k == 2 * int ( k / 2 ) .and. k .ne. 0 ) then
                 if ( y0 <= 1.0D+00 ) then
                   cbs = cbs + 2.0D+00 * cf
                 else
                   cbs = cbs + ( -1.0D+00 ) ** ( k / 2 ) * 2.0D+00 * cf
                 end if
                 csu = csu + ( -1.0D+00 ) ** ( k / 2 ) * cf / k
               else if ( 1 < k ) then
                 csv = csv + ( -1.0D+00 ) ** ( k / 2 ) * k / ( k * k - 1.0D+00 ) * cf
               end if
               cf2 = cf1
               cf1 = cf
             end do
             if ( y0 <= 1.0D+00 ) then
               cs0 = cbs + cf
             else
               cs0 = ( cbs + cf ) / cos ( z )
             end if
             do k = 0, nm
               cbj(k) = cbj(k) / cs0
             end do
             ce = log ( z / 2.0D+00 ) + el
             cby(0) = r2p * ( ce * cbj(0) - 4.0D+00 * csu / cs0 )
             cby(1) = r2p * ( - cbj(0) / z + ( ce - 1.0D+00 ) * cbj(1) &
               - 4.0D+00 * csv / cs0 )
           else
             ct1 = z - 0.25D+00 * pi
             cp0 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
             do k = 1, 4
               cp0 = cp0 + a(k) * z ** ( - 2 * k )
             end do
             cq0 = -0.125D+00 / z
             do k = 1, 4
               cq0 = cq0 + b(k) * z ** ( - 2 * k - 1 )
             end do
             cu = sqrt ( r2p / z )
             cbj0 = cu * ( cp0 * cos ( ct1 ) - cq0 * sin ( ct1 ) )
             cby0 = cu * ( cp0 * sin ( ct1 ) + cq0 * cos ( ct1 ) )
             cbj(0) = cbj0
             cby(0) = cby0
             ct2 = z - 0.75D+00 * pi
             cp1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
             do k = 1, 4
               cp1 = cp1 + a1(k) * z ** ( - 2 * k )
             end do
             cq1 = 0.375D+00 / z
             do k = 1, 4
               cq1 = cq1 + b1(k) * z ** ( - 2 * k - 1 )
             end do
             cbj1 = cu * ( cp1 * cos ( ct2 ) - cq1 * sin ( ct2 ) )
             cby1 = cu * ( cp1 * sin ( ct2 ) + cq1 * cos ( ct2 ) )
             cbj(1) = cbj1
             cby(1) = cby1
             do k = 2, nm
               cbjk = 2.0D+00 * ( k - 1.0D+00 ) / z * cbj1 - cbj0
               cbj(k) = cbjk
               cbj0 = cbj1
               cbj1 = cbjk
             end do
           end if
           cdj(0) = -cbj(1)
           do k = 1, nm
             cdj(k) = cbj(k-1) - k / z * cbj(k)
           end do
           if ( 1.0D+00 < abs ( cbj(0) ) ) then
             cby(1) = ( cbj(1) * cby(0) - 2.0D+00 / ( pi * z ) ) / cbj(0)
           end if
           do k = 2, nm
             if ( abs ( cbj(k-2) ) <= abs ( cbj(k-1) ) ) then
               cyy = ( cbj(k) * cby(k-1) - 2.0D+00 / ( pi * z ) ) / cbj(k-1)
             else
               cyy = ( cbj(k) * cby(k-2) - 4.0D+00 * ( k - 1.0D+00 ) &
                 / ( pi * z * z ) ) / cbj(k-2)
             end if
             cby(k) = cyy
           end do
           cdy(0) = -cby(1)
           do k = 1, nm
             cdy(k) = cby(k-1) - k / z * cby(k)
           end do
           return
         end subroutine cjynb
!
!     ===================================================
!     Purpose: Determine the starting point for backward
!              recurrence such that the magnitude of
!              Jn(x) at that point is about 10^(-MP)
!     Input :  x     --- Argument of Jn(x)
!              MP    --- Value of magnitude
!     Output:  MSTA1 --- Starting point
!     ===================================================
!
!
!  last revised: 15 January 2011
!
         integer function msta1(x,mp)
         implicit none
         integer :: mp,n0,n1,it,nn
         real(8) :: x, a0,f1,f,f0
         a0=dabs(x)
         n0=int(1.1*a0)+1
         f0=envj(n0,a0)-mp
         n1=n0+5
         f1=envj(n1,a0)-mp
         do it=1,20
            nn=n1-(n1-n0)/(1.0d0-f0/f1)
            f=envj(nn,a0)-mp
            if(abs(nn-n1).lt.1) exit
            n0=n1
            f0=f1
            n1=nn
            f1=f
         enddo
         msta1=nn
         end function msta1
!
!     ===================================================
!     Purpose: Determine the starting point for backward
!              recurrence such that all Jn(x) has MP
!              significant digits
!     Input :  x  --- Argument of Jn(x)
!              n  --- Order of Jn(x)
!              MP --- Significant digit
!     Output:  MSTA2 --- Starting point
!     ===================================================
!
!
!  last revised: 15 January 2011
!
         integer function msta2(x,n,mp)
         implicit none
         integer :: n,mp,n0,n1,it,nn
         real(8) :: x,a0,hmp,ejn,obj,f0,f1,f
         a0=dabs(x)
         hmp=0.5d0*dble(mp)
         ejn=envj(n,a0)
         if (ejn.le.hmp) then
            obj=mp
            n0=int(1.1*a0)
         else
            obj=hmp+ejn
            n0=n
         endif
         f0=envj(n0,a0)-obj
         n1=n0+5
         f1=envj(n1,a0)-obj
         do it=1,20
            nn=n1-(n1-n0)/(1.0d0-f0/f1)
            f=envj(nn,a0)-obj
            if (abs(nn-n1).lt.1) exit
            n0=n1
            f0=f1
            n1=nn
            f1=f
         enddo
         msta2=nn+10
         end function msta2

         real(8) function envj(n,x)
         implicit none
         integer :: n
         real(8) :: x
         n=max(1,abs(n))
         envj=0.5d0*dlog10(6.28d0*n)-n*dlog10(1.36d0*x/n)
         end function envj
!
!    vector coupling coefficients vc(w) = C(m,n|k,l|m+k,w), w = |n-l|,... n+l
!    uses downwards and upwards recurrence
!
!
!  last revised: 15 January 2011
!
         subroutine vcfunc(m,n,k,l,vcn)
         use numconstants
         implicit none
         integer :: m,n,k,l,wmax,wmin,w,mk
         real(8) :: vcn(0:n+l),t1,t2,t3,vcmax,vctest,rat
         vcn=0.d0
         wmax=n+l
         wmin=max(abs(n-l),abs(m+k))
         vcn(wmax)=bcof(n+m,l+k)*bcof(n-m,l-k)/bcof(n+n,l+l)
         if(wmin.eq.wmax) return
         vcn(wmax-1)=vcn(wmax)*(l*m-k*n)*fnr(2*(l+n)-1)/fnr(l)/fnr(n)&
        &  /fnr(n+l+m+k)/fnr(n+l-m-k)
         if(wmin.eq.wmax-1) return
         mk=m+k
         vcmax=abs(vcn(wmax))+abs(vcn(wmax-1))
!
!  a downwards recurrence is used initially
!
         do w=wmax,wmin+2,-1
            t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk)&
        &     *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
            t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1))&
        &    /dble(2*w*(w-1))
            t3=fnr(w-mk-1)*fnr(w+mk-1)*fnr(l-n+w-1)*fnr(n-l+w-1)&
        &     *fnr(n+l-w+2)*fnr(n+l+w)/(dble(2*(w-1))*fnr(2*w-3)&
        &     *fnr(2*w-1))
            vcn(w-2)=(t2*vcn(w-1)-vcn(w)/t1)/t3
            if(mod(wmax-w,2).eq.1) then
               vctest=abs(vcn(w-2))+abs(vcn(w-1))
               vcmax=max(vcmax,vctest)
               rat=vctest/vcmax
!
!  if/when the coefficients start to decrease in magnitude, an upwards recurrence takes over
!
               if(rat.lt.0.01d0) exit
            endif
         enddo
         if(w-2.gt.wmin) then
            wmax=w-3
            call vcfuncuprec(m,n,k,l,wmax,vcn)
         endif
         end subroutine vcfunc
!
!  upwards VC coefficient recurrence
!
!
!  last revised: 15 January 2011
!
         subroutine vcfuncuprec(m,n,k,l,wmax,vcn)
         use numconstants
         implicit none
         integer :: m,n,k,l,wmax,w,mk,nl,m1,n1,l1,k1,w1,w2
         real(8) :: vcn(0:n+l),t1,t2,t3,vc1
         mk=abs(m+k)
         nl=abs(n-l)
         if(nl.ge.mk) then
            w=nl
            if(n.ge.l) then
               m1=m
               n1=n
               l1=l
               k1=k
            else
               m1=k
               n1=l
               k1=m
               l1=n
            endif
            vc1=(-1)**(k1+l1)*bcof(l1+k1,w-m1-k1) &
               *bcof(l1-k1,w+m1+k1)/bcof(l1+l1,w+w+1)
         else
            w=mk
            if(m+k.ge.0) then
               vc1=(-1)**(n+m)*bcof(n-l+w,l-k)*bcof(l-n+w,n-m) &
                  /bcof(w+w+1,n+l-w)
            else
               vc1=(-1)**(l+k)*bcof(n-l+w,l+k)*bcof(l-n+w,n+m) &
                 /bcof(w+w+1,n+l-w)
            endif
         endif
         w1=w
         vcn(w)=vc1
         w=w1+1
         mk=m+k
         w2=min(wmax,n+l)
         if(w2.gt.w1) then
            t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk) &
              *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
            if(w1.eq.0) then
               t2=.5*dble(m-k)
            else
               t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1)) &
                 /dble(2*w*(w-1))
            endif
            vcn(w)=t1*t2*vcn(w1)
         endif
         do w=w1+2,w2
            t1=2*w*fnr(w+w+1)*fnr(w+w-1)/(fnr(w+mk)*fnr(w-mk) &
              *fnr(n-l+w)*fnr(l-n+w)*fnr(n+l-w+1)*fnr(n+l+w+1))
            t2=dble((m-k)*w*(w-1)-mk*n*(n+1)+mk*l*(l+1)) &
             /dble(2*w*(w-1))
            t3=fnr(w-mk-1)*fnr(w+mk-1)*fnr(l-n+w-1)*fnr(n-l+w-1) &
              *fnr(n+l-w+2)*fnr(n+l+w)/(dble(2*(w-1))*fnr(2*w-3) &
              *fnr(2*w-1))
            vcn(w)=t1*(t2*vcn(w-1)-t3*vcn(w-2))
         enddo
         end subroutine vcfuncuprec
!
!  Normalized associated legendre functions
!
!
!  last revised: 15 January 2011
!
         subroutine normalizedlegendre(cbe,mmax,nmax,dc)
         use numconstants
         implicit none
         integer :: nmax,mmax,m,n,im
         real(8) :: dc(-mmax:mmax,0:nmax),cbe,sbe
         sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
         dc=0.d0
         do m=0,mmax
            dc(m,m)=(-1)**m*(0.5d0*sbe)**m*bcof(m,m)
            if(m.eq.nmax) exit
            dc(m,m+1)=fnr(m+m+1)*cbe*dc(m,m)
            do n=m+1,nmax-1
               dc(m,n+1)=(-fnr(n-m)*fnr(n+m)*dc(m,n-1)+dble(n+n+1)*cbe*dc(m,n)) &
                         /(fnr(n+1-m)*fnr(n+1+m))
            enddo
         enddo
         do m=1,mmax
            im=(-1)**m
            do n=m,nmax
               dc(-m,n)=im*dc(m,n)
            enddo
         enddo
         end subroutine normalizedlegendre

!!
!  Generalized spherical functions
!
!  dc(m,n*(n+1)+k)=(-1)^(m + k)((n - k)!(n + k)!/(n - m)!/(n + m)!)^(1/2)
!  ((1 + x)/2)^((m + k)/2)((1 - x)/2)^((k - m)/2)JacobiP[n - k, k - m, k + m, x]
!
!  for |m| <= kmax, n=0,1,...nmax, |k| <= n
!
!
!  last revised: 15 January 2011
!
         subroutine rotcoef(cbe,kmax,nmax,dc)
         use numconstants
         implicit none
         integer :: kmax,nmax,k,m,sin,n,knmax,nn1,kn,im,m1
         real(8) :: cbe,sbe,dc(-kmax:kmax,0:nmax*(nmax+2)),cbe2,sbe2,dk0(-nmax-1:nmax+1),&
                    dk01(-nmax-1:nmax+1),sben,dkt,fmn,dkm0,dkm1,dkn1
!if(light_up) then
!write(*,'('' rot1 '',3es13.5)') cbe
!call flush(6)
!endif
         dc=0.d0
         if(abs(cbe).ge.1.d0) then
            sbe=0.d0
         else
            sbe=dsqrt(abs((1.d0+cbe)*(1.d0-cbe)))
         endif
         cbe2=.5d0*(1.d0+cbe)
         sbe2=.5d0*(1.d0-cbe)
         sin=1
         dk0(0)=1.d0
         sben=1.d0
         dc(0,0)=1.d0
         dk01(0)=0.d0
!if(light_up) then
!write(*,'('' rot2 '',3es13.5)') cbe,sbe
!call flush(6)
!endif
         do n=1,nmax
            knmax=min(n,kmax)
            nn1=n*(n+1)
            sin=-sin
            sben=sben*sbe/2.d0
            if(sben.lt.1.d-30) sben=0.d0
            dk0(n)=sin*sben*bcof(n,n)
            dk0(-n)=sin*dk0(n)
            dk01(n)=0.d0
            dk01(-n)=0.d0
            dc(0,nn1+n)=dk0(n)
            dc(0,nn1-n)=dk0(-n)
!if(light_up.and.abs(cbe).gt.0.99999d0) then
!write(*,'('' rot2b '',i5,3es13.5)') n,sben
!call flush(6)
!endif
            do k=-n+1,n-1
               kn=nn1+k
               dkt=dk01(k)
               dk01(k)=dk0(k)
               dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt)&
                     /(fnr(n+k)*fnr(n-k))
               dc(0,kn)=dk0(k)
            enddo
            im=1
            do m=1,knmax
               im=-im
               fmn=1.d0/fnr(n-m+1)/fnr(n+m)
               m1=m-1
               dkm0=0.d0
               do k=-n,n
                  kn=nn1+k
                  dkm1=dkm0
                  dkm0=dc(m1,kn)
                  if(k.eq.n) then
                     dkn1=0.d0
                  else
                     dkn1=dc(m1,kn+1)
                  endif
                  dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1 &
                          -fnr(n-k)*fnr(n+k+1)*sbe2*dkn1  &
                          -dble(k)*sbe*dc(m1,kn))*fmn
                  dc(-m,nn1-k)=dc(m,kn)*(-1)**(k)*im
               enddo
            enddo
         enddo
!if(light_up) then
!write(*,'('' rot3 '',3es13.5)') cbe
!call flush(6)
!endif
         end subroutine rotcoef
!
!  Generalized spherical functions: complex argument
!
!  dc(m,n*(n+1)+k)=(-1)^(m + k)((n - k)!(n + k)!/(n - m)!/(n + m)!)^(1/2)
!  ((1 + x)/2)^((m + k)/2)((1 - x)/2)^((k - m)/2)JacobiP[n - k, k - m, k + m, x]
!
!  for |m| <= kmax, n=0,1,...nmax, |k| <= n
!
!
!  New: 08/25/2011
!
         subroutine crotcoef(cbe,kmax,nmax,dc,sin_beta)
         use numconstants
         implicit none
         integer :: kmax,nmax,k,m,in,n,knmax,nn1,kn,im,m1
         real(8) :: fmn
         complex(8) :: cbe,sbe,dc(-kmax:kmax,0:nmax*(nmax+2)),cbe2,sbe2,dk0(-nmax-1:nmax+1),&
                    dk01(-nmax-1:nmax+1),sben,dkt,dkm0,dkm1,dkn1
         complex(8), optional :: sin_beta
         if(present(sin_beta)) then
            sbe=sin_beta
         else
            sbe=cdsqrt((1.d0+cbe)*(1.d0-cbe))
         endif
         cbe2=.5d0*(1.d0+cbe)
         sbe2=.5d0*(1.d0-cbe)
         in=1
         dk0(0)=1.d0
         sben=1.d0
         dc(0,0)=1.d0
         dk01(0)=0.
         do n=1,nmax
            knmax=min(n,kmax)
            nn1=n*(n+1)
            in=-in
            sben=sben*sbe/2.d0
            dk0(n)=in*sben*bcof(n,n)
            dk0(-n)=in*dk0(n)
            dk01(n)=0.
            dk01(-n)=0.
            dc(0,nn1+n)=dk0(n)
            dc(0,nn1-n)=dk0(-n)
            do k=-n+1,n-1
               kn=nn1+k
               dkt=dk01(k)
               dk01(k)=dk0(k)
               dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt)&
                     /(fnr(n+k)*fnr(n-k))
               dc(0,kn)=dk0(k)
            enddo
            im=1
            do m=1,knmax
               im=-im
               fmn=1.d0/fnr(n-m+1)/fnr(n+m)
               m1=m-1
               dkm0=0.
               do k=-n,n
                  kn=nn1+k
                  dkm1=dkm0
                  dkm0=dc(m1,kn)
                  if(k.eq.n) then
                     dkn1=0.
                  else
                     dkn1=dc(m1,kn+1)
                  endif
                  dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1 &
                          -fnr(n-k)*fnr(n+k+1)*sbe2*dkn1  &
                          -dble(k)*sbe*dc(m1,kn))*fmn
                  dc(-m,nn1-k)=dc(m,kn)*(-1)**(k)*im
               enddo
            enddo
         enddo
         end subroutine crotcoef
!
! vector spherical harmonic function
! november 2011
! april 2012: lr formulation
! 2020 : back to nm formulation.   Complex cb
!
         subroutine complexpivec(cb,nodr,pivec,icon,lr_model,azimuth_angle,index_model)
         use numconstants
         implicit none
         logical :: lrmod
         logical, optional :: lr_model
         integer :: nodr,n,m,p,mn,i,nn1,imod,mnp,q
         integer, optional :: icon,index_model
         real(8) :: fnm,const,alpha
         real(8), optional :: azimuth_angle
         complex(8) :: cb,pivec(2*nodr*(nodr+2),2),drot(-1:1,0:nodr*(nodr+2)),tau(2),ci,cin,ephi
         data ci/(0.d0,1.d0)/

         if(present(index_model)) then
            imod=index_model
         else
            imod=1
         endif
         if(present(azimuth_angle)) then
            alpha=azimuth_angle
         else
            alpha=0.d0
         endif
         if(present(icon)) then
            i=icon
         else
            i=1
         endif
         if(present(lr_model)) then
            lrmod=lr_model
         else
            lrmod=.false.
         endif
         const=0.5d0/sqrt(2.d0*pi)
         call crotcoef(cb,1,nodr,drot)
         do n=1,nodr
            nn1=n*(n+1)
            fnm=sqrt(dble(n+n+1)/2.d0)/4.d0*const
            cin=4.d0*(-i*ci)**(n+1)
            do m=-n,n
               ephi=cdexp(i*ci*dble(m)*alpha)
               mn=nn1+m
               tau(1)=-fnm*(-drot(-1,mn)+drot(1,mn))*ephi
               tau(2)=-fnm*(drot(-1,mn)+drot(1,mn))*ephi
               do p=1,2
                  mnp=amnpaddress(m,n,p,nodr,imod)
                  pivec(mnp,1)=cin*tau(p)
                  pivec(mnp,2)=i*ci*cin*tau(3-p)
               enddo
               if(lrmod) then
                  do q=1,2
                     do p=1,2
                        mnp=amnpaddress(m,n,p,nodr,imod)
                        tau(p) =pivec(mnp,q)
                     enddo
                     mnp=amnpaddress(m,n,1,nodr,imod)
                     pivec(mnp,q)=tau(1)+tau(2)
                     mnp=amnpaddress(m,n,2,nodr,imod)
                     pivec(mnp,q)=tau(1)-tau(2)
                  enddo
               endif
            enddo
         enddo
         end subroutine complexpivec
!
!  tau are the vector spherical harmonic functions, normalized
!
!
!  last revised: 15 January 2011
!
         subroutine taufunc(cb,nmax,tau)
         use numconstants
         implicit none
         integer :: nmax,n,m,nn1,mn
         real(8) :: drot(-1:1,0:nmax*(nmax+2)),tau(0:nmax+1,nmax,2),cb,fnm
         call rotcoef(cb,1,nmax,drot)
         do n=1,nmax
            nn1=n*(n+1)
            fnm=sqrt(dble(n+n+1)/2.d0)/4.d0
            do m=-n,-1
               mn=nn1+m
               tau(n+1,-m,1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(n+1,-m,2)=-fnm*(drot(-1,mn)+drot(1,mn))
            enddo
            do m=0,n
               mn=nn1+m
               tau(m,n,1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(m,n,2)=-fnm*(drot(-1,mn)+drot(1,mn))
            enddo
         enddo
         end subroutine taufunc
!
!  rotation of expansion coefficients amn by euler angles alpha,beta,gamma
!  idir=1: forward rotation, idir=-1, reverse rotation.
!
!
!  last revised: 15 January 2011
!
         subroutine rotvec(alpha,beta,gamma,nmax,mmax,amn,idir)
         use numconstants
         implicit none
         integer :: nmax,mmax,idir,k,n,m,in,kmax,ka,na,im,m1
         real(8) :: dc(-nmax-1:nmax+1,-nmax-1:nmax+1),dk0(-nmax-1:nmax+1), &
                    dk01(-nmax-1:nmax+1),sbe,cbe,sbe2,cbe2,sben,dkt, &
                    fmn,dkm0,dkm1,alpha,beta,gamma
         complex(8) :: ealpha,amn(0:nmax+1,nmax,2),ealpham(-nmax:nmax), &
                       amnt(2,-nmax:nmax),a,b,ci,egamma,egammam(-nmax:nmax)
         data ci/(0.d0,1.d0)/
         call init(nmax)
         dc=0.d0
         dk01=0.d0
         dk0=0.d0
         ealpha=cdexp(ci*alpha)
         egamma=cdexp(ci*gamma)
         cbe=cos(beta)
         sbe=sqrt((1.d0+cbe)*(1.d0-cbe))
         cbe2=.5d0*(1.d0+cbe)
         sbe2=.5d0*(1.d0-cbe)
         call ephicoef(ealpha,nmax,ealpham)
         call ephicoef(egamma,nmax,egammam)
         in=1
         dk0(0)=1.d0
         sben=1.d0
         dk01(0)=0.d0
         do n=1,nmax
            kmax=min(n,mmax)
            do k=-kmax,kmax
               if(k.le.-1) then
                  ka=n+1
                  na=-k
               else
                  ka=k
                  na=n
               endif
               if(idir.eq.1) then
                  amnt(1,k)=amn(ka,na,1)*ealpham(k)
                  amnt(2,k)=amn(ka,na,2)*ealpham(k)
               else
                  amnt(1,-k)=amn(ka,na,1)*egammam(k)
                  amnt(2,-k)=amn(ka,na,2)*egammam(k)
               endif
            enddo
            in=-in
            sben=sben*sbe/2.d0
            dk0(n)=in*sben*bcof(n,n)
            dk0(-n)=in*dk0(n)
            dk01(n)=0.d0
            dk01(-n)=0.d0
            dc(0,n)=dk0(n)
            dc(0,-n)=dk0(-n)
            do k=-n+1,n-1
               dkt=dk01(k)
               dk01(k)=dk0(k)
               dk0(k)=(cbe*(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt) &
                     /(fnr(n+k)*fnr(n-k))
               dc(0,k)=dk0(k)
            enddo
            im=1
            do m=1,kmax
               im=-im
               fmn=1./fnr(n-m+1)/fnr(n+m)
               m1=m-1
               dkm0=0.
               do k=-n,n
                  dkm1=dkm0
                  dkm0=dc(m1,k)
                  dc(m,k)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1 &
                         -fnr(n-k)*fnr(n+k+1)*sbe2*dc(m1,k+1) &
                         -k*sbe*dc(m1,k))*fmn
                  dc(-m,-k)=dc(m,k)*(-1)**(k)*im
               enddo
            enddo
            do m=-n,n
               if(m.le.-1) then
                  ka=n+1
                  na=-m
               else
                  ka=m
                  na=n
               endif
               a=0.
               b=0.
               do k=-kmax,kmax
                  a=a+dc(-k,-m)*amnt(1,k)
                  b=b+dc(-k,-m)*amnt(2,k)
               enddo
               if(idir.eq.1) then
                  amn(ka,na,1)=a*egammam(m)
                  amn(ka,na,2)=b*egammam(m)
               else
                  amn(ka,na,1)=a*ealpham(m)
                  amn(ka,na,2)=b*ealpham(m)
               endif
            enddo
         enddo
         end subroutine rotvec
!
!  regular vswf expansion coefficients for a plane wave: general case, complex cos beta
!
         subroutine genplanewavecoef(alpha,cb,nodr,pmnp0,lr_tran)
         use numconstants
         implicit none
         logical :: lrtran
         logical, optional :: lr_tran
         integer :: nodr,m,n,p,sp,nn1,mn
         real(8) :: alpha,fnm,ca,sa
         complex(8) :: drot(-1:1,0:nodr*(nodr+2)),tau(0:nodr+1,nodr,2), &
              taulr(0:nodr+1,nodr,2),cb,ealpha,ci,cin, &
              pmnp0(0:nodr+1,nodr,2,2),ealpham(-nodr:nodr)
         data ci/(0.d0,1.d0)/
         if(present(lr_tran)) then
            lrtran=lr_tran
         else
            lrtran=.true.
         endif
         call crotcoef(cb,1,nodr,drot)
         do n=1,nodr
            nn1=n*(n+1)
            fnm=sqrt(dble(n+n+1)/2.d0)/4.d0
            do m=-n,-1
               mn=nn1+m
               tau(n+1,-m,1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(n+1,-m,2)=-fnm*(drot(-1,mn)+drot(1,mn))
            enddo
            do m=0,n
               mn=nn1+m
               tau(m,n,1)=-fnm*(-drot(-1,mn)+drot(1,mn))
               tau(m,n,2)=-fnm*(drot(-1,mn)+drot(1,mn))
            enddo
         enddo
         ca=cos(alpha)
         sa=sin(alpha)
         ealpha=dcmplx(ca,sa)
         call ephicoef(ealpha,nodr,ealpham)
         if(lrtran) then
            taulr(:,:,1)=(tau(:,:,1)+tau(:,:,2))*.5d0
            taulr(:,:,2)=(tau(:,:,1)-tau(:,:,2))*.5d0
            do n=1,nodr
               cin=4.d0*ci**(n+1)
               do p=1,2
                  sp=-(-1)**p
                  do m=-n,-1
                     pmnp0(n+1,-m,p,1)=-cin*taulr(n+1,-m,p)*ealpham(-m)
                     pmnp0(n+1,-m,p,2)=sp*ci*cin*taulr(n+1,-m,p)*ealpham(-m)
                  enddo
                  do m=0,n
                     pmnp0(m,n,p,1)=-cin*taulr(m,n,p)*ealpham(-m)
                     pmnp0(m,n,p,2)=sp*ci*cin*taulr(m,n,p)*ealpham(-m)
                  enddo
               enddo
            enddo
         else
            do n=1,nodr
               cin=4.d0*ci**(n+1)
               do p=1,2
                  do m=-n,-1
                     pmnp0(n+1,-m,p,1)=-cin*tau(n+1,-m,p)*ealpham(-m)
                     pmnp0(n+1,-m,p,2)=ci*cin*tau(n+1,-m,3-p)*ealpham(-m)
                  enddo
                  do m=0,n
                     pmnp0(m,n,p,1)=-cin*tau(m,n,p)*ealpham(-m)
                     pmnp0(m,n,p,2)=ci*cin*tau(m,n,3-p)*ealpham(-m)
                  enddo
               enddo
            enddo
         endif
         end subroutine genplanewavecoef

         subroutine gaussianbeamcoef(alpha,cbeta,cbeam,nodr,pmnp0,lr_tran)
         use numconstants
         implicit none
         logical :: lrtran
         logical, optional :: lr_tran
         integer :: nodr,m,n,p,k
         real(8) :: alpha,cbeta,cbeam,gbn
         complex(8) :: ccb,pmnp0(0:nodr+1,nodr,2,2)
         if(present(lr_tran)) then
            lrtran=lr_tran
         else
            lrtran=.true.
         endif
         ccb=cbeta
         call genplanewavecoef(alpha,ccb,nodr,pmnp0,lr_tran=lrtran)
         do n=1,nodr
            gbn=dexp(-((dble(n)+.5d0)*cbeam)**2.)
            do p=1,2
               do k=1,2
                  do m=-n,-1
                     pmnp0(n+1,-m,p,k)=pmnp0(n+1,-m,p,k)*gbn
                  enddo
                  do m=0,n
                     pmnp0(m,n,p,k)=pmnp0(m,n,p,k)*gbn
                  enddo
               enddo
            enddo
         enddo
         end subroutine gaussianbeamcoef
!
!  axial translation coefficients calculated by the diamond recurrence formula
!  new: 10 october 2011
!  april 2012: lr formulation
!  may 2012: new ordering scheme:
!  input:  itype : 1 or 3 (regular, outgoing)
!     r: axial translation distance (positive)
!     ri: rank 2 complex array: L and R refractive indices of medium
!     nmax, lmax: largest row and column orders.
!     ndim: dimension of ac
!  output:
!     ac:  rank 1 complex array, dimension ndim, containing the matrix elements.
!  storage scheme:   for each degree m, with ordering m=0, -1, 1, -2, 2, ..min(nmax,lmax),
!  the elements for degree m are stored
!
         subroutine axialtrancoefrecurrence(itype,r,ri,nmax,lmax,ndim,ac)
         use numconstants
         implicit none
         integer :: itype,nmax,lmax,n,l,m,p,nlmin, &
                    wmin,wmax,ml,m1,np1,nm1,lm1,lp1,sp
         integer :: iadd,nlmax,iadd0,iadd1,ndim
         integer :: ma,blockdim
         integer, save :: nlmax0
         real(8) :: r
         complex(8) :: ri(2),ci,z(2),xi(0:nmax+lmax,2)
         complex(8) :: ac(ndim),act(nmax,lmax,2),actt(2,2)
         data ci,nlmax0/(0.d0,1.d0),0/
         nlmax=max(nmax,lmax)
         nlmin=min(nmax,lmax)
         if(nlmax.gt.nlmax0) then
            nlmax0=nlmax
            call axialtrancoefinit(nlmax)
         endif

         if(r.le.1.d-12) then
            ac=(0.d0,0.d0)
            if(itype.ne.1) return
            iadd0=0
            do ma=0,nlmin
               m1=max(1,ma)
               do m=-ma,ma,2*m1
                  blockdim=(nmax-m1+1)*(lmax-m1+1)*2
                  iadd1=iadd0+blockdim
                  act=0.d0
                  do l=m1,nlmin
                     act(l,l,1)=1.d0
                     act(l,l,2)=1.d0
                  enddo
                  ac(iadd0+1:iadd1)=reshape(act(m1:nmax,m1:lmax,1:2),(/blockdim/))
                  iadd0=iadd1
               enddo
            enddo
            return
         endif
         z=r*ri
         do p=1,2
            if(itype.eq.1) then
               call cricbessel(nmax+lmax,z(p),xi(0:,p))
            else
               call crichankel(nmax+lmax,z(p),xi(0:,p))
            endif
            xi(0:,p)=xi(0:,p)/z(p)
            if(z(1).eq.z(2)) then
               xi(0:nmax+lmax,2)=xi(0:nmax+lmax,1)
               exit
            endif
         enddo
         lm1=lmax-1

         iadd0=0
         do ma=0,nlmin
            m1=max(1,ma)
            lp1=m1+1
            do m=-ma,ma,2*m1
               blockdim=2*(nmax-m1+1)*(lmax-m1+1)
               iadd1=iadd0+blockdim
               n=m1
               do l=m1,lmax
                  wmin=abs(n-l)
                  wmax=n+l
                  iadd=iadd+1
                  ml=l*(l+1)+m
                  do p=1,2
                     actt(1,p)=sum(vcc_const(n,ml,wmin:wmax:2)*xi(wmin:wmax:2,p))
                     actt(2,p)=ci*sum(vcc_const(n,ml,wmin+1:wmax-1:2)*xi(wmin+1:wmax-1:2,p))
                  enddo
                  act(n,l,1)=actt(1,1)+actt(2,1)
                  act(n,l,2)=actt(1,2)-actt(2,2)
               enddo
               l=lmax
               ml=l*(l+1)+m
               do n=m1+1,nmax
                  wmin=abs(n-l)
                  wmax=n+l
                  do p=1,2
                     actt(1,p)=sum(vcc_const(n,ml,wmin:wmax:2)*xi(wmin:wmax:2,p))
                     actt(2,p)=ci*sum(vcc_const(n,ml,wmin+1:wmax-1:2)*xi(wmin+1:wmax-1:2,p))
                  enddo
                  act(n,l,1)=actt(1,1)+actt(2,1)
                  act(n,l,2)=actt(1,2)-actt(2,2)
               enddo

               if(m1.lt.nlmin) then
                  do n=m1,nmax-1
                     np1=n+1
                     nm1=n-1
                     do p=1,2
                        sp=-(-1)**p
                        act(np1,m1:lmax-1,p)= &
                          - act(n,m1+1:lmax,p)*fnp1_const(m,m1:lm1) &
                          + sp*(fn_const(m,m1:lm1)-fn_const(m,n))*ci*act(n,m1:lm1,p)
                        act(np1,m1+1:lm1,p)=act(np1,m1+1:lm1,p) &
                          + act(n,m1:lmax-2,p)*fnm1_const(m,lp1:lm1)
                        if(n.gt.m1) then
                           act(np1,m1:lm1,p)=act(np1,m1:lm1,p) &
                             + act(nm1,m1:lm1,p)*fnm1_const(m,n)
                        endif
                        act(np1,m1:lm1,p)=act(np1,m1:lm1,p)/fnp1_const(m,n)
                     enddo
                  enddo
               endif
               ac(iadd0+1:iadd1)=reshape(act(m1:nmax,m1:lmax,1:2),(/blockdim/))
               iadd0=iadd1
            enddo
         enddo
         end subroutine axialtrancoefrecurrence
!
!  constants for translation coefficient calculation
!
         subroutine axialtrancoefinit(nmax)
         use numconstants
         implicit none
         integer :: nmax,m,n,l,w,n21,ml,ll1,wmin,wmax,nlmin,lp1,lm1
         real(8) :: c1,c2,vc1(0:2*nmax),vc2(0:2*nmax)
         complex(8) :: ci,inlw
         data ci/(0.d0,1.d0)/
         if(allocated(vcc_const)) deallocate(vcc_const,fnm1_const,fn_const,fnp1_const)
         allocate(vcc_const(nmax,nmax*(nmax+2),0:2*nmax),fnm1_const(-nmax:nmax,nmax), &
                  fn_const(-nmax:nmax,nmax),fnp1_const(-nmax:nmax,nmax))
         do n=1,nmax
            n21=n+n+1
            do l=1,nmax
               c1=fnr(n21)*fnr(l+l+1)
               ll1=l*(l+1)
               call vcfunc(-1,n,1,l,vc2)
               wmin=abs(n-l)
               wmax=n+l
               nlmin=min(l,n)
               do m=-nlmin,nlmin
                  ml=ll1+m
                  c2=-c1*(-1)**m
                  call vcfunc(-m,n,m,l,vc1)
                  do w=wmin,wmax
                     inlw=ci**(n-l+w)
                     vcc_const(n,ml,w)=c2*vc1(w)*vc2(w)*(dble(inlw)+dimag(inlw))
                  enddo
               enddo
            enddo
         enddo
         fnm1_const=0.
         fn_const=0.
         fnp1_const=0.
         do m=-nmax,nmax
            do l=max(1,abs(m)),nmax
               lp1=l+1
               lm1=l-1
               fnm1_const(m,l)=fnr(lm1)*fnr(lp1)*fnr(l-m)*fnr(l+m)/fnr(lm1+l)/fnr(l+lp1)/dble(l)
               fn_const(m,l)=dble(m)/dble(l)/dble(lp1)
               fnp1_const(m,l)=fnr(l)*fnr(l+2)*fnr(lp1-m)*fnr(lp1+m)/fnr(l+lp1)/fnr(l+l+3)/dble(lp1)
            enddo
         enddo
         end subroutine axialtrancoefinit

         subroutine gentrancoefconstants(nodrmax)
         use numconstants
         implicit none
         integer :: nodrmax,v,w,wmax,wmin,n,l,m,k,m1m,mn,kl
         real(8) :: vc1(0:2*nodrmax),vc2(0:2*nodrmax)
         complex(8) :: ci,c,a
         data ci/(0.d0,1.d0)/
         if(allocated(tran_coef)) deallocate(tran_coef)
         allocate(tran_coef(nodrmax*(nodrmax+2),nodrmax*(nodrmax+2),0:2*nodrmax))
         tran_coef=0.d0
         do l=1,nodrmax
            do n=1,nodrmax
               wmax=n+l
               call vcfunc(-1,n,1,l,vc2)
               c=-ci**(n-l)*fnr(n+n+1)*fnr(l+l+1)
               do k=-l,l
                  kl=l*(l+1)+k
                  do m=-n,n
                     mn=n*(n+1)+m
                     m1m=(-1)**m
                     v=k-m
                     call vcfunc(-m,n,k,l,vc1)
                     wmin=max(abs(v),abs(n-l))
                     do w=wmax,wmin,-1
                        a=ci**w*c*m1m*vc1(w)*vc2(w)
                        if(mod(wmax-w,2).eq.0) then
                           tran_coef(mn,kl,w)=dble(a)
                        else
                           tran_coef(mn,kl,w)=dimag(a)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
         return
         end subroutine gentrancoefconstants

         subroutine gentranmatrix(nodr_s,nodr_t,translation_vector, &
                             refractive_index,ac_matrix,vswf_type, &
                             mode_s,mode_t,index_model)
         use numconstants
         implicit none
         integer :: nodr_s,nodr_t,nodrmax,wmax,p,n,m,k,l,mn,kl,imodel, &
                    nblks,nblkt,w,v,wmin,itype,nmodes,nmodet,nmode,mna,kla
         integer, optional :: vswf_type,mode_s,mode_t,index_model
         integer, save :: setnodrmax
         real(8) :: r,ct,xp(3),ymn(-nodr_s-nodr_t:nodr_s+nodr_t,0:nodr_s+nodr_t)
         real(8) :: translation_vector(3)
         complex(8) :: ri(2),ephi,rri,ephim(-nodr_s-nodr_t:nodr_s+nodr_t), &
               hn(0:nodr_s+nodr_t,2),a1,a2,b1,b2
         complex(8) :: ac_matrix(nodr_t*(nodr_t+2),nodr_s*(nodr_s+2),1:2)
         complex(8), optional :: refractive_index(2)
         data setnodrmax/0/
         nblks=nodr_s*(nodr_s+2)
         nblkt=nodr_t*(nodr_t+2)
         if(present(mode_s)) then
            nmodes=mode_s
         else
            nmodes=2
         endif
         if(present(mode_t)) then
            nmodet=mode_t
         else
            nmodet=2
         endif
         if(present(index_model)) then
            imodel=index_model
         else
            imodel=2
         endif
         nmode=max(nmodet,nmodes)
         xp=translation_vector
         if(present(refractive_index)) then
            ri=refractive_index
         else
            ri=(1.d0,0.d0)
         endif
         if(present(vswf_type)) then
            itype=vswf_type
         else
            itype=3
         endif
         nodrmax=max(nodr_s,nodr_t)
         if(nodrmax.gt.setnodrmax) then
            setnodrmax=nodrmax
            call gentrancoefconstants(nodrmax)
         endif
         wmax=nodr_s+nodr_t
         r=xp(1)*xp(1)+xp(2)*xp(2)+xp(3)*xp(3)
         if(r.eq.0.d0) then
            ac_matrix=0.d0
            if(itype.eq.1) then
               do n=1,min(nblks,nblkt)
                  ac_matrix(n,n,1:nmode)=1.d0
               enddo
            endif
         else
            r=sqrt(r)
            ct=xp(3)/r
            if(xp(1).eq.0.d0.and.xp(2).eq.0.d0) then
               ephi=(1.d0,0.d0)
            else
               ephi=dcmplx(xp(1),xp(2))/sqrt(xp(1)*xp(1)+xp(2)*xp(2))
            endif
            ephim(0)=1.d0
            do m=1,wmax
               ephim(m)=ephi*ephim(m-1)
               ephim(-m)=conjg(ephim(m))
            enddo
            call normalizedlegendre(ct,wmax,wmax,ymn)
            do p=1,2
               rri=r*ri(p)
               if(itype.eq.3) then
                  hn(0,p)=-(0.d0,1.d0)*cdexp((0.d0,1.d0)*rri)/rri
                  hn(1,p)=-cdexp((0.d0,1.d0)*rri)*((0.d0,1.d0)+rri)/rri/rri
                  do n=2,wmax
                     hn(n,p)=dble(n+n-1)/rri*hn(n-1,p)-hn(n-2,p)
                  enddo
               else
                  call cricbessel(wmax,rri,hn(:,p))
                  hn(:,p)=hn(:,p)/rri
               endif
               if(ri(2).eq.ri(1)) then
                  hn(:,2)=hn(:,1)
                  exit
               endif
            enddo
            do n=1,nodr_s
               do m=-n,n
                  mna=amnaddress(m,n,nodr_s,imodel)
                  mn=n*(n+1)+m
                  do l=1,nodr_t
                     wmax=n+l
                     do k=-l,l
                        kla=amnaddress(k,l,nodr_t,imodel)
                        kl=l*(l+1)+k
                        v=m-k
                        wmin=max(abs(v),abs(n-l))
                        a1=0.
                        a2=0.
                        b1=0.
                        b2=0.
                        do w=wmax,wmin,-1
                           if(mod(wmax-w,2).eq.0) then
                              a1=a1+hn(w,1)*ymn(v,w)*tran_coef(kl,mn,w)
                              if(nmode.eq.1) cycle
                              a2=a2+hn(w,2)*ymn(v,w)*tran_coef(kl,mn,w)
                           else
                              b1=b1+hn(w,1)*ymn(v,w)*tran_coef(kl,mn,w)
                              b2=b2+hn(w,2)*ymn(v,w)*tran_coef(kl,mn,w)
                           endif
                        enddo
                        if(nmode.eq.1) then
                           ac_matrix(kla,mna,1)=ephim(v)*a1
                        else
                           ac_matrix(kla,mna,1)=ephim(v)*(a1+(0.d0,1.d0)*b1)
                           ac_matrix(kla,mna,2)=ephim(v)*(a2-(0.d0,1.d0)*b2)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         endif
         end subroutine gentranmatrix
!
!  test to determine convergence of regular vswf addition theorem for max. order lmax
!  and translation distance r w/ refractive index ri.
!
!
!  last revised: 15 January 2011
!
         subroutine tranordertest(r,ri,lmax,eps,nmax)
         use numconstants
         implicit none
         integer :: nmax,lmax,n,l,m,w,n21,wmin,wmax
         integer, parameter :: nlim=200
         real(8) :: r,alnw,sum,eps
         real(8) :: vc1(0:nlim+lmax)
         complex(8) :: ri,ci,z,a,b,c
         complex(8) :: xi(0:nlim+lmax)
         data ci/(0.d0,1.d0)/
         if(r.eq.0.d0) then
            nmax=lmax
            return
         endif
         z=r*dble(ri)
         sum=0.d0
         do n=1,nlim
            call init(n+lmax)
            call cricbessel(n+lmax,z,xi)
            do l=0,n+lmax
               xi(l)=xi(l)/z*ci**l
            enddo
            n21=n+n+1
            l=lmax
            c=fnr(n21)*fnr(l+l+1)*ci**(n-l)
            call vcfunc(-1,n,1,l,vc1)
            wmin=abs(n-l)
            wmax=n+l
            m=1
            a=0.
            b=0.
            do w=wmin,wmax
               alnw=vc1(w)*vc1(w)
               if(mod(n+l+w,2).eq.0) then
                  a=a+alnw*xi(w)
               else
                  b=b+alnw*xi(w)
               endif
            enddo
            a=c*a
            b=c*b
            sum=sum+a*conjg(a)+b*conjg(b)
            if(abs(1.d0-sum).lt.eps) exit
         enddo
         nmax=min(n,nlim)
         nmax=max(nmax,lmax)
         end subroutine tranordertest
!
!  address for axial translation coefficient
!
!
!  last revised: 15 January 2011
!

         integer function atcadd(m,n,ntot)
         implicit none
         integer :: m,n,ntot
         atcadd=n-ntot+(max(1,m)*(1+2*ntot-max(1,m)))/2+ntot*min(1,m)
         end function atcadd

         integer function atcdim(ntot,ltot)
         implicit none
         integer :: ntot,ltot,nmin,nmax
         nmin=min(ntot,ltot)
         nmax=max(ntot,ltot)
         atcdim=2*(nmin*(1- nmin*nmin + 3*nmax*(2 + nmin)))/3
         end function atcdim
!
! the offset (integer) for the ntot X ltot translation matrix for degree m
!
         integer function moffset(m,ntot,ltot)
         implicit none
         integer :: m,ntot,ltot
         if(m.eq.0) then
            moffset=0
         elseif(m.lt.0) then
            moffset=2*(-((1+m)*(2+m)*(3+2*m +3*ntot)) &
                   -3*ltot*(2+ntot+m*(3+m+2*ntot)))/3
         else
            moffset=2*(-3*ltot*(-1 + m)**2 +6*ltot*m*ntot &
                   +(-1+m)*(m*(-4+2*m-3*ntot)+3*(1 + ntot)))/3
         endif
         end function moffset
!
! cartosphere takes the cartesian point (x,y,z) = xp(1), xp(2), xp(3)
! and converts to polar form: r: radius, ct: cos(theta), ep = exp(i phi)
!
!
!  last revised: 15 January 2011
!
         subroutine cartosphere(xp,r,ct,ep)
         implicit none
         real(8) :: xp(3),r,ct
         complex(8) :: ep
         r=xp(1)*xp(1)+xp(2)*xp(2)+xp(3)*xp(3)
         if(r.eq.0.d0) then
            ct=1.d0
            ep=(1.d0,0.d0)
            return
         endif
         r=sqrt(r)
         ct=xp(3)/r
         if(xp(1).eq.0.d0.and.xp(2).eq.0.d0) then
            ep=(1.d0,0.d0)
         else
            ep=dcmplx(xp(1),xp(2))/sqrt(xp(1)*xp(1)+xp(2)*xp(2))
         endif
         return
         end subroutine cartosphere

         subroutine cartospherevec(nt,xp,xps)
         implicit none
         integer :: nt,i
         real(8) :: xp(3,nt),xps(3,nt),r,ct,phi
         do i=1,nt
            r=xp(1,i)*xp(1,i)+xp(2,i)*xp(2,i)+xp(3,i)*xp(3,i)
            if(r.eq.0.d0) then
               ct=1.d0
               phi=0.d0
            else
               r=sqrt(r)
               ct=xp(3,i)/r
               if(xp(1,i).eq.0.d0.and.xp(2,i).eq.0.d0) then
                  phi=0.d0
               else
                  phi=datan2(xp(2,i),xp(1,i))
               endif
            endif
            xps(:,i)=(/ct,phi,r/)
         enddo
         return
         end subroutine cartospherevec

!
! euler rotation of a point (x,y,z) = xp(1), xp(2), xp(3)
! November 2012
!
         subroutine eulerrotation(xp,eulerangf,dir,xprot,num)
         implicit none
         integer :: dir,n,i
         integer, optional :: num
         real(8) :: xp(3,*),eulerangf(3),eulerang(3),cang(3),sang(3), &
                    mat1(3,3),mat2(3,3),mat3(3,3),xprot(3,*),xpt(3)
         if(present(num)) then
            n=num
         else
            n=1
         endif
         if(dir.eq.1) then
            eulerang=eulerangf
         else
            eulerang(1:3)=-eulerangf(3:1:-1)
         endif
         cang=cos(eulerang)
         sang=sin(eulerang)
         mat1(1,:) = (/cang(1),sang(1),0.d0/)
         mat1(2,:) = (/-sang(1),cang(1),0.d0/)
         mat1(3,:) = (/0.d0,0.d0,1.d0/)
         mat2(1,:) = (/cang(2),0.d0,-sang(2)/)
         mat2(2,:) = (/0.d0,1.d0,0.d0/)
         mat2(3,:) = (/sang(2),0.d0,cang(2)/)
         mat3(1,:) = (/cang(3),sang(3),0.d0/)
         mat3(2,:) = (/-sang(3),cang(3),0.d0/)
         mat3(3,:) = (/0.d0,0.d0,1.d0/)
         do i=1,n
            xpt=xp(:,i)
            xpt=matmul(mat1,xpt)
            xpt=matmul(mat2,xpt)
            xpt=matmul(mat3,xpt)
            xprot(:,i)=xpt
         enddo
         end subroutine eulerrotation
!
! ephicoef returns the complex array epm(m) = exp(i m phi) for
! m=-nodr,nodr.   ep =exp(i phi), and epm is dimensioned epm(-nd:nd)
!
!
!  last revised: 15 January 2011
!
         subroutine ephicoef(ep,nodr,epm)
         implicit none
         integer :: nodr,m
         complex(8) :: ep,epm(-nodr:nodr)
         epm(0)=(1.d0,0.d0)
         do m=1,nodr
            epm(m)=ep*epm(m-1)
            epm(-m)=dconjg(epm(m))
         enddo
         return
         end subroutine ephicoef
!
!  test to determine max order of vswf expansion of a plane wave at distance r
!
!
!  last revised: 15 January 2011
!
         subroutine planewavetruncationorder(r,rimedium,eps,nodr)
         implicit none
         integer :: nodr,n1,n
         real(8) :: r,eps,err
         complex(8), allocatable :: jn(:)
         complex(8) :: sum, ci,eir,rimedium(2),rri,rib
         data ci/(0.d0,1.d0)/
         rib=2.d0/(1.d0/rimedium(1)+1.d0/rimedium(2))
         n1=max(10,int(3.*r+1))
         allocate(jn(0:n1))
         rri=r*rib
         call cricbessel(n1,rri,jn)
         jn(0:n1)=jn(0:n1)/rri
         eir=cdexp(-ci*rri)
         sum=jn(0)*eir
         do n=1,n1
            sum=sum+ci**n*dble(n+n+1)*jn(n)*eir
            err=cdabs(1.d0-sum)
            if(err.lt.eps) then
               nodr=n
               deallocate(jn)
               return
            endif
         enddo
         nodr=n1
         deallocate(jn)
         end subroutine planewavetruncationorder
!
!  calculates the cartesian components of the vswf at position rpos, in ref. index ri.
!
!
!  original: 15 January 2011
!  revised: 23 February 2011: multiplied by root 2
!  april 2012: lr formulation
!
         subroutine vwhcalc(rpos,ri,nodr,itype,vwh,index_model,lr_to_mode)
         use numconstants
         implicit none
         logical, optional :: lr_to_mode
         logical :: lrtomode
         integer :: nodr,itype,n,nodrp1,nodrm1,nn1,np1,nm1,p,sp,imod,m,iadd(-nodr:nodr)
         integer, save :: nodrmax
         integer, optional :: index_model
         real(8) ::  rpos(3),r,ct
         real(8) pmn(0:0,0:(nodr+1)*(nodr+3))
         complex(8) :: ci,vwh(3,2*nodr*(nodr+2)),ri(2),ephi,a(2),vtemp(3,2)
         complex(8)  :: a1vec(-nodr:nodr), &
                       b1vec(-nodr:nodr),z1vec(-nodr:nodr),a2vec(-nodr:nodr), &
                       b2vec(-nodr:nodr),z2vec(-nodr:nodr)
         complex(8) :: umn(-nodr-2:nodr+2,0:nodr+1,2), hn(0:nodr+1,2), ephim(-nodr-1:nodr+1)
         data ci,nodrmax/(0.d0,1.d0),0/
         if(nodr.gt.nodrmax) then
            nodrmax=nodr
            call init(nodr+2)
         endif
         if(present(index_model)) then
            imod=index_model
         else
            imod=1
         endif
         if(present(lr_to_mode)) then
            lrtomode=lr_to_mode
         else
            lrtomode=.false.
         endif
         call cartosphere(rpos,r,ct,ephi)
         if(r.le.1.d-5) then
            vwh(:,1:2*nodr*(nodr+2))=(0.d0,0.d0)
            if(itype.eq.3) return
            do p=1,2
               vwh(1,amnpaddress(-1,1,p,nodr,imod))=.5d0*fnr(2)/fnr(3)
               vwh(2,amnpaddress(-1,1,p,nodr,imod))=-.5d0*ci*fnr(2)/fnr(3)
               vwh(3,amnpaddress(0,1,p,nodr,imod))=1.d0*fnr(2)/fnr(6)
               vwh(1,amnpaddress(1,1,p,nodr,imod))=-.5d0*fnr(2)/fnr(3)
               vwh(2,amnpaddress(1,1,p,nodr,imod))=-.5d0*ci*fnr(2)/fnr(3)
               if(lrtomode) exit
            enddo
            return
         endif
         nodrp1=nodr+1
         nodrm1=nodr-1
!
! this is now a vector operation w/ l/r form
!
         a=ri*r
         do p=1,2
            if(itype.eq.1) then
               call cricbessel(nodrp1,a(p),hn(0,p))
            else
               call crichankel(nodrp1,a(p),hn(0,p))
            endif
            hn(0:nodrp1,p)=hn(0:nodrp1,p)/a(p)
            if(a(2).eq.a(1)) then
               hn(0:nodrp1,2)=hn(0:nodrp1,1)
               exit
            endif
         enddo
         call rotcoef(ct,0,nodrp1,pmn)
         call ephicoef(ephi,nodrp1,ephim)
         umn=0.d0
         do p=1,2
            umn(0,0,p)=hn(0,p)*fnr(2)
            do n=1,nodrp1
               nn1=n*(n+1)
               umn(-n:n,n,p)=fnr(2)*pmn(0,nn1-n:nn1+n)*ephim(-n:n)*hn(n,p)
               umn(-n-1,n,p)=0.d0
               umn(n+1,n,p)=0.d0
            enddo
         enddo
         do p=1,2
            sp=-(-1)**p
            do n=1,nodr
               do m=-n,n
                  iadd(m)=amnpaddress(m,n,p,nodr,imod)
               enddo
               nn1=n*(n+1)
               np1=n+1
               nm1=n-1
               a1vec(-n:n)=vwh_coef(-n:n,n,1,1)*umn(-nm1:np1,np1,p) &
                  +vwh_coef(-n:n,n,1,-1)*umn(-nm1:np1,nm1,p)
               b1vec(-n:n)=vwh_coef(-n:n,n,-1,1)*umn(-np1:nm1,np1,p) &
                  +vwh_coef(-n:n,n,-1,-1)*umn(-np1:nm1,nm1,p)
               z1vec(-n:n)=vwh_coef(-n:n,n,0,1)*umn(-n:n,np1,p) &
                  +vwh_coef(-n:n,n,0,-1)*umn(-n:n,nm1,p)
               a2vec(-n:n)=vwh_coef(-n:n,n,1,0)*umn(-nm1:np1,n,p)
               b2vec(-n:n)=vwh_coef(-n:n,n,-1,0)*umn(-np1:nm1,n,p)
               z2vec(-n:n)=vwh_coef(-n:n,n,0,0)*umn(-n:n,n,p)
               vwh(1,iadd(-n:n))=-0.5d0*(a1vec(-n:n)+b1vec(-n:n)) &
                      -sp*0.5d0*ci*(a2vec(-n:n)+b2vec(-n:n))
               vwh(2,iadd(-n:n))=-0.5d0*ci*(-a1vec(-n:n)+b1vec(-n:n)) &
                      -sp*0.5d0*(a2vec(-n:n)-b2vec(-n:n))
               vwh(3,iadd(-n:n))=-z1vec(-n:n) &
                      -sp*ci*z2vec(-n:n)
            enddo
         enddo
         if(lrtomode) then
            do n=1,nodr
               do m=-n,n
                  do p=1,2
                     vtemp(:,p)=vwh(:,amnpaddress(m,n,p,nodr,imod))
                  enddo
                  vwh(:,amnpaddress(m,n,1,nodr,imod))=(vtemp(:,1)+vtemp(:,2))*0.5d0
                  vwh(:,amnpaddress(m,n,2,nodr,imod))=(vtemp(:,1)-vtemp(:,2))*0.5d0
               enddo
            enddo
         endif
         end subroutine vwhcalc

         subroutine scalar_wave_function(nodr,itype,x,y,z,ri,swf)
         use numconstants
         implicit none
         integer :: nodr,itype,n,m,mn
         real(8) :: x,y,z,r,ct,rho,ymn(0:nodr*(nodr+2)),c,c0
         complex(8) :: ri,swf(0:nodr*(nodr+2)),ephi,rri,hn(0:nodr)

         r=sqrt(x*x+y*y+z*z)
         if(r.lt.1.d-10) then
            swf=0.d0
            if(itype.eq.1) swf(0)=1.d0
            return
         endif
         ct=z/r
         if(x.eq.0.d0.and.y.eq.0.d0) then
            ephi=1.d0
            rho=0.d0
         else
            rho=sqrt(x*x+y*y)
            ephi=dcmplx(x,y)/rho
         endif
         call rotcoef(ct,0,nodr,ymn)
         rri=r*ri
         if(itype.eq.3) then
            hn(0)=-(0.d0,1.d0)*cdexp((0.d0,1.d0)*rri)/rri
            hn(1)=-cdexp((0.d0,1.d0)*rri)*((0.d0,1.d0)+rri)/rri/rri
            do n=2,nodr
               hn(n)=dble(n+n-1)/rri*hn(n-1)-hn(n-2)
            enddo
         else
            call cricbessel(nodr,rri,hn)
            hn=hn/rri
         endif

         c0=sqrt(1.d0/4.d0/pi)
         do n=0,nodr
            c=c0*sqrt(dble(n+n+1))
            do m=-n,n
               mn=n*(n+1)+m
               swf(mn)=hn(n)*ymn(mn)*ephi**m*c
            enddo
         enddo
         end subroutine scalar_wave_function

         subroutine reciprocal_scalar_wave_function(nodr,kx,ky,x,y,z,ri,swf)
         use numconstants
         implicit none
         integer :: nodr,n,m,mn
         real(8) :: kx,ky,x,y,z,k
         complex(8) :: ri,swf(0:nodr*(nodr+2)),ephi,kz,ymn(0:nodr*(nodr+2)),skz,c,cr
         k=sqrt(kx*kx+ky*ky)
         kz=cdsqrt((1.d0,0.d0)-k*k/ri/ri)
         if(z.gt.0.d0) then
            skz=kz
         else
            skz=-kz
         endif
         if(k.eq.0.d0) then
            ephi=1.d0
         else
            ephi=dcmplx(kx,ky)/k
         endif
         call crotcoef(skz,0,nodr,ymn)
         c=cdexp((0.d0,1.d0)*(kx*x+ky*y+ri*skz*z))/ri/ri/kz/sqrt(4.d0*pi)
         do n=0,nodr
            cr=((0.d0,-1.d0))**n*sqrt(dble(n+n+1))
            do m=-n,n
               mn=n*(n+1)+m
               swf(mn)=cr*c*ymn(mn)*ephi**m
            enddo
         enddo
         end subroutine reciprocal_scalar_wave_function
!
! inverse of a 2 X 2 complex matrix.
! March 2013
!
         subroutine twobytwoinverse(mat,imat)
         implicit none
         integer :: s,t,ss,st
         complex(8) :: mat(2,2),imat(2,2),tmat(2,2),det
         tmat=mat
         det=mat(1,1)*mat(2,2)-mat(2,1)*mat(1,2)
         do s=1,2
            ss=(-1)**s
            do t=1,2
               st=(-1)**t
               imat(s,t)=ss*st*tmat(3-t,3-s)/det
            enddo
         enddo
         end subroutine twobytwoinverse
!
! move between unequal size matrices.
! March 2013
!
         subroutine mtransfer(nin,nout,cin,cout)
         implicit none
         integer :: nin,nout,nmin
         complex(8) :: cin(0:nin+1,nin,2),cout(0:nout+1,nout,2), &
            ct(0:max(nin,nout)+1,max(nin,nout),2)
         nmin=min(nin,nout)
         ct=0.d0
         ct(0:nin+1,1:nin,1:2)=cin(0:nin+1,1:nin,1:2)
         cout=0.d0
         cout(0:nout+1,1:nout,1:2)=ct(0:nout+1,1:nout,1:2)
         end subroutine mtransfer
!
! packed address value for m,n pair
! model 1: mn=n*(n+1)+m
! model 2: mn = which(m>=0, (n-1)(L+2)+m+1,
!                     m<0, -(m+1)(L+2)+n+2)
!
         function amnaddress(m,n,l,model)
         implicit none
         integer :: amnaddress,m,n,l,model
         if(model.eq.1) then
            amnaddress=n*(n+1)+m
         else
            if(m.ge.0) then
               amnaddress=(n-1)*(l+2)+m+1
            else
               amnaddress=-(m+1)*(l+2)+n+2
            endif
         endif
         end function amnaddress

         function amnpaddress(m,n,p,l,model)
         implicit none
         integer :: amnpaddress,m,n,p,l,model
         if(model.eq.1) then
            amnpaddress=2*(n*(n+1)+m-1)+p
         else
            if(m.ge.0) then
               amnpaddress=(n-1)*(l+2)+m+1+(p-1)*l*(l+2)
            else
               amnpaddress=-(m+1)*(l+2)+n+2+(p-1)*l*(l+2)
            endif
         endif
         end function amnpaddress

         subroutine lr_mode_transformation(nodr,alr,amode,lr_to_mode)
         implicit none
         logical :: lrtomode
         logical, optional :: lr_to_mode
         integer :: nodr
         complex(8) :: alr(nodr*(nodr+2),2),amode(nodr*(nodr+2),2),at(nodr*(nodr+2),2)
         if(present(lr_to_mode)) then
            lrtomode=lr_to_mode
         else
            lrtomode=.true.
         endif
         if(lrtomode) then
            at=alr(:,:)
            amode(:,1)=at(:,1)+at(:,2)
            amode(:,2)=at(:,1)-at(:,2)
         else
            at=amode(:,:)
            alr(:,1)=.5d0*(at(:,1)+at(:,2))
            alr(:,2)=.5d0*(at(:,1)-at(:,2))
         endif
         end subroutine lr_mode_transformation

         subroutine degree_transformation(nodr,ain,aout)
         implicit none
         integer :: nodr,m,n,p,mnp,mnp2,im,m1
         complex(8) :: ain(2*nodr*(nodr+2)),aout(2*nodr*(nodr+2))
         do m=-nodr,nodr
            m1=max(abs(m),1)
            im=(-1)**m
            do n=m1,nodr
               do p=1,2
                  mnp=amnpaddress(m,n,p,nodr,2)
                  mnp2=amnpaddress(-m,n,p,nodr,2)
                  aout(mnp2)=im*ain(mnp)
               enddo
            enddo
         enddo
         end subroutine degree_transformation

         subroutine groupfilename(firststring,number,laststring,newstring)
         implicit none
         integer :: number
         character*256 :: firststring,laststring,newstring,sform,intfile
         if(number.lt.10) then
            sform='(a,i1,a,a)'
         elseif(number.lt.100) then
            sform='(a,i2,a,a)'
         elseif(number.lt.1000) then
            sform='(a,i3,a,a)'
         else
            sform='(a,i4,a,a)'
         endif
         write(intfile,fmt=sform) trim(firststring),number,'_',trim(laststring)
         read(intfile,'(a)') newstring
         end subroutine groupfilename

!*****************************************************************************80
!
!! QNG estimates an integral, using non-adaptive integration.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!    The routine is a simple non-adaptive automatic integrator, based on
!    a sequence of rules with increasing degree of algebraic
!    precision (Patterson, 1968).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983

         subroutine qng (n,a,b,epsabs,epsrel,f,resultf,abserr,neval,ier)
         implicit none
         integer :: n,ier,k,l,neval,ipx
         real(8) a,absc,abserr,b,centr,dhlgth,epsabs,epsrel,hlgth,w10(5),w21a(5),w21b(6), &
            w43a(10),w43b(12),w87a(21),w87b(23),x1(5),x2(5),x3(11),x4(22)
         complex(8) :: fcentr(n),fval(n),fval1(n),fval2(n),fv1(5,n),fv2(5,n),fv3(5,n),fv4(5,n), &
            resultf(n),res10(n),res21(n),res43(n),res87(n),savfun(21,n)
         external :: f
         data x1(1),x1(2),x1(3),x1(4),x1(5)/ &
              9.739065285171717E-01,     8.650633666889845E-01, &
              6.794095682990244E-01,     4.333953941292472E-01, &
              1.488743389816312E-01/
         data x2(1),x2(2),x2(3),x2(4),x2(5)/ &
              9.956571630258081E-01,     9.301574913557082E-01, &
              7.808177265864169E-01,     5.627571346686047E-01, &
              2.943928627014602E-01/
         data x3(1),x3(2),x3(3),x3(4),x3(5),x3(6),x3(7),x3(8),x3(9),x3(10), &
           x3(11)/ &
              9.993333609019321E-01,     9.874334029080889E-01, &
              9.548079348142663E-01,     9.001486957483283E-01, &
              8.251983149831142E-01,     7.321483889893050E-01, &
              6.228479705377252E-01,     4.994795740710565E-01, &
              3.649016613465808E-01,     2.222549197766013E-01, &
              7.465061746138332E-02/
         data x4(1),x4(2),x4(3),x4(4),x4(5),x4(6),x4(7),x4(8),x4(9),x4(10), &
           x4(11),x4(12),x4(13),x4(14),x4(15),x4(16),x4(17),x4(18),x4(19), &
           x4(20),x4(21),x4(22)/         9.999029772627292E-01, &
              9.979898959866787E-01,     9.921754978606872E-01, &
              9.813581635727128E-01,     9.650576238583846E-01, &
              9.431676131336706E-01,     9.158064146855072E-01, &
              8.832216577713165E-01,     8.457107484624157E-01, &
              8.035576580352310E-01,     7.570057306854956E-01, &
              7.062732097873218E-01,     6.515894665011779E-01, &
              5.932233740579611E-01,     5.314936059708319E-01, &
              4.667636230420228E-01,     3.994248478592188E-01, &
              3.298748771061883E-01,     2.585035592021616E-01, &
              1.856953965683467E-01,     1.118422131799075E-01, &
              3.735212339461987E-02/
         data w10(1),w10(2),w10(3),w10(4),w10(5)/ &
              6.667134430868814E-02,     1.494513491505806E-01, &
              2.190863625159820E-01,     2.692667193099964E-01, &
              2.955242247147529E-01/
         data w21a(1),w21a(2),w21a(3),w21a(4),w21a(5)/ &
              3.255816230796473E-02,     7.503967481091995E-02, &
              1.093871588022976E-01,     1.347092173114733E-01, &
              1.477391049013385E-01/
         data w21b(1),w21b(2),w21b(3),w21b(4),w21b(5),w21b(6)/ &
              1.169463886737187E-02,     5.475589657435200E-02, &
              9.312545458369761E-02,     1.234919762620659E-01, &
              1.427759385770601E-01,     1.494455540029169E-01/
         data w43a(1),w43a(2),w43a(3),w43a(4),w43a(5),w43a(6),w43a(7), &
           w43a(8),w43a(9),w43a(10)/     1.629673428966656E-02, &
              3.752287612086950E-02,     5.469490205825544E-02, &
              6.735541460947809E-02,     7.387019963239395E-02, &
              5.768556059769796E-03,     2.737189059324884E-02, &
              4.656082691042883E-02,     6.174499520144256E-02, &
              7.138726726869340E-02/
         data w43b(1),w43b(2),w43b(3),w43b(4),w43b(5),w43b(6),w43b(7), &
           w43b(8),w43b(9),w43b(10),w43b(11),w43b(12)/ &
              1.844477640212414E-03,     1.079868958589165E-02, &
              2.189536386779543E-02,     3.259746397534569E-02, &
              4.216313793519181E-02,     5.074193960018458E-02, &
              5.837939554261925E-02,     6.474640495144589E-02, &
              6.956619791235648E-02,     7.282444147183321E-02, &
              7.450775101417512E-02,     7.472214751740301E-02/
         data w87a(1),w87a(2),w87a(3),w87a(4),w87a(5),w87a(6),w87a(7), &
           w87a(8),w87a(9),w87a(10),w87a(11),w87a(12),w87a(13),w87a(14), &
           w87a(15),w87a(16),w87a(17),w87a(18),w87a(19),w87a(20),w87a(21)/ &
              8.148377384149173E-03,     1.876143820156282E-02, &
              2.734745105005229E-02,     3.367770731163793E-02, &
              3.693509982042791E-02,     2.884872430211531E-03, &
              1.368594602271270E-02,     2.328041350288831E-02, &
              3.087249761171336E-02,     3.569363363941877E-02, &
              9.152833452022414E-04,     5.399280219300471E-03, &
              1.094767960111893E-02,     1.629873169678734E-02, &
              2.108156888920384E-02,     2.537096976925383E-02, &
              2.918969775647575E-02,     3.237320246720279E-02, &
              3.478309895036514E-02,     3.641222073135179E-02, &
              3.725387550304771E-02/
         data w87b(1),w87b(2),w87b(3),w87b(4),w87b(5),w87b(6),w87b(7), &
           w87b(8),w87b(9),w87b(10),w87b(11),w87b(12),w87b(13),w87b(14), &
           w87b(15),w87b(16),w87b(17),w87b(18),w87b(19),w87b(20),w87b(21), &
           w87b(22),w87b(23)/            2.741455637620724E-04, &
              1.807124155057943E-03,     4.096869282759165E-03, &
              6.758290051847379E-03,     9.549957672201647E-03, &
              1.232944765224485E-02,     1.501044734638895E-02, &
              1.754896798624319E-02,     1.993803778644089E-02, &
              2.219493596101229E-02,     2.433914712600081E-02, &
              2.637450541483921E-02,     2.828691078877120E-02, &
              3.005258112809270E-02,     3.164675137143993E-02, &
              3.305041341997850E-02,     3.425509970422606E-02, &
              3.526241266015668E-02,     3.607698962288870E-02, &
              3.669860449845609E-02,     3.712054926983258E-02, &
              3.733422875193504E-02,     3.736107376267902E-02/
!
!  Test on validity of parameters.
!
         resultf = 0.0D+00
         abserr = 0.0D+00
         neval = 0

         hlgth = 5.0D-01 * ( b - a )
         dhlgth = abs ( hlgth )
         centr = 5.0D-01 * ( b + a )
         call f(n,centr,fcentr)
         neval = 21
         ier = 1

         do l = 1, 3
            if ( l == 1 ) then
               res10 = 0.0D+00
               res21 = w21b(6) * fcentr
               do k = 1, 5
                  absc = hlgth * x1(k)
                  call f(n,centr+absc,fval1)
                  call f(n,centr-absc,fval2)
                  fval = fval1 + fval2
                  res10 = res10 + w10(k)*fval
                  res21 = res21 + w21a(k)*fval
                  savfun(k,:) = fval(:)
                  fv1(k,:) = fval1(:)
                  fv2(k,:) = fval2(:)
               enddo
               ipx = 5
               do k = 1, 5
                  ipx = ipx + 1
                  absc = hlgth * x2(k)
                  call f(n,centr+absc,fval1)
                  call f(n,centr-absc,fval2)
                  fval = fval1 + fval2
                  res21 = res21 + w21b(k)*fval
                  savfun(ipx,:) = fval(:)
                  fv3(k,:) = fval1(:)
                  fv4(k,:) = fval2(:)
               enddo
               resultf = res21 * hlgth
               abserr = maxval(cdabs((res21-res10)*hlgth))
            elseif ( l == 2 ) then
               res43 = w43b(12)*fcentr
               neval = 43
               do k = 1, 10
                  res43 = res43 + savfun(k,:) * w43a(k)
               enddo
               do k = 1, 11
                  ipx = ipx + 1
                  absc = hlgth * x3(k)
                  call f(n,centr+absc,fval1)
                  call f(n,centr-absc,fval2)
                  fval=fval1+fval2
                  res43 = res43 + fval * w43b(k)
                  savfun(ipx,:) = fval(:)
               enddo
               resultf = res43 * hlgth
               abserr = maxval(cdabs((res43-res21)*hlgth))
            elseif ( l == 3 ) then
               res87 = w87b(23) * fcentr
               neval = 87
               do k = 1, 21
                  res87 = res87+savfun(k,:)*w87a(k)
               enddo
               do k = 1, 22
                  absc = hlgth * x4(k)
                  call f(n,centr+absc,fval1)
                  call f(n,centr-absc,fval2)
                  res87 = res87 + w87b(k)*(fval1+fval2)
               enddo
               resultf = res87 * hlgth
               abserr = maxval(cdabs((res87-res43)*hlgth))
            endif
            if ( abserr <= max(epsabs,epsrel*maxval(cdabs(resultf)))) then
               ier = 0
            endif
            if ( ier == 0 ) then
               exit
            endif
         enddo
         end subroutine qng

         recursive subroutine gkintegrate(ntot,t0,t1,qsub,qint,subdiv, &
            errorcodes,inteps,mindiv,maxnumdiv)
         implicit none
         integer, intent(in) :: ntot
         integer, intent(inout) :: subdiv,errorcodes
         integer :: nsteps,ier,subdiv1,subdiv2,maxnumdiv,ec1,ec2
         real(8), intent(in) :: t0,t1
         real(8) :: t00,tmid,t11,errstep,inteps,mindiv
         complex(8), intent(out) :: qint(ntot)
         complex(8) :: qint1(ntot),qint2(ntot)
         external :: qsub

         errorcodes=0
         call qng (ntot,t0,t1,inteps,inteps,qsub,qint,errstep,nsteps,ier)
         if(abs(t1-t0).lt.mindiv) then
            errorcodes=2
!            write(*,'('' min delta: subdiv,t0,t1,err,steps:'',i5,5es12.4,i4)') &
!                subdiv,t0,t1,t1-t0,maxval(cdabs(qint)),errstep,nsteps
            return
         endif
         if(ier.ne.0) then
            if(subdiv.ge.maxnumdiv) then
               errorcodes=1
!               write(*,'('' max sub: subdiv,t0,t1,err,steps:'',i5,5es12.4,i4)') &
!                  subdiv,t0,t1,t1-t0,maxval(cdabs(qint)),errstep,nsteps
            else
               subdiv=subdiv+1
               subdiv1=subdiv
               subdiv2=subdiv
               t00=t0
               tmid=(t0+t1)*0.5d0
               t11=t1
               call gkintegrate(ntot,t00,tmid,qsub,qint1,subdiv1,ec1,inteps,mindiv,maxnumdiv)
               call gkintegrate(ntot,tmid,t11,qsub,qint2,subdiv2,ec2,inteps,mindiv,maxnumdiv)
               subdiv=max(subdiv1,subdiv2)
               errorcodes=max(ec1,ec2)
               qint=qint1+qint2
            endif
         endif
         return
         end subroutine gkintegrate

         subroutine gauleg(x1,x2,x,w,n)
         implicit none
         integer :: n,m,j,i
         real(8) :: x1,x2,x(n),w(n),xm,xl,z,p1,p2,p3,pp,z1,dj
         real(8), parameter :: eps=3.d-14
         m=(n+1)/2
         xm=0.5d0*(x2+x1)
         xl=0.5d0*(x2-x1)
         do i=1,m
            z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
            z1=z-1.d0
            do while(abs(z-z1).gt.eps)
               p1=1.d0
               p2=0.d0
               do j=1,n
                  dj=dble(j)
                  p3=p2
                  p2=p1
                  p1=((2.d0*dj-1.d0)*z*p2-(dj-1.d0)*p3)/dj
               enddo
               pp=n*(z*p1-p2)/(z*z-1.d0)
               z1=z
               z=z1-p1/pp
           enddo
           x(i)=xm-xl*z
           x(n+1-i)=xm+xl*z
           w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
           w(n+1-i)=w(i)
         enddo
         return
         end subroutine gauleg

         subroutine realsort(nlimits0,limits,eps,nlimits)
         implicit none
         integer :: nlimits0,nlimits,imin(1),n
         real(8) :: limits(1:nlimits0),rtemp(1:nlimits0),eps
         rtemp(1:nlimits0)=limits(1:nlimits0)
         imin=minloc(rtemp(1:nlimits0))
         nlimits=1
         limits(nlimits)=rtemp(imin(1))
         rtemp(imin(1))=1.d10
         do n=2,nlimits0
            imin=minloc(rtemp(1:nlimits0))
            if(rtemp(imin(1))-limits(nlimits).gt.eps) then
               nlimits=nlimits+1
               limits(nlimits)=rtemp(imin(1))
            endif
            rtemp(imin(1))=1.d10
         enddo
         end subroutine realsort

      end module specialfuncs

      module surface_subroutines
      use numconstants
      use specialfuncs
      implicit none
      logical :: source_sum,include_direct_source, &
         pole_integration,plane_surface_present
      integer, parameter :: max_number_plane_boundaries=10,max_singular_points=100
      integer :: source_order,target_order,max_azimuth_mode,number_layers,source_layer,target_layer, &
         source_layer2,number_limits,max_gf_iterations,number_gf_iterations,source_order2, &
         error_codes(4),energy_kernel_region,number_singular_points,singular_point_polarization(max_singular_points)
      integer, target :: number_plane_boundaries,maximum_integration_subdivisions
      real(8) :: source_z,target_z,azimuth_angle,incident_field_boundary, &
                 radial_distance,max_s,gf_error_epsilon, &
                 s_sc1,s_sc2,max_gf,max_bf,max_pi,max_picon, &
                 top_boundary, bot_boundary,source_z2,integration_error,singular_points(max_singular_points), &
                 incident_field_scale(2),pole_integration_radius,singular_gf_value(max_singular_points), &
                 incident_lateral_vector(2),g_cut,g_sing_mag
      real(8), target :: layer_thickness(max_number_plane_boundaries),integration_limit_epsilon, &
         integration_error_epsilon,real_axis_integration_limit,gf_switch_factor,s_scale_constant, &
         minimum_integration_spacing,minimum_initial_segment_size
      real(8), allocatable :: plane_boundary_position(:),real_axis_limits(:)
      complex(8) :: source_ri,target_ri,pole_integration_s
      complex(8), target :: layer_ref_index(0:max_number_plane_boundaries)
      complex(8), allocatable :: source_coefficient(:,:),source_coefficient_1(:,:,:),source_coefficient_2(:,:,:)
      data real_axis_integration_limit,minimum_integration_spacing,minimum_initial_segment_size/100.d0,1.d-5,0.d0/
      data include_direct_source/.false./
      data max_gf_iterations,gf_switch_factor,number_gf_iterations,gf_error_epsilon/1000,.5d0,0,1.d-7/
      data s_scale_constant/0.01d0/
      data layer_thickness/max_number_plane_boundaries*1.d0/
      data layer_ref_index(0)/(1.d0,0.d0)/
      data number_plane_boundaries/0/
      data integration_limit_epsilon,integration_error_epsilon,maximum_integration_subdivisions/1.d-15,1.d-6,25/
      data layer_ref_index(1:max_number_plane_boundaries)/max_number_plane_boundaries*(1.d0,0.d0)/
      data pole_integration_radius/1.d-6/
      data g_cut,g_sing_mag/10.d0,1.d5/

      contains
         subroutine plane_boundary_initialization()
         implicit none
         integer :: i
         real(8) :: smax
         if(allocated(plane_boundary_position)) deallocate(plane_boundary_position)
         allocate(plane_boundary_position(1:max(number_plane_boundaries,1)))
         plane_boundary_position(1)=0.d0
         do i=1,number_plane_boundaries-1
            plane_boundary_position(i+1)=plane_boundary_position(i)+layer_thickness(i)
         enddo
         smax=maxval(dble(layer_ref_index(0:number_plane_boundaries)))
         top_boundary=plane_boundary_position(max(1,number_plane_boundaries))+1.d-8
         bot_boundary=-1.d-8
         if(number_plane_boundaries.gt.1) then
            call gfunc_sing_points(bot_boundary,top_boundary,g_cut,smax, &
               number_singular_points,singular_points,singular_point_polarization, &
               singular_gf_value)
         else
            number_singular_points=0
            singular_gf_value=1.d0
         endif
         end subroutine plane_boundary_initialization

         subroutine gfunc_sing_points(sourcez,targetz,gcut,smax,nsingpoints,singpoints,singpol,singval)
         implicit none
         integer :: p,i,nbrack,nsingpoints,singpol(*)
         real(8) :: sourcez,targetz,gcut,smax,singpoints(*),sbrack(2,max_singular_points),s0,s1,s2,fmax,sfmax,singval(*)
         nsingpoints=0
         do p=1,2
            call sing_point_bracket(sourcez,targetz,p,gcut,smax,nbrack,sbrack)
            do i=1,nbrack
               s0=sbrack(1,i)
               s2=sbrack(2,i)
               s1=0.5d0*(s0+s2)
               call maxgfunc(sourcez,targetz,p,s0,s1,s2,1.d-9,100,fmax,sfmax)
               if(fmax.lt.g_sing_mag) cycle
               nsingpoints=nsingpoints+1
               singpoints(nsingpoints)=sfmax
               singpol(nsingpoints)=p
               singval(nsingpoints)=fmax
            enddo
         enddo
         end subroutine gfunc_sing_points

         subroutine sing_point_bracket(sz,tz,p,gcut,smax,nbrack,sbrack)
         implicit none
         logical :: inbrack
         integer :: nbrack,p
         real(8) :: sz,tz,smax,sbrack(2,*),dels,fm,gcut
         complex(8) :: gf(2,2,2),skz,tkz,s
         data dels/1.d-3/
         inbrack=.false.
         nbrack=0
         s=dels*0.5d0
         do while (dble(s).lt.smax)
            call layer_gf(s,sz,tz,gf,skz,tkz)
            fm=dble(sum(cdabs(gf(:,:,p))))
            if(fm.gt.gcut) then
               if(.not.inbrack) then
                  nbrack=nbrack+1
                  sbrack(1,nbrack)=dble(s)
                  inbrack=.true.
               endif
            else
               if(inbrack) then
                  inbrack=.false.
                  sbrack(2,nbrack)=dble(s)
                  if(nbrack.eq.max_singular_points) then
                     write(*,'('' max number GF singular points exceeded'')')
                     exit
                  endif
               endif
            endif
            s=s+dels
         enddo
         end subroutine sing_point_bracket

         subroutine maxgfunc(sz,tz,p,ax,bx,cx,tol,maxsteps,gmax,xmax)
         implicit none
         integer :: n,maxsteps,p
         real(8) :: ax,bx,cx,tol,xmax,r,c,x0,x3,x1,x2,f1,f2,f3,f0,gmax,sz,tz
         complex(8) :: gf(2,2,2),skz,tkz,s
         data r,c/.61803399d0,.38196602d0/
         x0=ax
         x3=cx
         if(abs(cx-bx).gt.abs(bx-ax))then
            x1=bx
            x2=bx+c*(cx-bx)
         else
            x2=bx
            x1=bx-c*(bx-ax)
         endif
         s=x1
         call layer_gf(s,sz,tz,gf,skz,tkz)
         f1=dble(sum(cdabs(gf(:,:,p))))
         s=x2
         call layer_gf(s,sz,tz,gf,skz,tkz)
         f2=dble(sum(cdabs(gf(:,:,p))))
         n=1
         do while(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)).and.n.le.maxsteps)
            n=n+1
            if(f2.gt.f1)then
               x0=x1
               x1=x2
               x2=r*x1+c*x3
               f0=f1
               f1=f2
               s=x2
               call layer_gf(s,sz,tz,gf,skz,tkz)
               f2=dble(sum(cdabs(gf(:,:,p))))
            else
               x3=x2
               x2=x1
               x1=r*x2+c*x0
               f3=f2
               f2=f1
               s=x1
               call layer_gf(s,sz,tz,gf,skz,tkz)
               f1=dble(sum(cdabs(gf(:,:,p))))
            endif
         enddo
         if(f1.gt.f2)then
            gmax=f1
            xmax=x1
         else
            gmax=f2
            xmax=x2
         endif
         end subroutine maxgfunc

         subroutine layer_gfos(sourcelayer,targetlayer,sourcez,targetz,kz,omega,tm,gfs)
         implicit none
         integer :: n,sourcelayer,targetlayer,p,sourcedir,iter
         real(8) :: scale,scale0,sourcez,targetz
         complex(8) :: gf(2,number_plane_boundaries,2),kz(0:number_plane_boundaries),omega(0:number_plane_boundaries), &
            tm(2,2,0:number_plane_boundaries+1,2),gfs(2,2,2),gp(2,number_plane_boundaries,2),gftot(2,number_plane_boundaries,2), &
            tfup,tfdn
         if(targetlayer.gt.0) then
            tfup=cdexp((0.d0,1.d0)*layer_ref_index(targetlayer)*kz(targetlayer) &
               *(targetz-plane_boundary_position(targetlayer)))
         else
            tfup=0.d0
         endif
         if(targetlayer.lt.number_plane_boundaries) then
            tfdn=cdexp((0.d0,1.d0)*layer_ref_index(targetlayer)*kz(targetlayer) &
               *(plane_boundary_position(targetlayer+1)-targetz))
         else
            tfdn=0.d0
         endif
         do sourcedir=1,2
            gfs(:,sourcedir,:)=0.d0
            gf=0.d0
            if(sourcedir.eq.1) then
               if(sourcelayer.lt.number_plane_boundaries) then
                  gf(:,sourcelayer+1,:)=tm(:,1,sourcelayer+1,:) &
                     *cdexp((0.d0,1.d0)*layer_ref_index(sourcelayer)*kz(sourcelayer) &
                     *(plane_boundary_position(sourcelayer+1)-sourcez))
               else
                  cycle
               endif
            endif
            if(sourcedir.eq.2) then
               if(sourcelayer.gt.0) then
                  gf(:,sourcelayer,:)=tm(:,2,sourcelayer,:) &
                     *cdexp((0.d0,1.d0)*layer_ref_index(sourcelayer)*kz(sourcelayer) &
                     *(sourcez-plane_boundary_position(sourcelayer)))
               else
                  cycle
               endif
            endif
            gftot=gf
            if(targetlayer.gt.0) then
               gfs(1,sourcedir,:)=gftot(1,targetlayer,:)*tfup
            endif
            if(targetlayer.lt.number_plane_boundaries) then
               gfs(2,sourcedir,:)=gftot(2,targetlayer+1,:)*tfdn
            endif
            scale0=sqrt(sum(cdabs(gftot(:,:,:))**2))
            do iter=1,max_gf_iterations
               gp=0
               do p=1,2
                  do n=2,number_plane_boundaries
                     gp(:,n,p)=omega(n-1)*tm(:,1,n,p)*gf(1,n-1,p)
                  enddo
                  do n=1,number_plane_boundaries-1
                     gp(:,n,p)=gp(:,n,p)+omega(n)*tm(:,2,n,p)*gf(2,n+1,p)
                  enddo
               enddo
               gftot(:,:,:)=gftot(:,:,:)+gp(:,:,:)
               gf(:,:,:)=gp(:,:,:)
               if(targetlayer.gt.0) then
                  gfs(1,sourcedir,:)=gftot(1,targetlayer,:)*tfup
               endif
               if(targetlayer.lt.number_plane_boundaries) then
                  gfs(2,sourcedir,:)=gftot(2,targetlayer+1,:)*tfdn
               endif
               scale=sqrt(sum(cdabs(gftot(:,:,:))**2))
               if(abs(scale-scale0)/max(scale0,1.d-12).lt.gf_error_epsilon) then
                  number_gf_iterations=max(number_gf_iterations,iter)
                  exit
               endif
               scale0=scale
            enddo
            number_gf_iterations=max(number_gf_iterations,iter)
         enddo
         if(number_gf_iterations.gt.max_gf_iterations) then
            error_codes(1)=1
         endif
         end subroutine layer_gfos

         subroutine layer_gfrec(sourcelayer,targetlayer,sourcez,targetz,kz,omega,tm,gfs)
         implicit none
         integer :: n,sourcelayer,targetlayer,p,sourcedir
         real(8) :: sourcez,targetz
         complex(8) :: gf(2,0:number_plane_boundaries+1),kz(0:number_plane_boundaries),omega(0:number_plane_boundaries), &
            tm(2,2,0:number_plane_boundaries+1,2),gfs(2,2,2),tfup,tfdn
!         complex(16) :: um(2,2,0:number_plane_boundaries),sm(2,2,0:number_plane_boundaries), &
!            sv(2,0:number_plane_boundaries+1),tempv(2),umt(2,2),bv(2,0:number_plane_boundaries+1),umt0(2,2)
         complex(8) :: um(2,2,0:number_plane_boundaries),sm(2,2,0:number_plane_boundaries), &
            sv(2,0:number_plane_boundaries+1),tempv(2),umt(2,2),bv(2,0:number_plane_boundaries+1),umt0(2,2)
         if(targetlayer.gt.0) then
            tfup=cdexp((0.d0,1.d0)*layer_ref_index(targetlayer)*kz(targetlayer) &
               *(targetz-plane_boundary_position(targetlayer)))
         else
            tfup=0.d0
         endif
         if(targetlayer.lt.number_plane_boundaries) then
            tfdn=cdexp((0.d0,1.d0)*layer_ref_index(targetlayer)*kz(targetlayer) &
               *(plane_boundary_position(targetlayer+1)-targetz))
         else
            tfdn=0.d0
         endif
         do p=1,2
            do n=0,number_plane_boundaries
               um(1,1,n)=omega(n)*tm(1,1,n,p)
               um(1,2,n)=omega(n)*tm(1,2,n,p)
               um(2,1,n)=-omega(n)*tm(2,1,n+1,p)*tm(1,1,n,p)/tm(2,2,n+1,p)
               um(2,2,n)=(1.d0-omega(n)*omega(n)*tm(2,1,n+1,p)*tm(1,2,n,p))/(omega(n)*tm(2,2,n+1,p))
               sm(1,1,n)=1.d0
               sm(1,2,n)=0.d0
               sm(2,1,n)=-tm(2,1,n+1,p)/tm(2,2,n+1,p)
               sm(2,2,n)=-1.d0/(omega(n)*tm(2,2,n+1,p))
            enddo
            do sourcedir=1,2
               sv=0.d0
               gf(:,:)=0.d0
               gfs(:,sourcedir,p)=0.d0
               if(sourcedir.eq.1) then
                  if(sourcelayer.eq.number_plane_boundaries) then
                     cycle
                  else
                     sv(1,sourcelayer+1)=cdexp((0.d0,1.d0)*layer_ref_index(sourcelayer)*kz(sourcelayer) &
                        *(plane_boundary_position(sourcelayer+1)-sourcez))
                  endif
               elseif(sourcedir.eq.2) then
                  if(sourcelayer.eq.0) then
                     cycle
                  else
                     sv(2,sourcelayer)=cdexp((0.d0,1.d0)*layer_ref_index(sourcelayer)*kz(sourcelayer) &
                        *(sourcez-plane_boundary_position(sourcelayer)))
                  endif
               endif
               tempv(1)=sm(1,1,sourcelayer)*sv(1,sourcelayer+1)
               tempv(2)=sm(2,1,sourcelayer)*sv(1,sourcelayer+1)+sm(2,2,sourcelayer)*sv(2,sourcelayer)
               bv(:,sourcelayer+1)=tempv(:)
               do n=sourcelayer+1,number_plane_boundaries
                  bv(1,n+1)=um(1,1,n)*bv(1,n)+um(1,2,n)*bv(2,n)
                  bv(2,n+1)=um(2,1,n)*bv(1,n)+um(2,2,n)*bv(2,n)
               enddo
               umt(:,:)=um(:,:,0)
               do n=1,number_plane_boundaries
                  umt0=umt
                  umt(1,1)=um(1,1,n)*umt0(1,1)+um(1,2,n)*umt0(2,1)
                  umt(2,1)=um(2,1,n)*umt0(1,1)+um(2,2,n)*umt0(2,1)
                  umt(1,2)=um(1,1,n)*umt0(1,2)+um(1,2,n)*umt0(2,2)
                  umt(2,2)=um(2,1,n)*umt0(1,2)+um(2,2,n)*umt0(2,2)
               enddo
               gf(2,0)=-bv(2,number_plane_boundaries+1)/umt(2,2)
               do n=0,number_plane_boundaries
                  gf(1,n+1)=um(1,1,n)*gf(1,n)+um(1,2,n)*gf(2,n)+sm(1,1,n)*sv(1,n+1)
                  gf(2,n+1)=um(2,1,n)*gf(1,n)+um(2,2,n)*gf(2,n)+sm(2,1,n)*sv(1,n+1)+sm(2,2,n)*sv(2,n)
               enddo
               if(targetlayer.gt.0) then
                  gfs(1,sourcedir,p)=tfup*sum(tm(1,:,targetlayer,p)*gf(:,targetlayer))
               endif
               if(targetlayer.lt.number_plane_boundaries) then
                  gfs(2,sourcedir,p)=tfdn*sum(tm(2,:,targetlayer+1,p)*gf(:,targetlayer+1))
               endif
            enddo
         enddo
         end subroutine layer_gfrec

         subroutine layer_gf(s,sourcez,targetz,gfs,sourcekz,targetkz,include_direct)
         implicit none
         logical :: incdir,prop
         logical, optional :: include_direct
         integer :: n,sourcelayer,targetlayer,i
         real(8) :: sourcez,targetz,c
         complex(8) :: kz(0:number_plane_boundaries),omega(0:number_plane_boundaries), &
            den,tm(2,2,0:number_plane_boundaries+1,2), &
            gfs(2,2,2), &
            sourcekz,targetkz,s
         if(present(include_direct)) then
            incdir=include_direct
         else
            incdir=.false.
         endif
         sourcelayer=layer_id(sourcez)
         targetlayer=layer_id(targetz)
         do n=0,number_plane_boundaries
            kz(n)=cdsqrt((layer_ref_index(n)-s)*(layer_ref_index(n)+s))/layer_ref_index(n)
         enddo
         if(incdir) then
            if(sourcez.eq.targetz) then
               c=0.5d0
            else
               c=1.d0
            endif
         endif
         prop=dble(s).le.dble(layer_ref_index(targetlayer))
         if(number_plane_boundaries.eq.0) then
            sourcekz=kz(0)
            targetkz=kz(0)
            gfs=0.d0
            if(incdir) then
               if(targetz.le.sourcez) then
                  gfs(2,2,:)=c*cdexp((0.d0,1.d0)*layer_ref_index(targetlayer)*sourcekz*(abs(sourcez-targetz)))
               endif
               if(targetz.ge.sourcez) then
                  gfs(1,1,:)=c*cdexp((0.d0,1.d0)*layer_ref_index(targetlayer)*sourcekz*(abs(sourcez-targetz)))
               endif
            endif
            return
         endif
         omega(0)=1.d0
         omega(number_plane_boundaries)=1.d0
         do n=1,number_plane_boundaries-1
            omega(n)=cdexp((0.d0,1.d0)*layer_ref_index(n)*kz(n)*layer_thickness(n))
         enddo
         do n=1,number_plane_boundaries
            den=kz(n)*layer_ref_index(n-1)+kz(n-1)*layer_ref_index(n)
            tm(1,1,n,1)=2.d0*kz(n-1)*layer_ref_index(n-1)/den
            tm(1,2,n,1)=(kz(n)*layer_ref_index(n-1)-kz(n-1)*layer_ref_index(n))/den
            tm(2,1,n,1)=-tm(1,2,n,1)
            tm(2,2,n,1)=2.d0*kz(n)*layer_ref_index(n)/den
            den=kz(n-1)*layer_ref_index(n-1)+kz(n)*layer_ref_index(n)
            tm(1,1,n,2)=2.d0*kz(n-1)*layer_ref_index(n-1)/den
            tm(1,2,n,2)=-(kz(n-1)*layer_ref_index(n-1)-kz(n)*layer_ref_index(n))/den
            tm(2,1,n,2)=-tm(1,2,n,2)
            tm(2,2,n,2)=2.d0*kz(n)*layer_ref_index(n)/den
         enddo
         sourcekz=kz(sourcelayer)
         targetkz=kz(targetlayer)
         gfs=0.d0
         if(number_plane_boundaries.eq.1) then
            gfs=0.d0
            if(sourcelayer.eq.0.and.targetlayer.eq.0) then
               gfs(2,1,:)=cdexp((0.d0,1.d0)*layer_ref_index(0)*sourcekz*(abs(sourcez)+abs(targetz))) &
                  *tm(2,1,1,:)
            elseif(sourcelayer.eq.0.and.targetlayer.eq.1) then
               gfs(1,1,:)=cdexp((0.d0,1.d0)*(layer_ref_index(0)*sourcekz*abs(sourcez) &
                  +layer_ref_index(1)*targetkz*abs(targetz)))*tm(1,1,1,:)
            elseif(sourcelayer.eq.1.and.targetlayer.eq.1) then
               gfs(1,2,:)=cdexp((0.d0,1.d0)*layer_ref_index(1)*sourcekz*(abs(sourcez)+abs(targetz))) &
                  *tm(1,2,1,:)
            elseif(sourcelayer.eq.1.and.targetlayer.eq.0) then
               gfs(2,2,:)=cdexp((0.d0,1.d0)*(layer_ref_index(1)*sourcekz*abs(sourcez) &
                  +layer_ref_index(0)*targetkz*abs(targetz)))*tm(2,2,1,:)
            endif
         else
            do i=1,2
               tm(i,i,0,:)=1.d0
               tm(i,3-i,0,:)=0.d0
               tm(i,i,number_plane_boundaries+1,:)=1.d0
               tm(i,3-i,number_plane_boundaries+1,:)=0.d0
            enddo
            if(maxval(cdabs(omega(1:number_plane_boundaries-1))).lt.gf_switch_factor) then
               number_gf_iterations=0
               call layer_gfos(sourcelayer,targetlayer,sourcez,targetz,kz,omega,tm,gfs)
            else
               call layer_gfrec(sourcelayer,targetlayer,sourcez,targetz,kz,omega,tm,gfs)
            endif
         endif
         if(incdir.and.sourcelayer.eq.targetlayer) then
            if(targetz.le.sourcez) then
               gfs(2,2,:)=gfs(2,2,:)+c*cdexp((0.d0,1.d0)*layer_ref_index(targetlayer)*sourcekz*(abs(sourcez-targetz)))
            endif
            if(targetz.ge.sourcez) then
               gfs(1,1,:)=gfs(1,1,:) &
                  +c*cdexp((0.d0,1.d0)*layer_ref_index(targetlayer)*sourcekz*(abs(sourcez-targetz)))
            endif
         endif
         end subroutine layer_gf

         subroutine realaxiskernel(ntot,t,kernmat)
         implicit none
         integer :: m,m1,k,k1,mk,n,l,p,q,mn,kl,i,pol,dir,tsign,ssign,mmax,ntot
         integer, save :: count
         real(8) :: t
         complex(8) :: kernmat(ntot),temp,bfunc(0:source_order+target_order), &
            s,const,const2,sr,dsdt, &
            pivec(2,source_order*(source_order+2),2),picvec(2,target_order*(target_order+2),2), &
            bvec(2,2),sourcekz,targetkz,gfunc(2,2,2),sourceri,targetri,sources,targets
         data count/0/
         if(pole_integration) then
            dsdt=(0.d0,1.d0)*pole_integration_radius*cdexp((0.d0,1.d0)*t)
            s=pole_integration_s+pole_integration_radius*cdexp((0.d0,1.d0)*t)
         else
            s=t*dcmplx(1.d0,-s_sc1)+dcmplx(0.d0,-1.d0)*s_sc2
            dsdt=dcmplx(1.d0,-s_sc1)
         endif
         call layer_gf(s,source_z,target_z,gfunc,sourcekz,targetkz,include_direct_source)
         sourceri=layer_ref_index(source_layer)
         targetri=layer_ref_index(target_layer)
         sr=s*radial_distance
         sources=s/sourceri
         targets=s/targetri
         call complexpivec(sourcekz,source_order,pivec,1)
         call complexpivec(targetkz,target_order,picvec,-1)
         bfunc=0.d0
         mmax=target_order+source_order
         if(radial_distance.eq.0.d0) then
            bfunc(0)=1.d0
         else
            call bessel_integer_complex(target_order+source_order,sr,mmax,bfunc)
         endif
!         max_gf=max(max_gf,maxval(cdabs(gfunc)))
!         max_pi=max(max_pi,maxval(cdabs(pivec)))
!         max_picon=max(max_picon,maxval(cdabs(picvec)))
!         max_bf=max(max_bf,maxval(cdabs(bfunc)))
         const=4.d0*pi*s/sourcekz/sourceri/sourceri*dsdt
         i=0
         if(source_sum) then
            do m=-target_order,target_order
               bvec=0.d0
               m1=max(1,abs(m))
               do k=-source_order,source_order
                  k1=max(1,abs(k))
                  mk=abs(k-m)
                  if((radial_distance.eq.0.d0).and.mk.ne.0) cycle
                  if(mk.gt.mmax) cycle
                  const2=cdexp((0.d0,1.d0)*(k-m)*azimuth_angle)*((0.d0,1.d0)**mk)*bfunc(mk)*const
                  do l=k1,source_order
                     kl=l*(l+1)+k
                     do q=1,2
                        ssign=(-1)**(k+l+q-1)
                        do pol=1,2
                           do dir=1,2
                              bvec(dir,pol)=bvec(dir,pol)+const2*pivec(q,kl,pol)*(gfunc(dir,1,pol) &
                                 +ssign*(-1)**pol*gfunc(dir,2,pol))*source_coefficient(q,kl)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               do n=m1,target_order
                  mn=n*(n+1)+m
                  do p=1,2
                     tsign=(-1)**(m+n+p-1)
                     i=i+1
                     kernmat(i)=picvec(p,mn,1)*(bvec(1,1)-tsign*bvec(2,1)) &
                        +picvec(p,mn,2)*(bvec(1,2)+tsign*bvec(2,2))
                  enddo
               enddo
            enddo
         else
            do m=0,target_order
               m1=max(1,abs(m))
               do k=-source_order,source_order
                  k1=max(1,abs(k))
                  mk=abs(k-m)
                  if((radial_distance.eq.0.d0).and.mk.ne.0) then
                     cycle
                  endif
                  const2=((0.d0,1.d0)**mk)*bfunc(mk)*const
                  do n=m1,target_order
                     mn=m+n*(n+1)
                     do l=k1,source_order
                        kl=k+l*(l+1)
                        do p=1,2
                           tsign=(-1)**(m+n+p-1)
                           do q=1,2
                              ssign=(-1)**(k+l+q-1)
                              i=i+1
                              temp=0.d0
                              do pol=1,2
                                 temp=temp+const2*picvec(p,mn,pol)*pivec(q,kl,pol) &
                                    *(gfunc(1,1,pol)+(-1)**pol*(tsign*gfunc(2,1,pol)+ssign*gfunc(1,2,pol)) &
                                    +tsign*ssign*gfunc(2,2,pol))
                              enddo
                              kernmat(i)=temp
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         endif
         end subroutine realaxiskernel

         subroutine energykernel(ntot,t,kernmat)
         implicit none
         integer :: m,m1,k,k1,mk,n,l,p,q,mn,kl,pol,tsign,ssign,mmax,sourceorder(2), &
            nmax,spol,tlay,ntot
         integer, save :: count
         real(8) :: t,sourcez(2),targetz,const3(3,2)
         complex(8) :: kernmat(ntot),bfunc(0:source_order+source_order2), &
            s,const,const2,sr,cscale, &
            pivec1(2,source_order*(source_order+2),2), &
            pivec2(2,source_order2*(source_order2+2),2),bvec(2,2), &
            bvec1(2),bvec2(2),sourcekz(2),targetkz,gfunc1(2,2,2),gfunc2(2,2,2),sourceri(2),targetri
         data count/0/
         sourcez(1)=source_z
         sourcez(2)=source_z2
         targetz=target_z
         sourceri(1)=layer_ref_index(source_layer)
         sourceri(2)=layer_ref_index(source_layer2)
         sourceorder(1)=source_order
         sourceorder(2)=source_order2
         tlay=target_layer
         targetri=layer_ref_index(tlay)
         nmax=max(source_order,source_order2)
         if(energy_kernel_region.eq.0) then
            s=cdsqrt(((1.d0,0.d0)-t)*((1.d0,0.d0)+t))*dble(targetri)
         else
            s=cdsqrt((1.d0,0.d0)+t*t)*dble(targetri)
         endif
         call layer_gf(s,sourcez(1),targetz,gfunc1,sourcekz(1),targetkz,include_direct=.true.)
         call layer_gf(s,sourcez(2),targetz,gfunc2,sourcekz(2),targetkz,include_direct=.true.)
         cscale=t*dble(targetri)**2
         call complexpivec(sourcekz(1),sourceorder(1),pivec1,1)
         call complexpivec(sourcekz(2),sourceorder(2),pivec2,1)
         bfunc=0.d0
         sr=s*radial_distance
         mmax=sum(sourceorder(:))
         if(radial_distance.eq.0.d0) then
            bfunc(0)=1.d0
         else
            call bessel_integer_complex(sum(sourceorder(:)),sr,mmax,bfunc)
         endif
         kernmat=0.d0
         const=4.d0*pi*cscale/sourcekz(1)/sourceri(1)**2/dconjg(sourcekz(2)*sourceri(2)**2)
         const3(1,1)=(cdabs(targetkz)**2+s*s/cdabs(targetri)**2)*dble(targetri*targetkz)
         const3(2,1)=-const3(1,1)
         const3(3,1)=(cdabs(targetkz)**2-s*s/cdabs(targetri)**2)*dimag(targetri*targetkz)
         const3(1,2)=dble(targetri*targetkz)
         const3(2,2)=-const3(1,2)
         const3(3,2)=-dimag(targetri*targetkz)
         do spol=1,2
            do pol=1,2
               bvec=0.d0
               do m=-sourceorder(2),sourceorder(2)
                  bvec1=0.d0
                  bvec2=0.d0
                  m1=max(1,abs(m))
                  do k=-sourceorder(1),sourceorder(1)
                     k1=max(1,abs(k))
                     mk=abs(k-m)
                     if((radial_distance.eq.0.d0).and.mk.ne.0) cycle
                     if(mk.gt.mmax) cycle
                     const2=cdexp((0.d0,1.d0)*(k-m)*azimuth_angle)*((0.d0,1.d0)**mk)*bfunc(mk)*const
                     do l=k1,sourceorder(1)
                        kl=l*(l+1)+k
                        do q=1,2
                           ssign=(-1)**(k+l+q-1+pol)
                           bvec1(:)=bvec1(:)+const2*pivec1(q,kl,pol)*(gfunc1(:,1,pol) &
                              +ssign*gfunc1(:,2,pol))*source_coefficient_1(q,kl,spol)
                        enddo
                     enddo
                  enddo
                  do n=m1,sourceorder(2)
                     mn=n*(n+1)+m
                     do p=1,2
                        tsign=(-1)**(m+n+p-1+pol)
                        bvec2(:)=bvec2(:)+dconjg(pivec2(p,mn,pol)*(gfunc2(:,1,pol) &
                           +tsign*gfunc2(:,2,pol)) &
                           *source_coefficient_2(p,mn,spol))
                     enddo
                  enddo
                  do p=1,2
                     do q=1,2
                        bvec(p,q)=bvec(p,q)+bvec1(p)*bvec2(q)
                     enddo
                  enddo
               enddo
               kernmat(spol)=kernmat(spol)+const3(1,pol)*dble(bvec(1,1))+const3(2,pol)*dble(bvec(2,2)) &
                  +const3(3,pol)*dimag(bvec(1,2)-bvec(2,1))
            enddo
         enddo
         end subroutine energykernel

         real(8) function mnorm(n,m)
         implicit none
         integer :: n
         complex(8) :: m(n)
         mnorm=sqrt(dble(dot_product(m,m)))
         end function mnorm

         integer function layer_id(z)
         implicit none
         integer :: n
         real(8) :: z
         layer_id=0
         do n=1,number_plane_boundaries
            if(z.ge.plane_boundary_position(n)) then
               layer_id=n
            else
               return
            endif
         enddo
         end function layer_id

         subroutine plane_interaction(ntot,ltot,x,y,sourcez,targetz, &
            interactionmatrix,source_vector, &
            index_model,lr_transformation,make_symmetric,propagating_directions_only)
         implicit none
         logical :: lrtran,makesymmetric,propdir
         logical, optional :: lr_transformation,make_symmetric,propagating_directions_only
         integer :: ntot,ltot,l,k,n,m,p,q,mn,kl,tranmat(2,2),indexmodel, &
            mnp,qtot,m1,k1,i,rmatdim,rmataddress
         integer, optional :: index_model
         real(8) :: x,y,rho,sourcez,targetz,ssc,rail
         complex(8) :: ephi,ephim,interactionmatrix(*), &
            ctemp(2),c2temp(2,2),c2tempm(2,2)
         complex(8), allocatable :: rmat(:)
         complex(8), optional :: source_vector(2*ltot*(ltot+2))
         tranmat=reshape((/1,1,1,-1/),(/2,2/))
         target_order=ntot
         source_order=ltot
         if(present(index_model)) then
            indexmodel=index_model
         else
            indexmodel=1
         endif
         if(present(lr_transformation)) then
            lrtran=lr_transformation
         else
            lrtran=.true.
         endif
         if(present(make_symmetric)) then
            makesymmetric=make_symmetric
         else
            makesymmetric=.false.
         endif
         if(present(propagating_directions_only)) then
            propdir=propagating_directions_only
         else
            propdir=.false.
         endif
         if(present(source_vector)) then
            source_sum=.true.
            allocate(source_coefficient(2,ltot*(ltot+2)))
            do n=1,ltot
               do m=-n,n
                  mn=n*(n+1)+m
                  do p=1,2
                     mnp=amnpaddress(m,n,p,ltot,indexmodel)
                     ctemp(p)=source_vector(mnp)
                  enddo
                  if(lrtran) ctemp=matmul(tranmat,ctemp)
                  source_coefficient(:,mn)=ctemp(:)
               enddo
            enddo
         else
            source_sum=.false.
         endif
         if(x.eq.0.d0.and.y.eq.0.d0) then
            rho=0.d0
            ephi=1.d0
            azimuth_angle=0.d0
         else
            rho=sqrt(x*x+y*y)
            ephi=dcmplx(x,y)/rho
            azimuth_angle=datan2(y,x)
         endif
         if(rho.le.0.0001d0) then
            max_azimuth_mode=0
            rho=0.d0
         else
            max_azimuth_mode=source_order+target_order
         endif
         source_z=sourcez
         target_z=targetz
         radial_distance=rho
         if(.not.source_sum) then
            if(makesymmetric) then
               max_azimuth_mode=0
               rmatdim=2*atcdim(ntot,ltot)
            else
               rmatdim=4*ntot*(ntot+2)*ltot*(ltot+2)
            endif
         endif
         source_layer=layer_id(sourcez)
         target_layer=layer_id(targetz)
         allocate(real_axis_limits(number_plane_boundaries+1))
         real_axis_limits(1:number_plane_boundaries+1)=cdabs(layer_ref_index(0:number_plane_boundaries))
         call realsort(number_plane_boundaries+1,real_axis_limits,1.d-10,number_limits)
         if(source_sum) then
            qtot=2*ntot*(ntot+2)
         else
            qtot=0
            do m=0,ntot
               m1=max(1,m)
               do k=-ltot,ltot
                  if(abs(m-k).gt.max_azimuth_mode) cycle
                  k1=max(abs(k),1)
                  do n=m1,ntot
                     do l=k1,ltot
                        qtot=qtot+4
                     enddo
                  enddo
               enddo
            enddo
         endif
         allocate(rmat(qtot))
         max_gf=0.d0
         max_bf=0.d0
         max_pi=0.d0
         max_picon=0.d0

         if(propdir) then
            ssc=s_scale_constant
            rail=real_axis_integration_limit
            s_scale_constant=0.d0
            real_axis_integration_limit=dble(layer_ref_index(source_layer))-1.d-6
         endif

         call refmatrixrealaxis(qtot,rmat)

         if(propdir) then
            s_scale_constant=ssc
            real_axis_integration_limit=rail
         endif

         if(source_sum) then
            interactionmatrix(1:2*ntot*(ntot+2))=0.d0
            i=0
            do m=-ntot,ntot
               m1=max(1,abs(m))
               do n=m1,ntot
                  mn=n*(n+1)+m
                  do p=1,2
                     i=i+1
                     ctemp(p)=rmat(i)
                  enddo
                  if(lrtran) ctemp=matmul(tranmat,ctemp)/2.d0
                  do p=1,2
                     mnp=amnpaddress(m,n,p,ntot,indexmodel)
                     interactionmatrix(mnp)=ctemp(p)
                  enddo
               enddo
            enddo
            deallocate(rmat,source_coefficient)
         else
            interactionmatrix(1:rmatdim)=0.d0
            i=0
            do m=0,ntot
               m1=max(1,m)
               do k=-ltot,ltot
                  if(abs(m-k).gt.max_azimuth_mode) cycle
                  k1=max(abs(k),1)
                  ephim=ephi**(k-m)
                  do n=m1,ntot
                     do l=k1,ltot
                        do p=1,2
                           do q=1,2
                              i=i+1
                              c2temp(p,q)=rmat(i)*ephim
                              if(m.ne.0) c2tempm(p,q)=rmat(i)/ephim*((-1)**(abs(m-k)+p+q))
                           enddo
                        enddo
                        if(lrtran) then
                           c2temp=matmul(tranmat,matmul(c2temp,tranmat))/2.d0
                           if(m.ne.0) c2tempm=matmul(tranmat,matmul(c2tempm,tranmat))/2.d0
                        endif
                        do q=1,2
                           do p=1,2
                              if(makesymmetric) then
                                 rmataddress=2*moffset(m,ntot,ltot)+p+2*(n-m1)+2*(ntot-m1+1)*(q-1+2*(l-m1))
                              else
                                 mn=amnpaddress(m,n,p,ntot,indexmodel)
                                 kl=amnpaddress(k,l,q,ltot,indexmodel)
                                 rmataddress=mn+(kl-1)*2*ntot*(ntot+2)
                              endif
                              interactionmatrix(rmataddress)=c2temp(p,q)
                              if(m.ne.0) then
                                 if(makesymmetric) then
                                    rmataddress=2*moffset(-m,ntot,ltot)+p+2*(n-m1)+2*(ntot-m1+1)*(q-1+2*(l-m1))
                                 else
                                    mn=amnpaddress(-m,n,p,ntot,indexmodel)
                                    kl=amnpaddress(-k,l,q,ltot,indexmodel)
                                    rmataddress=mn+(kl-1)*2*ntot*(ntot+2)
                                 endif
                                 interactionmatrix(rmataddress)=c2tempm(p,q)
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            deallocate(rmat)
         endif
         deallocate(real_axis_limits)
         end subroutine plane_interaction

         subroutine sphere_boundary_scattering(ntot1,rpos1,scoef1,ntot2,rpos2,scoef2, &
             targetz,qsca,lr_to_mode)
         implicit none
         logical :: lr2mode
         logical, optional :: lr_to_mode
         integer :: ntot1,ntot2,n,m,p,mnp,mn,tranmat(2,2),nlimits,subdiv,nlimits0,ec
         real(8) :: rpos1(3),rpos2(3),qsca(2), &
            xyvec(2),limits(1:number_plane_boundaries+max_singular_points),s0,s1,errstep,targetz,riscale
         complex(8) :: scoef1(2*ntot1*(ntot1+2),2),scoef2(2*ntot2*(ntot2+2),2),ctemp(2,2),qmat(2)
         tranmat=reshape((/1,1,1,-1/),(/2,2/))

         if(present(lr_to_mode)) then
            lr2mode=lr_to_mode
         else
            lr2mode=.true.
         endif
         allocate(source_coefficient_1(2,ntot1*(ntot1+2),2))
         allocate(source_coefficient_2(2,ntot2*(ntot2+2),2))
         do n=1,ntot1
            do m=-n,n
               mn=n*(n+1)+m
               do p=1,2
                  mnp=amnpaddress(m,n,p,ntot1,2)
                  ctemp(p,:)=scoef1(mnp,:)
               enddo
               if(lr2mode) ctemp=matmul(tranmat,ctemp)
               source_coefficient_1(:,mn,:)=ctemp(:,:)
            enddo
         enddo
         do n=1,ntot2
            do m=-n,n
               mn=n*(n+1)+m
               do p=1,2
                  mnp=amnpaddress(m,n,p,ntot2,2)
                  ctemp(p,:)=scoef2(mnp,:)
               enddo
               if(lr2mode) ctemp=matmul(tranmat,ctemp)
               source_coefficient_2(:,mn,:)=ctemp(:,:)
            enddo
         enddo
         source_layer=layer_id(rpos1(3))
         source_layer2=layer_id(rpos2(3))
         target_layer=layer_id(targetz)
         source_order=ntot1
         source_order2=ntot2
         source_z=rpos1(3)
         source_z2=rpos2(3)
         target_z=targetz
         xyvec=rpos2(1:2)-rpos1(1:2)
         radial_distance=sqrt(sum(xyvec*xyvec))
         if(radial_distance.eq.0.d0) then
            azimuth_angle=0.d0
         else
            azimuth_angle=datan2(xyvec(2),xyvec(1))
         endif
         riscale=dble(layer_ref_index(target_layer))
         nlimits0=1
         limits(1)=1.d0
         do n=1,number_plane_boundaries+1
            if(n-1.ne.target_layer.and.dble(layer_ref_index(n-1)).lt.riscale) then
               nlimits0=nlimits0+1
               limits(nlimits0)=sqrt((1.d0-dble(layer_ref_index(n-1))/riscale) &
                  *(1.d0+dble(layer_ref_index(n-1))/riscale))
            endif
         enddo
         do n=1,number_singular_points
            if(singular_points(n).lt.riscale) then
               nlimits0=nlimits0+1
               limits(nlimits0) = sqrt((1.d0-singular_points(n)/riscale) &
                  *(1.d0+singular_points(n)/riscale))
            endif
         enddo
         call realsort(nlimits0,limits,1.d-10,nlimits)
         qsca=0.d0
         integration_error=0.d0
         s1=0.d0
         energy_kernel_region=0
         do n=1,nlimits
            qmat=0.d0
            s0=s1
            s1=limits(n)
            subdiv=0
            ec=0
            call gkintegrate(2,s0,s1,energykernel,qmat,subdiv,ec, &
               integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
            if(ec.eq.1) error_codes(4)=1
            if(ec.eq.2) error_codes(3)=1
            qsca=qsca+qmat
         enddo
         nlimits0=0
         do n=1,number_plane_boundaries+1
            if(n-1.ne.target_layer.and.dble(layer_ref_index(n-1)).gt.riscale) then
               nlimits0=nlimits0+1
               limits(nlimits0)=sqrt(dble(layer_ref_index(n-1))**2/riscale**2+1.d0)
            endif
         enddo
         do n=1,number_singular_points
            if(singular_points(n).gt.riscale) then
               nlimits0=nlimits0+1
               limits(nlimits0) = sqrt(1.d0+singular_points(n)**2/riscale**2)
            endif
         enddo
         s1=0.d0
         energy_kernel_region=1
         if(nlimits0.gt.0) then
            call realsort(nlimits0,limits,1.d-10,nlimits)
            do n=1,nlimits
               qmat=0.d0
               s0=s1
               s1=limits(n)
               subdiv=0
               ec=0
               call gkintegrate(2,s0,s1,energykernel,qmat,subdiv,ec, &
                  integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
               if(ec.eq.1) error_codes(4)=1
               if(ec.eq.2) error_codes(3)=1
               qsca=qsca+qmat
            enddo
         endif
         errstep=1.d0
         do while(errstep.gt.integration_limit_epsilon.and.s1.lt.real_axis_integration_limit)
            s0=s1
            s1=s1+0.5d0
            qmat=0.d0
            subdiv=0
            ec=0
            call gkintegrate(2,s0,s1,energykernel,qmat,subdiv,ec, &
               integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
            if(ec.eq.1) error_codes(4)=1
            if(ec.eq.2) error_codes(3)=1
            qsca=qsca+qmat
            errstep=cdabs(sum(qmat))/abs(sum(qsca))
         enddo
         deallocate(source_coefficient_1,source_coefficient_2)
         end subroutine sphere_boundary_scattering

         subroutine refmatrixrealaxis(qtot,rmat)
         implicit none
         integer :: qtot,n,limit,nseg,seg,n0,subdiv,ec
         real(8) :: t1,t2,delt,dnorm,norm0,r, &
            deltseg,dt1,dt2,errlim,t1t,t2t
         complex(8) :: rmat(qtot),drmat(qtot)
         if(pole_integration) then
            t1=0.d0
            t2=2.d0*pi
            subdiv=0.d0
            rmat=0.d0
            ec=0
            call gkintegrate(qtot,t1,t2,realaxiskernel,rmat,subdiv,ec, &
                  integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
            return
         endif
         r=sqrt(sqrt(radial_distance**2+(abs(source_z)+abs(target_z))**2))
         delt=.5d0/r
         do n=1,number_limits
            if(abs(cdabs(layer_ref_index(source_layer))-real_axis_limits(n)).le.1.d-8) then
               n0=n
               exit
            endif
         enddo
         delt=min(delt,.5d0)
         delt=max(delt,minimum_initial_segment_size)
         s_sc1=s_scale_constant/real_axis_limits(n0)
         s_sc2=0.d0
         rmat=0.d0
         t2=0.d0
         do limit=1,n0
            t1=t2
            t2=real_axis_limits(limit)
            nseg=ceiling((t2-t1)/delt)
            deltseg=(t2-t1)/dble(nseg)
            dt2=t1
            do seg=1,nseg
               dt1=dt2
               dt2=min(dt2+deltseg,real_axis_integration_limit)
               drmat=0.d0
               subdiv=0
               ec=0
               call gkintegrate(qtot,dt1,dt2,realaxiskernel,drmat,subdiv,ec, &
                  integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
!if(mstm_global_rank.eq.0) write(*,'('' p1 '',2es12.4,2i6)') dt1,dt2,ec,subdiv
               if(ec.eq.1) error_codes(4)=1
               if(ec.eq.2) error_codes(3)=1
               rmat=rmat+drmat
               if(dt2.ge.real_axis_integration_limit) return
            enddo
         enddo
         s_sc1=0.d0
         s_sc2=s_scale_constant
         do limit=n0,number_limits-1
            t1=t2
            t2=real_axis_limits(limit+1)
            nseg=ceiling((t2-t1)/delt)
            deltseg=(t2-t1)/dble(nseg)
            dt2=t1
            do seg=1,nseg
               dt1=dt2
               dt2=min(dt2+deltseg,real_axis_integration_limit)
               drmat=0.d0
               subdiv=0
               ec=0
               call gkintegrate(qtot,dt1,dt2,realaxiskernel,drmat,subdiv,ec, &
                  integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
!if(mstm_global_rank.eq.0) write(*,'('' p2 '',2es12.4,2i6)') dt1,dt2,ec,subdiv
               if(ec.eq.1) error_codes(4)=1
               if(ec.eq.2) error_codes(3)=1
               rmat=rmat+drmat
               if(dt2.ge.real_axis_integration_limit) return
            enddo
         enddo
         delt=max(1.d0,minimum_initial_segment_size)
         do while(t2.lt.real_axis_integration_limit)
            t1=t2
            t2=t2+delt
            subdiv=0
            t1t=t1
            t2t=t2
            ec=0
            call gkintegrate(qtot,t1t,t2t,realaxiskernel,drmat,subdiv,ec, &
               integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
!if(mstm_global_rank.eq.0) write(*,'('' p3 '',2es12.4,2i6)') t1t,t2t,ec,subdiv
            if(ec.eq.1) error_codes(4)=1
            if(ec.eq.2) error_codes(3)=1
            rmat=rmat+drmat
            dnorm=mnorm(qtot,drmat)
            norm0=max(mnorm(qtot,rmat),0.01*integration_limit_epsilon)
            errlim=dnorm/norm0
            if(errlim.lt.integration_limit_epsilon) return
         enddo
         error_codes(2)=1
         end subroutine refmatrixrealaxis

         subroutine boundary_energy_transfer(sinc,sdir,r,t,a,fres_r,fres_t)
         implicit none
         integer :: sdir,rdir,tdir
         real(8) :: sinc,zref,ztra,r(2),t(2),a(2)
         complex(8) :: riref,s,rkz,tkz,gfr(2,2,2),gft(2,2,2),ritra
         complex(8), optional :: fres_r(2),fres_t(2)
         if(sdir.eq.1) then
            riref=layer_ref_index(0)
            ritra=layer_ref_index(number_plane_boundaries)
            zref=-1.d-8
            ztra=plane_boundary_position(max(1,number_plane_boundaries))+1.d-8
            rdir=2
            tdir=1
         else
            riref=layer_ref_index(number_plane_boundaries)
            ritra=layer_ref_index(0)
            ztra=-1.d-8
            zref=plane_boundary_position(max(1,number_plane_boundaries))+1.d-8
            rdir=1
            tdir=2
         endif
         s=sinc
         if(number_plane_boundaries.eq.0) then
            gfr=0.d0
            gft=1.d0
            rkz=1.d0
            tkz=1.d0
         else
            call layer_gf(s,zref,zref,gfr,rkz,tkz)
            call layer_gf(s,zref,ztra,gft,rkz,tkz)
         endif
         r(:)=cdabs(gfr(rdir,sdir,:))**2
         t(1)=cdabs(gft(tdir,sdir,1))**2*dble(tkz*dconjg(ritra/riref)/rkz)
         t(2)=cdabs(gft(tdir,sdir,2))**2*dble(dconjg(tkz*ritra/rkz/riref))
         a(:)=1.d0-r(:)-t(:)
         if(present(fres_r)) then
            fres_r(:)=gfr(rdir,sdir,:)
         endif
         if(present(fres_t)) then
            fres_t(:)=gft(tdir,sdir,:)
         endif
         end subroutine boundary_energy_transfer

         subroutine incident_field_initialization(alpha,sinc,sdir)
         implicit none
         integer :: sdir,ssign,p,k,q,klq
         real(8) :: alpha,sinc,targetz,sourcez
         complex(8) s,pmnp(6,2),riinc,skz,tkz,gfs(2,2,2)
         incident_lateral_vector=(/sinc*cos(alpha),sinc*sin(alpha)/)
         if(number_plane_boundaries.eq.0) then
            incident_field_scale=1.d0
            incident_field_boundary=0.d0
            return
         endif
         if(sdir.eq.1) then
            riinc=layer_ref_index(0)
         else
            riinc=layer_ref_index(number_plane_boundaries)
         endif
         if(sinc.gt.dble(riinc)) then
            s=sinc
            if(sdir.eq.1) then
               sourcez=-1.d-8
            else
               sourcez=plane_boundary_position(number_plane_boundaries)+1.d-8
            endif
            targetz=0.5d0*plane_boundary_position(number_plane_boundaries)
            call layer_gf(s,sourcez,targetz,gfs,skz,tkz)
            call genplanewavecoef(alpha,tkz,1,pmnp)
            do p=1,2
               do k=-1,1
                  do q=1,2
                     klq=amnpaddress(k,1,q,1,2)
                     ssign=(-1)**(k+q+p)
                     pmnp(klq,p)=pmnp(klq,p)*(gfs(1,sdir,p)+ssign*gfs(2,sdir,p))
                  enddo
               enddo
               incident_field_scale(p)=sqrt(sum(cdabs(pmnp(1:3,p)**2)))
            enddo
         else
            if(sdir.eq.1) then
               sourcez=bot_boundary
            else
               sourcez=top_boundary
            endif
            incident_field_scale(:)=1.d0
         endif
         incident_field_scale(:)=maxval(incident_field_scale(1:2))
         incident_field_boundary=sourcez
         end subroutine incident_field_initialization

         subroutine layerplanewavecoef(alpha,sinc,sdir,rpos,nodr,pmnp,include_direct)
         implicit none
         logical, optional :: include_direct
         logical :: incdir,evanescent
         integer :: p,incregion,sdir,nodr,nblk,layer
         real(8) :: alpha,ca,sa,rpos(3),sourcez,targetz,sinc
         complex(8) :: pmnp(2*nodr*(nodr+2),2),riinc,cbinc, &
                       s,phasefaclat,skz,tkz,gfs(2,2,2)
         complex(8), allocatable :: pmnpinc(:,:),pmnpup(:,:), &
                pmnpdn(:,:),pmnptot(:,:)
         if(present(include_direct)) then
            incdir=include_direct
         else
            incdir=.true.
         endif
         nblk=2*nodr*(nodr+2)
         ca=cos(alpha)
         sa=sin(alpha)
         layer=layer_id(rpos(3))
         if(sdir.eq.1) then
            riinc=layer_ref_index(0)
            incregion=0
         else
            riinc=layer_ref_index(number_plane_boundaries)
            incregion=number_plane_boundaries
         endif
         if(sinc.gt.dble(riinc)) then
            evanescent=.true.
            incdir=.false.
         else
            evanescent=.false.
         endif
         sourcez=incident_field_boundary
         s=sinc
         cbinc=cdsqrt((1.d0-sinc/riinc)*(1.d0+sinc/riinc))*(3-2*sdir)
         pmnp=0.d0
         targetz=rpos(3)
         phasefaclat=cdexp((0.d0,1.d0)*s*(ca*rpos(1)+sa*rpos(2)))
         allocate(pmnptot(nblk,2))
         pmnptot=0.d0
         if(layer.eq.incregion.and.incdir) then
            allocate(pmnpinc(nblk,2))
            call genplanewavecoef(alpha,cbinc,nodr,pmnpinc)
            pmnptot=pmnpinc*cdexp((0.d0,1.d0)*riinc*cbinc*(rpos(3)-sourcez))
            deallocate(pmnpinc)
         endif
         if(number_plane_boundaries.gt.0) then
            call layer_gf(s,sourcez,targetz,gfs,skz,tkz)
            allocate(pmnpup(nblk,2),pmnpdn(nblk,2))
            call genplanewavecoef(alpha,tkz,nodr,pmnpup)
            call genplanewavecoef(alpha,-tkz,nodr,pmnpdn)
            do p=1,2
               pmnptot(:,p)=pmnptot(:,p)+pmnpup(:,p)*gfs(1,sdir,p)  &
                  +pmnpdn(:,p)*gfs(2,sdir,p)
            enddo
            deallocate(pmnpup,pmnpdn)
         else
            tkz=abs(cbinc)
            skz=abs(cbinc)
         endif
         do p=1,2
            pmnp(:,p)=pmnptot(:,p)*phasefaclat/incident_field_scale(p)
         enddo
         deallocate(pmnptot)
         end subroutine layerplanewavecoef

         subroutine layervsh(s,alpha,targetz,tdir,rpos,nodr,pmnp)
         implicit none
         integer :: p,tdir,nodr,nblk,slayer,tlayer,k,l,q,klq,ssign
         real(8) :: alpha,ca,sa,rpos(3),sourcez,targetz
         complex(8) :: pmnp(2*nodr*(nodr+2),2),targetri,sourceri, &
                       s,phasefaclat,skz,tkz,gfs(2,2,2)

         nblk=2*nodr*(nodr+2)
         ca=cos(alpha)
         sa=sin(alpha)
         slayer=layer_id(rpos(3))
         tlayer=layer_id(targetz)
         targetri=layer_ref_index(tlayer)
         sourceri=layer_ref_index(slayer)
         sourcez=rpos(3)
         phasefaclat=cdexp(-(0.d0,1.d0)*s*(ca*rpos(1)+sa*rpos(2)))
         call layer_gf(s,sourcez,targetz,gfs,skz,tkz,include_direct=.true.)
         call genplanewavecoef(alpha,dconjg(skz),nodr,pmnp,lr_tran=.false.)
         do p=1,2
            do l=1,nodr
               do k=-l,l
                  do q=1,2
                     klq=amnpaddress(k,l,q,nodr,2)
                     ssign=(-1)**(k+l+q-1+p)
                     pmnp(klq,p)=dconjg(pmnp(klq,p))*(gfs(tdir,1,p)+ssign*gfs(tdir,2,p))
                  enddo
               enddo
            enddo
         enddo
         pmnp=pmnp*phasefaclat/4.d0/sourceri/sourceri/skz
         end subroutine layervsh

      end module surface_subroutines

      module periodic_lattice_subroutines
      use numconstants
      use specialfuncs
      use surface_subroutines
      use mpidefs
      implicit none
      logical, target :: periodic_lattice,time_it,phase_shift_form,finite_lattice
      integer :: pl_max_subdivs,pl_rs_nmax,pl_error_codes(6),pl_fs_method,pl_rs_imax, &
         q1d_number_segments,q2d_number_segments,s_max_q2
      integer, target :: pl_integration_method
      real(8) :: lattice_integration_segment,pl_rs_eps,time_count(4),time_0
      real(8), target :: cell_width(2),rs_dz_min,pl_integration_error_epsilon, &
         pl_integration_limit_epsilon
      data lattice_integration_segment/1.d0/
      data pl_rs_nmax,pl_rs_eps,pl_rs_imax/200,1.d-7,0/
      data time_it,rs_dz_min,s_max_q2/.true.,100.d0,100/
      data pl_integration_error_epsilon,pl_integration_limit_epsilon,pl_integration_method/1.d-7,1.d-9,1/
      data time_count/0.d0,0.d0,0.d0,0.d0/
      data phase_shift_form,finite_lattice/.false.,.false./

      contains

         subroutine plane_boundary_lattice_interaction(nodrt,nodrs,x0,y0,zt,zs, &
            matrix,source_vector,include_source,lr_transformation,index_model)
         implicit none
         logical :: incsrc,lrtran,rsincsrc,fsincsrc
         logical, optional :: include_source,lr_transformation
         integer :: nodrt,nodrs,i,ix,iy,nmax,nterms,n,p,q,imodl,m,l,k,mnp,klq,mn,kl,tranmat(2,2), &
            slay,tlay,np(2)
         integer, optional :: index_model
         real(8) :: x0,y0,zt,zs,w(2),wx,wy,k0x,k0y,kconst,kx,ky,eps,cerr,asum,asum0,x,y,wcrit
         complex(8) :: matrix(*),csum(2,2),ri
         complex(8), optional :: source_vector(2*nodrs*(nodrs+2))
         complex(8), allocatable :: kernel(:,:,:,:),dkernel(:,:,:,:),fsmat(:,:,:),rsmat(:,:)

         tranmat=reshape((/1,1,1,-1/),(/2,2/))
         if(present(include_source)) then
            incsrc=include_source
         else
            incsrc=.false.
         endif
         if(present(lr_transformation)) then
            lrtran=lr_transformation
         else
            lrtran=.true.
         endif
         if(present(index_model)) then
            imodl=index_model
         else
            imodl=2
         endif
         slay=layer_id(zs)
         tlay=layer_id(zt)
         if(.not.plane_surface_present) then
            if(present(source_vector)) then
               call free_space_lattice_translation_matrix(nodrt,nodrs,(/x0,y0,zt-zs/), &
                  layer_ref_index(0),matrix,source_vector=source_vector, &
                  include_source=incsrc,lr_transformation=lrtran,index_model=imodl)
            else
               call free_space_lattice_translation_matrix(nodrt,nodrs,(/x0,y0,zt-zs/), &
                  layer_ref_index(0),matrix, &
                  include_source=incsrc,lr_transformation=lrtran,index_model=imodl)
            endif
            return
         elseif(slay.eq.tlay) then
            if(present(source_vector)) then
               call common_layer_lattice_translation_matrix(nodrt,nodrs,x0,y0,zt,zs, &
                  matrix,source_vector=source_vector, &
                  include_source=incsrc,lr_transformation=lrtran,index_model=imodl)
            else
               call common_layer_lattice_translation_matrix(nodrt,nodrs,x0,y0,zt,zs, &
                  matrix,include_source=incsrc,lr_transformation=lrtran,index_model=imodl)
            endif
            return
         endif
         w=cell_width
         wcrit=min(rs_dz_min,minval(cell_width)/2.d0)
         k0x=incident_lateral_vector(1)
         k0y=incident_lateral_vector(2)
         np(1)=floor((x0+cell_width(1)/2.d0)/cell_width(1))
         np(2)=floor((y0+cell_width(2)/2.d0)/cell_width(2))
         x=x0-cell_width(1)*dble(np(1))
         y=y0-cell_width(2)*dble(np(2))
         nmax=pl_rs_nmax
         eps=pl_rs_eps
         if(incsrc) then
            ri=layer_ref_index(slay)
            if(slay.eq.tlay) then
               pl_fs_method=0
               rsincsrc=.false.
               fsincsrc=.true.
            else
               pl_fs_method=1
               rsincsrc=.true.
               fsincsrc=.false.
            endif
         else
            rsincsrc=.false.
            fsincsrc=.false.
         endif
         allocate(kernel(2,nodrt*(nodrt+2),2,nodrs*(nodrs+2)),dkernel(2,nodrt*(nodrt+2),2,nodrs*(nodrs+2)))
         wx=w(1)
         wy=w(2)
         kconst=8.d0*pi*pi/cell_width(1)/cell_width(2)
         kernel=0.d0
         kx=incident_lateral_vector(1)
         ky=incident_lateral_vector(2)
         call plane_boundary_lattice_kernel(nodrt,nodrs,kx,ky,x,y,zt,zs,kernel,include_source=rsincsrc)
         do n=1,pl_rs_nmax
            dkernel=0.d0
            do i=0,8*n-1
               q=i/(2*n)
               p=i-2*q*n
               if(q.eq.0) then
                  ix=n
                  iy=-n+p
               elseif(q.eq.1) then
                  ix=n-p
                  iy=n
               elseif(q.eq.2) then
                  ix=-n
                  iy=n-p
               else
                  ix=-n+p
                  iy=-n
               endif
               kx=2.d0*pi*dble(ix)/cell_width(1)+incident_lateral_vector(1)
               ky=2.d0*pi*dble(iy)/cell_width(2)+incident_lateral_vector(2)
               call plane_boundary_lattice_kernel(nodrt,nodrs,kx,ky,x,y,zt,zs,dkernel,include_source=rsincsrc)
            enddo
            kernel=kernel+dkernel
            asum0=sum(cdabs(kernel))
            asum=sum(cdabs(dkernel))
            cerr=asum/asum0
            if(cerr.lt.pl_rs_eps) exit
         enddo
         nterms=n
         pl_rs_imax=nterms
         if(nterms.ge.pl_rs_nmax) pl_error_codes(3)=1
         deallocate(dkernel)
         allocate(rsmat(2*nodrt*(nodrt+2),2*nodrs*(nodrs+2)))
         do l=1,nodrs
            do k=-l,l
               kl=l*(l+1)+k
               do n=1,nodrt
                  do m=-n,n
                     mn=n*(n+1)+m
                     do q=1,2
                        do p=1,2
                           csum(p,q)=kernel(p,mn,q,kl)
                        enddo
                     enddo
                     if(lrtran) then
                        csum=matmul(tranmat,matmul(csum,tranmat))/2.d0
                     endif
                     do q=1,2
                        klq=amnpaddress(k,l,q,nodrs,imodl)
                        do p=1,2
                           mnp=amnpaddress(m,n,p,nodrt,imodl)
                           rsmat(mnp,klq)=csum(p,q)*kconst
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         deallocate(kernel)
         if(slay.eq.tlay) then
            allocate(fsmat(nodrt*(nodrt+2),nodrs*(nodrs+2),2))
            ri=layer_ref_index(slay)
            call free_space_lattice_translation_matrix(nodrt,nodrs,(/x,y,zt-zs/), &
              ri,fsmat,include_source=fsincsrc,lr_transformation=lrtran,index_model=imodl)
            do l=1,nodrs
               do k=-l,l
                  kl=amnaddress(k,l,nodrs,imodl)
                  do n=1,nodrt
                     do m=-n,n
                        mn=amnaddress(m,n,nodrt,imodl)
                        do p=1,2
                           klq=amnpaddress(k,l,p,nodrs,imodl)
                           mnp=amnpaddress(m,n,p,nodrt,imodl)
                           rsmat(mnp,klq)=rsmat(mnp,klq)+fsmat(mn,kl,p)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            deallocate(fsmat)
         endif
         if(present(source_vector)) then
            matrix(1:2*nodrt*(nodrt+2)) &
               =matmul(rsmat(:,:),source_vector(:))*cdexp((0.d0,1.d0)*sum(incident_lateral_vector*dble(np)*cell_width))
         else
            matrix(1:4*nodrt*(nodrt+2)*nodrs*(nodrs+2)) &
               =reshape(rsmat,(/4*nodrt*(nodrt+2)*nodrs*(nodrs+2)/))*cdexp((0.d0,1.d0) &
                  *sum(incident_lateral_vector*dble(np)*cell_width))
         endif
         deallocate(rsmat)
         end subroutine plane_boundary_lattice_interaction

         subroutine free_space_lattice_translation_matrix(nodrt,nodrs,rpos0, &
             ri,matrix,source_vector,include_source,lr_transformation,index_model)
         implicit none
         logical :: incsrc,lrtran
         logical, optional :: include_source,lr_transformation
         integer :: nodrt,nodrs,wmax,l,m,mn,n,np(2), &
                    m1m,k,kl,nn1,ll1,v,w,wmin,vw,imodl,nterms,klp,p
         integer, optional :: index_model
         real(8) :: rpos(3),vc1(0:nodrs+nodrt),vc2(0:nodrs+nodrt),rpos0(3),wcrit
         complex(8) :: ci,c,a,b,ri,pshift, &
             ysum(0:(nodrt+nodrs)*(nodrt+nodrs+2)), &
             swf(0:(nodrt+nodrs)*(nodrt+nodrs+2)),matrix(*)
         complex(8), allocatable :: fsmat(:,:,:)
         complex(8), optional :: source_vector(nodrs*(nodrs+2)*2)
         data ci/(0.d0,1.d0)/
         if(present(include_source)) then
            incsrc=include_source
         else
            incsrc=.true.
         endif
         if(present(lr_transformation)) then
            lrtran=lr_transformation
         else
            lrtran=.true.
         endif
         if(present(index_model)) then
            imodl=index_model
         else
            imodl=2
         endif
         wcrit=min(rs_dz_min,minval(cell_width)/2.d0)
         rpos(3)=rpos0(3)
         np(1)=floor((rpos0(1)+cell_width(1)/2.d0)/cell_width(1))
         np(2)=floor((rpos0(2)+cell_width(2)/2.d0)/cell_width(2))
         rpos(1)=rpos0(1)-dble(np(1))*cell_width(1)
         rpos(2)=rpos0(2)-dble(np(2))*cell_width(2)
         pshift=cdexp((0.d0,1.d0)*sum(incident_lateral_vector*dble(np)*cell_width))
!
! phase shift form=.true. scales scattering coefficients with lateral phase shift.   Used only for testing purposes.
!
         if(phase_shift_form) pshift=cdexp((0.d0,-1.d0)*sum(incident_lateral_vector*rpos(1:2)))
         pl_max_subdivs=0
         allocate(fsmat(nodrt*(nodrt+2),nodrs*(nodrs+2),2))
!
! finite lattice uses cell-centered translation matrix H^i-j instead of L.   Included only for testing/expeimentation purposes
!
         if(finite_lattice) then
            call gentranmatrix(nodrs,nodrt,translation_vector=rpos, &
                             refractive_index=(/ri,ri/),ac_matrix=fsmat,vswf_type=3, &
                             index_model=2)
         else
            wmax=nodrt+nodrs
            ysum=0.d0
            if(abs(rpos(3)).lt.wcrit) then
               pl_fs_method=0
               if(time_it) time_0=mstm_mpi_wtime()
               call swf_lattice_sum(wmax,rpos(1),rpos(2),rpos(3),cell_width,incident_lateral_vector, &
                  ri,ysum,include_source=incsrc)
               if(time_it) time_count(1)=mstm_mpi_wtime()-time_0+time_count(1)
               if(pl_max_subdivs.ge.maximum_integration_subdivisions) pl_error_codes(1)=1
            else
               pl_fs_method=1
               if(time_it) time_0=mstm_mpi_wtime()
               call reciprocal_space_swf_lattice_sum(wmax,rpos(1),rpos(2),rpos(3),cell_width, &
                  incident_lateral_vector,ri,pl_rs_nmax,pl_rs_eps,nterms,ysum)
               if(time_it) time_count(2)=mstm_mpi_wtime()-time_0+time_count(2)
               if(.not.incsrc) then
                  call scalar_wave_function(wmax,3,rpos(1),rpos(2),rpos(3),ri,swf)
                  ysum=ysum-swf
               endif
               if(nterms.ge.pl_rs_nmax) pl_error_codes(2)=1
            endif
            do l=1,nodrs
               ll1=l*(l+1)
               do n=1,nodrt
                  nn1=n*(n+1)
                  wmax=n+l
                  call vcfunc(-1,n,1,l,vc2)
                  c=-sqrt(4.d0*pi)*ci**(n-l)*fnr(n+n+1)*fnr(l+l+1)
                  do k=-l,l
                     kl=amnaddress(k,l,nodrs,imodl)
                     do m=-n,n
                        m1m=(-1)**m
                        mn=amnaddress(m,n,nodrt,imodl)
                        v=k-m
                        call vcfunc(-m,n,k,l,vc1)
                        wmin=max(abs(v),abs(n-l))
                        a=0.
                        b=0.
                        do w=wmax,wmin,-1
                           vw=w*(w+1)+v
                           if(mod(wmax-w,2).eq.0) then
                              a=a+(ci**w)*vc1(w)*vc2(w)*ysum(vw)/sqrt(dble(w+w+1))
                           else
                              b=b+(ci**w)*vc1(w)*vc2(w)*ysum(vw)/sqrt(dble(w+w+1))
                           endif
                        enddo
                        if(lrtran) then
                           fsmat(mn,kl,1)=m1m*c*(a+b)
                           fsmat(mn,kl,2)=m1m*c*(a-b)
                        else
                           fsmat(mn,kl,1)=m1m*c*a
                           fsmat(mn,kl,2)=m1m*c*b
                        endif
                     enddo
                  enddo
               enddo
            enddo
         endif
         if(present(source_vector)) then
            do p=1,2
               do n=1,nodrt
                  do m=-n,n
                     mn=amnaddress(m,n,nodrt,imodl)
                     c=0.d0
                     do l=1,nodrs
                        do k=-l,l
                           kl=amnaddress(k,l,nodrs,imodl)
                           klp=amnpaddress(k,l,p,nodrs,imodl)
                           c=c+fsmat(mn,kl,p)*source_vector(klp)
                        enddo
                     enddo
                     mn=amnpaddress(m,n,p,nodrt,imodl)
                     matrix(mn)=c*pshift
                  enddo
               enddo
            enddo
         else
            matrix(1:2*nodrt*(nodrt+2)*nodrs*(nodrs+2))=reshape(fsmat,(/2*nodrt*(nodrt+2)*nodrs*(nodrs+2)/))*pshift
         endif
         deallocate(fsmat)
         end subroutine free_space_lattice_translation_matrix

         subroutine common_layer_lattice_translation_matrix(nodrt,nodrs,x0,y0,zt,zs, &
            matrix,source_vector,include_source,lr_transformation,index_model)
         implicit none
         logical :: incsrc,lrtran
         logical, optional :: include_source,lr_transformation
         integer :: nodrt,nodrs,slay,nodrw,nterms,nblk,sdirs(2),tdirs(2),n,i,q,p,ix,iy, &
            ll1,nn1,wmax,l,k,m,mn,kl,mnp,klq,tsign,ssign,iw,wmin,v,np(2),imodl,sdir,tdir, &
            tranmat(2,2),vw,w,nelem
         integer, optional :: index_model
         real(8) :: x,y,zt,zs,kx,ky,vc1(0:nodrt+nodrs), &
            vcp1m1(0:nodrt+nodrs),vcm1m1(0:nodrt+nodrs),x0,y0,a0mag
         real(8), allocatable :: asum(:),asum0(:),cerr(:)
         complex(8) :: ri,c,tmat(-1:2),csum(2,2),matrix(*)
         complex(8), allocatable :: qsum(:,:,:,:),dqsum(:,:,:,:),mat(:,:),fsmat(:,:,:)
         complex(8), optional :: source_vector(2*nodrs*(nodrs+2))
         tranmat=reshape((/1,1,1,-1/),(/2,2/))
         if(present(include_source)) then
            incsrc=include_source
         else
            incsrc=.false.
         endif
         if(present(lr_transformation)) then
            lrtran=lr_transformation
         else
            lrtran=.true.
         endif
         if(present(index_model)) then
            imodl=index_model
         else
            imodl=2
         endif
         np(1)=floor((x0+cell_width(1)/2.d0)/cell_width(1))
         np(2)=floor((y0+cell_width(2)/2.d0)/cell_width(2))
         x=x0-cell_width(1)*dble(np(1))
         y=y0-cell_width(2)*dble(np(2))
         slay=layer_id(zs)
         ri=layer_ref_index(slay)
         nodrw=nodrs+nodrt
         nblk=nodrw*(nodrw+2)
         if(slay.eq.0) then
            sdirs=(/1,1/)
            tdirs=(/2,2/)
         elseif(slay.eq.number_plane_boundaries) then
            sdirs=(/2,2/)
            tdirs=(/1,1/)
         else
            sdirs=(/1,2/)
            tdirs=(/1,2/)
         endif
         nelem=3*(nblk+1)*(tdirs(2)-tdirs(1)+1)*(sdirs(2)-sdirs(1)+1)
         allocate(qsum(-1:1,0:nblk,tdirs(1):tdirs(2),sdirs(1):sdirs(2)), &
            dqsum(-1:1,0:nblk,tdirs(1):tdirs(2),sdirs(1):sdirs(2)),asum(nelem),asum0(nelem),cerr(nelem))
         kx=incident_lateral_vector(1)
         ky=incident_lateral_vector(2)
         qsum=0.d0
         call common_layer_lattice_kernel(nodrw,kx,ky,x,y,zt,zs,tdirs,sdirs,qsum)
         do n=1,pl_rs_nmax
            dqsum=0.d0
            do i=0,8*n-1
               q=i/(2*n)
               p=i-2*q*n
               if(q.eq.0) then
                  ix=n
                  iy=-n+p
               elseif(q.eq.1) then
                  ix=n-p
                  iy=n
               elseif(q.eq.2) then
                  ix=-n
                  iy=n-p
               else
                  ix=-n+p
                  iy=-n
               endif
               kx=2.d0*pi*dble(ix)/cell_width(1)+incident_lateral_vector(1)
               ky=2.d0*pi*dble(iy)/cell_width(2)+incident_lateral_vector(2)
               call common_layer_lattice_kernel(nodrw,kx,ky,x,y,zt,zs,tdirs,sdirs,dqsum)
            enddo
            qsum=qsum+dqsum
            asum0=reshape(cdabs(qsum),(/nelem/))
            asum=reshape(cdabs(dqsum),(/nelem/))
            a0mag=maxval(asum0)
            cerr=0.d0
            do i=1,nelem
               if(asum0(i).gt.1.d-15*a0mag) cerr(i)=asum(i)/asum0(i)
            enddo
            if(maxval(cerr).lt.pl_rs_eps) exit
         enddo
         nterms=n
         if(nterms.ge.pl_rs_nmax) pl_error_codes(3)=1
         deallocate(dqsum)
         qsum=qsum/cell_width(1)/cell_width(2)*pi**1.5d0
         allocate(mat(2*nodrt*(nodrt+2),2*nodrs*(nodrs+2)))
         do l=1,nodrs
            ll1=l*(l+1)
            do n=1,nodrt
               nn1=n*(n+1)
               wmax=n+l
               call vcfunc(1,n,-1,l,vcp1m1)
               call vcfunc(-1,n,-1,l,vcm1m1)
               c=(0.d0,1.d0)**(n-l)*fnr(n+n+1)*fnr(l+l+1)
               do k=-l,l
                  kl=ll1+k
                  do m=-n,n
                     mn=nn1+m
                     v=k-m
                     call vcfunc(-m,n,k,l,vc1)
                     wmin=max(abs(v),abs(n-l))
                     csum=0.d0
                     do q=1,2
                        klq=amnpaddress(k,l,q,nodrs,imodl)
                        do p=1,2
                           mnp=amnpaddress(m,n,p,nodrt,imodl)
                           tmat=0.d0
                           do sdir=sdirs(1),sdirs(2)
                              if(sdir.eq.1) then
                                 ssign=1
                              else
                                 ssign=(-1)**(k+l+q-1)
                              endif
                              do tdir=tdirs(1),tdirs(2)
                                 if(tdir.eq.1) then
                                    tsign=ssign
                                 else
                                    tsign=(-1)**(m+n+p-1)*ssign
                                 endif
                                 do w=wmax,wmin,-1
                                    iw=(-1)**w
                                    vw=w*(w+1)+v
                                    if(w.ge.2) then
                                       tmat(-1)=tmat(-1)+vc1(w)*vcm1m1(w)*qsum(-1,vw,tdir,sdir)*tsign
                                       tmat(1)=tmat(1)+iw*vc1(w)*vcm1m1(w)*qsum(1,vw,tdir,sdir)*tsign
                                    endif
                                    tmat(0)=tmat(0)+vc1(w)*vcp1m1(w)*qsum(0,vw,tdir,sdir)*tsign
                                    tmat(2)=tmat(2)+iw*vc1(w)*vcp1m1(w)*qsum(0,vw,tdir,sdir)*tsign
                                 enddo
                              enddo
                           enddo
                           csum(p,q)=(-1)**m*c*((-1)**q*tmat(-1)+(-1)**(p+q)*tmat(0) &
                              +(-1)**(n+l)*tmat(2)+(-1)**(n+l+p)*tmat(1))
                        enddo
                     enddo
                     if(lrtran) then
                        csum=matmul(tranmat,matmul(csum,tranmat))/2.d0
                     endif
                     do q=1,2
                        klq=amnpaddress(k,l,q,nodrs,imodl)
                        do p=1,2
                           mnp=amnpaddress(m,n,p,nodrt,imodl)
                           mat(mnp,klq)=csum(p,q)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         deallocate(qsum)
         allocate(fsmat(nodrt*(nodrt+2),nodrs*(nodrs+2),2))
         call free_space_lattice_translation_matrix(nodrt,nodrs,(/x,y,zt-zs/), &
           ri,fsmat,include_source=incsrc,lr_transformation=lrtran,index_model=imodl)
         do l=1,nodrs
            do k=-l,l
               kl=amnaddress(k,l,nodrs,imodl)
               do n=1,nodrt
                  do m=-n,n
                     mn=amnaddress(m,n,nodrt,imodl)
                     do p=1,2
                        klq=amnpaddress(k,l,p,nodrs,imodl)
                        mnp=amnpaddress(m,n,p,nodrt,imodl)
                        mat(mnp,klq)=mat(mnp,klq)+fsmat(mn,kl,p)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         deallocate(fsmat)
         if(present(source_vector)) then
            matrix(1:2*nodrt*(nodrt+2)) &
               =matmul(mat(:,:),source_vector(:))*cdexp((0.d0,1.d0)*sum(incident_lateral_vector*dble(np)*cell_width))
         else
            matrix(1:4*nodrt*(nodrt+2)*nodrs*(nodrs+2)) &
               =reshape(mat,(/4*nodrt*(nodrt+2)*nodrs*(nodrs+2)/))*cdexp((0.d0,1.d0) &
                  *sum(incident_lateral_vector*dble(np)*cell_width))
         endif
         deallocate(mat)
         end subroutine common_layer_lattice_translation_matrix

         subroutine swf_lattice_sum(nodr,x,y,z,w,k0,ri,swfsum,include_source)
         implicit none
         logical :: addsource
         logical, optional :: include_source
         integer :: nodr,n,nn1,k,m,mn
         real(8) :: x,y,z,w(2),k0(2),drot(-nodr:nodr,0:nodr*(nodr+2))
         complex(8) :: ri,swfsum(0:nodr*(nodr+2)),swf(0:nodr*(nodr+2)),act(-nodr:nodr),csum

         if(present(include_source)) then
            if(x.ne.0.d0.or.y.ne.0.d0.or.z.ne.0.d0) then
               addsource=include_source
            else
               addsource=.false.
            endif
         else
            addsource=.false.
         endif
         call swfyzlatticesum(nodr,z,-y,x,(/w(2),w(1)/),-k0(2),k0(1),ri,swfsum)
         call rotcoef(0.d0,nodr,nodr,drot)
         do n=0,nodr
            nn1=n*(n+1)
            do k=-n,n
               act(k)=swfsum(nn1+k)
            enddo
            do m=-n,n
               mn=nn1+m
               csum=0.
               do k=-n,n
                  csum=csum+drot(m,nn1-k)*act(k)
               enddo
               swfsum(nn1+m)=((-1)**n)*csum
            enddo
         enddo
         if(addsource) then
            call scalar_wave_function(nodr,3,x,y,z,ri,swf)
            swfsum=swfsum+swf
         endif
         end subroutine swf_lattice_sum

         subroutine swfyzlatticesum(nodr,x,y,z,w,k0y,k0z,ri,swfyzsum)
         implicit none
         logical :: convrg
         integer :: nodr,m,n,s,ntz,l,smaxp,smaxn
         real(8) :: x,k0y,k0z,w(2),y,z,kz,mag,terr(0:nodr*(nodr+2)),derr(0:nodr*(nodr+2))
         complex(8) :: q1d(0:nodr*(nodr+2)),q2d(-nodr:nodr), &
             swfyzsum(0:nodr*(nodr+2)),ri,ci, &
             c,ct,ymn(0:nodr*(nodr+2))
         complex(8), allocatable :: sum1(:,:)
         data ci/(0.d0,1.d0)/
         ntz=max(s_max_q2,ceiling(w(2)/2.d0/pi))
         allocate(sum1(-nodr:nodr,0:ntz))
         smaxp=ntz
         smaxn=-ntz
         swfyzsum=0.d0

         call q1dbnosource(nodr,x,y,z,w(2),k0z,ri,q1d)
         do n=0,nodr
            do m=-n,n
               swfyzsum(m+n*(n+1))=q1d(m+n*(n+1))
            enddo
         enddo

         sum1=0.d0
         convrg=.false.
         do s=0,ntz
            kz=k0z+2.d0*pi*dble(s)/w(2)
            call q2db(nodr,x,y,w(1),k0y,kz,ri,q2d)
            sum1(:,s)=q2d
            ct=kz/ri
            call crotcoef(ct,0,nodr,ymn)
            terr=0.d0
            derr=cdabs(swfyzsum(:))
            do n=0,nodr
               c=-((0.d0,-1.d0)**n)/w(2)/ri/sqrt(4.d0*pi/dble(n+n+1))
               do m=-n,n
                  l=m+n*(n+1)
                  swfyzsum(l)=swfyzsum(l) &
                      +c*ymn(l)*sum1(m,s)*cdexp((0.d0,1.d0)*kz*z)
               enddo
            enddo
            terr=cdabs(swfyzsum(:))
            if(abs(kz).gt.1.d0) then
               do n=0,nodr*(nodr+2)
                  if(terr(n).ne.0.d0) then
                     terr(n)=abs(terr(n)-derr(n))/terr(n)
                  endif
               enddo
               if(maxval(terr).lt.pl_rs_eps) then
                  smaxp=s
                  convrg=.true.
                  exit
               endif
            endif
         enddo
         if(.not.convrg) pl_error_codes(6)=1
         if(k0z.eq.0.d0) then
            do s=1,smaxp
               kz=-2.d0*pi*dble(s)/w(2)
               ct=kz/ri
               call crotcoef(ct,0,nodr,ymn)
               do n=0,nodr
                  c=-((0.d0,-1.d0)**n)/w(2)/ri/sqrt(4.d0*pi/dble(n+n+1))
                  do m=-n,n
                     l=m+n*(n+1)
                     swfyzsum(l)=swfyzsum(l) &
                         +c*ymn(l)*sum1(m,s)*cdexp((0.d0,1.d0)*kz*z)
                  enddo
               enddo
            enddo
         else
            convrg=.false.
            do s=1,ntz
               kz=k0z-2.d0*pi*dble(s)/w(2)
               call q2db(nodr,x,y,w(1),k0y,kz,ri,q2d)
               sum1(:,s)=q2d
               ct=kz/ri
               call crotcoef(ct,0,nodr,ymn)
               terr=0.d0
               derr=cdabs(swfyzsum(:))
               do n=0,nodr
                  c=-((0.d0,-1.d0)**n)/w(2)/ri/sqrt(4.d0*pi/dble(n+n+1))
                  do m=-n,n
                     l=m+n*(n+1)
                     swfyzsum(l)=swfyzsum(l) &
                         +c*ymn(l)*sum1(m,s)*cdexp((0.d0,1.d0)*kz*z)
                  enddo
               enddo
               terr=cdabs(swfyzsum(:))
               if(abs(kz).gt.1.d0) then
                  do n=0,nodr*(nodr+2)
                     if(terr(n).ne.0.d0) then
                        terr(n)=abs(terr(n)-derr(n))/terr(n)
                     endif
                  enddo
                  if(maxval(terr).lt.pl_rs_eps) then
                     convrg=.true.
                     exit
                  endif
               endif
            enddo
            if(.not.convrg) pl_error_codes(6)=1
         endif
         deallocate(sum1)
         end subroutine swfyzlatticesum

         subroutine swfyzlatticesum0(nodr,x,y,z,w,k0y,k0z,ri,swfyzsum)
         implicit none
         logical :: convrg
         integer :: nodr,m,n,s,ntz,l,smaxp,smaxn
         real(8) :: x,k0y,k0z,w(2),y,z,kz,mag,tmag
         complex(8) :: q1d(0:nodr*(nodr+2)),q2d(-nodr:nodr), &
             swfyzsum(0:nodr*(nodr+2)),ri,ci, &
             c,ct,ymn(0:nodr*(nodr+2))
         complex(8), allocatable :: sum1(:,:)
         data ci/(0.d0,1.d0)/
         ntz=max(20,ceiling(w(2)/2.d0/pi))
         allocate(sum1(-nodr:nodr,-ntz:ntz))
         smaxp=ntz
         smaxn=-ntz
         swfyzsum=0.d0
         sum1=0.d0
         convrg=.false.
         do s=0,ntz
            kz=k0z+2.d0*pi*dble(s)/w(2)
            call q2db(nodr,x,y,w(1),k0y,kz,ri,q2d)
            sum1(:,s)=q2d
            mag=sum(cdabs(q2d))
            if((abs(kz).gt.1.d0).and.(mag.lt.1.d-8)) then
               smaxp=s
               convrg=.true.
               exit
            endif
         enddo
         if(.not.convrg) pl_error_codes(6)=1
         if(k0z.eq.0.) then
            smaxn=-smaxp
            do s=1,smaxp
               sum1(:,-s)=sum1(:,s)
            enddo
         else
            convrg=.false.
            do s=1,ntz
               kz=k0z-2.d0*pi*dble(s)/w(2)
               call q2db(nodr,x,y,w(1),k0y,kz,ri,q2d)
               sum1(:,-s)=q2d
               mag=sum(cdabs(q2d))
               if((abs(kz).gt.1.d0).and.(mag.lt.1.d-8)) then
                  smaxn=-s
                  convrg=.true.
                  exit
               endif
            enddo
            if(.not.convrg) pl_error_codes(6)=1
         endif
         call q1dbnosource(nodr,x,y,z,w(2),k0z,ri,q1d)
         do n=0,nodr
            do m=-n,n
               swfyzsum(m+n*(n+1))=q1d(m+n*(n+1))
            enddo
         enddo
         do s=smaxn,smaxp
            kz=k0z+2.d0*pi*dble(s)/w(2)
            ct=kz/ri
            call crotcoef(ct,0,nodr,ymn)
            do n=0,nodr
               c=-((0.d0,-1.d0)**n)/w(2)/ri/sqrt(4.d0*pi/dble(n+n+1))
               do m=-n,n
                  l=m+n*(n+1)
                  swfyzsum(l)=swfyzsum(l) &
                      +c*ymn(l)*sum1(m,s)*cdexp((0.d0,1.d0)*kz*z)
               enddo
            enddo
         enddo
         end subroutine swfyzlatticesum0

         subroutine qkernel2d(ntot,t,qfunc)
         implicit none
         integer :: ntot,integrationmodel,nodr,m
         real(8) :: t,x,y,w,kz,tt,dt,k0y,z
         complex(8) :: qfunc(ntot),ci,u,efunc1,efunc2,pfunc1,pfunc2,c,du,v,rkz,ri
         common/qkernelcommon/nodr,integrationmodel,x,y,z,w,k0y,kz,ri
         data ci/(0.d0,1.d0)/
         if(integrationmodel.eq.0) then
            tt=t
            dt=1.d0
         else
            tt=1.d0/t
            dt=tt/t
         endif
         rkz=cdsqrt((1.d0,0.d0)-kz*kz/ri/ri)
         c=tt*tt-2.d0*ci*rkz
         u=tt*cdsqrt(c)
         du=(c+tt*tt)/cdsqrt(c)
         v=cdsqrt(1.d0-kz*kz/ri/ri-u*u)
         pfunc1=(u-ci*v)/rkz
         pfunc2=(u+ci*v)/rkz
         efunc1=cdexp(ci*(k0y*w+ri*v*(w-y)))/(cdexp(ci*(k0y/ri+v)*w*ri)-1.d0)
         efunc2=-cdexp(ci*v*ri*(w+y))/(cdexp(ci*k0y*w)-cdexp(ci*v*ri*w))
         qfunc=0.d0
         do m=-nodr,nodr
            qfunc(m+1+nodr)=efunc1*pfunc1**m+efunc2*pfunc2**m
         enddo
         qfunc=qfunc*cdexp(ci*u*x*ri)*du/v*dt
         end subroutine qkernel2d

         subroutine q2db(nodr,x,y,w,k0y,kz,ri,qint)
         implicit none
         integer :: nodr,ntot,nodrc,integrationmodel,subdiv,nseg,ec
         real(8) :: x,y,w,k0y,kz,xc,yc,wc,k0yc,kzc,zc,t0,t1,cerr,dseg
         complex(8) :: qint(-nodr:nodr),qintt(-nodr:nodr),qintp(-nodr:nodr),ri,ric
         common/qkernelcommon/nodrc,integrationmodel,xc,yc,zc,wc,k0yc,kzc,ric
         nodrc=nodr
         xc=x
         yc=y
         wc=w
         k0yc=k0y
         kzc=kz
         ric=ri
         qint=0.
         ntot=2*nodr+1
         integrationmodel=0
         cerr=1.d0
         t0=0.d0
         nseg=0
         dseg=lattice_integration_segment/wc
         do while(cerr.gt.pl_integration_limit_epsilon)
            t1=t0
            t0=t0-dseg
            nseg=nseg+1
            qintt=0.
            subdiv=0
            if(nseg.eq.2.and.pl_integration_method.eq.0) then
               integrationmodel=1
               t0=-1.d0/dseg
               t1=0.d0
            endif
            ec=0
            call gkintegrate(ntot,t0,t1,qkernel2d,qintt,subdiv,ec, &
               pl_integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
            if(ec.ne.0) pl_error_codes(4)=1
            cerr=sum(cdabs(qintt))
            qint=qint+qintt
            cerr=cerr/sum(cdabs(qint))
            pl_max_subdivs=max(pl_max_subdivs,subdiv)
            if(nseg.eq.2.and.pl_integration_method.eq.0) exit
         enddo
         q2d_number_segments=nseg
         integrationmodel=0
         qintp=0.d0
         cerr=1.d0
         t1=0.d0
         nseg=0
         do while(cerr.gt.pl_integration_limit_epsilon)
            t0=t1
            t1=t1+dseg
            nseg=nseg+1
            qintt=0.
            subdiv=0
            if(nseg.eq.2.and.pl_integration_method.eq.0) then
               integrationmodel=1
               t1=1.d0/dseg
               t0=0.d0
            endif
            ec=0
            call gkintegrate(ntot,t0,t1,qkernel2d,qintt,subdiv,ec, &
               pl_integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
            if(ec.ne.0) pl_error_codes(4)=1
            cerr=sum(cdabs(qintt))
            qintp=qintp+qintt
            cerr=cerr/sum(cdabs(qintp))
            pl_max_subdivs=max(pl_max_subdivs,subdiv)
         if(nseg.eq.2.and.pl_integration_method.eq.0) exit
         enddo
         q2d_number_segments=max(q2d_number_segments,nseg)
         qint=qint+qintp
         end subroutine q2db

         subroutine qkernel1d(ntot,t,qfunc)
         implicit none
         integer :: ntot,integrationmodel,nodr,i,n,m,nodrtemp,mlim
         real(8) :: t,x,y,z,w,k0z,tt,dt,rho,k0y
         complex(8) :: qfunc(ntot),ymn1(0:nodr*(nodr+2)),ci,u,st1,rhos, &
                       ephi,bfunc(-nodr:nodr),expzfunc1,expzfunc2,c,ri
         common/qkernelcommon/nodr,integrationmodel,x,y,z,w,k0y,k0z,ri
         data ci/(0.d0,1.d0)/
         if(integrationmodel.eq.0) then
            tt=t
            dt=1.d0
         else
            tt=1.d0/t
            dt=tt/t
         endif
         u=dcmplx(1.d0,tt)
         rho=sqrt(x*x+y*y)
         st1=cdsqrt((1.d0-u)*(1.d0+u))
         rhos=rho*st1*ri
         bfunc=0.d0
         if(rho.eq.0.) then
            ephi=1.d0
            bfunc(0)=1.d0
            nodrtemp=0
         else
            ephi=dcmplx(x,y)/rho
            nodrtemp=nodr
            call bessel_integer_complex(nodr,rhos,nodrtemp,bfunc(0:nodr))
            do m=1,nodrtemp
               bfunc(-m)=(-1)**m*bfunc(m)
            enddo
         endif
         expzfunc1=cdexp(ci*u*ri*(w+z))/(cdexp(ci*k0z*w)-cdexp(ci*ri*u*w))
         expzfunc2=cdexp(ci*u*ri*(w-z))/(cdexp(-ci*k0z*w)-cdexp(ci*ri*u*w))
         call crotcoef(u,0,nodr,ymn1)
         qfunc=0.
         do n=0,nodr
            mlim=min(n,nodrtemp)
            do m=-mlim,mlim
               i=n*(n+1)+m+1
               c=-ci*((-ci)**(n-m))*(ephi**m)/2.d0/sqrt(pi/dble(n+n+1))
               qfunc(i)=dt*c*ymn1(n*(n+1)+m)*bfunc(m)*(expzfunc1+((-1)**(n+m))*expzfunc2)
            enddo
         enddo
         end subroutine qkernel1d

         subroutine q1dbnosource(nodr,x,y,z,w,kz,ri,qint)
         implicit none
         integer :: nodr,ntot,nodrc,integrationmodel,subdiv,nseg,ec
         real(8) :: x,y,z,w,kz,xc,yc,zc,wc,kzc,k0yc,t0,t1,cerr,dseg
         complex(8) :: qint(0:nodr*(nodr+2)),qintt(0:nodr*(nodr+2)),ri,ric
         common/qkernelcommon/nodrc,integrationmodel,xc,yc,zc,wc,k0yc,kzc,ric
         nodrc=nodr
         xc=x
         yc=y
         zc=z
         wc=w
         kzc=kz
         ric=ri
         ntot=1+nodr*(nodr+2)
         integrationmodel=0
         qint=0.
         t1=0.d0
         cerr=1.d0
         nseg=0
         dseg=lattice_integration_segment/wc
         do while(cerr.gt.pl_integration_limit_epsilon)
            t0=t1
            t1=t1+dseg
            nseg=nseg+1
            subdiv=0
            qintt=0.
            if(nseg.eq.2.and.pl_integration_method.eq.0) then
               integrationmodel=1
               t1=1.d0/dseg
               t0=0.d0
            endif
            ec=0
            call gkintegrate(ntot,t0,t1,qkernel1d,qintt,subdiv,ec, &
               pl_integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
            if(ec.ne.0) pl_error_codes(5)=1
            cerr=sum(cdabs(qintt))
            qint=qint+qintt
            cerr=cerr/sum(cdabs(qint))
            pl_max_subdivs=max(pl_max_subdivs,subdiv)
            if(nseg.eq.2.and.pl_integration_method.eq.0) exit
         enddo
         q1d_number_segments=nseg
         end subroutine q1dbnosource

         subroutine reciprocal_space_swf_lattice_sum(nodr,x0,y0,z0,w,k0,ri,nmax,eps,nterms,wf)
         implicit none
         integer :: nodr,i,ix,iy,nmax,nterms,n,p,q
         real(8) :: x0,y0,z0,w(2),wx,wy,k0(2),kconst,kx,ky,eps,cerr(0:nodr*(nodr+2)), &
            asum(0:nodr*(nodr+2)),asum0(0:nodr*(nodr+2)),a0mag
         complex(8) :: ri,wf(0:nodr*(nodr+2)),kswf(0:nodr*(nodr+2)),dwf(0:nodr*(nodr+2))
         wx=w(1)
         wy=w(2)
         kconst=2.d0*pi/wx/wy
         kx=k0(1)
         ky=k0(2)
         call reciprocal_scalar_wave_function(nodr,kx,ky,x0,y0,z0,ri,kswf)
         wf=kswf*kconst
         do n=1,nmax
            dwf=0.d0
            do i=0,8*n-1
               q=i/(2*n)
               p=i-2*q*n
               if(q.eq.0) then
                  ix=n
                  iy=-n+p
               elseif(q.eq.1) then
                  ix=n-p
                  iy=n
               elseif(q.eq.2) then
                  ix=-n
                  iy=n-p
               else
                  ix=-n+p
                  iy=-n
               endif
               kx=2.d0*pi*dble(ix)/wx+k0(1)
               ky=2.d0*pi*dble(iy)/wy+k0(2)
               call reciprocal_scalar_wave_function(nodr,kx,ky,x0,y0,z0,ri,kswf)
               dwf=dwf+kswf*kconst
            enddo
            wf=wf+dwf
            asum0=cdabs(wf)
            asum=cdabs(dwf)
            a0mag=maxval(asum0)
            cerr=0.d0
            do i=0,nodr*(nodr+2)
               if(asum0(i).gt.1.d-15*a0mag) cerr(i)=asum(i)/asum0(i)
            enddo
            if(maxval(cerr).lt.eps) exit
         enddo
         nterms=n
         end subroutine reciprocal_space_swf_lattice_sum

         subroutine common_layer_lattice_kernel(nodr,kx,ky,x,y,zt,zs,tdirs,sdirs,kernel)
         implicit none
         integer :: nodr,n,m,k,mn,sdirs(2),tdirs(2),sdir,tdir,slay,i,is,id,it
         real(8) :: kx,ky,x,y,zs,zt,kr
         complex(8) :: s,gfunc(2,2,2),skz,tkz, &
            drot(-2:2,0:nodr*(nodr+2)),ealpha,ri,c,c2, &
            kernel(-1:1,0:nodr*(nodr+2),tdirs(1):tdirs(2),sdirs(1):sdirs(2))
         if(time_it) time_0=mstm_mpi_wtime()
         slay=layer_id(zs)
         ri=layer_ref_index(slay)
         kr=sqrt(kx*kx+ky*ky)
         if(kr.eq.0.d0) then
            ealpha=1.d0
         else
            ealpha=dcmplx(kx,ky)/kr
         endif
         s=kr
         call layer_gf(s,zs,zt,gfunc,skz,tkz)
         call crotcoef(skz,2,nodr,drot)
         c=cdexp((0.d0,1.d0)*(kx*x+ky*y))/ri/ri/skz/sqrt(4.d0*pi)
         do k=-2,2,2
            i=k/2
            is=(-1)**i
            do n=abs(k),nodr
               do m=-n,n
                  mn=n*(n+1)+m
                  c2=c*ealpha**m*drot(k,mn)
                  do sdir=sdirs(1),sdirs(2)
                     do tdir=tdirs(1),tdirs(2)
                        if(sdir.eq.tdir) then
                           id=1
                           it=-1
                        else
                           id=-1
                           it=1
                        endif
                        kernel(i,mn,tdir,sdir)=kernel(i,mn,tdir,sdir)+c2*it*(gfunc(tdir,sdir,1)+is*id*gfunc(tdir,sdir,2))
                     enddo
                  enddo
               enddo
            enddo
         enddo
         if(time_it) time_count(4)=mstm_mpi_wtime()-time_0+time_count(4)
         end subroutine common_layer_lattice_kernel

         subroutine plane_boundary_lattice_kernel(nodrt,nodrs,kx,ky,x,y,zt,zs,kernel, &
            include_source)
         implicit none
         logical :: incsrc
         logical, optional :: include_source
         integer :: nodrs,nodrt,n,m,p,k,l,q,mn,kl,ssign,tsign,pol
         real(8) :: kx,ky,x,y,zs,zt,kr
         complex(8) :: kernel(2,nodrt*(nodrt+2),2,nodrs*(nodrs+2)),s,gfunc(2,2,2),skz,tkz, &
            pivec(2,nodrs*(nodrs+2),2),picvec(2,nodrt*(nodrt+2),2),ealpha,ri,ekm,csum(2,2),c
         if(present(include_source)) then
            incsrc=include_source
         else
            incsrc=.false.
         endif
         if(zt.eq.zs) incsrc=.false.
         if(time_it) time_0=mstm_mpi_wtime()
         ri=layer_ref_index(layer_id(zs))
         kr=sqrt(kx*kx+ky*ky)
         if(kr.eq.0.d0) then
            ealpha=1.d0
         else
            ealpha=dcmplx(kx,ky)/kr
         endif
         s=kr
         call layer_gf(s,zs,zt,gfunc,skz,tkz,incsrc)
         call complexpivec(skz,nodrs,pivec,1)
         call complexpivec(tkz,nodrt,picvec,-1)
         c=cdexp((0.d0,1.d0)*(kx*x+ky*y))/ri/ri/skz
         do n=1,nodrt
            do m=-n,n
               mn=n*(n+1)+m
               do l=1,nodrs
                  do k=-l,l
                     kl=l*(l+1)+k
                     ekm=ealpha**(k-m)
                     csum=0.d0
                     do p=1,2
                        tsign=(-1)**(m+n+p-1)
                        do q=1,2
                           ssign=(-1)**(k+l+q-1)
                           do pol=1,2
                              csum(p,q)=csum(p,q)+picvec(p,mn,pol)*pivec(q,kl,pol) &
                                    *(gfunc(1,1,pol)+(-1)**pol*(tsign*gfunc(2,1,pol)+ssign*gfunc(1,2,pol)) &
                                    +tsign*ssign*gfunc(2,2,pol))
                           enddo
                        enddo
                     enddo
                     kernel(:,mn,:,kl)=kernel(:,mn,:,kl)+csum(:,:)*c*ekm
                  enddo
               enddo
            enddo
         enddo
         if(time_it) time_count(3)=mstm_mpi_wtime()-time_0+time_count(3)
         end subroutine plane_boundary_lattice_kernel

      end module periodic_lattice_subroutines
!
! module spheredata: used to 1) input sphere data, 2) dimension sphere data
! arrays, and 3) provide common access to the data in other subroutines.
!
!
!  last revised: 15 January 2011
!
!  30 March 2011: added optical activity
!  March 2013: a few more options
!
      module spheredata
      use specialfuncs
      use numconstants
      use surface_subroutines
      use periodic_lattice_subroutines
      type linked_sphere_list
         integer :: sphere
         type(linked_sphere_list), pointer :: next
      end type linked_sphere_list
      type host_list
         integer :: number
         type(linked_sphere_list), pointer :: sphere_list
      end type host_list
      type(host_list), allocatable :: sphere_links(:,:)
      logical :: one_side_only,recalculate_surface_matrix,any_optically_active
      logical, target :: store_translation_matrix,effective_medium_simulation, &
          store_surface_matrix,fft_translation_option,normalize_solution_error
      logical, allocatable :: optically_active(:)
      integer, target :: number_spheres,run_print_unit,translation_switch_order, &
         max_t_matrix_order
      integer :: max_mie_order,number_host_spheres,number_eqns,t_matrix_order,max_sphere_depth
      integer, allocatable :: host_sphere(:),sphere_order(:), &
                 sphere_block(:),sphere_offset(:),translation_order(:), &
                 number_field_expansions(:),mie_offset(:), &
                 mie_block_offset(:),sphere_layer(:),sphere_depth(:)
      real(8) :: cluster_origin(3),vol_radius,sphere_mean_position(3),area_mean_radius, &
         sphere_min_position(3),sphere_max_position(3),mean_qext_mie,mean_qabs_mie, &
         circumscribing_radius,cross_section_radius,effective_cluster_radius
      real(8), target :: gaussian_beam_constant,gaussian_beam_focal_point(3)
      real(8), allocatable :: qext_mie(:),qabs_mie(:)
      real(8), allocatable, target :: sphere_radius(:),sphere_position(:,:)
      complex(8) :: effective_ref_index
      complex(8), allocatable :: an_mie(:),cn_mie(:),un_mie(:),vn_mie(:),dn_mie(:),an_inv_mie(:)
      complex(8), allocatable, target :: sphere_ref_index(:,:)

      data run_print_unit/6/
      data store_translation_matrix,store_surface_matrix/.false.,.true./
      data recalculate_surface_matrix/.true./
      data translation_switch_order/3/
      data fft_translation_option/.false./
      data max_t_matrix_order/100/
      data normalize_solution_error/.true./
      data gaussian_beam_constant,gaussian_beam_focal_point/0.d0,0.d0,0.d0,0.d0/
      data effective_medium_simulation,effective_ref_index/.false.,(1.d0,0.d0)/

      contains

         subroutine sphere_layer_initialization()
         implicit none
         type(linked_sphere_list), pointer ::slist
         integer :: i,j,l

         if(allocated(sphere_layer)) deallocate(sphere_layer)
         allocate(sphere_layer(number_spheres))
         sphere_layer=0
         if(number_plane_boundaries.gt.0) then
            do i=1,number_spheres
               do j=1,number_plane_boundaries
                  if(sphere_position(3,i).gt.plane_boundary_position(j)) then
                     sphere_layer(i)=j
                  else
                     exit
                  endif
               enddo
            enddo
         endif

         if(number_plane_boundaries.gt.0) then
            top_boundary=max(plane_boundary_position(max(1,number_plane_boundaries)),sphere_max_position(3))+1.d-5
            bot_boundary=min(0.d0,sphere_min_position(3))-1.d-5
         else
            top_boundary=sphere_max_position(3)+1.d-5
            bot_boundary=sphere_min_position(3)-1.d-5
         endif
            top_boundary=max(plane_boundary_position(max(1,number_plane_boundaries)),sphere_max_position(3))+1.d-5
            bot_boundary=min(0.d0,sphere_min_position(3))-1.d-5

         if(allocated(sphere_links)) then
            do l=0,ubound(sphere_links,2)-1
               do i=0,ubound(sphere_links,1)-1
                  call clear_host_list(sphere_links(i,l))
               enddo
            enddo
            deallocate(sphere_links)
         endif
         allocate(sphere_links(0:number_spheres,0:number_plane_boundaries))
         do l=0,number_plane_boundaries
            do i=0,number_spheres
               sphere_links(i,l)%number=0
               allocate(sphere_links(i,l)%sphere_list)
            enddo
         enddo
         do l=0,number_plane_boundaries
            do j=0,number_spheres
               slist=>sphere_links(j,l)%sphere_list
               do i=1,number_spheres
                  if(host_sphere(i).eq.j.and.sphere_layer(i).eq.l) then
                     sphere_links(j,l)%number=sphere_links(j,l)%number+1
                     slist%sphere=i
                     allocate(slist%next)
                     slist=>slist%next
                  endif
               enddo
            enddo
         enddo

         if(allocated(sphere_depth)) deallocate(sphere_depth)
         allocate(sphere_depth(number_spheres))
         max_sphere_depth=0
         do i=1,number_spheres
            sphere_depth(i)=0
            j=host_sphere(i)
            do while(j.ne.0)
               sphere_depth(i)=sphere_depth(i)+1
               j=host_sphere(j)
            enddo
            max_sphere_depth=max(max_sphere_depth,sphere_depth(i))
         enddo

         end subroutine sphere_layer_initialization

         subroutine clear_host_list(hlist)
         implicit none
         type(host_list) :: hlist
         type(linked_sphere_list), pointer :: slist,slist2
         integer :: n,i
         if(.not.associated(hlist%sphere_list)) return
         n=hlist%number
         slist=>hlist%sphere_list
         do i=1,n
            slist2=>slist%next
            deallocate(slist)
            slist=>slist2
         enddo
         end subroutine clear_host_list
!
!  findhostspheres finds the host sphere of each sphere in the set.   host=0 for
!  external sphere.
!
!  december 2011
!  march 2013: something changed
!
         subroutine findhostspheres()
         implicit none
         integer :: i,j,m
         real(8) :: xij(3),rij,xspmin,zmax,zmin

         host_sphere=0
         sphere_mean_position=0.d0
         sphere_min_position=1.d10
         sphere_max_position=-1.d10
         do i=1,number_spheres
            sphere_mean_position=sphere_mean_position+sphere_position(:,i)
            do m=1,3
               sphere_min_position(m)=min(sphere_min_position(m),sphere_position(m,i))
               sphere_max_position(m)=max(sphere_max_position(m),sphere_position(m,i))
            enddo
            xspmin=1.d6
            do j=1,number_spheres
               if(sphere_radius(j).gt.sphere_radius(i)) then
                  xij(:)=sphere_position(:,j)-sphere_position(:,i)
                  rij=sqrt(dot_product(xij,xij))
!                  if(rij.le.sphere_radius(j)-sphere_radius(i).and.sphere_radius(j).lt.xspmin) then
                  if(rij.le.sphere_radius(j).and.sphere_radius(j).lt.xspmin) then
                     host_sphere(i)=j
                     xspmin=sphere_radius(j)
                  endif
               endif
            enddo
         enddo
         sphere_mean_position=sphere_mean_position/dble(number_spheres)
         number_field_expansions=1
         number_host_spheres=0
         circumscribing_radius=0.d0
         do i=1,number_spheres
            rij=sqrt(sum((sphere_position(:,i)-cluster_origin(:))**2))+sphere_radius(i)
!            rij=sqrt(sum((sphere_position(:,i)-sphere_mean_position(:))**2))+sphere_radius(i)
            circumscribing_radius=max(circumscribing_radius,rij)
            j=host_sphere(i)
            if(j.ne.0) then
               number_field_expansions(j)=2
               number_host_spheres=number_host_spheres+1
            endif
         enddo
         zmax=maxval(sphere_position(3,1:number_spheres))
         zmin=minval(sphere_position(3,1:number_spheres))
         one_side_only=(zmax*zmin.gt.0.d0)

         end subroutine findhostspheres

      end module spheredata
!
! module miecoefdata: used to 1) calculate single sphere mie coefficient values,
! 2) store values in an allocated array, 3) provide common access to values, and
! 4) perform multiplication of coefficient values with vectors containing VWH scattering
! coefficients.
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!
      module mie

      contains
!
!  calculation of the max order of sphere expansions and storage of mie coefficients
!
!
!  last revised: 15 January 2011
!  30 March 2011: added optical activity
!  april 2012: all spheres are assumed OA: l/r formulation
!  february 2013: tmatrix file option.
!
         subroutine miecoefcalc(qeps)
         use mpidefs
         use spheredata
         implicit none
         integer :: i,nodrn,nsphere,ntermstot,nblktot,nterms, &
                    n1,n2,j,n
         real(8) :: qext,qsca,qeps,qabs
         complex(8) :: rihost(2)
         complex(8), allocatable :: anp(:,:,:),cnp(:,:,:),unp(:,:,:), &
                     vnp(:,:,:),dnp(:,:,:),anpinv(:,:,:)
         nsphere=number_spheres
         if(allocated(sphere_order)) deallocate(sphere_order,mie_offset,sphere_block, &
                      qext_mie,qabs_mie,optically_active,sphere_offset)
         allocate(sphere_order(nsphere),mie_offset(nsphere+1),sphere_block(nsphere), &
                      qext_mie(nsphere),qabs_mie(nsphere),optically_active(nsphere), &
                      sphere_offset(nsphere+1))
         ntermstot=0
         nblktot=0
         max_mie_order=0
         mean_qext_mie=0.d0
         mean_qabs_mie=0.d0
         area_mean_radius=0.d0
!
!
! march 2013: calculates orders, and
!             forces host to have at least same order as constituents.
!
         do i=1,number_spheres
            call exteriorrefindex(i,rihost)
            call mieoa(sphere_radius(i),sphere_ref_index(:,i),nodrn,qeps,qext,qsca,qabs, &
                    ri_medium=rihost)
            sphere_order(i)=nodrn
         enddo
         do i=1,number_spheres
            j=host_sphere(i)
            do while(j.ne.0)
               sphere_order(j)=max(sphere_order(j),sphere_order(i))
               j=host_sphere(j)
            enddo
         enddo
!
!  calculate the order limits and efficiencies
!
         n=0
         any_optically_active=.false.
         do i=1,number_spheres
            call exteriorrefindex(i,rihost)
            if(sphere_ref_index(1,i).eq.sphere_ref_index(2,i)) then
               optically_active(i)=.false.
            else
               optically_active(i)=.true.
               any_optically_active=.true.
            endif
            nodrn=sphere_order(i)
            sphere_block(i)=2*nodrn*(nodrn+2)
            call mieoa(sphere_radius(i),sphere_ref_index(:,i),nodrn,0.d0,qext,qsca,qabs, &
                    ri_medium=rihost)
            nterms=4*nodrn
            mie_offset(i)=ntermstot
            ntermstot=ntermstot+nterms
            qext_mie(i)=qext
            qabs_mie(i)=qabs
            if(host_sphere(i).eq.0) then
               mean_qext_mie=mean_qext_mie+qext*sphere_radius(i)**2
               mean_qabs_mie=mean_qabs_mie+qabs*sphere_radius(i)**2
               area_mean_radius=area_mean_radius+sphere_radius(i)**2
               n=n+1
            endif
            sphere_offset(i)=nblktot
            nblktot=nblktot+sphere_block(i)*number_field_expansions(i)
            max_mie_order=max(max_mie_order,sphere_order(i))
         enddo
         n=max(n,1)
         mie_offset(nsphere+1)=ntermstot
         sphere_offset(nsphere+1)=nblktot
         number_eqns=nblktot
         area_mean_radius=sqrt(area_mean_radius/dble(n))
         mean_qext_mie=mean_qext_mie/area_mean_radius**2/dble(n)
         mean_qabs_mie=mean_qabs_mie/area_mean_radius**2/dble(n)
!
! calculate the mie coefficients, and store in memory
!
         if(allocated(an_mie)) deallocate(an_mie,cn_mie,un_mie,vn_mie,dn_mie, &
                               an_inv_mie)
         allocate(an_mie(ntermstot),cn_mie(ntermstot),un_mie(ntermstot), &
                  vn_mie(ntermstot),dn_mie(ntermstot),an_inv_mie(ntermstot))
         do i=1,number_spheres
            nodrn=sphere_order(i)
            call exteriorrefindex(i,rihost)
            allocate(anp(2,2,nodrn),cnp(2,2,nodrn),unp(2,2,nodrn), &
               vnp(2,2,nodrn),dnp(2,2,nodrn),anpinv(2,2,nodrn))
            call mieoa(sphere_radius(i),sphere_ref_index(:,i),nodrn,0.d0,qext,qsca,qabs, &
                 anp_mie=anp,cnp_mie=cnp,dnp_mie=dnp, &
                 unp_mie=unp,vnp_mie=vnp,anp_inv_mie=anpinv, &
                 ri_medium=rihost)
            nterms=4*nodrn
            n1=mie_offset(i)+1
            n2=mie_offset(i)+nterms
            an_mie(n1:n2)=reshape(anp(1:2,1:2,1:nodrn),(/nterms/))
            cn_mie(n1:n2)=reshape(cnp(1:2,1:2,1:nodrn),(/nterms/))
            dn_mie(n1:n2)=reshape(dnp(1:2,1:2,1:nodrn),(/nterms/))
            un_mie(n1:n2)=reshape(unp(1:2,1:2,1:nodrn),(/nterms/))
            vn_mie(n1:n2)=reshape(vnp(1:2,1:2,1:nodrn),(/nterms/))
            an_inv_mie(n1:n2)=reshape(anpinv(1:2,1:2,1:nodrn),(/nterms/))
            deallocate(anp,cnp,unp,vnp,dnp,anpinv)
         enddo
         end subroutine miecoefcalc

         subroutine exteriorrefindex(i,rihost)
         use spheredata
         use surface_subroutines
         implicit none
         integer :: i
         complex(8) :: rihost(2)
         if(plane_surface_present.and.(host_sphere(i).eq.0)) then
            rihost=layer_ref_index(sphere_layer(i))
         else
            rihost=sphere_ref_index(:,host_sphere(i))
         endif
         end subroutine exteriorrefindex
!
! transformation between lr and te tm basis
!
         subroutine lrtomodetran(at,am)
         implicit none
         complex(8) :: a(2,2),am(2,2),at(2,2)
         a=at
         am(1,1)=(a(1,1) + a(1,2) + a(2,1) + a(2,2))/2.
         am(1,2)=(a(1,1) - a(1,2) + a(2,1) - a(2,2))/2.
         am(2,1)=(a(1,1) + a(1,2) - a(2,1) - a(2,2))/2.
         am(2,2)=(a(1,1) - a(1,2) - a(2,1) + a(2,2))/2.
         end subroutine lrtomodetran
!
! optically active lorenz/mie coefficients
! original 30 March 2011
! April 2012: generalized LR formulation, generalized mie coefficients
!
         subroutine mieoa(x,ri,nodr0,qeps,qext,qsca,qabs,anp_mie,dnp_mie, &
           unp_mie,vnp_mie,cnp_mie,ri_medium,anp_inv_mie,dnp_eff_mie,anp_eff_mie)
         use specialfuncs
         implicit none
         integer :: nstop,n,i,p,q,nodr0,s,t,ss,st
         real(8) :: x,qeps,qext,qsca,fn1,err,qextt,qabs
         real(8), allocatable :: qext1(:)
         complex(8) :: ci,ri(2),xri(2,2),rii(2,2),ribulk(2),psip(2,2), &
                     xip(2,2),gmatinv(2,2),bmatinv(2,2),gmat(2,2),bmat(2,2),&
                     amat(2,2),dmat(2,2),umat(2,2),vmat(2,2),amatinv(2,2), &
                     psipn(2,2),xipn(2,2)
         complex(8), optional :: anp_mie(2,2,*),dnp_mie(2,2,*),ri_medium(2), &
                                 unp_mie(2,2,*),vnp_mie(2,2,*),cnp_mie(2,2,*), &
                                 anp_inv_mie(2,2,*),dnp_eff_mie(2,2,*),anp_eff_mie(2,2,*)
         complex(8), allocatable :: psi(:,:,:),xi(:,:,:)
         data ci/(0.d0,1.d0)/

         if(present(ri_medium)) then
            rii(:,1)=ri_medium
         else
            rii(:,1)=(/(1.d0,0.d0),(1.d0,0.d0)/)
         endif
         rii(:,2)=ri
         ribulk(:)=2.d0/(1.d0/rii(1,:)+1.d0/rii(2,:))
         xri=x*rii
         if(qeps.gt.0.) then
            nstop=nint(x+4.*x**(1./3.))+5.
         elseif(qeps.lt.0.) then
            nstop=ceiling(-qeps)
            nodr0=nstop
         else
            nstop=nodr0
         endif
         allocate(psi(0:nstop+1,2,2),xi(0:nstop+1,2,2),qext1(nstop))

         do i=1,2
            do p=1,2
               call cricbessel(nstop+1,xri(p,i),psi(0,p,i))
               call crichankel(nstop+1,xri(p,i),xi(0,p,i))
            enddo
         enddo
         qabs=0.d0
         qsca=0.d0
         qext=0.d0
! 2/6/2019: scaled psip,xip
         do n=1,nstop
            psip(:,:) = psi(n-1,:,:)-dble(n)*psi(n,:,:)/xri(:,:)
            xip(:,:) = xi(n-1,:,:)-dble(n)*xi(n,:,:)/xri(:,:)
            psipn(:,:) = psi(n-1,:,:)/psi(n,:,:)-dble(n)/xri(:,:)
            xipn(:,:) = xi(n-1,:,:)/xi(n,:,:)-dble(n)/xri(:,:)
            do s=1,2
               ss=(-1)**s
               do t=1,2
                  st=(-1)**t
                  gmat(s,t)=(ss*ribulk(1)+st*ribulk(2))/(rii(s,2)*rii(t,1)) &
                    *(psipn(s,2)*xi(n,t,1)-ss*st*xip(t,1))
                  amat(s,t)=(ss*ribulk(1)+st*ribulk(2))/(rii(s,2)*rii(t,1)) &
                    *(psipn(s,2)*psi(n,t,1)-ss*st*psip(t,1))
!                  umat(s,t)=(ss*ribulk(2)+st*ribulk(2))/(rii(s,2)*rii(t,2)) &
!                    *(xip(s,2)*psi(n,t,2)-ss*st*xi(n,s,2)*psip(t,2))
                  umat(s,t)=(ss*ribulk(2)+st*ribulk(2))/(rii(s,2)*rii(t,2))*ci/psi(n,s,2)
                  bmat(s,t)=(ss*ribulk(2)+st*ribulk(1))/(rii(s,1)*rii(t,2)) &
                    *(xipn(s,1)*psi(n,t,2)-ss*st*psip(t,2))
!                  dmat(s,t)=(ss*ribulk(1)+st*ribulk(1))/(rii(s,1)*rii(t,1)) &
!                    *(psip(s,1)*xi(n,t,1)-ss*st*psi(n,s,1)*xip(t,1))
                  dmat(s,t)=-(ss*ribulk(1)+st*ribulk(1))/(rii(s,1)*rii(t,1))*ci/xi(n,s,1)
                  vmat(s,t)=(ss*ribulk(2)+st*ribulk(1))/(rii(s,1)*rii(t,2)) &
                    *(xipn(s,1)*xi(n,t,2)-ss*st*xip(t,2))
               enddo
            enddo

!write(*,'(8e12.4)') gmat
!call flush(6)

            call twobytwoinverse(gmat,gmatinv)
            call twobytwoinverse(bmat,bmatinv)

            amat=-matmul(gmatinv,amat)
            umat=-matmul(gmatinv,umat)
            dmat=-matmul(bmatinv,dmat)
            vmat=-matmul(bmatinv,vmat)

            if(present(anp_mie)) then
               anp_mie(:,:,n)=amat(:,:)
            endif
            if(present(dnp_mie)) then
               dnp_mie(:,:,n)=dmat(:,:)
            endif
            if(present(unp_mie)) then
               unp_mie(:,:,n)=umat(:,:)
            endif
            if(present(vnp_mie)) then
               vnp_mie(:,:,n)=vmat(:,:)
            endif
            if(present(cnp_mie)) then
               call twobytwoinverse(amat,amatinv)
               amatinv=matmul(dmat,amatinv)
               cnp_mie(:,:,n)=amatinv(:,:)
            endif
            if(present(anp_inv_mie)) then
               call twobytwoinverse(amat,amatinv)
               anp_inv_mie(:,:,n)=amatinv(:,:)
            endif
            if(present(dnp_eff_mie)) then
               call twobytwoinverse(umat,amatinv)
               amatinv=matmul(amatinv,amat)
               amatinv=matmul(vmat,amatinv)
               dnp_eff_mie(:,:,n)=dmat-amatinv
            endif
            if(present(anp_eff_mie)) then
               call twobytwoinverse(umat,amatinv)
               amatinv=matmul(amatinv,amat)
               anp_eff_mie(:,:,n)=-amatinv
            endif
            qext1(n)=0.d0
            fn1=n+n+1
            do p=1,2
               do q=1,2
                  qsca=qsca+fn1*cdabs(amat(p,q))*cdabs(amat(p,q))
               enddo
               qext1(n)=qext1(n)-fn1*dble(amat(p,p))
            enddo
            qext=qext+qext1(n)
         enddo

         if(qeps.gt.0.d0) then
            qextt=qext
            qext=0.
            do n=1,nstop
               qext=qext+qext1(n)
               err=abs(1.d0-qext/qextt)
               if(err.lt.qeps.or.n.eq.nstop) exit
            enddo
            nodr0=n
         endif
         qsca=2./x/x*qsca
         qext=2./x/x*qext
         qabs=qext-qsca
         nstop=min(n,nstop)
         return
         end subroutine mieoa

!
!  multiplies coefficients for sphere i by appropriate lm coefficient.
!  lr, oa model
!  april 2012
!
         subroutine onemiecoeffmult(i,nodr,cx,cy,mie_coefficient)
         use spheredata
         implicit none
         integer :: i,n,p,nodr,n1,n2,nterms
         complex(8) :: cx(0:nodr+1,nodr,2),cy(0:nodr+1,nodr,2)
         complex(8), allocatable :: an1(:,:,:)
         character*1, optional :: mie_coefficient
         character*1 :: miecoefficient

         if(present(mie_coefficient)) then
            miecoefficient=mie_coefficient
         else
            miecoefficient='a'
         endif
         nterms=4*nodr
         n1=mie_offset(i)+1
         n2=mie_offset(i)+nterms
         allocate(an1(2,2,nodr))
         if(miecoefficient.eq.'a') then
            an1=reshape(an_mie(n1:n2),(/2,2,nodr/))
         elseif(miecoefficient.eq.'c') then
            an1=reshape(cn_mie(n1:n2),(/2,2,nodr/))
         elseif(miecoefficient.eq.'d') then
            an1=reshape(dn_mie(n1:n2),(/2,2,nodr/))
         elseif(miecoefficient.eq.'u') then
            an1=reshape(un_mie(n1:n2),(/2,2,nodr/))
         elseif(miecoefficient.eq.'v') then
            an1=reshape(vn_mie(n1:n2),(/2,2,nodr/))
         elseif(miecoefficient.eq.'i') then
            an1=reshape(an_inv_mie(n1:n2),(/2,2,nodr/))
         endif
         do n=1,nodr
            do p=1,2
               cy(n+1,n:1:-1,p)=an1(p,1,n)*cx(n+1,n:1:-1,1) &
                 +an1(p,2,n)*cx(n+1,n:1:-1,2)
               cy(0:n,n,p)=an1(p,1,n)*cx(0:n,n,1) &
                 +an1(p,2,n)*cx(0:n,n,2)
            enddo
         enddo
         deallocate(an1)
         end subroutine onemiecoeffmult
!
! generalized mie coefficient mult:
!  (a,f) = (generalized mie matrix)*(g,b)
! idir not = 1 does the transpose.
! aout is written over in this one.
! february 2013: tmatrix file option
!
         subroutine multmiecoeffmult(neqns,nrhs,idir,ain,aout,rhs_list)
         use spheredata
         implicit none
         integer :: neqns,i,n,p,q,nodri,nblki,n1,n2,b11,b12,b21,b22,idir, &
                    nterms,j,nrhs
         complex(8) :: ain(neqns,nrhs),aout(neqns,nrhs)
         complex(8), allocatable :: gin_t(:,:,:),aout_t(:,:,:), &
                     an1(:,:,:),dn1(:,:,:),un1(:,:,:),vn1(:,:,:), &
                     bin_t(:,:,:),fout_t(:,:,:)
         logical, optional :: rhs_list(nrhs)
         logical :: rhslist(nrhs)
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif

         do j=1,nrhs
            do i=1,number_spheres
               nodri=sphere_order(i)
               nblki=2*nodri*(nodri+2)
               allocate(gin_t(0:nodri+1,nodri,2),aout_t(0:nodri+1,nodri,2), &
                        an1(2,2,nodri))
               b11=sphere_offset(i)+1
               b12=sphere_offset(i)+nblki
               gin_t=reshape(ain(b11:b12,j),(/nodri+2,nodri,2/))
               nterms=4*nodri
               n1=mie_offset(i)+1
               n2=mie_offset(i)+nterms
               an1=reshape(an_mie(n1:n2),(/2,2,nodri/))
               if(number_field_expansions(i).eq.1) then
                  aout_t=0.d0
                  if(rhslist(j).and.idir.eq.1) then
                     do n=1,nodri
                        do p=1,2
                           do q=1,2
                              aout_t(n+1,n:1:-1,p)=aout_t(n+1,n:1:-1,p) &
                                  +an1(p,q,n)*gin_t(n+1,n:1:-1,q)
                              aout_t(0:n,n,p)=aout_t(0:n,n,p)&
                                  +an1(p,q,n)*gin_t(0:n,n,q)
                            enddo
                        enddo
                     enddo
                  elseif(rhslist(j).and.idir.ne.1) then
                     do n=1,nodri
                        do p=1,2
                           do q=1,2
                              aout_t(n+1,n:1:-1,p)=aout_t(n+1,n:1:-1,p) &
                                  +an1(q,p,n)*gin_t(n+1,n:1:-1,q)
                              aout_t(0:n,n,p)=aout_t(0:n,n,p)&
                                  +an1(q,p,n)*gin_t(0:n,n,q)
                           enddo
                        enddo
                     enddo
                  endif
                  aout(b11:b12,j)&
                     =reshape(aout_t(0:nodri+1,1:nodri,1:2),(/nblki/))
               else
                  allocate(bin_t(0:nodri+1,nodri,2),&
                      fout_t(0:nodri+1,nodri,2), &
                      dn1(2,2,nodri),un1(2,2,nodri),vn1(2,2,nodri))
                  b21=sphere_offset(i)+nblki+1
                  b22=sphere_offset(i)+2*nblki
                  bin_t=reshape(ain(b21:b22,j),(/nodri+2,nodri,2/))
                  dn1=reshape(dn_mie(n1:n2),(/2,2,nodri/))
                  un1=reshape(un_mie(n1:n2),(/2,2,nodri/))
                  vn1=reshape(vn_mie(n1:n2),(/2,2,nodri/))
                  aout_t=0.d0
                  fout_t=0.d0
                  if(rhslist(j).and.idir.eq.1) then
                     do n=1,nodri
                        do p=1,2
                           do q=1,2
                              aout_t(n+1,n:1:-1,p)=aout_t(n+1,n:1:-1,p) &
                                + an1(p,q,n)*gin_t(n+1,n:1:-1,q) &
                                + un1(p,q,n)*bin_t(n+1,n:1:-1,q)
                              aout_t(0:n,n,p)=aout_t(0:n,n,p) &
                                + an1(p,q,n)*gin_t(0:n,n,q) &
                                + un1(p,q,n)*bin_t(0:n,n,q)
                              fout_t(n+1,n:1:-1,p)=fout_t(n+1,n:1:-1,p) &
                                + dn1(p,q,n)*gin_t(n+1,n:1:-1,q) &
                                + vn1(p,q,n)*bin_t(n+1,n:1:-1,q)
                              fout_t(0:n,n,p)=fout_t(0:n,n,p) &
                                + dn1(p,q,n)*gin_t(0:n,n,q) &
                                + vn1(p,q,n)*bin_t(0:n,n,q)
                           enddo
                        enddo
                     enddo
                  elseif(rhslist(j).and.idir.ne.1) then
                     do n=1,nodri
                        do p=1,2
                           do q=1,2
                              aout_t(n+1,n:1:-1,p)=aout_t(n+1,n:1:-1,p) &
                                + an1(q,p,n)*gin_t(n+1,n:1:-1,q) &
                                + dn1(q,p,n)*bin_t(n+1,n:1:-1,q)
                              aout_t(0:n,n,p)=aout_t(0:n,n,p) &
                                + an1(q,p,n)*gin_t(0:n,n,q) &
                                + dn1(q,p,n)*bin_t(0:n,n,q)
                              fout_t(n+1,n:1:-1,p)=fout_t(n+1,n:1:-1,p) &
                                + un1(q,p,n)*gin_t(n+1,n:1:-1,q) &
                                + vn1(q,p,n)*bin_t(n+1,n:1:-1,q)
                              fout_t(0:n,n,p)=fout_t(0:n,n,p) &
                                + un1(q,p,n)*gin_t(0:n,n,q) &
                                + vn1(q,p,n)*bin_t(0:n,n,q)

                           enddo
                        enddo
                     enddo
                  endif
                  aout(b11:b12,j)&
                     =reshape(aout_t(0:nodri+1,1:nodri,1:2),(/nblki/))
                  aout(b21:b22,j)&
                     =reshape(fout_t(0:nodri+1,1:nodri,1:2),(/nblki/))
                  deallocate(bin_t,fout_t,un1,vn1,dn1)
               endif
               deallocate(gin_t,aout_t,an1)
            enddo
         enddo
         end subroutine multmiecoeffmult

      end module mie
!
!  module translation contains subroutines for VSWF translation and rotation
!
!
!  last revised: 15 January 2011
!  february 2013: generalized sphere configurations
!  march 2013: new mpi group convention
!
      module translation
      use mpidefs
      use intrinsics
      use numconstants
      use specialfuncs
      use surface_subroutines
      use periodic_lattice_subroutines
      use spheredata
      use mie

      implicit none
      type translation_data
         logical :: matrix_calculated,rot_op,zero_translation
         integer :: vswf_type
         real(8) :: translation_vector(3)
         real(8), pointer :: rot_mat(:,:)
         complex(8) :: refractive_index(2)
         complex(8), pointer :: phi_mat(:),z_mat(:),gen_mat(:,:,:)
      end type translation_data
      type surface_ref_data
         logical :: symmetrical
         integer :: row_order,col_order
         complex(8), pointer :: matrix(:)
      end type surface_ref_data
      type pl_translation_data
         integer :: row_order,col_order
         complex(8), pointer :: matrix(:)
      end type pl_translation_data
      real(8), target :: interaction_radius
      type(translation_data), target,allocatable :: stored_trans_mat(:)
      type(surface_ref_data), target,allocatable :: stored_ref(:)
      type(pl_translation_data), target,allocatable :: stored_plmat(:)
      data interaction_radius/1.d10/

      contains

         subroutine clear_stored_trans_mat(mat)
         implicit none
         integer :: n,i
         type(translation_data), target, allocatable :: mat(:)
         if(.not.allocated(mat)) return
         n=size(mat)
if(light_up) then
write(*,'('' cstm 1 '',3i10)') mstm_global_rank,n
call flush(6)
endif
         do i=1,n
            if(.not.mat(i)%zero_translation) then
               if(mat(i)%rot_op) then
                  if(associated(mat(i)%rot_mat)) deallocate(mat(i)%rot_mat)
                  nullify(mat(i)%rot_mat)
                  if(associated(mat(i)%phi_mat)) deallocate(mat(i)%phi_mat)
                  nullify(mat(i)%phi_mat)
                  if(associated(mat(i)%z_mat)) deallocate(mat(i)%z_mat)
                  nullify(mat(i)%z_mat)
               else
                  if(associated(mat(i)%gen_mat)) deallocate(mat(i)%gen_mat)
                  nullify(mat(i)%gen_mat)
               endif
            endif
         enddo
if(light_up) then
write(*,'('' cstm 2 '',3i10)') mstm_global_rank,n
call flush(6)
endif
         deallocate(mat)
         end subroutine clear_stored_trans_mat

         subroutine clear_stored_ref_mat(mat)
         implicit none
         integer :: n,i
         type(surface_ref_data), allocatable :: mat(:)
         if(.not.allocated(mat)) return
         n=size(mat)
         do i=1,n
            if(associated(mat(i)%matrix)) deallocate(mat(i)%matrix)
         enddo
         deallocate(mat)
         end subroutine clear_stored_ref_mat

         subroutine clear_stored_pl_mat(mat)
         implicit none
         integer :: n,i
         type(pl_translation_data), allocatable :: mat(:)
         if(.not.allocated(mat)) return
         n=size(mat)
         do i=1,n
            if(associated(mat(i)%matrix)) deallocate(mat(i)%matrix)
         enddo
         deallocate(mat)
         end subroutine clear_stored_pl_mat

         subroutine general_interaction_matrix(matrix,mie_mult,mpi_comm)
         implicit none
         logical :: miemult
         logical, optional :: mie_mult
         integer :: j,i,nbi,nbj,ilay,jlay,vtype,i1,i2,j1,j2,mpicomm,rank,numprocs,task,nsend
         integer, optional :: mpi_comm
         real(8) :: rp(3),rdist
         complex(8) :: matrix(number_eqns,number_eqns),rimedium(2)
         complex(8), allocatable :: fsmat(:,:,:),acmat(:,:),anp(:)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(mie_mult)) then
            miemult=mie_mult
         else
            miemult=.true.
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)

         matrix=0.d0
         task=0
         do i=1,number_spheres
            nbi=sphere_order(i)*(sphere_order(i)+2)
            ilay=layer_id(sphere_position(3,i))
            do j=1,number_spheres
               if(host_sphere(i).ne.host_sphere(j).and.host_sphere(j).ne.i.and.host_sphere(i).ne.j) cycle
               task=task+1
               if(mod(task,numprocs).ne.rank) cycle
               rp=sphere_position(:,i)-sphere_position(:,j)
               rdist=sqrt(sum(rp**2))
               nbj=sphere_order(j)*(sphere_order(j)+2)
               jlay=layer_id(sphere_position(3,j))
               allocate(acmat(2*nbi,2*nbj))
               acmat=0.d0
               if(host_sphere(i).ne.0.or.host_sphere(j).ne.0) then
                  if(host_sphere(j).eq.i) then
                     rimedium=sphere_ref_index(:,i)
                     vtype=1
                     i1=sphere_offset(i)+sphere_block(i)+1
                     i2=sphere_offset(i)+2*sphere_block(i)
                     j1=sphere_offset(j)+1
                     j2=sphere_offset(j)+sphere_block(j)
                  elseif(host_sphere(i).eq.j) then
                     rimedium=sphere_ref_index(:,j)
                     vtype=1
                     j1=sphere_offset(j)+sphere_block(j)+1
                     j2=sphere_offset(j)+2*sphere_block(j)
                     i1=sphere_offset(i)+1
                     i2=sphere_offset(i)+sphere_block(i)
                  elseif(host_sphere(i).eq.host_sphere(j)) then
                     if(rdist.gt.interaction_radius) cycle
                     rimedium=sphere_ref_index(:,host_sphere(i))
                     i1=sphere_offset(i)+1
                     i2=sphere_offset(i)+sphere_block(i)
                     j1=sphere_offset(j)+1
                     j2=sphere_offset(j)+sphere_block(j)
                     vtype=3
                  else
                     cycle
                  endif
                  allocate(fsmat(nbi,nbj,2))
                  call gentranmatrix(sphere_order(j),sphere_order(i),translation_vector=rp, &
                     refractive_index=rimedium,ac_matrix=fsmat,vswf_type=vtype, &
                     mode_s=2,mode_t=2,index_model=2)
                  acmat(1:nbi,1:nbj)=fsmat(1:nbi,1:nbj,1)
                  acmat(nbi+1:2*nbi,nbj+1:2*nbj)=fsmat(1:nbi,1:nbj,2)
                  deallocate(fsmat)
               else
                  if(rdist.gt.interaction_radius) cycle
                  i1=sphere_offset(i)+1
                  i2=sphere_offset(i)+sphere_block(i)
                  j1=sphere_offset(j)+1
                  j2=sphere_offset(j)+sphere_block(j)
                  if(periodic_lattice) then
                     if(plane_surface_present) then
                        call plane_boundary_lattice_interaction(sphere_order(i),sphere_order(j), &
                           rp(1),rp(2),sphere_position(3,i),sphere_position(3,j), &
                           acmat,include_source=.true.,lr_transformation=.true.,index_model=2)
                     else
                        allocate(fsmat(nbi,nbj,2))
                        call plane_boundary_lattice_interaction(sphere_order(i),sphere_order(j), &
                           rp(1),rp(2),sphere_position(3,i),sphere_position(3,j), &
                           fsmat,include_source=.true.,lr_transformation=.true.,index_model=2)
                        acmat(1:nbi,1:nbj)=acmat(1:nbi,1:nbj)+fsmat(1:nbi,1:nbj,1)
                        acmat(nbi+1:2*nbi,nbj+1:2*nbj)=acmat(nbi+1:2*nbi,nbj+1:2*nbj)+fsmat(1:nbi,1:nbj,2)
                        deallocate(fsmat)
                     endif
                  else
                     if(plane_surface_present) then
                        call plane_interaction(sphere_order(i),sphere_order(j), &
                           rp(1),rp(2),sphere_position(3,j),sphere_position(3,i), &
                           acmat,index_model=2,lr_transformation=.true., &
                           make_symmetric=.false.)
                     endif
                     if(ilay.eq.jlay) then
                        allocate(fsmat(nbi,nbj,2))
                        call exteriorrefindex(i,rimedium)
                        call gentranmatrix(sphere_order(j),sphere_order(i),translation_vector=rp, &
                                refractive_index=rimedium,ac_matrix=fsmat,vswf_type=3, &
                                mode_s=2,mode_t=2,index_model=2)
                        acmat(1:nbi,1:nbj)=acmat(1:nbi,1:nbj)+fsmat(1:nbi,1:nbj,1)
                        acmat(nbi+1:2*nbi,nbj+1:2*nbj)=acmat(nbi+1:2*nbi,nbj+1:2*nbj)+fsmat(1:nbi,1:nbj,2)
                        deallocate(fsmat)
                     endif
                  endif
               endif
               matrix(i1:i2,j1:j2)=acmat(1:2*nbi,1:2*nbj)
               deallocate(acmat)
            enddo
         enddo
         if(numprocs.gt.1) then
            nsend=number_eqns*number_eqns
            call mstm_mpi(mpi_command='reduce',mpi_operation=mstm_mpi_sum,mpi_rank=0, &
               mpi_recv_buf_dc=matrix,mpi_number=nsend,mpi_comm=mpicomm)
         endif
         if(rank.eq.0) then
            if(miemult) then
               allocate(anp(number_eqns))
               do j=1,number_eqns
                  if(any(cdabs(matrix(1:number_eqns,j)).ne.0.d0)) then
                    call multmiecoeffmult(number_eqns,1,1,matrix(1:number_eqns,j),anp)
                    matrix(1:number_eqns,j)=-anp(1:number_eqns)
                  endif
                  matrix(j,j)=matrix(j,j)+1.d0
               enddo
               deallocate(anp)
            endif
         endif
         end subroutine general_interaction_matrix

         subroutine periodic_lattice_sphere_interaction(neqns,nrhs,ain,aout, &
            initial_run,rhs_list,mpi_comm,con_tran,store_matrix_option)
         implicit none
         logical :: initrun,rhslist(nrhs),calcmat,contran(nrhs),smopt
         logical, optional :: initial_run,rhs_list(nrhs),con_tran(nrhs),store_matrix_option
         integer :: neqns,rank,numprocs,nrhs,nmat,mpicomm,proc,i,j,rhs, &
                    i1,i2,j1,j2,task,rank0,nbi,nbj,rmatdim
         integer, save :: nmat_tot
         integer, optional :: mpi_comm
         real(8) :: rp(3),time1,time2
         complex(8) :: ain(neqns,nrhs),aout(neqns,nrhs),ri
         type(pl_translation_data), target :: rmat
         type(pl_translation_data), pointer :: loc_rmat

         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(initial_run)) then
            initrun=initial_run
         else
            initrun=.false.
         endif
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif
         if(present(con_tran)) then
            contran=con_tran
         else
            contran=.false.
         endif
         if(present(store_matrix_option)) then
            smopt=store_matrix_option
         else
            smopt=.true.
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank0)
         calcmat=initrun
         if(.not.(smopt.and.store_translation_matrix)) calcmat=.true.

         if(calcmat) then
            task=0
            nmat=0
            do i=1,number_spheres
               do j=1,number_spheres
                  if((host_sphere(i).eq.0).and.(host_sphere(j).eq.0)) then
                     task=task+1
                     proc=mod(task,numprocs)
                     if(proc.eq.rank) then
                        nmat=nmat+1
                     endif
                  endif
               enddo
            enddo
            nmat_tot=nmat
            if(smopt.and.store_translation_matrix) then
               call clear_stored_pl_mat(stored_plmat)
               allocate(stored_plmat(nmat))
            endif
         endif

         task=0
         nmat=0
         aout=0.d0
         time1=mstm_mpi_wtime()
         do i=1,number_spheres
            i1=sphere_offset(i)+1
            i2=sphere_offset(i)+sphere_block(i)
            nbi=sphere_order(i)*(sphere_order(i)+2)
            do j=1,number_spheres
               if((host_sphere(i).eq.0).and.(host_sphere(j).eq.0)) then
                  task=task+1
                  proc=mod(task,numprocs)
                  if(proc.eq.rank) then
                     j1=sphere_offset(j)+1
                     j2=sphere_offset(j)+sphere_block(j)
                     nbj=sphere_order(j)*(sphere_order(j)+2)
                     nmat=nmat+1
                     if(calcmat) then
                        rp=sphere_position(:,i)-sphere_position(:,j)
                        ri=layer_ref_index(sphere_layer(i))
                        if(plane_surface_present) then
                           rmatdim=4*nbi*nbj
                        else
                           rmatdim=2*nbi*nbj
                        endif
                        if(smopt.and.store_translation_matrix) then
                           stored_plmat(nmat)%row_order=sphere_order(i)
                           stored_plmat(nmat)%col_order=sphere_order(j)
                           allocate(stored_plmat(nmat)%matrix(rmatdim))
                           loc_rmat=>stored_plmat(nmat)
                        else
                           rmat%row_order=sphere_order(i)
                           rmat%col_order=sphere_order(j)
                           allocate(rmat%matrix(rmatdim))
                           loc_rmat=>rmat
                        endif
                        call plane_boundary_lattice_interaction(sphere_order(i),sphere_order(j), &
                           rp(1),rp(2),sphere_position(3,i),sphere_position(3,j), &
                           loc_rmat%matrix,include_source=.true.,lr_transformation=.true.,index_model=2)
                        if(smopt.and.store_translation_matrix.and.rank0.eq.0) then
                           time2=mstm_mpi_wtime()
                           if(time2-time1.ge.15.d0) then
                              write(run_print_unit,'('' assembling pl matrix '',i5,''/'',i5)') &
                                 nmat,nmat_tot
                              call flush(run_print_unit)
                              time1=time2
                           endif
                        endif
!if(mstm_global_rank.le.3) then
! write(*,'('' rank,comm:'',4i12,l2)') rank,mstm_global_rank,i,j,contran(1)
!call flush(6)
!endif
                    else
                        loc_rmat=>stored_plmat(nmat)
                     endif
                     do rhs=1,nrhs
                        if(.not.rhslist(rhs)) cycle
                        if(plane_surface_present) then
                           if(.not.contran(rhs)) then
                              call pl_matrix_mult(sphere_order(i),sphere_order(j),ain(j1:j2,rhs),aout(i1:i2,rhs), &
                                 .false.,pb_mat=loc_rmat%matrix)
                           else
                              call pl_matrix_mult(sphere_order(i),sphere_order(j),aout(j1:j2,rhs),ain(i1:i2,rhs), &
                                 .true.,pb_mat=loc_rmat%matrix)
                           endif
                        else
                           if(.not.contran(rhs)) then
                              call pl_matrix_mult(sphere_order(i),sphere_order(j),ain(j1:j2,rhs),aout(i1:i2,rhs), &
                                 .false.,fs_mat=loc_rmat%matrix)
                           else
                              call pl_matrix_mult(sphere_order(i),sphere_order(j),aout(j1:j2,rhs),ain(i1:i2,rhs), &
                                 .true.,fs_mat=loc_rmat%matrix)
                           endif
                        endif
                     enddo
                     if(.not.(smopt.and.store_translation_matrix)) deallocate(rmat%matrix)
                  endif
               endif
            enddo
         enddo
         if(store_surface_matrix) recalculate_surface_matrix=.false.
         end subroutine periodic_lattice_sphere_interaction

         subroutine pl_matrix_mult(nodrt,nodrs,as,at,tran,fs_mat,pb_mat)
         implicit none
         logical :: tran
         integer :: nodrt,nodrs,p
         complex(8) :: as(nodrs*(nodrs+2),2),at(nodrt*(nodrt+2),2)
         complex(8), optional :: fs_mat(nodrt*(nodrt+2),nodrs*(nodrs+2),2),pb_mat(nodrt*(nodrt+2),2,nodrs*(nodrs+2),2)

         if(.not.tran) then
            if(present(fs_mat)) then
               do p=1,2
                  at(:,p)=at(:,p)+matmul(fs_mat(:,:,p),as(:,p))
               enddo
            else
               do p=1,2
                  at(:,p)=at(:,p)+matmul(pb_mat(:,p,:,1),as(:,1))+matmul(pb_mat(:,p,:,2),as(:,2))
               enddo
            endif
         else
            if(present(fs_mat)) then
               do p=1,2
                  as(:,p)=as(:,p)+matmul(at(:,p),fs_mat(:,:,p))
               enddo
            else
               do p=1,2
                  as(:,p)=as(:,p)+matmul(at(:,1),pb_mat(:,1,:,p))+matmul(at(:,2),pb_mat(:,2,:,p))
               enddo
            endif
         endif
         end subroutine pl_matrix_mult

         subroutine spheresurfaceinteraction(neqns,nrhs,ain,aout, &
            initial_run,rhs_list,mpi_comm,con_tran)
         implicit none
         logical :: initrun,rhslist(nrhs),calcmat,contran(nrhs),rmatsymm
         logical, optional :: initial_run,rhs_list(nrhs),con_tran(nrhs)
         integer :: neqns,rank,numprocs,nrhs,nmat,mpicomm,proc,i,j,rhs, &
                    i1,i2,j1,j2,task,jstart,rmatdim,rank0
         integer, save :: nmat_tot
         integer, optional :: mpi_comm
         real(8) :: rp(3),time1,time2
         complex(8) :: ain(neqns,nrhs),aout(neqns,nrhs)
         complex(8), allocatable :: atempi(:),atempj(:), &
            atempi2(:),atempj2(:)
         type(surface_ref_data), target :: rmat
         type(surface_ref_data), pointer :: loc_rmat

         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(initial_run)) then
            initrun=initial_run
         else
            initrun=.false.
         endif
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif
         if(present(con_tran)) then
            contran=con_tran
         else
            contran=.false.
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank0)
         calcmat=initrun.and.recalculate_surface_matrix
         if(.not.store_surface_matrix) calcmat=.true.

         if(calcmat) then
            task=0
            nmat=0
            do i=1,number_spheres
               if(one_side_only) then
                  jstart=i
               else
                  jstart=1
               endif
               do j=jstart,number_spheres
                  if((host_sphere(i).eq.0).and.(host_sphere(j).eq.0)) then
                     task=task+1
                     proc=mod(task,numprocs)
                     if(proc.eq.rank) then
                        nmat=nmat+1
                     endif
                  endif
               enddo
            enddo
            nmat_tot=nmat
            if(store_surface_matrix) then
               call clear_stored_ref_mat(stored_ref)
               allocate(stored_ref(nmat))
            endif
         endif

         task=0
         nmat=0
         aout=0.d0
         time1=mstm_mpi_wtime()
         do i=1,number_spheres
            i1=sphere_offset(i)+1
            i2=sphere_offset(i)+sphere_block(i)
            if(one_side_only) then
               jstart=i
            else
               jstart=1
            endif
            do j=jstart,number_spheres
               if((host_sphere(i).eq.0).and.(host_sphere(j).eq.0)) then
                  task=task+1
                  proc=mod(task,numprocs)
                  if(proc.eq.rank) then
                     j1=sphere_offset(j)+1
                     j2=sphere_offset(j)+sphere_block(j)
                     nmat=nmat+1
                     if(calcmat) then
!if(mstm_global_rank.eq.0) then
!write(*,'(2i6,l2)') i,j,calcmat
!call flush(6)
!endif
                        rp=sphere_position(:,i)-sphere_position(:,j)
                        rmatsymm=(abs(rp(1)).lt.1.d-4.and.abs(rp(2)).lt.1.d-4 &
                        .and.(.not.periodic_lattice))
                        if(rmatsymm) then
                           rmatdim=2*atcdim(sphere_order(i),sphere_order(j))
                        else
                           rmatdim=sphere_block(i)*sphere_block(j)
                        endif
                        if(store_surface_matrix) then
                           stored_ref(nmat)%row_order=sphere_order(i)
                           stored_ref(nmat)%col_order=sphere_order(j)
                           stored_ref(nmat)%symmetrical=rmatsymm
                           allocate(stored_ref(nmat)%matrix(rmatdim))
                           loc_rmat=>stored_ref(nmat)
                        else
                           rmat%row_order=sphere_order(i)
                           rmat%col_order=sphere_order(j)
                           rmat%symmetrical=rmatsymm
                           allocate(rmat%matrix(rmatdim))
                           loc_rmat=>rmat
                        endif
                        call plane_interaction(sphere_order(i),sphere_order(j), &
                           rp(1),rp(2),sphere_position(3,j),sphere_position(3,i), &
                           loc_rmat%matrix, &
                           index_model=2,lr_transformation=.true., &
                           make_symmetric=rmatsymm)
                        if(store_surface_matrix.and.rank0.eq.0) then
                           time2=mstm_mpi_wtime()
                           if(time2-time1.ge.15.d0) then
                              write(run_print_unit,'('' assembling surf matrix '',i5,''/'',i5)') &
                                 nmat,nmat_tot
                              call flush(run_print_unit)
                              time1=time2
                           endif
                        endif
                     else
                        loc_rmat=>stored_ref(nmat)
                     endif
                     allocate(atempj(sphere_block(j)),atempi(sphere_block(i)))
                     if((j.ne.i).and.one_side_only) allocate(atempj2(sphere_block(j)),atempi2(sphere_block(i)))
                     do rhs=1,nrhs
                        if(.not.rhslist(rhs)) cycle
                        if(.not.contran(rhs)) then
                           call surface_interaction_matrix_mult(sphere_order(j),sphere_order(i),ain(j1:j2,rhs),atempi,loc_rmat,1)
                           aout(i1:i2,rhs)=aout(i1:i2,rhs)+atempi(1:sphere_block(i))
                        else
                           call surface_interaction_matrix_mult(sphere_order(i),sphere_order(j),ain(i1:i2,rhs),atempj,loc_rmat,2)
                           aout(j1:j2,rhs)=aout(j1:j2,rhs)+atempj(1:sphere_block(j))
                        endif
                        if((j.ne.i).and.one_side_only) then
                           if(.not.contran(rhs)) then
                              call degree_transformation(sphere_order(i), &
                                 ain(i1:i2,rhs),atempi)
                                 call surface_interaction_matrix_mult(sphere_order(i),sphere_order(j),atempi,atempj,loc_rmat,2)
                              call degree_transformation(sphere_order(j), &
                                 atempj,atempj2)
                              aout(j1:j2,rhs)=aout(j1:j2,rhs) &
                                 +atempj2(1:sphere_block(j))
                           else
                              call degree_transformation(sphere_order(j), &
                                 ain(j1:j2,rhs),atempj)
                              call surface_interaction_matrix_mult(sphere_order(j),sphere_order(i),atempj,atempi,loc_rmat,1)
                              call degree_transformation(sphere_order(i), &
                                 atempi,atempi2)
                              aout(i1:i2,rhs)=aout(i1:i2,rhs)+atempi2(1:sphere_block(i))
                           endif
                        endif
                     enddo
                     deallocate(atempi,atempj)
                     if((j.ne.i).and.one_side_only) deallocate(atempi2,atempj2)
                     if(.not.store_surface_matrix) deallocate(rmat%matrix)
                  endif
               endif
            enddo
         enddo
         if(store_surface_matrix) recalculate_surface_matrix=.false.
         end subroutine spheresurfaceinteraction

         subroutine surface_interaction_matrix_mult(nin,nout,ain,aout,rmat,dir)
         implicit none
         integer :: nin,nout,bin,bout,m,dir,m1,n,l,p,q,mnp,klq,i,moff
         complex(8) :: ain(2*nin*(nin+2)),aout(2*nout*(nout+2))
         type(surface_ref_data), pointer :: rmat
         bin=2*nin*(nin+2)
         bout=2*nout*(nout+2)
         aout=0.d0
         if(rmat%symmetrical) then
            do m=-min(nin,nout),min(nin,nout)
               m1=max(abs(m),1)
               moff=2*moffset(m,nin,nout)
               do n=m1,nout
                  do p=1,2
                     mnp=amnpaddress(m,n,p,nout,2)
                     do l=m1,nin
                        do q=1,2
                           klq=amnpaddress(m,l,q,nin,2)
                           if(dir.eq.1) then
                              i=moff+p+2*(n-m1)+2*(nout-m1+1)*(q-1+2*(l-m1))
                           else
                              i=moff+q+2*(l-m1)+2*(nin-m1+1)*(p-1+2*(n-m1))
                           endif
                           aout(mnp)=aout(mnp)+rmat%matrix(i)*ain(klq)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         else
            if(dir.eq.1) then
               do n=1,bout
                  do l=1,bin
                     i=n+bout*(l-1)
                     aout(n)=aout(n)+rmat%matrix(i)*ain(l)
                  enddo
               enddo
            else
               do n=1,bout
                  do l=1,bin
                     i=l+bin*(n-1)
                     aout(n)=aout(n)+rmat%matrix(i)*ain(l)
                  enddo
               enddo
            endif
         endif
         end subroutine surface_interaction_matrix_mult
!
! outgoing translation operation:  a(i) = H(i-j) a(j).
! February 2013: number of rhs is a required argument.   mpi comm option added.
! this does not perform an allgather on output arrays.  that operation will be needed
! to use the results
!
         subroutine external_to_external_expansion(neqns,nrhs,ain,gout, &
                    store_matrix_option,initial_run,rhs_list, &
                    mpi_comm,con_tran)
         implicit none
         logical :: smopt,rhslist(nrhs),contran(nrhs),rot
         logical, save :: calcmat,firstrun
         logical, optional :: store_matrix_option,initial_run, &
                              rhs_list(nrhs),con_tran(nrhs)
         integer :: neqns,rank,numprocs,nsphere,nrhs,mpicomm, &
                    i,j,npi1,npi2,npj1,npj2,noj,noi,task,proc, &
                    nmin,nmax,ndim,tdim,idim,rhs
         integer, optional :: mpi_comm
         real(8) :: rdist
         complex(8)  :: ain(neqns,nrhs),gout(neqns,nrhs),rimedium(2)
         type(translation_data), pointer :: loc_tranmat
         type(translation_data), target :: tranmat
         data firstrun/.true./
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nsphere=number_spheres
         gout=0.
         if(present(store_matrix_option)) then
            smopt=store_matrix_option
         else
            smopt=.true.
         endif
         if(present(initial_run)) then
            firstrun=initial_run
         endif
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif
         if(present(con_tran)) then
            contran=con_tran
         else
            contran=.false.
         endif
!
!  compute offsets for scattering coefficients
!
         if(firstrun) then
            if(smopt.and.store_translation_matrix) then
               task=0
               ndim=0
               do i=1,nsphere-1
                  do j=i+1,nsphere
                     if(host_sphere(j).eq.host_sphere(i) &
                      .and.sphere_layer(j).eq.sphere_layer(i)) then
                        rdist=sqrt(sum((sphere_position(:,i)-sphere_position(:,j))**2))
                        if(rdist.gt.interaction_radius) cycle
                        task=task+1
                        proc=mod(task,numprocs)
                        if(proc.eq.rank) ndim=ndim+1
                     endif
                  enddo
               enddo
               if(allocated(stored_trans_mat)) then
                  call clear_stored_trans_mat(stored_trans_mat)
               endif
               allocate(stored_trans_mat(ndim))
            endif
            calcmat=.true.
            firstrun=.false.
         else
            calcmat=(.not.(smopt.and.store_translation_matrix))
         endif

         idim=0
         task=0
         do i=1,nsphere-1
            do j=i+1,nsphere
               if(host_sphere(j).eq.host_sphere(i) &
                .and.sphere_layer(j).eq.sphere_layer(i)) then
                  rdist=sqrt(sum((sphere_position(:,i)-sphere_position(:,j))**2))
                  if(rdist.gt.interaction_radius) cycle
                  task=task+1
                  proc=mod(task,numprocs)
                  if(proc.eq.rank) then
                     idim=idim+1
                     call exteriorrefindex(i,rimedium)
                     noi=sphere_order(i)
                     npi1=sphere_offset(i)+1
                     npi2=sphere_offset(i)+sphere_block(i)
                     noj=sphere_order(j)
                     npj1=sphere_offset(j)+1
                     npj2=sphere_offset(j)+sphere_block(j)
                     if(calcmat) then
                        nmin=min(noi,noj)
                        nmax=max(noi,noj)
                        rot=(nmin.ge.translation_switch_order)
                        tdim=atcdim(noi,noj)
                        if(smopt.and.store_translation_matrix) then
                           stored_trans_mat(idim)%matrix_calculated=.false.
                           stored_trans_mat(idim)%vswf_type=3
                           stored_trans_mat(idim)%translation_vector=sphere_position(:,i)-sphere_position(:,j)
                           stored_trans_mat(idim)%refractive_index=rimedium
                           stored_trans_mat(idim)%rot_op=rot
                           loc_tranmat=>stored_trans_mat(idim)
                        else
                           tranmat%matrix_calculated=.false.
                           tranmat%vswf_type=3
                           tranmat%translation_vector=sphere_position(:,i)-sphere_position(:,j)
                           tranmat%refractive_index=rimedium
                           tranmat%rot_op=rot
                           loc_tranmat=>tranmat
                        endif
                     else
                        loc_tranmat=>stored_trans_mat(idim)
                     endif
                     do rhs=1,nrhs
                        if(.not.rhslist(rhs)) cycle
                        call coefficient_translation(noj,2,noi,2,ain(npj1:npj2,rhs),gout(npi1:npi2,rhs), &
                           loc_tranmat,shift_op=.false.,tran_op=contran(rhs))
                        call coefficient_translation(noi,2,noj,2,ain(npi1:npi2,rhs),gout(npj1:npj2,rhs), &
                           loc_tranmat,shift_op=.true.,tran_op=contran(rhs))
                     enddo
                     if(calcmat.and.(.not.(smopt.and.store_translation_matrix))) then
                        if(rot) then
                           deallocate(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                        else
                           deallocate(tranmat%gen_mat)
                        endif
                     endif
                  endif
               endif
            enddo
         enddo

         end subroutine external_to_external_expansion
!
!  calculation of bmnp(i) = J(i-j) amnp(j) for i internal, host j=i, and
!  gmnp(i) = J(i-j) f(j), for host i = j.    This is the regular translation operation.   g and b
!  are returned ordered as (a,f) in the output array.
!  february 2013: number of rhs option
!  the routine does not perform an mpi reduce.
!
         subroutine external_to_internal_expansion(neqns,nrhs,ain,bout, &
            rhs_list,mpi_comm,con_tran)
         implicit none
         logical :: rhslist(nrhs),contran(nrhs)
         logical, optional :: rhs_list(nrhs),con_tran(nrhs)
         integer :: neqns,rank,numprocs,nsphere,nrhs,mpicomm, &
                    i,j,task,proc,extsurf,intsurf,ext1,ext2,&
                    int1,int2,noext,noint,rhs
         integer :: count
         integer, optional :: mpi_comm
         complex(8)  :: ain(neqns,nrhs),bout(neqns,nrhs),rimedium(2)
         type(translation_data), pointer :: loc_tranmat
         type(translation_data), target :: tranmat
         data count/0/
         count=count+1
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(rhs_list)) then
            rhslist=rhs_list
         else
            rhslist=.true.
         endif
         if(present(con_tran)) then
            contran=con_tran
         else
            contran=.false.
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         bout=0.d0
         nsphere=number_spheres
         task=0

         do i=1,nsphere-1
            do j=i+1,nsphere
               if(host_sphere(j).eq.i.or.host_sphere(i).eq.j) then
                  task=task+1
                  proc=mod(task,numprocs)
                  if(proc.eq.rank) then
                     if(host_sphere(j).eq.i) then
                        extsurf=j
                        intsurf=i
                     else
                        extsurf=i
                        intsurf=j
                     endif
                     noext=sphere_order(extsurf)
                     noint=sphere_order(intsurf)
                     rimedium=sphere_ref_index(:,intsurf)
                     tranmat%matrix_calculated=.false.
                     tranmat%vswf_type=1
                     tranmat%translation_vector=sphere_position(:,intsurf)-sphere_position(:,extsurf)
                     tranmat%refractive_index=sphere_ref_index(:,intsurf)
                     tranmat%rot_op=.true.
                     loc_tranmat=>tranmat
                     ext1=sphere_offset(extsurf)+1
                     ext2=ext1-1+sphere_block(extsurf)
                     int1=sphere_offset(intsurf)+1+sphere_block(intsurf)
                     int2=int1-1+sphere_block(intsurf)
                     do rhs=1,nrhs
                        call coefficient_translation(noext,2,noint,2,ain(ext1:ext2,rhs), &
                           bout(int1:int2,rhs),loc_tranmat, &
                           shift_op=.false.,tran_op=contran(rhs))
                        call coefficient_translation(noint,2,noext,2,ain(int1:int2,rhs), &
                           bout(ext1:ext2,rhs),loc_tranmat, &
                           shift_op=.true.,tran_op=contran(rhs))
                     enddo
                     if(.not.tranmat%zero_translation) deallocate(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                  endif
               endif
            enddo
         enddo
         end subroutine external_to_internal_expansion

         subroutine coefficient_translation(nodra,nmodea,nodrg,nmodeg, &
             acoef,gcoef,tranmat,shift_op,tran_op)
         implicit none
         logical :: sop,top,rot
         logical, optional :: shift_op,tran_op
         integer :: nodra,nodrg,nblka,nblkg,lengtha, &
                    lengthg,nmodea,nmodeg,nmode,shiftvec(2),n,m,im,nn1,nn2,n1, &
                    p,nmin,m1,offset,blocksize,nmax,tdim,vtype,nodrs,nodrt
         real(8) :: r,ct,rtran(3)
         complex(8) :: acoef(*),gcoef(*), &
            a_t(0:nodra+1,nodra,nmodea),g_t(0:nodrg+1,nodrg,nmodeg), &
            a_tt(-nodra:nodra,nodra,2),g_tt(-nodrg:nodrg,nodrg,2), &
            atc(max(nodra,nodrg),max(nodra,nodrg),2), &
            a_t2(nodra*(nodra+2),2),g_t2(nodrg*(nodrg+2),2),rimed(2),ephi
         type(translation_data) :: tranmat
         nblka=nodra*(nodra+2)
         nblkg=nodrg*(nodrg+2)
         nmin=min(nodra,nodrg)
         nmax=max(nodra,nodrg)
         nmode=max(nmodea,nmodeg)
         if(present(shift_op)) then
            sop=shift_op
         else
            sop=.false.
         endif
         if(present(tran_op)) then
            top=tran_op
         else
            top=.false.
         endif
         rot=tranmat%rot_op

         if(.not.tranmat%matrix_calculated) then
            vtype=tranmat%vswf_type
            rimed=tranmat%refractive_index
            rtran=tranmat%translation_vector
            r=dot_product(rtran,rtran)
            if(r.lt.1.d-12) then
               tranmat%zero_translation=.true.
            else
               tranmat%zero_translation=.false.
            endif
            if(sop) then
               nodrs=nodrg
               nodrt=nodra
            else
               nodrs=nodra
               nodrt=nodrg
            endif
            if(.not.tranmat%zero_translation) then
               if(rot) then
                  tdim=atcdim(nodrt,nodrs)
                  allocate(tranmat%rot_mat(-nmin:nmin,0:nmax*(nmax+2)))
                  allocate(tranmat%phi_mat(-nmax:nmax))
                  allocate(tranmat%z_mat(1:tdim))
                  call cartosphere(rtran,r,ct,ephi)
                  call rotcoef(ct,nmin,nmax,tranmat%rot_mat)
                  call axialtrancoefrecurrence(vtype,r,rimed,nodrt,nodrs, &
                     tdim,tranmat%z_mat)
                  call ephicoef(ephi,nmax,tranmat%phi_mat)
               else
                  allocate(tranmat%gen_mat(nodrt*(nodrt+2),nodrs*(nodrs+2),2))
                  call gentranmatrix(nodrs,nodrt,translation_vector=rtran, &
                       refractive_index=rimed,ac_matrix=tranmat%gen_mat,vswf_type=vtype, &
                       mode_s=2,mode_t=2)
               endif
            endif
            tranmat%matrix_calculated=.true.
         endif
         lengtha=nblka*nmodea
         lengthg=nblkg*nmodeg
         shiftvec=(/1,1/)
         if(sop.neqv.top) then
            shiftvec=-shiftvec
         endif
         if(sop) then
            im=-1
         else
            im=1
         endif
         if(tranmat%zero_translation) then
            if(tranmat%vswf_type.eq.1) then
               call mtransfer(nodra,nodrg,acoef(1:lengtha),g_t)
               gcoef(1:lengthg)=gcoef(1:lengthg) &
                  +reshape(g_t(0:nodrg+1,1:nodrg,1:nmodeg),(/lengthg/))
            endif
         else
            if(rot) then
               call shiftcoefficient(nodra,nmodea,shiftvec(1),shiftvec(2), &
                  acoef(1:lengtha),a_t(0:nodra+1,1:nodra,1:2))
               a_tt(0,1:nodra,1:2)=a_t(0,1:nodra,1:2)
               do m=1,nodra
                  a_tt(m,m:nodra,1:2)=a_t(m,m:nodra,1:2)*tranmat%phi_mat(im*m)
                  a_tt(-m,m:nodra,1:2)=a_t(m+1:nodra+1,m,1:2)*tranmat%phi_mat(-im*m)
               enddo
               do n=1,nodra
                  nn1=n*(n+1)-n
                  nn2=nn1+(2*n+1)-1
                  n1=min(n,nodrg)
                  a_tt(-n1:n1,n,1:2) &
                     = matmul(tranmat%rot_mat(-n1:n1,nn1:nn2),a_tt(-n:n,n,1:2))
               enddo
               do m=-nmin,nmin
                  m1=max(1,abs(m))
                  if(sop) then
                     offset=moffset(m,nodra,nodrg)
                     blocksize=(nodrg-m1+1)*(nodra-m1+1)*2
                     atc(m1:nodra,m1:nodrg,1:2)= &
                         reshape(tranmat%z_mat(offset+1:offset+blocksize), &
                         (/nodra-m1+1,nodrg-m1+1,2/))
                     do p=1,2
                        g_tt(m,m1:nodrg,p) &
                           =matmul(a_tt(m,m1:nodra,p),atc(m1:nodra,m1:nodrg,p))
                     enddo
                  else
                     offset=moffset(m,nodrg,nodra)
                     blocksize=(nodrg-m1+1)*(nodra-m1+1)*2
                     atc(m1:nodrg,m1:nodra,1:2)= &
                         reshape(tranmat%z_mat(offset+1:offset+blocksize), &
                         (/nodrg-m1+1,nodra-m1+1,2/))
                     do p=1,2
                        g_tt(m,m1:nodrg,p) &
                           =matmul(atc(m1:nodrg,m1:nodra,p),a_tt(m,m1:nodra,p))
                     enddo
                  endif
               enddo
               do n=1,nodrg
                  nn1=n*(n+1)-n
                  nn2=nn1+(2*n+1)-1
                  n1=min(n,nodra)
                  g_tt(-n:n,n,1)=matmul(g_tt(-n1:n1,n,1),tranmat%rot_mat(-n1:n1,nn1:nn2))
                  g_tt(-n:n,n,2)=matmul(g_tt(-n1:n1,n,2),tranmat%rot_mat(-n1:n1,nn1:nn2))
               enddo
               g_t(0,1:nodrg,1:2)=g_tt(0,1:nodrg,1:2)
               do m=1,nodrg
                  g_t(m,m:nodrg,1:2)=g_tt(m,m:nodrg,1:2)*tranmat%phi_mat(-im*m)
                  g_t(m+1:nodrg+1,m,1:2)=g_tt(-m,m:nodrg,1:2)*tranmat%phi_mat(im*m)
               enddo
               call shiftcoefficient(nodrg,nmodeg,shiftvec(1),shiftvec(2), &
                  g_t,g_t)
            else
               call shiftcoefficient(nodra,nmodea,shiftvec(1),shiftvec(2), &
                  acoef(1:lengtha), &
                  a_t2(1:nblka,1:2))
               if(sop) then
                  g_t2(:,1)=matmul(a_t2(1:nblka,1),tranmat%gen_mat(1:nblka,1:nblkg,1))
               else
                  g_t2(:,1)=matmul(tranmat%gen_mat(1:nblkg,1:nblka,1),a_t2(1:nblka,1))
               endif
               if(nmodeg.eq.2.and.nmodea.eq.1) then
                  if(sop) then
                     g_t2(:,2)=matmul(a_t2(1:nblka,1),tranmat%gen_mat(:,:,2))
                     g_t2=0.5d0*g_t2
                  else
                     g_t2(:,2)=matmul(tranmat%gen_mat(:,:,2),a_t2(1:nblka,1))
                     g_t2=0.5d0*g_t2
                  endif
               elseif(nmodeg.eq.1.and.nmodea.eq.2) then
                  if(sop) then
                     g_t2(:,1)=g_t2(:,1)+matmul(a_t2(1:nblka,2),tranmat%gen_mat(:,:,2))
                  else
                     g_t2(:,1)=g_t2(:,1)+matmul(tranmat%gen_mat(:,:,2),a_t2(1:nblka,2))
                  endif
               elseif(nmodeg.eq.2.and.nmodea.eq.2) then
                  if(sop) then
                     g_t2(:,2)=matmul(a_t2(1:nblka,2),tranmat%gen_mat(1:nblka,1:nblkg,2))
                  else
                     g_t2(:,2)=matmul(tranmat%gen_mat(1:nblkg,1:nblka,2),a_t2(1:nblka,2))
                  endif
               endif
               call shiftcoefficient(nodrg,nmodeg,shiftvec(1),shiftvec(2), &
                  g_t2,g_t)
            endif
            gcoef(1:lengthg) &
               =gcoef(1:lengthg) &
                 +reshape(g_t(0:nodrg+1,1:nodrg,1:nmodeg),(/lengthg/))
         endif
         end subroutine coefficient_translation

         subroutine shiftcoefficient(nodr,nmode,msign,mflip, &
             ain,aout)
         implicit none
         integer :: nodr,nmode,m,n,msign,mflip,im
         complex(8) :: ain(0:nodr+1,nodr,nmode),at(nmode), &
            aout(0:nodr+1,nodr,nmode)
         if(msign.eq.1.and.mflip.eq.1) then
            aout=ain
         else
            aout(0,1:nodr,1:nmode)=ain(0,1:nodr,1:nmode)
            if(mflip.eq.-1) then
               im=1
               do m=1,nodr
                  im=im*msign
                  do n=m,nodr
                     at=ain(n+1,m,1:nmode)
                     aout(n+1,m,1:nmode)=im*ain(m,n,1:nmode)
                     aout(m,n,1:nmode)=im*at
                  enddo
               enddo
            else
               im=1
               do m=1,nodr
                  im=im*msign
                  do n=m,nodr
                     aout(m,n,1:nmode)=im*ain(m,n,1:nmode)
                     aout(n+1,m,1:nmode)=im*ain(n+1,m,1:nmode)
                  enddo
               enddo
            endif
         endif
         end subroutine shiftcoefficient

      end module translation
