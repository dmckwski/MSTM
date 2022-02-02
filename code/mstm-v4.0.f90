!
!  numerical constants
!
!
!  last revised: 15 January 2011
!
      module numconstants
      use mpidefs
      implicit none
      logical :: light_up
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
         z=r*ri
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

         subroutine lu_decomposition(a,n,indx,d,ierr)
         implicit none
         integer :: n,indx(n),i,j,k,imax,ierr
         real(8) :: tiny,aamax,d,vv(n),dum
         complex(8) ::  a(n,n),sum,cdum
         data tiny/1.d-20/
         ierr=0
         d=1.d0
         do i=1,n
            aamax=0.d0
            do j=1,n
               if (cdabs(a(i,j)).gt.aamax) aamax=cdabs(a(i,j))
            enddo
            if (aamax.eq.0.d0) then
               ierr=1
               return
            endif
            vv(i)=1.d0/aamax
         enddo
         do j=1,n
            if (j.gt.1) then
               do i=1,j-1
                  sum=a(i,j)
                  if (i.gt.1) then
                     do k=1,i-1
                        sum=sum-a(i,k)*a(k,j)
                     enddo
                     a(i,j)=sum
                  endif
               enddo
            endif
            aamax=0.d0
            do i=j,n
               sum=a(i,j)
               if (j.gt.1) then
                  do k=1,j-1
                     sum=sum-a(i,k)*a(k,j)
                  enddo
                  a(i,j)=sum
               endif
               dum=vv(i)*cdabs(sum)
               if (dum.ge.aamax) then
                  imax=i
                  aamax=dum
               endif
            enddo
            if (j.ne.imax) then
               do k=1,n
                  cdum=a(imax,k)
                  a(imax,k)=a(j,k)
                  a(j,k)=cdum
               enddo
               d=-d
               vv(imax)=vv(j)
            endif
            indx(j)=imax
            if (j.ne.n) then
               if (cdabs(a(j,j)).eq.0.d0) a(j,j)=tiny
               cdum=1.d0/a(j,j)
               do i=j+1,n
                  a(i,j)=a(i,j)*cdum
               enddo
            endif
         enddo
         if (cdabs(a(n,n)).eq.0.d0) a(n,n)=tiny
         end subroutine lu_decomposition

         subroutine lu_backsubstitution(a,n,indx,b)
         implicit none
         integer :: n,indx(n),i,ii,j,ll
         complex(8) ::  a(n,n),b(n),sum
         ii=0
         do i=1,n
            ll=indx(i)
            sum=b(ll)
            b(ll)=b(i)
            if (ii.ne.0) then
               do j=ii,i-1
                  sum=sum-a(i,j)*b(j)
               enddo
            elseif (cdabs(sum).ne.0.d0) then
               ii=i
            endif
            b(i)=sum
         enddo
         do i=n,1,-1
            sum=b(i)
            if (i.lt.n) then
               do j=i+1,n
                  sum=sum-a(i,j)*b(j)
               enddo
            endif
            b(i)=sum/a(i,i)
         enddo
         end subroutine lu_backsubstitution

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
         pole_integration
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
         minimum_integration_spacing
      real(8), allocatable :: plane_boundary_position(:),real_axis_limits(:)
      complex(8) :: source_ri,target_ri,pole_integration_s
      complex(8), target :: layer_ref_index(0:max_number_plane_boundaries)
      complex(8), allocatable :: source_coefficient(:,:),source_coefficient_1(:,:,:),source_coefficient_2(:,:,:)
      data real_axis_integration_limit,minimum_integration_spacing/100.d0,1.d-5/
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
               if(ec.eq.1) error_codes(4)=1
               if(ec.eq.2) error_codes(3)=1
               rmat=rmat+drmat
               if(dt2.ge.real_axis_integration_limit) return
            enddo
         enddo
         delt=1.d0
         do while(t2.lt.real_axis_integration_limit)
            t1=t2
            t2=t2+delt
            subdiv=0
            t1t=t1
            t2t=t2
            ec=0
            call gkintegrate(qtot,t1t,t2t,realaxiskernel,drmat,subdiv,ec, &
               integration_error_epsilon,minimum_integration_spacing,maximum_integration_subdivisions)
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
         q1d_number_segments,q2d_number_segments
      integer, target :: pl_integration_method
      real(8) :: lattice_integration_segment,pl_rs_eps,time_count(4),time_0
      real(8), target :: cell_width(2),rs_dz_min,pl_integration_error_epsilon, &
         pl_integration_limit_epsilon
      data lattice_integration_segment/1.d0/
      data pl_rs_nmax,pl_rs_eps,pl_rs_imax/200,1.d-6,0/
      data time_it,rs_dz_min/.true.,100.d0/
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
         if(number_plane_boundaries.eq.0) then
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
         wcrit=min(rs_dz_min,minval(cell_width))
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
            tranmat(2,2),vw,w
         integer, optional :: index_model
         real(8) :: x,y,zt,zs,kx,ky,cerr,asum,asum0,vc1(0:nodrt+nodrs), &
            vcp1m1(0:nodrt+nodrs),vcm1m1(0:nodrt+nodrs),x0,y0
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
         allocate(qsum(-1:1,0:nblk,tdirs(1):tdirs(2),sdirs(1):sdirs(2)), &
            dqsum(-1:1,0:nblk,tdirs(1):tdirs(2),sdirs(1):sdirs(2)))
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
            asum0=sum(cdabs(qsum))
            asum=sum(cdabs(dqsum))
            cerr=asum/asum0
            if(cerr.lt.pl_rs_eps) exit
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
         real(8) :: x,k0y,k0z,w(2),y,z,kz,mag
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
         end subroutine swfyzlatticesum

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
         real(8) :: x0,y0,z0,w(2),wx,wy,k0(2),kconst,kx,ky,eps,cerr,asum,asum0
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
            asum0=sum(cdabs(wf))
            asum=sum(cdabs(dwf))
            cerr=asum/asum0
            if(cerr.lt.eps) exit
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
      logical, target :: store_translation_matrix,plane_surface_present, &
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
         circumscribing_radius,cross_section_radius
      real(8), target :: gaussian_beam_constant,gaussian_beam_focal_point(3)
      real(8), allocatable :: qext_mie(:),qabs_mie(:)
      real(8), allocatable, target :: sphere_radius(:),sphere_position(:,:)
      complex(8), target :: medium_ref_index
      complex(8), allocatable :: an_mie(:),cn_mie(:),un_mie(:),vn_mie(:),dn_mie(:),an_inv_mie(:)
      complex(8), allocatable, target :: sphere_ref_index(:,:)

      data run_print_unit/6/
      data medium_ref_index/(1.d0,0.d0)/
      data store_translation_matrix,store_surface_matrix/.false.,.true./
      data recalculate_surface_matrix/.true./
      data translation_switch_order/3/
      data fft_translation_option/.false./
      data max_t_matrix_order/100/
      data normalize_solution_error/.true./
      data gaussian_beam_constant,gaussian_beam_focal_point/0.d0,0.d0,0.d0,0.d0/

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

         top_boundary=max(plane_boundary_position(max(1,number_plane_boundaries)),sphere_max_position(3))+1.d2
         bot_boundary=min(0.d0,sphere_min_position(3))-1.d2

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
                  if(rij.le.sphere_radius(j)-sphere_radius(i).and.sphere_radius(j).lt.xspmin) then
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
           unp_mie,vnp_mie,cnp_mie,ri_medium,anp_inv_mie)
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
                                 anp_inv_mie(2,2,*)
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
      type(translation_data), target,allocatable :: stored_trans_mat(:)
      type(surface_ref_data), target,allocatable :: stored_ref(:)
      type(pl_translation_data), target,allocatable :: stored_plmat(:)

      contains

         subroutine clear_stored_trans_mat(mat)
         implicit none
         integer :: n,i,ierr
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
      module fft_translation
      use mpidefs
      use intrinsics
      use numconstants
      use specialfuncs
      use spheredata
      use translation
      use mie
      implicit none
      type node_data
         integer :: number_elements
         type(linked_ilist), pointer :: members
      end type node_data
      type linked_ilist
         integer :: index
         type(linked_ilist), pointer :: next=>null()
      end type linked_ilist
      type coefficient_list
         complex(8), pointer :: coefficient_vector(:,:)
      end type coefficient_list

      integer :: cell_dim(3),number_neighbor_nodes
      integer, private :: neighbor_node(3,0:26),fft_local_host
      integer, target :: node_order,neighbor_node_model
      integer, allocatable, private :: sphere_node(:,:)
      real(8) :: d_cell
      real(8), private :: cell_origin(3),cell_boundary(3)
      real(8), target :: cell_volume_fraction
      complex(8), private :: host_ref_index(2)
      complex(4), allocatable, private :: cell_translation_matrix(:,:,:,:,:,:)
      type(node_data), allocatable, private :: cell_list(:,:,:),sphere_local_interaction_list(:)
      type(translation_data), target, allocatable, private :: stored_local_j_mat(:),stored_local_h_mat(:)
      data fft_local_host/0/
      data neighbor_node_model/2/
      data node_order/3/
      data cell_volume_fraction/0.2d0/


      contains

         subroutine clear_fft_matrix(clear_h)
         implicit none
         logical :: clearh
         logical, optional :: clear_h
         if(present(clear_h)) then
            clearh=clear_h
         else
            clearh=.false.
         endif
if(light_up) then
write(*,'('' fft cfm 1'',2i10,l)') mstm_global_rank,size(stored_local_j_mat),allocated(stored_local_j_mat)
call flush(6)
endif
         call clear_stored_trans_mat(stored_local_j_mat)
if(light_up) then
write(*,'('' fft cfm 2'',2i10,l)') mstm_global_rank,size(stored_local_h_mat),allocated(stored_local_h_mat)
call flush(6)
endif
         call clear_stored_trans_mat(stored_local_h_mat)
         if(clearh) then
            if(allocated(cell_translation_matrix)) deallocate(cell_translation_matrix)
         endif
         if(allocated(sphere_node)) deallocate(sphere_node)
if(light_up) then
write(*,'('' fft cfm 3'',i3,l)') mstm_global_rank,allocated(cell_translation_matrix)
call flush(6)
endif

         end subroutine clear_fft_matrix

         subroutine fft_external_to_external_expansion(neqns,nrhs,ain,gout, &
                    store_matrix_option,initial_run,rhs_list, &
                    mpi_comm,con_tran)
         implicit none
         logical :: smopt,rhslist(nrhs),contran(nrhs)
         logical, save :: firstrun,inp1,inp2
         logical, optional :: store_matrix_option,initial_run, &
                              rhs_list(nrhs),con_tran(nrhs)
         integer :: neqns,rank,numprocs,nsphere,nrhs,mpicomm,p,nsend, &
            i,rhs,noff,groupsize,mpigroup,syncgroup,oddnumproc
         integer, save :: pgroup,pcomm,synccomm1,synccomm2,p1,p2,prank
         integer, allocatable :: grouplist(:)
         integer, optional :: mpi_comm
         complex(8) :: ain(neqns,nrhs),gout(neqns,nrhs)
         complex(8), allocatable :: anode(:,:,:,:,:),gnode(:,:,:,:,:), &
            gout_loc(:,:),ain_t(:,:),gout_t(:,:)
         data firstrun/.true./
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
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
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nsphere=number_spheres
         gout=0.

         if(firstrun) then
            if(numprocs.gt.1) then
               oddnumproc=mod(numprocs,2)
               pgroup=floor(dble(2*rank)/dble(numprocs))+1
               p1=pgroup
               p2=p1
               call mstm_mpi(mpi_command='split', &
                  mpi_color=pgroup,mpi_key=rank, &
                  mpi_new_comm=pcomm, &
                  mpi_comm=mpicomm)
               call mstm_mpi(mpi_command='rank',mpi_rank=prank,mpi_comm=pcomm)
               call mstm_mpi(mpi_command='group',mpi_group=mpigroup,mpi_comm=mpicomm)
               groupsize=numprocs/2+1
               allocate(grouplist(groupsize))
               grouplist(1)=0
               do i=1,groupsize-1
                  grouplist(i+1)=i+(numprocs/2)-1+oddnumproc
               enddo
               inp1=.false.
               do i=1,groupsize
                  if(rank.eq.grouplist(i)) then
                     inp1=.true.
                     exit
                  endif
               enddo
               call mstm_mpi(mpi_command='incl', &
                  mpi_group=mpigroup,&
                  mpi_size=groupsize,&
                  mpi_new_group_list=grouplist,&
                  mpi_new_group=syncgroup)
               call mstm_mpi(mpi_command='create',&
                  mpi_group=syncgroup,&
                  mpi_comm=mpicomm,&
                  mpi_new_comm=synccomm1)
               deallocate(grouplist)
               groupsize=groupsize+oddnumproc
               allocate(grouplist(groupsize))
               grouplist(1)=numprocs/2+oddnumproc
               do i=1,groupsize-1
                  grouplist(i+1)=i-1
               enddo
               inp2=.false.
               do i=1,groupsize
                  if(rank.eq.grouplist(i)) then
                     inp2=.true.
                     exit
                  endif
               enddo
               call mstm_mpi(mpi_command='incl', &
                  mpi_group=mpigroup,&
                  mpi_size=groupsize,&
                  mpi_new_group_list=grouplist,&
                  mpi_new_group=syncgroup)
               call mstm_mpi(mpi_command='create',&
                  mpi_group=syncgroup,&
                  mpi_comm=mpicomm,&
                  mpi_new_comm=synccomm2)
               deallocate(grouplist)
            else
               pgroup=1
               p1=1
               p2=2
               pcomm=mpicomm
            endif
!            call node_selection(cell_volume_fraction)
if(light_up) then
write(*,'('' fft1 '',i3)') mstm_global_rank
call flush(6)
endif
            call fft_translation_initialization(host_ref_index,node_order,p1,p2)
         endif

         allocate(anode(cell_dim(1),cell_dim(2),cell_dim(3),node_order*(node_order+2)*2,nrhs), &
            gnode(cell_dim(1),cell_dim(2),cell_dim(3),node_order*(node_order+2)*2,nrhs), &
            gout_loc(number_eqns,nrhs),ain_t(neqns,nrhs),gout_t(neqns,nrhs))

if(light_up) then
write(*,'('' fft2 '',i3)') mstm_global_rank
call flush(6)
endif

         do rhs=1,nrhs
            noff=0
            do i=1,number_spheres
               if(contran(rhs)) then
                  call shiftcoefficient(sphere_order(i),2,-1,-1, &
                        ain(noff+1:noff+2*sphere_order(i)*(sphere_order(i)+2),rhs), &
                        ain_t(noff+1:noff+2*sphere_order(i)*(sphere_order(i)+2),rhs))
               else
                  call shiftcoefficient(sphere_order(i),2,1,1, &
                        ain(noff+1:noff+2*sphere_order(i)*(sphere_order(i)+2),rhs), &
                        ain_t(noff+1:noff+2*sphere_order(i)*(sphere_order(i)+2),rhs))
               endif
               noff=noff+2*sphere_order(i)*(sphere_order(i)+2)*number_field_expansions(i)
            enddo
         enddo

         anode=0.d0
         gnode=0.d0
         gout_t=0.d0
         gout_loc=0.d0
if(light_up) then
write(*,'('' fft3 '',i3)') mstm_global_rank
call flush(6)
endif

         call local_sphere_to_node_translation(nrhs,ain_t,anode, &
            store_matrix_option=smopt,initial_run=firstrun, &
            mpi_comm=mpicomm,local_host=fft_local_host,sphere_to_node=.true., &
            merge_procs=.true.)
if(light_up) then
write(*,'('' fft4 '',i3)') mstm_global_rank
call flush(6)
endif

         do p=p1,p2
            do rhs=1,nrhs
               call fft_node_to_node_translation(anode(:,:,:,:,rhs), &
                  cell_translation_matrix(:,:,:,:,:,p), &
                  gnode(:,:,:,:,rhs),p,mpi_comm=pcomm)
            enddo
         enddo

         call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
         if(numprocs.gt.1) then
            nsend=cell_dim(1)*cell_dim(2)*cell_dim(3)*node_order*(node_order+2)
            do rhs=1,nrhs
               if(inp1) then
                  call mstm_mpi(mpi_command='bcast', &
                     mpi_send_buf_dc=gnode(:,:,:,1:node_order*(node_order+2),rhs), &
                     mpi_number=nsend, &
                     mpi_rank=0, &
                     mpi_comm=synccomm1)
               endif
               if(inp2) then
                  call mstm_mpi(mpi_command='bcast', &
                     mpi_send_buf_dc=gnode(:,:,:, &
                        node_order*(node_order+2)+1:node_order*(node_order+2)*2,rhs), &
                     mpi_number=nsend, &
                     mpi_rank=0, &
                     mpi_comm=synccomm2)
               endif
            enddo
         endif
         call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)

!         if(numprocs/2.gt.1) then
!            nsend=cell_dim(1)*cell_dim(2)*cell_dim(3)*node_order*(node_order+2)*2*nrhs
!            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=gnode, &
!                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
!         endif
if(light_up) then
write(*,'('' fft5 '',i3)') mstm_global_rank
call flush(6)
endif

         call local_sphere_to_node_translation(nrhs,gout_t,gnode, &
                 store_matrix_option=smopt, &
                 mpi_comm=mpicomm,local_host=fft_local_host,sphere_to_node=.false., &
                 merge_procs=.false.)

         call local_sphere_to_sphere_expansion(nrhs,ain_t,gout_loc, &
                 store_matrix_option=smopt,initial_run=firstrun, &
                 mpi_comm=mpicomm,merge_procs=.false., &
                 local_host=fft_local_host)

         gout_t=gout_t+gout_loc

if(light_up) then
write(*,'('' fft6 '',i3)') mstm_global_rank
call flush(6)
endif

         do rhs=1,nrhs
            noff=0
            do i=1,number_spheres
               if(contran(rhs)) then
                  call shiftcoefficient(sphere_order(i),2,-1,-1, &
                        gout_t(noff+1:noff+2*sphere_order(i)*(sphere_order(i)+2),rhs), &
                        gout(noff+1:noff+2*sphere_order(i)*(sphere_order(i)+2),rhs))
               else
                  call shiftcoefficient(sphere_order(i),2,1,1, &
                        gout_t(noff+1:noff+2*sphere_order(i)*(sphere_order(i)+2),rhs), &
                        gout(noff+1:noff+2*sphere_order(i)*(sphere_order(i)+2),rhs))
               endif
               noff=noff+2*sphere_order(i)*(sphere_order(i)+2)*number_field_expansions(i)
            enddo
         enddo

if(light_up) then
write(*,'('' fft7 '',i3)') mstm_global_rank
call flush(6)
endif
         call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)

         deallocate(anode,gnode,gout_loc,ain_t,gout_t)
         firstrun=.false.

         end subroutine fft_external_to_external_expansion

         subroutine local_sphere_to_sphere_expansion(nrhs,ain,gout, &
                    store_matrix_option,initial_run,rhs_list, &
                    mpi_comm,con_tran,merge_procs,local_host)
         implicit none
         logical :: smopt,rhslist(nrhs),contran(nrhs),mergeprocs
         logical, save :: calcmat,firstrun
         logical, optional :: store_matrix_option,initial_run, &
                              rhs_list(nrhs),con_tran(nrhs),merge_procs
         integer :: rank,numprocs,nsphere,nrhs,mpicomm,rhs, &
                    i,j,npi1,npi2,npj1,npj2,noj,noi,task,proc, &
                    ndim,idim,localhost,npairs,n,nsend
         integer, optional :: mpi_comm,local_host
         complex(8)  :: ain(number_eqns,nrhs),gout(number_eqns,nrhs),rimedium(2)
         type(translation_data), pointer :: loc_tranmat
         type(translation_data), target :: tranmat
         type(linked_ilist), pointer :: llist
         data firstrun/.true./
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
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
         if(present(local_host)) then
            localhost=local_host
         else
            localhost=0
         endif
         if(present(merge_procs)) then
            mergeprocs=merge_procs
         else
            mergeprocs=.true.
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nsphere=number_spheres
         gout=0.
!
!  compute offsets for scattering coefficients
!
         if(firstrun) then
            if(smopt.and.store_translation_matrix) then
               task=0
               ndim=0
               do i=1,nsphere-1
                  if(host_sphere(i).ne.localhost) cycle
                  task=task+1
                  proc=mod(task,numprocs)
                  if(proc.eq.rank) then
                     ndim=ndim+sphere_local_interaction_list(i)%number_elements
                  endif
               enddo
               if(allocated(stored_local_h_mat)) deallocate(stored_local_h_mat)
               allocate(stored_local_h_mat(ndim))
            endif
            calcmat=.true.
            firstrun=.false.
         else
            calcmat=(.not.(smopt.and.store_translation_matrix))
         endif

         idim=0
         task=0
         do i=1,nsphere-1
            if(host_sphere(i).ne.localhost) cycle
            task=task+1
            proc=mod(task,numprocs)
            if(proc.eq.rank) then
               call exteriorrefindex(i,rimedium)
               npairs=sphere_local_interaction_list(i)%number_elements
               llist=>sphere_local_interaction_list(i)%members
               do n=1,npairs
                  idim=idim+1
                  j=llist%index
                  if(n.lt.npairs) llist=>llist%next
                  noi=sphere_order(i)
                  npi1=sphere_offset(i)+1
                  npi2=sphere_offset(i)+sphere_block(i)
                  noj=sphere_order(j)
                  npj1=sphere_offset(j)+1
                  npj2=sphere_offset(j)+sphere_block(j)
                  if(calcmat) then
                     if(smopt.and.store_translation_matrix) then
                        stored_local_h_mat(idim)%matrix_calculated=.false.
                        stored_local_h_mat(idim)%vswf_type=3
                        stored_local_h_mat(idim)%translation_vector=sphere_position(:,i)-sphere_position(:,j)
                        stored_local_h_mat(idim)%zero_translation=.false.
                        stored_local_h_mat(idim)%refractive_index=rimedium
                        stored_local_h_mat(idim)%rot_op=(min(noi,noj).ge.translation_switch_order)
                        loc_tranmat=>stored_local_h_mat(idim)
                     else
                        tranmat%matrix_calculated=.false.
                        tranmat%vswf_type=3
                        tranmat%translation_vector=sphere_position(:,i)-sphere_position(:,j)
                        tranmat%zero_translation=.false.
                        tranmat%refractive_index=rimedium
                        tranmat%rot_op=(min(noi,noj).ge.translation_switch_order)
                        loc_tranmat=>tranmat
                     endif
                  else
                     loc_tranmat=>stored_local_h_mat(idim)
                  endif
                  do rhs=1,nrhs
                     if(rhslist(rhs)) then
                        call coefficient_translation(noj,2,noi,2,ain(npj1:npj2,rhs),gout(npi1:npi2,rhs), &
                           loc_tranmat,shift_op=.false.,tran_op=contran(rhs))
                        call coefficient_translation(noi,2,noj,2,ain(npi1:npi2,rhs),gout(npj1:npj2,rhs), &
                           loc_tranmat,shift_op=.true.,tran_op=contran(rhs))
                     endif
                  enddo
                  if(calcmat.and.(.not.(smopt.and.store_translation_matrix))) then
                     if(tranmat%rot_op) then
                        deallocate(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                     else
                        deallocate(tranmat%gen_mat)
                     endif
                  endif
               enddo
            endif
         enddo

         if(numprocs.gt.1.and.mergeprocs) then
            nsend=number_eqns*nrhs
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=gout, &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         endif
         end subroutine local_sphere_to_sphere_expansion

         subroutine local_sphere_to_node_translation(nrhs,asphere,anode, &
                    store_matrix_option,initial_run,rhs_list, &
                    mpi_comm,con_tran,local_host,sphere_to_node,merge_procs)
         implicit none
         logical :: smopt,rhslist(nrhs),contran(nrhs),spheretonode,mergeprocs
         logical, save :: calcmat,firstrun
         logical, optional :: store_matrix_option,initial_run,merge_procs, &
                              rhs_list(nrhs),con_tran(nrhs),sphere_to_node
         integer :: neqns,rank,numprocs,nsphere,nrhs,mpicomm, &
                    i,npi1,npi2,noi,task,proc,nsend,rhs, &
                    ndim,idim,localhost,nodei(3)
         integer, optional :: mpi_comm,local_host
         real(8) :: rtran(3),r
         complex(8)  :: asphere(number_eqns,nrhs),anode(cell_dim(1),cell_dim(2),cell_dim(3),node_order*(node_order+2)*2,nrhs), &
            rimedium(2),anodet(node_order*(node_order+2)*2)
         type(translation_data), pointer :: loc_tranmat
         type(translation_data), target :: tranmat
         data firstrun/.true./
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
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
         if(present(local_host)) then
            localhost=local_host
         else
            localhost=0
         endif
         if(present(sphere_to_node)) then
            spheretonode=sphere_to_node
         else
            spheretonode=.true.
         endif
         if(present(merge_procs)) then
            mergeprocs=merge_procs
         else
            mergeprocs=.true.
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nsphere=number_spheres
         neqns=number_eqns
!
!  compute offsets for scattering coefficients
!
         if(firstrun) then
            if(smopt.and.store_translation_matrix) then
               task=0
               ndim=0
               do i=1,nsphere
                  if(host_sphere(i).eq.localhost) then
                     task=task+1
                     proc=mod(task,numprocs)
                     if(proc.eq.rank) ndim=ndim+1
                  endif
               enddo
               if(allocated(stored_local_j_mat)) deallocate(stored_local_j_mat)
               allocate(stored_local_j_mat(ndim))
            endif
            calcmat=.true.
            firstrun=.false.
         else
            calcmat=(.not.(smopt.and.store_translation_matrix))
         endif

         if(spheretonode) then
            anode=0.d0
         else
            asphere=0.d0
         endif

         idim=0
         task=0
         do i=1,nsphere
            if(host_sphere(i).eq.localhost) then
               task=task+1
               proc=mod(task,numprocs)
               if(proc.eq.rank) then
                  idim=idim+1
                  call exteriorrefindex(i,rimedium)
                  nodei(:)=sphere_node(:,i)
                  noi=sphere_order(i)
                  npi1=sphere_offset(i)+1
                  npi2=sphere_offset(i)+sphere_block(i)
                  rtran=d_cell*(dble(nodei)-.5d0)-sphere_position(:,i)+cell_boundary(:)
                  r=sum(rtran**2)
                  if(calcmat) then
                     if(smopt.and.store_translation_matrix) then
                        stored_local_j_mat(idim)%matrix_calculated=.false.
                        stored_local_j_mat(idim)%vswf_type=1
                        stored_local_j_mat(idim)%translation_vector=rtran
                        stored_local_j_mat(idim)%zero_translation=r.lt.1.d-6
                        stored_local_j_mat(idim)%refractive_index=rimedium
                        stored_local_j_mat(idim)%rot_op=(min(noi,node_order).ge.translation_switch_order)
                        loc_tranmat=>stored_local_j_mat(idim)
                     else
                        tranmat%matrix_calculated=.false.
                        tranmat%vswf_type=1
                        tranmat%translation_vector=rtran
                        tranmat%zero_translation=r.lt.1.d-6
                        tranmat%refractive_index=rimedium
                        tranmat%rot_op=(min(noi,node_order).ge.translation_switch_order)
                        loc_tranmat=>tranmat
                     endif
                  else
                     loc_tranmat=>stored_local_j_mat(idim)
                  endif
                  if(spheretonode) then
                     do rhs=1,nrhs
                        if(rhslist(rhs)) then
                           anodet=0.d0
                           call coefficient_translation(noi,2,node_order,2,asphere(npi1:npi2,rhs), &
                              anodet,loc_tranmat,shift_op=.false.,tran_op=contran(rhs))
                           anode(nodei(1),nodei(2),nodei(3),:,rhs) &
                              =anode(nodei(1),nodei(2),nodei(3),:,rhs)+anodet(:)
                        endif
                     enddo
                  else
                     do rhs=1,nrhs
                        if(rhslist(rhs)) then
                           anodet(:)=anode(nodei(1),nodei(2),nodei(3),:,rhs)
                           call coefficient_translation(node_order,2,noi,2,anodet, &
                              asphere(npi1:npi2,rhs),loc_tranmat, &
                              shift_op=.true.,tran_op=contran(rhs))
                        endif
                     enddo
                  endif

                  if(calcmat.and.(.not.(smopt.and.store_translation_matrix))) then
                     if(.not.tranmat%zero_translation) then
                        if(tranmat%rot_op) then
                           deallocate(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                        else
                           deallocate(tranmat%gen_mat)
                        endif
                     endif
                  endif
               endif
            endif
         enddo

         if(numprocs.gt.1.and.mergeprocs) then
            if(spheretonode) then
               nsend=cell_dim(1)*cell_dim(2)*cell_dim(3)*node_order*(node_order+2)*2*nrhs
               call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=anode, &
                    mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
            else
               nsend=number_eqns*nrhs
               call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=asphere, &
                    mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
            endif
         endif
         end subroutine local_sphere_to_node_translation

         subroutine node_selection(fva,target_min,target_max)
         implicit none
         integer :: nsphere,m,spherenode(3),n,i,node,ix,iy,iz,ir,j,icell,ncell,cell
         real(8) :: fva,r,fv,amean,tvol,dd,targetmin(3),targetmax(3)
         real(8), optional :: target_min(3),target_max(3)
         type(linked_ilist), pointer :: ilist,ilist2

         if(present(target_min)) then
            targetmin=target_min
         else
            targetmin(:)=sphere_min_position
         endif
         if(present(target_max)) then
            targetmax=target_max
         else
            targetmax(:)=sphere_max_position
         endif
         cell_boundary=targetmin
         if(fft_local_host.eq.0) then
            host_ref_index=medium_ref_index
         else
            host_ref_index=sphere_ref_index(:,fft_local_host)
         endif

         if(fva.le.0.d0) then
            amean=sum(sphere_radius(:))/dble(number_spheres)
            tvol=1.d0
            do i=1,3
               dd=targetmax(i)-targetmin(i)
               dd=max(dd,amean)
               tvol=tvol*dd
            enddo
            fv=4.d0*pi/3.d0*vol_radius**3/tvol
            fv=min(fv,1.d0)
            fv=max(fv,.02d0)
         else
            fv=fva
         endif
         fva=fv

!write(*,'('' host ri:'',2e13.5)') host_ref_index
!call flush(6)

         nsphere=number_spheres
         d_cell=(4.d0*pi/3.d0/fv/dble(number_spheres))**(1.d0/3.d0)*vol_radius
         cell_dim=ceiling((targetmax(:)-targetmin(:))/d_cell)
         do m=1,3
            cell_dim(m)=correctn235(cell_dim(m))
         enddo
         d_cell=maxval((targetmax(:)-targetmin(:))/dble(cell_dim(:)))

!write(*,'('' dcell:'',e13.5,3i5)') d_cell,cell_dim
!call flush(6)
         if(allocated(cell_list)) deallocate(cell_list)
         if(allocated(sphere_node)) deallocate(sphere_node)
         allocate(cell_list(cell_dim(1),cell_dim(2),cell_dim(3)),&
            sphere_node(3,number_spheres))
         cell_list(:,:,:)%number_elements=0
         do iz=1,cell_dim(3)
            do iy=1,cell_dim(2)
               do ix=1,cell_dim(1)
                  allocate(cell_list(ix,iy,iz)%members)
               enddo
            enddo
         enddo

!write(*,'('' d,n:'',e13.5,3i8)') d_cell,cell_dim
!call flush(6)

         cell_origin(:)=d_cell*dble(cell_dim(:))/2.d0
         do i=1,number_spheres
            do m=1,3
               r=sphere_position(m,i)-targetmin(m)
               node=floor(r/d_cell)+1
               spherenode(m)=node
               spherenode(m)=min(spherenode(m),cell_dim(m))
               spherenode(m)=max(spherenode(m),1)
            enddo
            sphere_node(:,i)=spherenode
            n=cell_list(spherenode(1),spherenode(2),spherenode(3))%number_elements
!write(*,'('' i:'',5i5)') i,spherenode,n
!call flush(6)
            ilist=>cell_list(spherenode(1),spherenode(2),spherenode(3))%members
            do m=1,n
               if(m.eq.n) allocate(ilist%next)
               ilist=>ilist%next
            enddo
            ilist%index=i
            ilist%next=>null()
            cell_list(spherenode(1),spherenode(2),spherenode(3))%number_elements=n+1
         enddo

         n=-1
         do iz=-1,1
            do iy=-1,1
               do ix=-1,1
                  ir=ix*ix+iy*iy+iz*iz
                  if(neighbor_node_model.eq.0) then
                     if(ir.gt.0) cycle
                  elseif(neighbor_node_model.eq.1) then
                     if(ir.gt.1) cycle
                  elseif(neighbor_node_model.eq.2) then
                     if(ir.gt.2) cycle
                  endif
                  n=n+1
                  neighbor_node(:,n)=(/ix,iy,iz/)
               enddo
            enddo
         enddo
         number_neighbor_nodes=n

!write(*,'('' nn:'',i8)') number_neighbor_nodes
!call flush(6)

         if(allocated(sphere_local_interaction_list)) deallocate(sphere_local_interaction_list)
         allocate(sphere_local_interaction_list(number_spheres))
         sphere_local_interaction_list(:)%number_elements=0

         do i=1,number_spheres
            allocate(sphere_local_interaction_list(i)%members)
            ilist2=>sphere_local_interaction_list(i)%members
            n=0
            do cell=0,number_neighbor_nodes
               spherenode=sphere_node(:,i)-neighbor_node(:,cell)
               if(any(spherenode.lt.(/1,1,1/)).or.any(spherenode.gt.cell_dim)) cycle
               ncell=cell_list(spherenode(1),spherenode(2),spherenode(3))%number_elements
               ilist=>cell_list(spherenode(1),spherenode(2),spherenode(3))%members
               do icell=1,ncell
                  j=ilist%index
                  if(j.gt.i) then
                     ilist2%index=j
                     n=n+1
                     allocate(ilist2%next)
                     ilist2=>ilist2%next
                  endif
                  if(icell.lt.ncell) ilist=>ilist%next
               enddo
            enddo
            sphere_local_interaction_list(i)%number_elements=n
         enddo

         end subroutine node_selection

         subroutine fft_node_to_node_translation(acoef,tranmat,gcoef,pmode,mpi_comm,tran_op)
         implicit none
         logical :: tranop
         logical, optional :: tran_op
         integer :: nblk,celldim2(3),n,l,mpicomm,numprocs,rank,ncells,task,proc,nsend,pmode,ncells8
         integer, optional :: mpi_comm
         complex(4) :: tranmat(8*cell_dim(1)*cell_dim(2)*cell_dim(3),node_order*(node_order+2),node_order*(node_order+2))
         complex(8) :: acoef(cell_dim(1)*cell_dim(2)*cell_dim(3),node_order*(node_order+2),2), &
            aft(8*cell_dim(1)*cell_dim(2)*cell_dim(3)), &
            gcoef(cell_dim(1)*cell_dim(2)*cell_dim(3),node_order*(node_order+2),2), &
            gft(8*cell_dim(1)*cell_dim(2)*cell_dim(3),node_order*(node_order+2))

         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(tran_op)) then
            tranop=tran_op
         else
            tranop=.false.
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         nblk=node_order*(node_order+2)
         celldim2=2*cell_dim
         ncells=cell_dim(1)*cell_dim(2)*cell_dim(3)
         ncells8=ncells*8
         gft=0.d0
         task=0
         do n=1,nblk
            task=task+1
            proc=mod(task,numprocs)
            if(proc.eq.rank) then
               aft=0.d0
               call fftmtx(acoef(1:ncells,n,pmode),aft(1:ncells8),1,cell_dim,celldim2,1)
               if(tranop) then
                  do l=1,nblk
                     gft(1:ncells8,l)=gft(1:ncells8,l)+tranmat(1:ncells8,n,l)*aft(1:ncells8)
                  enddo
               else
                  do l=1,nblk
                     gft(1:ncells8,l)=gft(1:ncells8,l)+tranmat(1:ncells8,l,n)*aft(1:ncells8)
                  enddo
               endif
            endif
         enddo

         if(numprocs.gt.1) then
            nsend=nblk*ncells8
            call mstm_mpi(mpi_command='allreduce', &
               mpi_recv_buf_dc=gft(1:ncells8,1:nblk), &
               mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         endif

         task=0
         do n=1,nblk
            task=task+1
            proc=mod(task,numprocs)
            if(proc.eq.rank) then
               gcoef(1:ncells,n,pmode)=0.d0
               call fftmtx(gcoef(1:ncells,n,pmode),gft(1:ncells8,n),1,cell_dim,celldim2,-1)
            endif
         enddo
         if(numprocs.gt.1) then
            nsend=nblk*ncells
            call mstm_mpi(mpi_command='allreduce', &
               mpi_recv_buf_dc=gcoef(1:ncells,1:nblk,pmode), &
               mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         endif
         gcoef(1:ncells,1:nblk,pmode)=gcoef(1:ncells,1:nblk,pmode)/dble(ncells8)
         end subroutine fft_node_to_node_translation

         subroutine fft_translation_initialization(ri,nodr,p1,p2)
         implicit none
         logical :: inhole
         integer :: nodr,nx,ny,nz,is,nblk,ncells2(3),isx,isy,isz,nx1,ny1,nz1, &
            n,i,node(3),p1,p2,p,l
         real(8) :: x,y,z,xp(3),r
         complex(8) :: hij(nodr*(nodr+2),nodr*(nodr+2),2),ri(2), &
            htemp(1:2*cell_dim(1),1:2*cell_dim(2),1:2*cell_dim(3))

         nblk=nodr*(nodr+2)
         ncells2=2*cell_dim

if(light_up) then
write(*,'('' fft ti '',i3,l)') mstm_global_rank,allocated(cell_translation_matrix)
call flush(6)
endif
if(allocated(cell_translation_matrix)) return
         if(allocated(cell_translation_matrix)) deallocate(cell_translation_matrix)
         allocate(cell_translation_matrix(0:2*cell_dim(1)-1, &
            0:2*cell_dim(2)-1,0:2*cell_dim(3)-1,1:nblk,1:nblk,p1:p2))
         cell_translation_matrix=0.d0

         do nz=0,cell_dim(3)
            z=d_cell*dble(nz)
            do ny=0,cell_dim(2)
               y=d_cell*dble(ny)
               do nx=0,cell_dim(1)
                  x=d_cell*dble(nx)
                  inhole=.false.
                  do i=0,number_neighbor_nodes
                     node(:)=(/nx,ny,nz/)-neighbor_node(:,i)
                     if(all(node.eq.(/0,0,0/))) then
                        inhole=.true.
                        exit
                     endif
                  enddo
                  if(inhole) cycle

                  r=sqrt(x*x+y*y+z*z)
                  do is=0,7
                     isx=1-2*mod(is,2)
                     isy=1-2*mod(int(is/2),2)
                     isz=1-2*mod(int(is/4),2)
                     nx1=isx*nx-(isx-1)*cell_dim(1)
                     ny1=isy*ny-(isy-1)*cell_dim(2)
                     nz1=isz*nz-(isz-1)*cell_dim(3)
                     if(nx1.lt.ncells2(1).and.ny1.lt.ncells2(2) &
                          .and.nz1.lt.ncells2(3)) then
                        xp(1)=isx*x
                        xp(2)=isy*y
                        xp(3)=isz*z
                        call gentranmatrix(nodr,nodr,translation_vector=xp, &
                           refractive_index=ri,ac_matrix=hij,vswf_type=3, &
                           mode_s=2,mode_t=2)
                        cell_translation_matrix(nx1,ny1,nz1,1:nblk,1:nblk,p1:p2) &
                           =cell_translation_matrix(nx1,ny1,nz1,1:nblk,1:nblk,p1:p2) &
                           +hij(1:nblk,1:nblk,p1:p2)
                     endif
                  enddo
               enddo
            enddo
         enddo
         n=2*nblk*nblk
         do p=p1,p2
            do n=1,nblk
               do l=1,nblk
                  htemp(:,:,:)=cell_translation_matrix(:,:,:,l,n,p)
                  call fftmtx(htemp,htemp,1,ncells2,ncells2,1)
                  cell_translation_matrix(:,:,:,l,n,p)=htemp(:,:,:)
               enddo
            enddo
         enddo

         end subroutine fft_translation_initialization

         subroutine fft1don3d(ain,aout,nblk,ntot1,ntot2,ntot3in,ntot3out,ndimin, &
                    ndimout,isign,looporder,trig)
         implicit none
         integer :: nblk,ntot1,ntot2,ntot3in,ntot3out,l,m,n,isign, &
                    looporder(3), &
                    triplet(3),ndimin(3),ndimout(3),ntot3,i1,i2,i3
         integer, parameter :: mxtrig=1000
         real(8) :: trig(mxtrig),ar_temp(nblk,max(ntot3in,ntot3out)), &
                    ai_temp(nblk,max(ntot3in,ntot3out))
         complex(8) :: aout(nblk,ndimout(1),ndimout(2),ndimout(3)), &
                    ain(nblk,ndimin(1),ndimin(2),ndimin(3))
         i1=looporder(1)
         i2=looporder(2)
         i3=looporder(3)
         ntot3=max(ntot3in,ntot3out)
         do l=1,ntot1
            triplet(i1)=l
            do m=1,ntot2
               triplet(i2)=m
               do n=1,ntot3in
                  triplet(i3)=n
                  ar_temp(1:nblk,n)=dble(ain(1:nblk,triplet(1),triplet(2),triplet(3)))
                  ai_temp(1:nblk,n)=dimag(ain(1:nblk,triplet(1),triplet(2),triplet(3)))
               enddo
               do n=ntot3in+1,ntot3out
                  ar_temp(1:nblk,n)=0.d0
                  ai_temp(1:nblk,n)=0.d0
               enddo
               call cgpfa(ar_temp(1:nblk,1:ntot3),ai_temp(1:nblk,1:ntot3), &
                    trig,nblk,ntot3,isign)
               do n=1,ntot3out
                  triplet(i3)=n
                  aout(1:nblk,triplet(1),triplet(2),triplet(3)) &
                    =dcmplx(ar_temp(1:nblk,n),ai_temp(1:nblk,n))
               enddo
            enddo
         enddo
         end subroutine fft1don3d

         subroutine fftmtx(am,amf,nblk,ntot,ntot2,isign)
         implicit none
         integer :: nblk,ntot(3),ntot2(3),isign, &
                    ntotxold,ntotyold,ntotzold, &
                    nblkold
         integer, parameter :: mxtrig=1000
         real(8) :: trig(mxtrig,3)
         save :: trig,ntotxold,ntotyold,ntotzold,nblkold
         complex(8) :: amf(nblk,ntot2(1),ntot2(2),ntot2(3)), &
                    am(nblk,ntot(1),ntot(2),ntot(3))
         data ntotxold,ntotyold,ntotzold,nblkold/0,0,0,0/
         if(ntot2(1).ne.ntotxold .or. ntot2(2).ne.ntotyold &
             .or. ntot2(3).ne.ntotzold) then
           ntotxold=ntot2(1)
           ntotyold=ntot2(2)
           ntotzold=ntot2(3)
           call setgpfa(trig(:,1),ntot2(1))
           call setgpfa(trig(:,2),ntot2(2))
           call setgpfa(trig(:,3),ntot2(3))
         endif
         if(isign.eq.1) then
            call fft1don3d(am,amf,nblk,ntot(1),ntot(2),ntot(3),ntot2(3),ntot, &
                    ntot2,1,(/1,2,3/),trig(:,3))
            call fft1don3d(amf,amf,nblk,ntot2(3),ntot(1),ntot(2),ntot2(2),ntot2, &
                    ntot2,1,(/3,1,2/),trig(:,2))
            call fft1don3d(amf,amf,nblk,ntot2(3),ntot2(2),ntot(1),ntot2(1),ntot2, &
                    ntot2,1,(/3,2,1/),trig(:,1))
         else
            call fft1don3d(amf,amf,nblk,ntot2(1),ntot2(2),ntot2(3),ntot(3),ntot2, &
                    ntot2,-1,(/1,2,3/),trig(:,3))
            call fft1don3d(amf,amf,nblk,ntot(3),ntot2(1),ntot2(2),ntot(2),ntot2, &
                    ntot2,-1,(/3,1,2/),trig(:,2))
            call fft1don3d(amf,am,nblk,ntot(3),ntot(2),ntot2(1),ntot(1),ntot2, &
                    ntot,-1,(/3,2,1/),trig(:,1))
         endif
         end subroutine fftmtx

         logical function checkn235(n)
         implicit none
         integer :: n,nn,ifac,ll,kk
         nn = n
         ifac = 2
         do ll = 1 , 3
            kk = 0
            do while (mod(nn,ifac).eq.0)
               kk = kk + 1
               nn = nn / ifac
            enddo
            ifac = ifac + ll
         enddo
         checkn235=nn.eq.1
         end function checkn235

         integer function correctn235(n)
         implicit none
         integer :: n,n1
         n1=n
         do while(.not.checkn235(n1))
            n1=n1+1
         enddo
         correctn235=n1
         end function correctn235


      subroutine cgpfa(cr,ci,trig,nblk,m,isign)
!      use iso_c_binding
      implicit none
      integer :: m,isign,nblk,n,i,inc
      real(8) :: trig(*),cr(nblk*m),ci(nblk*m)
!      real(8), pointer :: cr(:)
!      real(8) :: cr(2*nblk*m)
!      complex(8), target :: c(nblk*m)
!      call C_F_POINTER(C_LOC(c), cr, [2*nblk*m])
!      cr(1:2*nblk*m-1:2)=dble(c(1:nblk*m))
!      cr(2:2*nblk*m:2)=dimag(c(1:nblk*m))
!      inc=2*nblk
      inc=nblk
      do n=1,nblk
!         i=2*n-1
         i=n
         call gpfa(cr(i:),ci(i:),trig,inc,1,m,1,isign)
      enddo
!      c(1:nblk*m)=dcmplx(cr(1:2*nblk*m-1:2),cr(2:2*nblk*m:2))
      end subroutine cgpfa


!*********************************************************************
!                                                                    *
!     gpfapack - fortran implementation of the self-sorting          *
!     in-place generalized prime factor (complex) fft [gpfa]         *
!                                                                    *
!     written by clive temperton                                     *
!     recherche en prevision numerique / ecmwf                       *
!                                                                    *
!     the package consists of the setup routine setgpfa, together    *
!     with the routines gpfa, gpfa2f, gpfa3f, gpfa5f                 *
!                                                                    *
!*********************************************************************
!
!        subroutine 'setgpfa'
!        setup routine for self-sorting in-place
!            generalized prime factor (complex) fft [gpfa]
!
!        call setgpfa(trigs,n)
!
!        input :
!        -----
!        n is the length of the transforms. n must be of the form:
!          -----------------------------------
!            n = (2**ip) * (3**iq) * (5**ir)
!          -----------------------------------
!
!        output:
!        ------
!        trigs is a table of twiddle factors,
!          of length 2*ipqr (real) words, where:
!          --------------------------------------
!            ipqr = (2**ip) + (3**iq) + (5**ir)
!          --------------------------------------
!
!        written by clive temperton 1990
!
!----------------------------------------------------------------------
!
      subroutine setgpfa(trigs,n)
      implicit none
      integer :: n,nn,ifac,ll,kk,nj(3),ip,iq,ir,ni,irot,kink,k,i
      real(8) :: trigs(*),twopi,del,angle
!
!     decompose n into factors 2,3,5
!     ------------------------------
      nn = n
      ifac = 2
!
      do ll = 1 , 3
         kk = 0
         do while (mod(nn,ifac).eq.0)
            kk = kk + 1
            nn = nn / ifac
         enddo
         nj(ll) = kk
         ifac = ifac + ll
      enddo
!
      if (nn.ne.1) then
         write(6,40) n
   40    format(' *** warning!!!',i10,' is not a legal value of n ***')
         return
      endif
!
      ip = nj(1)
      iq = nj(2)
      ir = nj(3)
!
!     compute list of rotated twiddle factors
!     ---------------------------------------
      nj(1) = 2**ip
      nj(2) = 3**iq
      nj(3) = 5**ir
!
      twopi = 8.0d0 *datan(1.0d0)
      i = 1
!
      do ll = 1 , 3
         ni = nj(ll)
         if (ni.eq.1) cycle
!
         del = twopi / dble(ni)
         irot = n / ni
         kink = mod(irot,ni)
         kk = 0
!
         do k = 1 , ni
            angle = dble(kk) * del
            trigs(i) = cos(angle)
            trigs(i+1) = sin(angle)
            i = i + 2
            kk = kk + kink
            if (kk.gt.ni) kk = kk - ni
         enddo
      enddo
      end subroutine setgpfa
!        subroutine 'gpfa'
!        self-sorting in-place generalized prime factor (complex) fft
!
!        *** this is the all-fortran version ***
!            -------------------------------
!
!        call gpfa(a,b,trigs,inc,jump,n,lot,isign)
!
!        a is first real input/output vector
!        b is first imaginary input/output vector
!        trigs is a table of twiddle factors, precalculated
!              by calling subroutine 'setgpfa'
!        inc is the increment within each data vector
!        jump is the increment between data vectors
!        n is the length of the transforms:
!          -----------------------------------
!            n = (2**ip) * (3**iq) * (5**ir)
!          -----------------------------------
!        lot is the number of transforms
!        isign = +1 for forward transform
!              = -1 for inverse transform
!
!        written by clive temperton
!        recherche en prevision numerique
!        atmospheric environment service, canada
!
!----------------------------------------------------------------------
!
!        definition of transform
!        -----------------------
!
!        x(j) = sum(k=0,...,n-1)(c(k)*exp(isign*2*i*j*k*pi/n))
!
!---------------------------------------------------------------------
!
!        for a mathematical development of the algorithm used,
!        see:
!
!        c temperton : "a generalized prime factor fft algorithm
!          for any n = (2**p)(3**q)(5**r)",
!          siam j. sci. stat. comp., may 1992.
!
!----------------------------------------------------------------------
!
      subroutine gpfa(a,b,trigs,inc,jump,n,lot,isign)
      implicit none
      integer :: inc,jump,n,lot,isign,nn,ifac,ll,kk,nj(3),ip,iq,ir,i
      real(8) :: a(*),b(*),trigs(*)
!
!     decompose n into factors 2,3,5
!     ------------------------------
      nn = n
      ifac = 2
!
      do ll = 1 , 3
         kk = 0
         do while (mod(nn,ifac).eq.0)
            kk = kk + 1
            nn = nn / ifac
         enddo
         nj(ll) = kk
         ifac = ifac + ll
      enddo
!
      if (nn.ne.1) then
         write(6,40) n
   40    format(' *** warning!!!',i10,' is not a legal value of n ***')
         return
      endif
!
      ip = nj(1)
      iq = nj(2)
      ir = nj(3)
!
!     compute the transform
!     ---------------------
      i = 1
      if (ip.gt.0) then
         call gpfa2f(a,b,trigs,inc,jump,n,ip,lot,isign)
         i = i + 2 * ( 2**ip)
      endif
      if (iq.gt.0) then
         call gpfa3f(a,b,trigs(i),inc,jump,n,iq,lot,isign)
         i = i + 2 * (3**iq)
      endif
      if (ir.gt.0) then
         call gpfa5f(a,b,trigs(i),inc,jump,n,ir,lot,isign)
      endif
!
      end subroutine gpfa
!     fortran version of *gpfa2* -
!     radix-2 section of self-sorting, in-place, generalized pfa
!     central radix-2 and radix-8 passes included
!      so that transform length can be any power of 2
!
!-------------------------------------------------------------------
!
      subroutine gpfa2f(a,b,trigs,inc,jump,n,mm,lot,isign)
      implicit none
      integer :: inc,jump,n,mm,lot,isign,lvr,n2,inq,jstepx,ninc,ink, &
                 m2,m8,m,mh,nblox,left,nb,nvex,la,mu,ipass,jstep,jstepl, &
                 jjj,ja,nu,jb,jc,jd,j,l,kk,k,je,jf,jg,jh,laincl,ji,jj,jk,&
                 jl,jm,jn,jo,jp,istart,ll
      real(8) :: a(*),b(*),trigs(*),s,ss,t0,t2,t1,t3,u0,u2,u1,u3,co1,si1, &
                 co2,si2,co3,si3,c1,c2,c3,co4,si4,co5,si5,co6,si6,co7,si7
      data lvr/64/
!
!     ***************************************************************
!     *                                                             *
!     *  n.b. lvr = length of vector registers, set to 128 for c90. *
!     *  reset to 64 for other cray machines, or to any large value *
!     *  (greater than or equal to lot) for a scalar computer.      *
!     *                                                             *
!     ***************************************************************
!
      n2 = 2**mm
      inq = n/n2
      jstepx = (n2-n) * inc
      ninc = n * inc
      ink = inc * inq
!
      m2 = 0
      m8 = 0
      if (mod(mm,2).eq.0) then
         m = mm/2
      else if (mod(mm,4).eq.1) then
         m = (mm-1)/2
         m2 = 1
      else if (mod(mm,4).eq.3) then
         m = (mm-3)/2
         m8 = 1
      endif
      mh = (m+1)/2
!
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = dble(isign)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do 500 nb = 1 , nblox
!
         if (left.le.lvr) then
            nvex = left
         else if (left.lt.(2*lvr)) then
            nvex = left/2
            nvex = nvex + mod(nvex,2)
         else
            nvex = lvr
         endif
         left = left - nvex
!
         la = 1
!
!  loop on type i radix-4 passes
!  -----------------------------
         mu = mod(inq,4)
         if (isign.eq.-1) mu = 4 - mu
         ss = 1.0
         if (mu.eq.3) ss = -1.0
!
         if (mh.eq.0) go to 200
!
         do 160 ipass = 1 , mh
            jstep = (n*inc) / (4*la)
            jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
            do 120 jjj = 0 , (n-1)*inc , 4*jstep
               ja = istart + jjj
!
!     "transverse" loop
!     -----------------
               do 115 nu = 1 , inq
                  jb = ja + jstepl
                  if (jb.lt.istart) jb = jb + ninc
                  jc = jb + jstepl
                  if (jc.lt.istart) jc = jc + ninc
                  jd = jc + jstepl
                  if (jd.lt.istart) jd = jd + ninc
                  j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
                  do 110 l = 1 , nvex
                     t0 = a(ja+j) + a(jc+j)
                     t2 = a(ja+j) - a(jc+j)
                     t1 = a(jb+j) + a(jd+j)
                     t3 = ss * ( a(jb+j) - a(jd+j) )
                     u0 = b(ja+j) + b(jc+j)
                     u2 = b(ja+j) - b(jc+j)
                     u1 = b(jb+j) + b(jd+j)
                     u3 = ss * ( b(jb+j) - b(jd+j) )
                     a(ja+j) = t0 + t1
                     a(jc+j) = t0 - t1
                     b(ja+j) = u0 + u1
                     b(jc+j) = u0 - u1
                     a(jb+j) = t2 - u3
                     a(jd+j) = t2 + u3
                     b(jb+j) = u2 + t3
                     b(jd+j) = u2 - t3
                     j = j + jump
  110             continue
                  ja = ja + jstepx
                  if (ja.lt.istart) ja = ja + ninc
  115          continue
  120       continue
!
!  finished if n2 = 4
!  ------------------

            if (n2.eq.4) go to 490
            kk = 2 * la
!
!  loop on nonzero k
!  -----------------
            do 150 k = ink , jstep-ink , ink
               co1 = trigs(kk+1)
               si1 = s*trigs(kk+2)
               co2 = trigs(2*kk+1)
               si2 = s*trigs(2*kk+2)
               co3 = trigs(3*kk+1)
               si3 = s*trigs(3*kk+2)
!
!  loop along transform
!  --------------------
               do 140 jjj = k , (n-1)*inc , 4*jstep
                  ja = istart + jjj
!
!     "transverse" loop
!     -----------------
                  do 135 nu = 1 , inq
                     jb = ja + jstepl
                     if (jb.lt.istart) jb = jb + ninc
                     jc = jb + jstepl
                     if (jc.lt.istart) jc = jc + ninc
                     jd = jc + jstepl
                     if (jd.lt.istart) jd = jd + ninc
                     j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep,shortloop
                     do 130 l = 1 , nvex
                        t0 = a(ja+j) + a(jc+j)
                        t2 = a(ja+j) - a(jc+j)
                        t1 = a(jb+j) + a(jd+j)
                        t3 = ss * ( a(jb+j) - a(jd+j ) )
                        u0 = b(ja+j) + b(jc+j)
                        u2 = b(ja+j) - b(jc+j)
                        u1 = b(jb+j) + b(jd+j)
                        u3 = ss * ( b(jb+j) - b(jd+j) )
                        a(ja+j) = t0 + t1
                        b(ja+j) = u0 + u1
                        a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
                        b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
                        a(jc+j) = co2*(t0-t1) - si2*(u0-u1)
                        b(jc+j) = si2*(t0-t1) + co2*(u0-u1)
                        a(jd+j) = co3*(t2+u3) - si3*(u2-t3)
                        b(jd+j) = si3*(t2+u3) + co3*(u2-t3)
                        j = j + jump
  130                continue
!-----( end of loop across transforms )
                     ja = ja + jstepx
                     if (ja.lt.istart) ja = ja + ninc
  135             continue
  140          continue
!-----( end of loop along transforms )
               kk = kk + 2*la
  150       continue
!-----( end of loop on nonzero k )
            la = 4*la
  160    continue
!-----( end of loop on type i radix-4 passes)
!
!  central radix-2 pass
!  --------------------
  200 continue
      if (m2.eq.0) go to 300
!
      jstep = (n*inc) / (2*la)
      jstepl = jstep - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 220 jjj = 0 , (n-1)*inc , 2*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 215 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 210 l = 1 , nvex
      t0 = a(ja+j) - a(jb+j)
      a(ja+j) = a(ja+j) + a(jb+j)
      a(jb+j) = t0
      u0 = b(ja+j) - b(jb+j)
      b(ja+j) = b(ja+j) + b(jb+j)
      b(jb+j) = u0
      j = j + jump
  210 continue
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  215 continue
  220 continue
!
!  finished if n2=2
!  ----------------
      if (n2.eq.2) go to 490
!
      kk = 2 * la
!
!  loop on nonzero k
!  -----------------
      do 260 k = ink , jstep - ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
!
!  loop along transforms
!  ---------------------
      do 250 jjj = k , (n-1)*inc , 2*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 245 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      if (kk.eq.n2/2) then
!cdir$ ivdep, shortloop
      do 230 l = 1 , nvex
      t0 = ss * ( a(ja+j) - a(jb+j) )
      a(ja+j) = a(ja+j) + a(jb+j)
      a(jb+j) = ss * ( b(jb+j) - b(ja+j) )
      b(ja+j) = b(ja+j) + b(jb+j)
      b(jb+j) = t0
      j = j + jump
  230 continue
!
      else
!
!cdir$ ivdep, shortloop
      do 240 l = 1 , nvex
      t0 = a(ja+j) - a(jb+j)
      a(ja+j) = a(ja+j) + a(jb+j)
      u0 = b(ja+j) - b(jb+j)
      b(ja+j) = b(ja+j) + b(jb+j)
      a(jb+j) = co1*t0 - si1*u0
      b(jb+j) = si1*t0 + co1*u0
      j = j + jump
  240 continue
!
      endif
!
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  245 continue
  250 continue
!-----(end of loop along transforms)
      kk = kk + 2 * la
  260 continue
!-----(end of loop on nonzero k)
!-----(end of radix-2 pass)
!
      la = 2 * la
      go to 400
!
!  central radix-8 pass
!  --------------------

  300 continue
      if (m8.eq.0) go to 400
      jstep = (n*inc) / (8*la)
      jstepl = jstep - ninc
      mu = mod(inq,8)
      if (isign.eq.-1) mu = 8 - mu
      c1 = 1.0
      if (mu.eq.3.or.mu.eq.7) c1 = -1.0
      c2 = sqrt(0.5)
      if (mu.eq.3.or.mu.eq.5) c2 = -c2
      c3 = c1 * c2
!
!  stage 1
!  -------
      do 320 k = 0 , jstep - ink , ink
      do 315 jjj = k , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 312 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
!cdir$ ivdep, shortloop
      do 310 l = 1 , nvex
      t0 = a(ja+j) - a(je+j)
      a(ja+j) = a(ja+j) + a(je+j)
      t1 = c1 * ( a(jc+j) - a(jg+j) )
      a(je+j) = a(jc+j) + a(jg+j)
      t2 = a(jb+j) - a(jf+j)
      a(jc+j) = a(jb+j) + a(jf+j)
      t3 = a(jd+j) - a(jh+j)
      a(jg+j) = a(jd+j) + a(jh+j)
      a(jb+j) = t0
      a(jf+j) = t1
      a(jd+j) = c2 * ( t2 - t3 )
      a(jh+j) = c3 * ( t2 + t3 )
      u0 = b(ja+j) - b(je+j)
      b(ja+j) = b(ja+j) + b(je+j)
      u1 = c1 * ( b(jc+j) - b(jg+j) )
      b(je+j) = b(jc+j) + b(jg+j)
      u2 = b(jb+j) - b(jf+j)
      b(jc+j) = b(jb+j) + b(jf+j)
      u3 = b(jd+j) - b(jh+j)
      b(jg+j) = b(jd+j) + b(jh+j)
      b(jb+j) = u0
      b(jf+j) = u1
      b(jd+j) = c2 * ( u2 - u3 )
      b(jh+j) = c3 * ( u2 + u3 )
      j = j + jump
  310 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  312 continue
  315 continue
  320 continue
!
!  stage 2
!  -------
!
!  k=0 (no twiddle factors)
!  ------------------------
      do 330 jjj = 0 , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 328 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
!cdir$ ivdep, shortloop
      do 325 l = 1 , nvex
      t0 = a(ja+j) + a(je+j)
      t2 = a(ja+j) - a(je+j)
      t1 = a(jc+j) + a(jg+j)
      t3 = c1 * ( a(jc+j) - a(jg+j) )
      u0 = b(ja+j) + b(je+j)
      u2 = b(ja+j) - b(je+j)
      u1 = b(jc+j) + b(jg+j)
      u3 = c1 * ( b(jc+j) - b(jg+j ) )
      a(ja+j) = t0 + t1
      a(je+j) = t0 - t1
      b(ja+j) = u0 + u1
      b(je+j) = u0 - u1
      a(jc+j) = t2 - u3
      a(jg+j) = t2 + u3
      b(jc+j) = u2 + t3
      b(jg+j) = u2 - t3
      t0 = a(jb+j) + a(jd+j)
      t2 = a(jb+j) - a(jd+j)
      t1 = a(jf+j) - a(jh+j)
      t3 = a(jf+j) + a(jh+j)
      u0 = b(jb+j) + b(jd+j)
      u2 = b(jb+j) - b(jd+j)
      u1 = b(jf+j) - b(jh+j)
      u3 = b(jf+j) + b(jh+j)
      a(jb+j) = t0 - u3
      a(jh+j) = t0 + u3
      b(jb+j) = u0 + t3
      b(jh+j) = u0 - t3
      a(jd+j) = t2 + u1
      a(jf+j) = t2 - u1
      b(jd+j) = u2 - t1
      b(jf+j) = u2 + t1
      j = j + jump
  325 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  328 continue
  330 continue
!
      if (n2.eq.8) go to 490
!
!  loop on nonzero k
!  -----------------
      kk = 2 * la
!
      do 350 k = ink , jstep - ink , ink
!
      co1 = trigs(kk+1)
      si1 = s * trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s * trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s * trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s * trigs(4*kk+2)
      co5 = trigs(5*kk+1)
      si5 = s * trigs(5*kk+2)
      co6 = trigs(6*kk+1)
      si6 = s * trigs(6*kk+2)
      co7 = trigs(7*kk+1)
      si7 = s * trigs(7*kk+2)
!
      do 345 jjj = k , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 342 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
!cdir$ ivdep, shortloop
      do 340 l = 1 , nvex
      t0 = a(ja+j) + a(je+j)
      t2 = a(ja+j) - a(je+j)
      t1 = a(jc+j) + a(jg+j)
      t3 = c1 * ( a(jc+j) - a(jg+j) )
      u0 = b(ja+j) + b(je+j)
      u2 = b(ja+j) - b(je+j)
      u1 = b(jc+j) + b(jg+j)
      u3 = c1 * ( b(jc+j) - b(jg+j ) )
      a(ja+j) = t0 + t1
      b(ja+j) = u0 + u1
      a(je+j) = co4*(t0-t1) - si4*(u0-u1)
      b(je+j) = si4*(t0-t1) + co4*(u0-u1)
      a(jc+j) = co2*(t2-u3) - si2*(u2+t3)
      b(jc+j) = si2*(t2-u3) + co2*(u2+t3)
      a(jg+j) = co6*(t2+u3) - si6*(u2-t3)
      b(jg+j) = si6*(t2+u3) + co6*(u2-t3)
      t0 = a(jb+j) + a(jd+j)
      t2 = a(jb+j) - a(jd+j)
      t1 = a(jf+j) - a(jh+j)
      t3 = a(jf+j) + a(jh+j)
      u0 = b(jb+j) + b(jd+j)
      u2 = b(jb+j) - b(jd+j)
      u1 = b(jf+j) - b(jh+j)
      u3 = b(jf+j) + b(jh+j)
      a(jb+j) = co1*(t0-u3) - si1*(u0+t3)
      b(jb+j) = si1*(t0-u3) + co1*(u0+t3)
      a(jh+j) = co7*(t0+u3) - si7*(u0-t3)
      b(jh+j) = si7*(t0+u3) + co7*(u0-t3)
      a(jd+j) = co3*(t2+u1) - si3*(u2-t1)
      b(jd+j) = si3*(t2+u1) + co3*(u2-t1)
      a(jf+j) = co5*(t2-u1) - si5*(u2+t1)
      b(jf+j) = si5*(t2-u1) + co5*(u2+t1)
      j = j + jump
  340 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  342 continue
  345 continue
      kk = kk + 2 * la
  350 continue
!
      la = 8 * la
!
!  loop on type ii radix-4 passes
!  ------------------------------
  400 continue
      mu = mod(inq,4)
      if (isign.eq.-1) mu = 4 - mu
      ss = 1.0
      if (mu.eq.3) ss = -1.0
!
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (4*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 4*jstep
!
      do 420 jjj = ll , (n-1)*inc , 4*la*ink
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 415 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = ja + laincl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = je + laincl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jj + jstepl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = ji + laincl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jo + jstepl
      if (jp.lt.istart) jp = jp + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 410 l = 1 , nvex
      t0 = a(ja+j) + a(jc+j)
      t2 = a(ja+j) - a(jc+j)
      t1 = a(jb+j) + a(jd+j)
      t3 = ss * ( a(jb+j) - a(jd+j) )
      a(jc+j) = a(ji+j)
      u0 = b(ja+j) + b(jc+j)
      u2 = b(ja+j) - b(jc+j)
      u1 = b(jb+j) + b(jd+j)
      u3 = ss * ( b(jb+j) - b(jd+j) )
      a(jb+j) = a(je+j)
      a(ja+j) = t0 + t1
      a(ji+j) = t0 - t1
      b(ja+j) = u0 + u1
      b(jc+j) = u0 - u1
      b(jd+j) = b(jm+j)
      a(je+j) = t2 - u3
      a(jd+j) = t2 + u3
      b(jb+j) = u2 + t3
      b(jm+j) = u2 - t3
!----------------------
      t0 = a(jb+j) + a(jg+j)
      t2 = a(jb+j) - a(jg+j)
      t1 = a(jf+j) + a(jh+j)
      t3 = ss * ( a(jf+j) - a(jh+j) )
      a(jg+j) = a(jj+j)
      u0 = b(je+j) + b(jg+j)
      u2 = b(je+j) - b(jg+j)
      u1 = b(jf+j) + b(jh+j)
      u3 = ss * ( b(jf+j) - b(jh+j) )
      b(je+j) = b(jb+j)
      a(jb+j) = t0 + t1
      a(jj+j) = t0 - t1
      b(jg+j) = b(jj+j)
      b(jb+j) = u0 + u1
      b(jj+j) = u0 - u1
      a(jf+j) = t2 - u3
      a(jh+j) = t2 + u3
      b(jf+j) = u2 + t3
      b(jh+j) = u2 - t3
!----------------------
      t0 = a(jc+j) + a(jk+j)
      t2 = a(jc+j) - a(jk+j)
      t1 = a(jg+j) + a(jl+j)
      t3 = ss * ( a(jg+j) - a(jl+j) )
      u0 = b(ji+j) + b(jk+j)
      u2 = b(ji+j) - b(jk+j)
      a(jl+j) = a(jo+j)
      u1 = b(jg+j) + b(jl+j)
      u3 = ss * ( b(jg+j) - b(jl+j) )
      b(ji+j) = b(jc+j)
      a(jc+j) = t0 + t1
      a(jk+j) = t0 - t1
      b(jl+j) = b(jo+j)
      b(jc+j) = u0 + u1
      b(jk+j) = u0 - u1
      a(jg+j) = t2 - u3
      a(jo+j) = t2 + u3
      b(jg+j) = u2 + t3
      b(jo+j) = u2 - t3
!----------------------
      t0 = a(jm+j) + a(jl+j)
      t2 = a(jm+j) - a(jl+j)
      t1 = a(jn+j) + a(jp+j)
      t3 = ss * ( a(jn+j) - a(jp+j) )
      a(jm+j) = a(jd+j)
      u0 = b(jd+j) + b(jl+j)
      u2 = b(jd+j) - b(jl+j)
      u1 = b(jn+j) + b(jp+j)
      u3 = ss * ( b(jn+j) - b(jp+j) )
      a(jn+j) = a(jh+j)
      a(jd+j) = t0 + t1
      a(jl+j) = t0 - t1
      b(jd+j) = u0 + u1
      b(jl+j) = u0 - u1
      b(jn+j) = b(jh+j)
      a(jh+j) = t2 - u3
      a(jp+j) = t2 + u3
      b(jh+j) = u2 + t3
      b(jp+j) = u2 - t3
      j = j + jump
  410 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  415 continue
  420 continue
  430 continue
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
      if (ipass.eq.m) go to 490
!
      kk = 2*la
!
!     loop on nonzero k
!     -----------------
      do 470 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 4*jstep
!
      do 450 jjj = ll , (n-1)*inc , 4*la*ink
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = ja + laincl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = je + laincl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jj + jstepl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = ji + laincl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jo + jstepl
      if (jp.lt.istart) jp = jp + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 440 l = 1 , nvex
      t0 = a(ja+j) + a(jc+j)
      t2 = a(ja+j) - a(jc+j)
      t1 = a(jb+j) + a(jd+j)
      t3 = ss * ( a(jb+j) - a(jd+j) )
      a(jc+j) = a(ji+j)
      u0 = b(ja+j) + b(jc+j)
      u2 = b(ja+j) - b(jc+j)
      u1 = b(jb+j) + b(jd+j)
      u3 = ss * ( b(jb+j) - b(jd+j) )
      a(jb+j) = a(je+j)
      a(ja+j) = t0 + t1
      b(ja+j) = u0 + u1
      a(je+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
      b(jd+j) = b(jm+j)
      a(ji+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jc+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jd+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jm+j) = si3*(t2+u3) + co3*(u2-t3)
!----------------------------------------
      t0 = a(jb+j) + a(jg+j)
      t2 = a(jb+j) - a(jg+j)
      t1 = a(jf+j) + a(jh+j)
      t3 = ss * ( a(jf+j) - a(jh+j) )
      a(jg+j) = a(jj+j)
      u0 = b(je+j) + b(jg+j)
      u2 = b(je+j) - b(jg+j)
      u1 = b(jf+j) + b(jh+j)
      u3 = ss * ( b(jf+j) - b(jh+j) )
      b(je+j) = b(jb+j)
      a(jb+j) = t0 + t1
      b(jb+j) = u0 + u1
      b(jg+j) = b(jj+j)
      a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jj+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jj+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jh+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jh+j) = si3*(t2+u3) + co3*(u2-t3)
!----------------------------------------
      t0 = a(jc+j) + a(jk+j)
      t2 = a(jc+j) - a(jk+j)
      t1 = a(jg+j) + a(jl+j)
      t3 = ss * ( a(jg+j) - a(jl+j) )
      u0 = b(ji+j) + b(jk+j)
      u2 = b(ji+j) - b(jk+j)
      a(jl+j) = a(jo+j)
      u1 = b(jg+j) + b(jl+j)
      u3 = ss * ( b(jg+j) - b(jl+j) )
      b(ji+j) = b(jc+j)
      a(jc+j) = t0 + t1
      b(jc+j) = u0 + u1
      b(jl+j) = b(jo+j)
      a(jg+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jg+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jk+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jk+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jo+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jo+j) = si3*(t2+u3) + co3*(u2-t3)
!----------------------------------------
      t0 = a(jm+j) + a(jl+j)
      t2 = a(jm+j) - a(jl+j)
      t1 = a(jn+j) + a(jp+j)
      t3 = ss * ( a(jn+j) - a(jp+j) )
      a(jm+j) = a(jd+j)
      u0 = b(jd+j) + b(jl+j)
      u2 = b(jd+j) - b(jl+j)
      a(jn+j) = a(jh+j)
      u1 = b(jn+j) + b(jp+j)
      u3 = ss * ( b(jn+j) - b(jp+j) )
      b(jn+j) = b(jh+j)
      a(jd+j) = t0 + t1
      b(jd+j) = u0 + u1
      a(jh+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jh+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jl+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jl+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jp+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jp+j) = si3*(t2+u3) + co3*(u2-t3)
      j = j + jump
  440 continue
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 4*la
  480 continue
!-----( end of loop on type ii radix-4 passes )
!-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
      return
      end subroutine gpfa2f
!     fortran version of *gpfa3* -
!     radix-3 section of self-sorting, in-place
!        generalized pfa
!
!-------------------------------------------------------------------
!
      subroutine gpfa3f(a,b,trigs,inc,jump,n,mm,lot,isign)
      implicit none
      integer :: inc,jump,n,mm,lot,isign,lvr,inq,jstepx,ninc,ink, &
                 m,mh,nblox,left,nb,nvex,la,mu,ipass,jstep,jstepl, &
                 jjj,ja,nu,jb,jc,jd,j,l,kk,k,je,jf,jg,jh,laincl,ji, &
                 n3,istart,ll
      real(8) :: a(*),b(*),trigs(*),s,t2,t1,t3,u2,u1,u3,co1,si1, &
                 co2,si2,c1, &
                 sin60
      data sin60/0.866025403784437/
      data lvr/64/
!
!     ***************************************************************
!     *                                                             *
!     *  n.b. lvr = length of vector registers, set to 128 for c90. *
!     *  reset to 64 for other cray machines, or to any large value *
!     *  (greater than or equal to lot) for a scalar computer.      *
!     *                                                             *
!     ***************************************************************
!
      n3 = 3**mm
      inq = n/n3
      jstepx = (n3-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,3)
      if (isign.eq.-1) mu = 3-mu
      m = mm
      mh = (m+1)/2
      s = float(isign)
      c1 = sin60
      if (mu.eq.2) c1 = -c1
!
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = float(isign)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do 500 nb = 1 , nblox
!
      if (left.le.lvr) then
         nvex = left
      else if (left.lt.(2*lvr)) then
         nvex = left/2
         nvex = nvex + mod(nvex,2)
      else
         nvex = lvr
      endif
      left = left - nvex
!
      la = 1
!
!  loop on type i radix-3 passes
!  -----------------------------
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
      do 120 jjj = 0 , (n-1)*inc , 3*jstep
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 115 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 110 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jb+j) = t2 - u3
      b(jb+j) = u2 + t3
      a(jc+j) = t2 + u3
      b(jc+j) = u2 - t3
      j = j + jump
  110 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  115 continue
  120 continue
!
!  finished if n3 = 3
!  ------------------
      if (n3.eq.3) go to 490
      kk = 2 * la
!
!  loop on nonzero k
!  -----------------
      do 150 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
!
!  loop along transform
!  --------------------
      do 140 jjj = k , (n-1)*inc , 3*jstep
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep,shortloop
      do 130 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jc+j) = co2*(t2+u3) - si2*(u2-t3)
      b(jc+j) = si2*(t2+u3) + co2*(u2-t3)
      j = j + jump
  130 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
!-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
!-----( end of loop on nonzero k )
      la = 3*la
  160 continue
!-----( end of loop on type i radix-3 passes)
!
!  loop on type ii radix-3 passes
!  ------------------------------
!
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
      laincl = la*ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 3*jstep
!
      do 420 jjj = ll , (n-1)*inc , 3*la*ink
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 415 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = ja + laincl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jd + laincl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 410 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      a(jb+j) = a(jd+j)
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      b(jb+j) = b(jd+j)
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jd+j) = t2 - u3
      b(jd+j) = u2 + t3
      a(jc+j) = t2 + u3
      b(jc+j) = u2 - t3
!----------------------
      t1 = a(je+j) + a(jf+j)
      t2 = a(jb+j) - 0.5 * t1
      t3 = c1 * ( a(je+j) - a(jf+j) )
      a(jf+j) = a(jh+j)
      u1 = b(je+j) + b(jf+j)
      u2 = b(jb+j) - 0.5 * u1
      u3 = c1 * ( b(je+j) - b(jf+j) )
      b(jf+j) = b(jh+j)
      a(jb+j) = a(jb+j) + t1
      b(jb+j) = b(jb+j) + u1
      a(je+j) = t2 - u3
      b(je+j) = u2 + t3
      a(jh+j) = t2 + u3
      b(jh+j) = u2 - t3
!----------------------
      t1 = a(jf+j) + a(ji+j)
      t2 = a(jg+j) - 0.5 * t1
      t3 = c1 * ( a(jf+j) - a(ji+j) )
      t1 = a(jg+j) + t1
      a(jg+j) = a(jc+j)
      u1 = b(jf+j) + b(ji+j)
      u2 = b(jg+j) - 0.5 * u1
      u3 = c1 * ( b(jf+j) - b(ji+j) )
      u1 = b(jg+j) + u1
      b(jg+j) = b(jc+j)
      a(jc+j) = t1
      b(jc+j) = u1
      a(jf+j) = t2 - u3
      b(jf+j) = u2 + t3
      a(ji+j) = t2 + u3
      b(ji+j) = u2 - t3
      j = j + jump
  410 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  415 continue
  420 continue
  430 continue
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
      if (ipass.eq.m) go to 490
!
      kk = 2*la
!
!     loop on nonzero k
!     -----------------
      do 470 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 3*jstep
!
      do 450 jjj = ll , (n-1)*inc , 3*la*ink
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = ja + laincl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jd + laincl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
!cdir$ ivdep, shortloop
      do 440 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      a(jb+j) = a(jd+j)
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      b(jb+j) = b(jd+j)
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jd+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jd+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jc+j) = co2*(t2+u3) - si2*(u2-t3)
      b(jc+j) = si2*(t2+u3) + co2*(u2-t3)
!----------------------
      t1 = a(je+j) + a(jf+j)
      t2 = a(jb+j) - 0.5 * t1
      t3 = c1 * ( a(je+j) - a(jf+j) )
      a(jf+j) = a(jh+j)
      u1 = b(je+j) + b(jf+j)
      u2 = b(jb+j) - 0.5 * u1
      u3 = c1 * ( b(je+j) - b(jf+j) )
      b(jf+j) = b(jh+j)
      a(jb+j) = a(jb+j) + t1
      b(jb+j) = b(jb+j) + u1
      a(je+j) = co1*(t2-u3) - si1*(u2+t3)
      b(je+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jh+j) = co2*(t2+u3) - si2*(u2-t3)
      b(jh+j) = si2*(t2+u3) + co2*(u2-t3)
!----------------------
      t1 = a(jf+j) + a(ji+j)
      t2 = a(jg+j) - 0.5 * t1
      t3 = c1 * ( a(jf+j) - a(ji+j) )
      t1 = a(jg+j) + t1
      a(jg+j) = a(jc+j)
      u1 = b(jf+j) + b(ji+j)
      u2 = b(jg+j) - 0.5 * u1
      u3 = c1 * ( b(jf+j) - b(ji+j) )
      u1 = b(jg+j) + u1
      b(jg+j) = b(jc+j)
      a(jc+j) = t1
      b(jc+j) = u1
      a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
      a(ji+j) = co2*(t2+u3) - si2*(u2-t3)
      b(ji+j) = si2*(t2+u3) + co2*(u2-t3)
      j = j + jump
  440 continue
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 3*la
  480 continue
!-----( end of loop on type ii radix-3 passes )
!-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
      return
      end subroutine gpfa3f
!     fortran version of *gpfa5* -
!     radix-5 section of self-sorting, in-place,
!        generalized pfa
!
!-------------------------------------------------------------------
!
      subroutine gpfa5f(a,b,trigs,inc,jump,n,mm,lot,isign)
      implicit none
      integer :: inc,jump,n,mm,lot,isign,lvr,inq,jstepx,ninc,ink, &
                 m,mh,nblox,left,nb,nvex,la,mu,ipass,jstep,jstepl, &
                 jjj,ja,nu,jb,jc,jd,j,l,kk,k,je,jf,jg,jh,laincl,ji,jj,jk,&
                 jl,jm,jn,jo,jp,n5,jq,jr,js,jt,ju,jv,jw,jx,jy,istart,ll
      real(8) :: a(*),b(*),trigs(*),s,t2,t1,t3,u2,u1,u3,co1,si1, &
                 co2,si2,co3,si3,c1,c2,c3,co4,si4, &
                 sin36,sin72,qrt5,t4,t5,t6,t7,t8,t9,t10,t11, &
                 ax,bx,u4,u5,u6,u7,u8,u9,u10,u11
      data sin36/0.587785252292473/, sin72/0.951056516295154/, &
           qrt5/0.559016994374947/
      data lvr/64/
!
!     ***************************************************************
!     *                                                             *
!     *  n.b. lvr = length of vector registers, set to 128 for c90. *
!     *  reset to 64 for other cray machines, or to any large value *
!     *  (greater than or equal to lot) for a scalar computer.      *
!     *                                                             *
!     ***************************************************************
!
      n5 = 5 ** mm
      inq = n / n5
      jstepx = (n5-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,5)
      if (isign.eq.-1) mu = 5 - mu
!
      m = mm
      mh = (m+1)/2
      s = float(isign)
      c1 = qrt5
      c2 = sin72
      c3 = sin36
      if (mu.eq.2.or.mu.eq.3) then
         c1 = -c1
         c2 = sin36
         c3 = sin72
      endif
      if (mu.eq.3.or.mu.eq.4) c2 = -c2
      if (mu.eq.2.or.mu.eq.4) c3 = -c3
!
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = float(isign)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do 500 nb = 1 , nblox
!
      if (left.le.lvr) then
         nvex = left
      else if (left.lt.(2*lvr)) then
         nvex = left/2
         nvex = nvex + mod(nvex,2)
      else
         nvex = lvr
      endif
      left = left - nvex
!
      la = 1
!
!  loop on type i radix-5 passes
!  -----------------------------
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      kk = 0
!
!  loop on k
!  ---------
      do 150 k = 0 , jstep-ink , ink
!
      if (k.gt.0) then
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s*trigs(4*kk+2)
      endif
!
!  loop along transform
!  --------------------
      do 140 jjj = k , (n-1)*inc , 5*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      if (k.eq.0) then
!
!cdir$ ivdep, shortloop
      do 110 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jb+j) = t8 - u11
      b(jb+j) = u8 + t11
      a(je+j) = t8 + u11
      b(je+j) = u8 - t11
      a(jc+j) = t9 - u10
      b(jc+j) = u9 + t10
      a(jd+j) = t9 + u10
      b(jd+j) = u9 - t10
      j = j + jump
  110 continue
!
      else
!
!cdir$ ivdep,shortloop
      do 130 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jb+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jb+j) = si1*(t8-u11) + co1*(u8+t11)
      a(je+j) = co4*(t8+u11) - si4*(u8-t11)
      b(je+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jc+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jc+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jd+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jd+j) = si3*(t9+u10) + co3*(u9-t10)
      j = j + jump
  130 continue
!
      endif
!
!-----( end of loop across transforms )
!
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
!-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
!-----( end of loop on nonzero k )
      la = 5*la
  160 continue
!-----( end of loop on type i radix-5 passes)
!
      if (n.eq.5) go to 490
!
!  loop on type ii radix-5 passes
!  ------------------------------
!
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc
      kk = 0
!
!     loop on k
!     ---------
      do 470 k = 0 , jstep-ink , ink
!
      if (k.gt.0) then
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s*trigs(4*kk+2)
      endif
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 5*jstep
!
      do 450 jjj = ll , (n-1)*inc , 5*la*ink
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = ja + laincl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jf + laincl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = jl + jstepl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jk + laincl
      if (jp.lt.istart) jp = jp + ninc
      jq = jp + jstepl
      if (jq.lt.istart) jq = jq + ninc
      jr = jq + jstepl
      if (jr.lt.istart) jr = jr + ninc
      js = jr + jstepl
      if (js.lt.istart) js = js + ninc
      jt = js + jstepl
      if (jt.lt.istart) jt = jt + ninc
      ju = jp + laincl
      if (ju.lt.istart) ju = ju + ninc
      jv = ju + jstepl
      if (jv.lt.istart) jv = jv + ninc
      jw = jv + jstepl
      if (jw.lt.istart) jw = jw + ninc
      jx = jw + jstepl
      if (jx.lt.istart) jx = jx + ninc
      jy = jx + jstepl
      if (jy.lt.istart) jy = jy + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      if (k.eq.0) then
!
!cdir$ ivdep, shortloop
      do 410 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      a(jb+j) = a(jf+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jc+j) = a(jk+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      b(jb+j) = b(jf+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jc+j) = b(jk+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jf+j) = t8 - u11
      b(jf+j) = u8 + t11
      a(je+j) = t8 + u11
      b(je+j) = u8 - t11
      a(jk+j) = t9 - u10
      b(jk+j) = u9 + t10
      a(jd+j) = t9 + u10
      b(jd+j) = u9 - t10
!----------------------
      t1 = a(jg+j) + a(jj+j)
      t2 = a(jh+j) + a(ji+j)
      t3 = a(jg+j) - a(jj+j)
      t4 = a(jh+j) - a(ji+j)
      a(jh+j) = a(jl+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jb+j) - 0.25 * t5
      a(jb+j) = a(jb+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ji+j) = a(jq+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jg+j) + b(jj+j)
      u2 = b(jh+j) + b(ji+j)
      u3 = b(jg+j) - b(jj+j)
      u4 = b(jh+j) - b(ji+j)
      b(jh+j) = b(jl+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jb+j) - 0.25 * u5
      b(jb+j) = b(jb+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ji+j) = b(jq+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jg+j) = t8 - u11
      b(jg+j) = u8 + t11
      a(jj+j) = t8 + u11
      b(jj+j) = u8 - t11
      a(jl+j) = t9 - u10
      b(jl+j) = u9 + t10
      a(jq+j) = t9 + u10
      b(jq+j) = u9 - t10
!----------------------
      t1 = a(jh+j) + a(jo+j)
      t2 = a(jm+j) + a(jn+j)
      t3 = a(jh+j) - a(jo+j)
      t4 = a(jm+j) - a(jn+j)
      a(jn+j) = a(jr+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jc+j) - 0.25 * t5
      a(jc+j) = a(jc+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jo+j) = a(jw+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jh+j) + b(jo+j)
      u2 = b(jm+j) + b(jn+j)
      u3 = b(jh+j) - b(jo+j)
      u4 = b(jm+j) - b(jn+j)
      b(jn+j) = b(jr+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jc+j) - 0.25 * u5
      b(jc+j) = b(jc+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jo+j) = b(jw+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jh+j) = t8 - u11
      b(jh+j) = u8 + t11
      a(jw+j) = t8 + u11
      b(jw+j) = u8 - t11
      a(jm+j) = t9 - u10
      b(jm+j) = u9 + t10
      a(jr+j) = t9 + u10
      b(jr+j) = u9 - t10
!----------------------
      t1 = a(ji+j) + a(jt+j)
      t2 = a(jn+j) + a(js+j)
      t3 = a(ji+j) - a(jt+j)
      t4 = a(jn+j) - a(js+j)
      a(jt+j) = a(jx+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jp+j) - 0.25 * t5
      ax = a(jp+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jp+j) = a(jd+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(jd+j) = ax
      u1 = b(ji+j) + b(jt+j)
      u2 = b(jn+j) + b(js+j)
      u3 = b(ji+j) - b(jt+j)
      u4 = b(jn+j) - b(js+j)
      b(jt+j) = b(jx+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jp+j) - 0.25 * u5
      bx = b(jp+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jp+j) = b(jd+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(jd+j) = bx
      a(ji+j) = t8 - u11
      b(ji+j) = u8 + t11
      a(jx+j) = t8 + u11
      b(jx+j) = u8 - t11
      a(jn+j) = t9 - u10
      b(jn+j) = u9 + t10
      a(js+j) = t9 + u10
      b(js+j) = u9 - t10
!----------------------
      t1 = a(jv+j) + a(jy+j)
      t2 = a(jo+j) + a(jt+j)
      t3 = a(jv+j) - a(jy+j)
      t4 = a(jo+j) - a(jt+j)
      a(jv+j) = a(jj+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ju+j) - 0.25 * t5
      ax = a(ju+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ju+j) = a(je+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(je+j) = ax
      u1 = b(jv+j) + b(jy+j)
      u2 = b(jo+j) + b(jt+j)
      u3 = b(jv+j) - b(jy+j)
      u4 = b(jo+j) - b(jt+j)
      b(jv+j) = b(jj+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ju+j) - 0.25 * u5
      bx = b(ju+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ju+j) = b(je+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(je+j) = bx
      a(jj+j) = t8 - u11
      b(jj+j) = u8 + t11
      a(jy+j) = t8 + u11
      b(jy+j) = u8 - t11
      a(jo+j) = t9 - u10
      b(jo+j) = u9 + t10
      a(jt+j) = t9 + u10
      b(jt+j) = u9 - t10
      j = j + jump
  410 continue
!
      else
!
!cdir$ ivdep, shortloop
      do 440 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      a(jb+j) = a(jf+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jc+j) = a(jk+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      b(jb+j) = b(jf+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jc+j) = b(jk+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jf+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jf+j) = si1*(t8-u11) + co1*(u8+t11)
      a(je+j) = co4*(t8+u11) - si4*(u8-t11)
      b(je+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jk+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jk+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jd+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jd+j) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      t1 = a(jg+j) + a(jj+j)
      t2 = a(jh+j) + a(ji+j)
      t3 = a(jg+j) - a(jj+j)
      t4 = a(jh+j) - a(ji+j)
      a(jh+j) = a(jl+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jb+j) - 0.25 * t5
      a(jb+j) = a(jb+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ji+j) = a(jq+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jg+j) + b(jj+j)
      u2 = b(jh+j) + b(ji+j)
      u3 = b(jg+j) - b(jj+j)
      u4 = b(jh+j) - b(ji+j)
      b(jh+j) = b(jl+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jb+j) - 0.25 * u5
      b(jb+j) = b(jb+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ji+j) = b(jq+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jg+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jg+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jj+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jj+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jl+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jl+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jq+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jq+j) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      t1 = a(jh+j) + a(jo+j)
      t2 = a(jm+j) + a(jn+j)
      t3 = a(jh+j) - a(jo+j)
      t4 = a(jm+j) - a(jn+j)
      a(jn+j) = a(jr+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jc+j) - 0.25 * t5
      a(jc+j) = a(jc+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jo+j) = a(jw+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jh+j) + b(jo+j)
      u2 = b(jm+j) + b(jn+j)
      u3 = b(jh+j) - b(jo+j)
      u4 = b(jm+j) - b(jn+j)
      b(jn+j) = b(jr+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jc+j) - 0.25 * u5
      b(jc+j) = b(jc+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jo+j) = b(jw+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jh+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jh+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jw+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jw+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jm+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jm+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jr+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jr+j) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      t1 = a(ji+j) + a(jt+j)
      t2 = a(jn+j) + a(js+j)
      t3 = a(ji+j) - a(jt+j)
      t4 = a(jn+j) - a(js+j)
      a(jt+j) = a(jx+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jp+j) - 0.25 * t5
      ax = a(jp+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jp+j) = a(jd+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(jd+j) = ax
      u1 = b(ji+j) + b(jt+j)
      u2 = b(jn+j) + b(js+j)
      u3 = b(ji+j) - b(jt+j)
      u4 = b(jn+j) - b(js+j)
      b(jt+j) = b(jx+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jp+j) - 0.25 * u5
      bx = b(jp+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jp+j) = b(jd+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(jd+j) = bx
      a(ji+j) = co1*(t8-u11) - si1*(u8+t11)
      b(ji+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jx+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jx+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jn+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jn+j) = si2*(t9-u10) + co2*(u9+t10)
      a(js+j) = co3*(t9+u10) - si3*(u9-t10)
      b(js+j) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      t1 = a(jv+j) + a(jy+j)
      t2 = a(jo+j) + a(jt+j)
      t3 = a(jv+j) - a(jy+j)
      t4 = a(jo+j) - a(jt+j)
      a(jv+j) = a(jj+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ju+j) - 0.25 * t5
      ax = a(ju+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ju+j) = a(je+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(je+j) = ax
      u1 = b(jv+j) + b(jy+j)
      u2 = b(jo+j) + b(jt+j)
      u3 = b(jv+j) - b(jy+j)
      u4 = b(jo+j) - b(jt+j)
      b(jv+j) = b(jj+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ju+j) - 0.25 * u5
      bx = b(ju+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ju+j) = b(je+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(je+j) = bx
      a(jj+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jj+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jy+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jy+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jo+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jo+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jt+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jt+j) = si3*(t9+u10) + co3*(u9-t10)
      j = j + jump
  440 continue
!
      endif
!
!-----(end of loop across transforms)
!
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 5*la
  480 continue
!-----( end of loop on type ii radix-5 passes )
!-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
      return
      end subroutine gpfa5f

      end module fft_translation
!
! scatprops module: various subroutines for calculation of observables from the solution
!
!
!  last revised: 15 January 2011
!
      module scatprops
      use mpidefs
      use intrinsics
      use numconstants
      use specialfuncs
      use surface_subroutines
      use periodic_lattice_subroutines
      use spheredata
      use mie
      use translation
      use fft_translation
      implicit none
      contains

! The general sphere interaction driver
! LR formulation
! february 2013: number of rhs now a required argument list, and mpi comm is an option
!
         subroutine sphereinteraction(neqns,nrhs,ain,aout,initial_run, &
                    rhs_list,mpi_comm,con_tran,mie_mult,fft_option, &
                    store_matrix_option,skip_external_translation)
         implicit none
         logical :: initrun,rhslist(nrhs),contran(nrhs),miemult(nrhs), &
            fftopt,smopt,etopt
         logical, optional :: initial_run,rhs_list(nrhs),con_tran(nrhs), &
            mie_mult(nrhs),fft_option,store_matrix_option,skip_external_translation
         integer :: neqns,rank,numprocs,nrhs,nsend, &
                    mpicomm,rhs
         integer, optional :: mpi_comm
         complex(8)  :: aout_t(neqns,nrhs),ain_t(neqns,nrhs), &
                        ain(neqns,nrhs),aout(neqns,nrhs), &
                        aout_t2(neqns,nrhs)
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
         if(present(mie_mult)) then
            miemult=mie_mult
         else
            miemult=.true.
         endif
         if(present(fft_option)) then
            fftopt=fft_option
         else
            fftopt=fft_translation_option
         endif
         if(present(store_matrix_option)) then
            smopt=store_matrix_option
         else
            smopt=.true.
         endif
         if(present(skip_external_translation)) then
            etopt=.not.skip_external_translation
         else
            etopt=.true.
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
!
! sphere-to-sphere H (i-j) interaction
!
if(light_up) then
write(*,'('' si1 '',i3)') mstm_global_rank
call flush(6)
endif
         do rhs=1,nrhs
            if(contran(rhs)) then
               ain_t(:,rhs)=conjg(ain(:,rhs))
               if(miemult(rhs)) then
                  call multmiecoeffmult(neqns,1,-1,ain_t(:,rhs),ain_t(:,rhs))
               endif
            else
               ain_t(:,rhs)=ain(:,rhs)
            endif
         enddo

         call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
         aout_t=0.
         if(etopt) then
            if(periodic_lattice) then
!write(*,*) 'step 3'
!call flush(6)
if(light_up) then
write(*,'('' si2a '',i3)') mstm_global_rank
call flush(6)
endif
               call periodic_lattice_sphere_interaction(neqns,nrhs,ain_t,aout_t, &
                  store_matrix_option=smopt,initial_run=initrun, &
                  rhs_list=rhslist,mpi_comm=mpicomm,con_tran=contran)
            else
if(light_up) then
write(*,'('' si2b '',i3)') mstm_global_rank
call flush(6)
endif
               if(fftopt) then
                  call fft_external_to_external_expansion(neqns,nrhs,ain_t,aout_t, &
                          store_matrix_option=smopt,initial_run=initrun, &
                          rhs_list=rhslist,mpi_comm=mpicomm,con_tran=contran)
               else
                  call external_to_external_expansion(neqns,nrhs,ain_t,aout_t, &
                          store_matrix_option=smopt,initial_run=initrun, &
                          rhs_list=rhslist,mpi_comm=mpicomm,con_tran=contran)
               endif
            endif
         endif

         call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
         if(number_host_spheres.gt.0) then
            aout_t2=0.
            call external_to_internal_expansion(neqns,nrhs,ain_t,aout_t2, &
                 rhs_list=rhslist,mpi_comm=mpicomm,con_tran=contran)
            aout_t=aout_t+aout_t2
         endif

         if(plane_surface_present.and.(.not.periodic_lattice)) then
            aout_t2=0.
            call spheresurfaceinteraction(neqns,nrhs,ain_t,aout_t2, &
               initial_run=initrun,rhs_list=rhslist, &
               mpi_comm=mpicomm,con_tran=contran)
            aout_t=aout_t+aout_t2
         endif

         do rhs=1,nrhs
            if(.not.contran(rhs)) then
               if(miemult(rhs)) then
                  call multmiecoeffmult(neqns,1,1,aout_t(:,rhs),aout(:,rhs))
               else
                  aout(:,rhs)=aout_t(:,rhs)
               endif
            else
               aout(:,rhs)=conjg(aout_t(:,rhs))
            endif
         enddo

         if(numprocs.gt.1) then
            nsend=neqns*nrhs
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=aout, &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         endif
if(light_up) then
write(*,'('' si3 '',i3)') mstm_global_rank
call flush(6)
endif

         end subroutine sphereinteraction

!
!  determination of maximum orders for target--based expansions
!
!
!  last revised: 15 January 2011
!
         subroutine tranorders(eps,ntran,nodrt)
         implicit none
         integer :: nodrt,ntran(*),i,host,nodrmax
         real(8) :: r,eps,rpos0(3),xi0(3)
         complex(8) :: ri0,riext(2)
         nodrt=0
         nodrmax=max_mie_order
         do i=1,number_spheres
            host=host_sphere(i)
            call exteriorrefindex(i,riext)
            ri0=2.d0/(1.d0/riext(1)+1.d0/riext(2))
            if(host.eq.0) then
               rpos0=cluster_origin
            else
               rpos0(:)=sphere_position(:,host)
            endif
            xi0(:)=sphere_position(:,i)-rpos0(:)
            r=sqrt(dot_product(xi0,xi0))
            call tranordertest(r,ri0,sphere_order(i),eps,ntran(i))
            if(host.eq.0) nodrt=max(nodrt,ntran(i),nodrmax)
         enddo
         end subroutine tranorders
!
!  plane wave expansion coefficients at sphere origins.  uses a phase shift.
!
!
!  last revised: 15 January 2011
!
! may 2019: k is now rightmost column
         subroutine sphereplanewavecoef(alpha,sinc,dir,pmnp,excited_spheres)
         implicit none
         logical :: exsphere(number_spheres)
         logical, optional :: excited_spheres(number_spheres)
         integer :: p,i,dir
         real(8) :: alpha,sinc
         complex(8) :: pmnp(number_eqns,2)
         complex(8), allocatable :: pmnptot(:,:)
         if(present(excited_spheres)) then
            exsphere=excited_spheres
         else
            exsphere=.true.
         endif
         pmnp=0.d0
         do i=1,number_spheres
            if(host_sphere(i).ne.0) cycle
            if(.not.exsphere(i)) cycle
            allocate(pmnptot(sphere_block(i),2))
            if(gaussian_beam_constant.eq.0.d0) then
               call layerplanewavecoef(alpha,sinc,dir,sphere_position(:,i),sphere_order(i), &
                   pmnptot)
            else
               call layergaussbeamcoef(alpha,sinc,dir,sphere_position(:,i),sphere_order(i), &
                   pmnptot)
            endif
            do p=1,2
               pmnp(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),p) &
                  =pmnptot(1:sphere_block(i),p)
            enddo
            deallocate(pmnptot)
         enddo
         end subroutine sphereplanewavecoef

         subroutine layergaussbeamcoef(alpha,sinc,sdir,rpos,nodr,pmnp,include_direct, &
            include_indirect)
         implicit none
         logical :: incdir,incindir,shiftgb
         logical, optional :: include_direct,include_indirect
         integer :: p,incregion,sdir,nodr,layer,nodrgb
         real(8) :: alpha,rpos(3),sinc,cbinc,rtran(3),r,rshft(3),zs
         complex(8) :: pmnp(2*nodr*(nodr+2),2),riinc
         complex(8), allocatable :: pmnp0(:,:),ptvec(:)
         type(translation_data), pointer :: loc_tranmat
         type(translation_data), target :: tranmat
         if(sdir.eq.1) then
            riinc=layer_ref_index(0)
            incregion=0
         else
            riinc=layer_ref_index(number_plane_boundaries)
            incregion=number_plane_boundaries
         endif
         if(present(include_direct)) then
            incdir=include_direct
         else
            incdir=.true.
         endif
         if(present(include_indirect)) then
            incindir=include_indirect
         else
            incindir=.true.
         endif
         layer=layer_id(rpos(3))
         cbinc=dble(cdsqrt((1.d0-sinc/riinc)*(1.d0+sinc/riinc))*(3-2*sdir))
         pmnp=0.d0
         rtran=rpos(:)-gaussian_beam_focal_point(:)
         r=sqrt(sum(rtran**2))
         call tranordertest(r,riinc,nodr,1.d-6,nodrgb)
         nodrgb=max(nodrgb,nodr)
!write(*,'('' gb order:'',i8)') nodrgb
         nodrgb=min(nodrgb,max_t_matrix_order)
         allocate(pmnp0(2*nodrgb*(nodrgb+2),2))
         pmnp=0.d0
!write(*,*) 's1'
         call gaussianbeamcoef(alpha,cbinc,gaussian_beam_constant,nodrgb,pmnp0)
         if(layer.eq.incregion.and.incdir) then
            tranmat%matrix_calculated=.false.
            tranmat%vswf_type=1
            tranmat%translation_vector=rtran
            tranmat%refractive_index=(/riinc,riinc/)
            tranmat%rot_op=(nodrgb.ge.translation_switch_order)
            loc_tranmat=>tranmat
            do p=1,2
!write(*,*) 's2', p
               call coefficient_translation(nodrgb,2,nodr,2, &
                  pmnp0(:,p),pmnp(:,p),loc_tranmat)
            enddo
            if(.not.tranmat%zero_translation) then
               if(tranmat%rot_op) then
                  if(associated(tranmat%rot_mat)) deallocate(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                  nullify(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                  nullify(loc_tranmat)
               else
                  if(associated(tranmat%gen_mat)) deallocate(tranmat%gen_mat)
                  nullify(tranmat%gen_mat)
                  nullify(loc_tranmat)
               endif
            endif
         endif
         if(number_plane_boundaries.gt.0.and.incindir) then
            shiftgb=.false.
            if(sdir.eq.1.and.gaussian_beam_focal_point(3).ge.0.d0) then
               shiftgb=.true.
               rshft=(/0.d0,0.d0,-1.d-5-gaussian_beam_focal_point(3)/)
               zs=-1.d-5
            elseif(sdir.ne.1.and. &
               gaussian_beam_focal_point(3).lt.plane_boundary_position(number_plane_boundaries)) then
               shiftgb=.true.
               rshft=(/0.d0,0.d0,plane_boundary_position(number_plane_boundaries)+1.d-5-gaussian_beam_focal_point(3)/)
               zs=plane_boundary_position(number_plane_boundaries)+1.d-5
            else
               zs=gaussian_beam_focal_point(3)
            endif
            if(shiftgb) then
               allocate(ptvec(2*nodrgb*(nodrgb+2)))
               tranmat%matrix_calculated=.false.
               tranmat%vswf_type=1
               tranmat%translation_vector=rshft
               tranmat%refractive_index=(/riinc,riinc/)
               tranmat%rot_op=(nodrgb.ge.translation_switch_order)
               loc_tranmat=>tranmat
               do p=1,2
                  ptvec(:)=pmnp0(:,p)
                  pmnp0(:,p)=0.d0
                  call coefficient_translation(nodrgb,2,nodrgb,2, &
                     ptvec(:),pmnp0(:,p),loc_tranmat)
               enddo
               if(.not.tranmat%zero_translation) then
                  if(tranmat%rot_op) then
                     if(associated(tranmat%rot_mat)) deallocate(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                     nullify(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                     nullify(loc_tranmat)
                  else
                     if(associated(tranmat%gen_mat)) deallocate(tranmat%gen_mat)
                     nullify(tranmat%gen_mat)
                     nullify(loc_tranmat)
                  endif
               endif
               deallocate(ptvec)
            endif

!write(*,*) 's3'
!            allocate(rmat(2*nodr*(nodr+2),2*nodrgb*(nodrgb+2)))
            allocate(ptvec(2*nodr*(nodr+2)))
!write(*,*) 's4',nodr,nodrgb

            do p=1,2
               ptvec=0.d0
               call plane_interaction(nodr,nodrgb, &
                  rtran(1),rtran(2),zs,rpos(3), &
                  ptvec,index_model=2,lr_transformation=.true., &
                  make_symmetric=.false.,propagating_directions_only=.true., &
                  source_vector=pmnp0(:,p))
!               pmnp(:,p)=pmnp(:,p)+0.5d0*matmul(rmat(:,:),pmnp0(:,p))
               pmnp(:,p)=pmnp(:,p)+0.5d0*ptvec(:)
            enddo
!write(*,*) 's5'
            deallocate(ptvec)
         endif
         deallocate(pmnp0)
         end subroutine layergaussbeamcoef

         subroutine phase_shift(amnp,dir)
         implicit none
         integer :: dir,i
         real(8), save :: ilv(2)
         complex(8) :: amnp(number_eqns,2)

         if(dir.eq.1) then
            do i=1,number_spheres
               if(host_sphere(i).ne.0) cycle
               amnp(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),1:2) &
                  =amnp(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),1:2) &
                  *cdexp((0.d0,-1.d0)*sum(incident_lateral_vector*sphere_position(1:2,i)))
            enddo
            ilv=incident_lateral_vector
!            incident_lateral_vector=0.d0
         else
!            incident_lateral_vector=ilv
            do i=1,number_spheres
               if(host_sphere(i).ne.0) cycle
               amnp(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),1:2) &
                  =amnp(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),1:2) &
                  *cdexp((0.d0,1.d0)*sum(incident_lateral_vector*sphere_position(1:2,i)))
            enddo
         endif
         end subroutine phase_shift
!
!  translation of sphere-based expansions to common target origin
!
!
!  last revised: 15 January 2011
!  april 2012: l/r formulation: amnp is in l/r basis, amnp0 is in e/m basis
!  february 2013: number rhs, mpi comm options.   This is for general sphere configurations.
!
         subroutine merge_to_common_origin(nodrt,amnp,amnp0,number_rhs, &
            single_sphere,origin_position,mpi_comm,merge_procs,merge_radius)
         implicit none
         logical :: mergeprocs
         logical, optional :: merge_procs
         integer :: nodrt,i,m,n,nblk,ntrani,nrhs,task, &
                    rank,numprocs,nsend,proc,mpicomm,noi, &
                    startsphere,endsphere,rhs
         integer, optional :: number_rhs,mpi_comm,single_sphere
         real(8) :: r0(3),mrad,ri0
         real(8), optional :: origin_position(3),merge_radius
         complex(8) :: amnp(number_eqns,*),amnp0(0:nodrt+1,nodrt,2,*),rimedium(2)
         complex(8), allocatable :: amnpt(:,:,:)
         type(translation_data), pointer :: loc_tranmat
         type(translation_data), target :: tranmat

         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(origin_position)) then
            r0=origin_position
         else
            r0=0.d0
         endif
         if(present(merge_radius)) then
            mrad=merge_radius
         else
            mrad=1.d10
         endif
         if(present(number_rhs)) then
            nrhs=number_rhs
         else
            nrhs=1
         endif
         if(present(merge_procs)) then
            mergeprocs=merge_procs
         else
            mergeprocs=.true.
         endif
         if(present(single_sphere)) then
            startsphere=single_sphere
            endsphere=single_sphere
         else
            startsphere=1
            endsphere=number_spheres
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         amnp0(:,:,:,1:nrhs)=(0.d0,0.d0)
         task=0
         do i=startsphere,endsphere
            nblk=sphere_order(i)*(sphere_order(i)+2)*2
            if(host_sphere(i).eq.0) then
               ri0=sqrt(sum(r0-sphere_position(:,i))**2)
               if(ri0.gt.mrad) cycle
               task=task+1
               proc=mod(task,numprocs)
               if(proc.eq.rank) then
                  call exteriorrefindex(i,rimedium)
                  noi=sphere_order(i)
                  ntrani=min(nodrt,translation_order(i))
!                  ntrani=max(ntrani,noi)
                  allocate(amnpt(0:ntrani+1,ntrani,2))
                  amnpt=(0.d0,0.d0)
                  tranmat%matrix_calculated=.false.
                  tranmat%vswf_type=1
                  tranmat%translation_vector=r0-sphere_position(:,i)
                  tranmat%refractive_index=rimedium
                  tranmat%rot_op=(max(noi,ntrani).ge.translation_switch_order)
                  loc_tranmat=>tranmat
                  do rhs=1,nrhs
                     call coefficient_translation(noi,2,ntrani,2, &
                        amnp(sphere_offset(i)+1:sphere_offset(i)+nblk,rhs),amnpt, &
                        loc_tranmat)
                     do n=1,ntrani
                        do m=0,ntrani+1
                           amnp0(m,n,1,rhs)=amnp0(m,n,1,rhs) &
                                +amnpt(m,n,1)+amnpt(m,n,2)
                           amnp0(m,n,2,rhs)=amnp0(m,n,2,rhs) &
                                +amnpt(m,n,1)-amnpt(m,n,2)
                        enddo
                     enddo
                  enddo
                  if(.not.tranmat%zero_translation) then
                     if(tranmat%rot_op) then
                        if(associated(tranmat%rot_mat)) deallocate(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                        nullify(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                        nullify(loc_tranmat)
                     else
                        if(associated(tranmat%gen_mat)) deallocate(tranmat%gen_mat)
                        nullify(tranmat%gen_mat)
                        nullify(loc_tranmat)
                     endif
                  endif
                  deallocate(amnpt)
               endif
            endif
         enddo
         if(numprocs.gt.1.and.mergeprocs) then
            nsend=2*nodrt*(nodrt+2)*nrhs
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=amnp0, &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         endif
         end subroutine merge_to_common_origin

         subroutine distribute_from_common_origin(nodrt,amnp0,amnp,number_rhs, &
            origin_position,origin_host,vswf_type,mpi_comm,merge_procs, &
            single_sphere)
         implicit none
         logical :: mergeprocs
         logical, optional :: merge_procs
         integer :: nodrt,i,nblk,ntrani,nrhs,task,rhs,startsphere,endsphere, &
                    rank,numprocs,nsend,proc,mpicomm,noi,vtype,ohost
         integer, optional :: number_rhs,mpi_comm,vswf_type,origin_host,single_sphere
         real(8) :: r0(3)
         real(8), optional :: origin_position(3)
         complex(8) :: amnp(number_eqns,*),amnp0(0:nodrt+1,nodrt,2,*),rimedium(2)
         type(translation_data), pointer :: loc_tranmat
         type(translation_data), target :: tranmat

         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(origin_position)) then
            r0=origin_position
         else
            r0=0.d0
         endif
         if(present(origin_host)) then
            ohost=origin_host
         else
            ohost=0
         endif
         if(present(number_rhs)) then
            nrhs=number_rhs
         else
            nrhs=1
         endif
         if(present(vswf_type)) then
            vtype=vswf_type
         else
            vtype=1
         endif
         if(present(merge_procs)) then
            mergeprocs=merge_procs
         else
            mergeprocs=.true.
         endif
         if(present(single_sphere)) then
            if(single_sphere.ne.0) then
               startsphere=single_sphere
               endsphere=single_sphere
            else
               startsphere=1
               endsphere=number_spheres
            endif
         else
            startsphere=1
            endsphere=number_spheres
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         amnp(1:number_eqns,1:nrhs)=(0.d0,0.d0)
         task=0
         do i=startsphere,endsphere
            nblk=sphere_order(i)*(sphere_order(i)+2)*2
            if(host_sphere(i).eq.ohost) then
               task=task+1
               proc=mod(task,numprocs)
               if(proc.eq.rank) then
                  call exteriorrefindex(i,rimedium)
                  noi=sphere_order(i)
                  ntrani=nodrt
                  tranmat%matrix_calculated=.false.
                  tranmat%vswf_type=vtype
                  tranmat%translation_vector=sphere_position(:,i)-r0
                  tranmat%refractive_index=rimedium
                  tranmat%rot_op=(max(noi,ntrani).ge.translation_switch_order)
                  loc_tranmat=>tranmat
                  do rhs=1,nrhs
                     call coefficient_translation(ntrani,2,noi,2, &
                        amnp0(0:nodrt+1,1:nodrt,1:2,rhs), &
                        amnp(sphere_offset(i)+1:sphere_offset(i)+nblk,rhs), &
                        loc_tranmat)
                  enddo
                  if(.not.tranmat%zero_translation) then
                     if(tranmat%rot_op) then
                        if(associated(tranmat%rot_mat)) deallocate(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                        nullify(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
                        nullify(loc_tranmat)
                     else
                        if(associated(tranmat%gen_mat)) deallocate(tranmat%gen_mat)
                        nullify(tranmat%gen_mat)
                        nullify(loc_tranmat)
                     endif
                  endif
               endif
            endif
         enddo
         if(numprocs.gt.1.and.mergeprocs) then
            nsend=number_eqns*nrhs
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=amnp, &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         endif
         end subroutine distribute_from_common_origin
!
!  general efficiency factor calcuation
!  L/R formulation
!  April 2012
!  march 2013: something is always changing in this one
!
         subroutine lrsphereqeff(ri0,nodr,npol,xsp,anp,gnpinc,gnp,qe,qs,qa)
         implicit none
         integer :: nodr,npol,m,n,p,p1,p2,k,ma,na,s,t,ss,st
         real(8) :: qe(2*npol-1),qa(2*npol-1),qs(2*npol-1),const,xsp,qi(2*npol-1)
         complex(8) :: anp(0:nodr+1,nodr,2,npol),gnp(0:nodr+1,nodr,2,npol), &
                       gnpinc(0:nodr+1,nodr,2,npol), &
                       psi(0:nodr,2),xi(0:nodr,2),psip(2),xip(2),xri(2),&
                       rib,aamat(2,2),ggmat(2,2),agmat(2,2), &
                       gamat(2,2),cterm,ri0(2)
         qe=0.d0
         qa=0.d0
         qs=0.d0
         qi=0.d0
         xri=xsp*ri0
         rib=2.d0/(1.d0/ri0(1)+1.d0/ri0(2))
         do p=1,2
            call cricbessel(nodr,xri(p),psi(0,p))
            call crichankel(nodr,xri(p),xi(0,p))
         enddo
         do n=1,nodr
            do s=1,2
               psip(s)=psi(n-1,s)-n*psi(n,s)/xri(s)
               xip(s)=xi(n-1,s)-n*xi(n,s)/xri(s)
            enddo
            do s=1,2
               ss=(-1)**s
               do t=1,2
                  st=(-1)**t
                  cterm=dcmplx(0.d0,1.d0)*conjg(1./ri0(t))/ri0(s)
                  aamat(s,t)=cterm*(xip(s)*conjg(xi(n,t)) &
                        -ss*st*xi(n,s)*conjg(xip(t)))
                  ggmat(s,t)=cterm*(psip(s)*conjg(psi(n,t)) &
                        -ss*st*psi(n,s)*conjg(psip(t)))
                  agmat(s,t)=cterm*(xip(s)*conjg(psi(n,t)) &
                        -ss*st*xi(n,s)*conjg(psip(t)))
                  gamat(s,t)=cterm*(psip(s)*conjg(xi(n,t)) &
                        -ss*st*psi(n,s)*conjg(xip(t)))
               enddo
            enddo
            do m=-n,n
               if(m.le.-1) then
                  ma=n+1
                  na=-m
               else
                  ma=m
                  na=n
               endif
               do p1=1,npol
                  do p2=1,npol
                     if(p1.eq.1.and.p2.eq.1) then
                        k=1
                        const=1.d0
                     elseif(p1.eq.2.and.p2.eq.2) then
                        k=2
                        const=1.d0
                     else
                        k=3
                        const=.5d0
                     endif
                     do s=1,2
                        do t=1,2
!                           qi(k)=qi(k)+const*ggmat(s,t)*(gnp(ma,na,s,p1)-gnpinc(ma,na,s,p1)) &
!                             *conjg(rib*(gnp(ma,na,t,p2)-gnpinc(ma,na,t,p2)))
                           qi(k)=qi(k)+const*agmat(s,t)*anp(ma,na,s,p1)*conjg(gnp(ma,na,t,p2)-gnpinc(ma,na,t,p2)) &
                               *(conjg(rib)+(-1)**(s+t)*rib)
                           qa(k)=qa(k)+const &
                              *(aamat(s,t)*anp(ma,na,s,p1)*conjg(rib*anp(ma,na,t,p2)) &
                              +ggmat(s,t)*gnp(ma,na,s,p1)*conjg(rib*gnp(ma,na,t,p2)) &
                              +agmat(s,t)*anp(ma,na,s,p1)*conjg(gnp(ma,na,t,p2)) &
                               *(conjg(rib)+(-1)**(s+t)*rib))
                           qs(k)=qs(k)-const &
                              *aamat(s,t)*anp(ma,na,s,p1)*conjg(rib*anp(ma,na,t,p2))
                           qe(k)=qe(k)+const &
                              *(agmat(s,t)*anp(ma,na,s,p1)*conjg(gnpinc(ma,na,t,p2)) &
                              *(conjg(rib)+(-1)**(s+t)*rib))
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         do k=1,2*npol-1
!            qi(k)=qi(k)*2./xsp/xsp
            qe(k)=qe(k)*2./xsp/xsp
            qa(k)=qa(k)*2./xsp/xsp
            qs(k)=(qs(k))*2./xsp/xsp
!write(*,'(i5,4e13.5)') k,qe(k),qa(k)+qs(k)-qi(k)
         enddo

         end subroutine lrsphereqeff
!
! calling routine for efficiency calculation
! april 2012: lr formulation
! february 2013:  number of rhs, mpi comm options added.
!
         subroutine qefficiencyfactors(nsphere,npol,amnp,gmnp0,qeff, &
                    mpi_comm)
         implicit none
         integer :: nsphere,i,nblk,noff,neqns,nodr, &
                    npol,mpicomm,p, &
                    b11,b12,rank,numprocs
         integer, optional :: mpi_comm
         real(8) :: qeff(3,2*npol-1,nsphere),qe(2*npol-1),qa(2*npol-1),qs(2*npol-1)
         real(8), allocatable :: qeffi(:,:)
         complex(8) :: amnp(number_eqns,npol),gmnp0(number_eqns,npol),ri0(2)
         complex(8), allocatable :: amnpi(:,:,:,:),gmnpi(:,:,:,:),fmnpi(:,:,:,:), &
                     gmnp(:,:)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         neqns=number_eqns
         allocate(gmnp(neqns,npol),qeffi(3,2*npol-1))
         gmnp=0.d0
if(light_up) then
write(*,'('' qe1 '',i3)') mstm_global_rank
call flush(6)
endif
         do p=1,npol
            if(number_host_spheres.eq.0) then
               noff=0
               do i=1,number_spheres
                  nodr=sphere_order(i)
                  nblk=2*nodr*(nodr+2)
                  b11=noff+1
                  b12=b11+nblk-1
                  call onemiecoeffmult(i,nodr,amnp(b11:b12,p),gmnp(b11:b12,p),'i')
                  noff=noff+nblk*number_field_expansions(i)
               enddo
            else
               call sphereinteraction(neqns,1,amnp(:,p),gmnp(:,p), &
                    mpi_comm=mpicomm,mie_mult=(/.false./), &
                    store_matrix_option=.false.)
               gmnp(:,p)=gmnp(:,p)+gmnp0(:,p)
            endif
         enddo
if(light_up) then
write(*,'('' qe2 '',i3)') mstm_global_rank
call flush(6)
endif
         qeff(1:3,1:2*npol-1,1:nsphere)=0.d0
         noff=0
         do i=1,number_spheres
            nodr=sphere_order(i)
            nblk=2*nodr*(nodr+2)
            allocate(amnpi(0:nodr+1,nodr,2,npol), &
                     gmnpi(0:nodr+1,nodr,2,npol), &
                     fmnpi(0:nodr+1,nodr,2,npol))
            call exteriorrefindex(i,ri0)
!            ri0=sphere_ref_index(:,host_sphere(i))
            b11=noff+1
            b12=b11+nblk-1
            do p=1,npol
               amnpi(0:nodr+1,1:nodr,1:2,p)=reshape(amnp(b11:b12,p), &
                     (/nodr+2,nodr,2/))
               fmnpi(0:nodr+1,1:nodr,1:2,p)=reshape(gmnp0(b11:b12,p), &
                     (/nodr+2,nodr,2/))
               gmnpi(0:nodr+1,1:nodr,1:2,p)=reshape(gmnp(b11:b12,p), &
                     (/nodr+2,nodr,2/))
            enddo
            call lrsphereqeff(ri0,nodr,npol,sphere_radius(i),amnpi,fmnpi,gmnpi, &
               qe,qs,qa)
            qeffi(1,:)=qe(:)
            qeffi(3,:)=qs(:)
            qeffi(2,:)=qa(:)
            noff=noff+nblk*number_field_expansions(i)
            deallocate(gmnpi,amnpi,fmnpi)
            if(npol.eq.1) then
               qeff(:,1,i)=qeffi(:,1)
            else
               qeff(:,1,i)=.5*(qeffi(:,1)+qeffi(:,2))
               qeff(:,2,i)=qeffi(:,1)
               qeff(:,3,i)=qeffi(:,2)
            endif
         enddo
         deallocate(gmnp,qeffi)
if(light_up) then
write(*,'('' qe3 '',i3)') mstm_global_rank
call flush(6)
endif
         end subroutine qefficiencyfactors

         subroutine qtotcalc(nsphere,nrow,xgeff,qeffp,qabsvol,qefftot)
         implicit none
         integer :: nsphere,nrow,i,j
         real(8) :: qeffp(3,nrow,nsphere),qefftot(3,nrow),qabsvol(nrow,nsphere),xgeff
         do i=1,nsphere
            qabsvol(:,i)=qeffp(2,:,i)
            do j=1,nsphere
               if(host_sphere(j).eq.i) then
                  qabsvol(:,i)=qabsvol(:,i)-qeffp(2,:,j) &
                     *sphere_radius(j)**2 &
                     /sphere_radius(i)**2
               endif
            enddo
         enddo
!         do i=1,nsphere
!            if(dimag(sphere_ref_index(1,i)).eq.0.d0 &
!             .and.dimag(sphere_ref_index(2,i)).eq.0.d0) then
!               qabsvol(:,i)=0.d0
!            endif
!         enddo
         qeffp(2,:,:)=0.d0
         do i=1,nsphere
            qeffp(2,:,i)=qabsvol(:,i)
            do j=1,nsphere
               if(host_sphere(j).eq.i) then
                  qeffp(2,:,i)=qeffp(2,:,i)+qabsvol(:,j) &
                     *sphere_radius(j)**2 &
                     /sphere_radius(i)**2
               endif
            enddo
         enddo
         qefftot=0.
         do i=1,nsphere
            if(host_sphere(i).eq.0) then
               qefftot(:,:)=qefftot(:,:)+qeffp(:,:,i)*sphere_radius(i)**2/xgeff**2
            endif
         enddo
         end subroutine qtotcalc

         subroutine waveguide_mode_scattering(amnp,qsevan,mpi_comm)
         implicit none
         integer :: i,p,pole,mpicomm
         integer, optional :: mpi_comm
         real(8) :: qsevan(2)
         complex(8) :: amnp(number_eqns,2)
         complex(8), allocatable :: amn(:),gmn(:)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         allocate(amn(number_eqns),gmn(number_eqns))
         pole_integration=.true.
         qsevan=0.
         do pole=1,number_singular_points
            pole_integration_s=singular_points(pole)
            recalculate_surface_matrix=.true.
            do p=1,2
               amn(:)=amnp(:,p)
               gmn=0.d0
               call sphereinteraction(number_eqns,1,amn,gmn, &
                  mie_mult=(/.false./),initial_run=.true.,skip_external_translation=.true., &
                  mpi_comm=mpicomm)
               do i=1,number_spheres
                  if(host_sphere(i).ne.0) cycle
                  call lr_mode_transformation(sphere_order(i),amn(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i)), &
                    amn(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i)))
                  call lr_mode_transformation(sphere_order(i),gmn(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i)), &
                    gmn(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i)))
                  qsevan(p)=qsevan(p)+sum(dble(dconjg(amn(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i))) &
                     *gmn(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i)))) &
                     /layer_ref_index(sphere_layer(i))
               enddo
            enddo
         enddo
         qsevan=qsevan*2.d0/cross_section_radius**2
         deallocate(gmn,amn)
         pole_integration=.false.
         end subroutine waveguide_mode_scattering

         subroutine boundary_scattering(amn,qbsca,mpi_comm)
         implicit none
         integer :: i,j,rank,numprocs,mpicomm,task,nsend,k
         integer, optional :: mpi_comm
         real(8) :: qbsca(2,0:number_plane_boundaries+1),qt(2),targetz, &
            qt2(2,0:number_plane_boundaries+1)
         complex(8) :: amn(number_eqns,2)
         complex(8), allocatable :: amn1(:,:),amn2(:,:)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='size',mpi_size=numprocs, &
              mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank, &
              mpi_comm=mpicomm)

         qbsca=0.d0
         task=0
         do k=0,number_plane_boundaries+1
            if(k.eq.0) then
               targetz=bot_boundary
            elseif(k.eq.number_plane_boundaries+1) then
               targetz=top_boundary
            elseif(k.eq.1) then
               targetz=plane_boundary_position(k)+1.d-8
            else
               targetz=plane_boundary_position(k)-1.d-8
            endif
            if(k.eq.0.and.dimag(layer_ref_index(0)).gt.1.d-6) then
               qbsca(:,k)=0.d0
            elseif(k.eq.number_plane_boundaries+1.and. &
             dimag(layer_ref_index(number_plane_boundaries)).gt.1.d-6) then
               qbsca(:,k)=0.d0
            else
               do i=1,number_spheres
                  if(host_sphere(i).ne.0) cycle
                  allocate(amn1(sphere_block(i),2))
                  amn1=amn(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),1:2)
                  do j=1,number_spheres
                     if(host_sphere(j).ne.0) cycle
                     task=task+1
                     if(mod(task,numprocs).ne.rank) cycle
                     allocate(amn2(sphere_block(j),2))
                     amn2=amn(sphere_offset(j)+1:sphere_offset(j)+sphere_block(j),1:2)
                     call sphere_boundary_scattering(sphere_order(i),sphere_position(:,i), &
                        amn1,sphere_order(j),sphere_position(:,j),amn2,targetz,qt)
                     qbsca(:,k)=qbsca(:,k)+qt
                     deallocate(amn2)
                  enddo
                  deallocate(amn1)
               enddo
            endif
         enddo
         qbsca=qbsca/cross_section_radius**2
         if(numprocs.gt.1) then
            nsend=2*(number_plane_boundaries+2)
            qt2=qbsca
            call mstm_mpi(mpi_command='allreduce',mpi_number=nsend, &
                 mpi_send_buf_dp=qt2,mpi_recv_buf_dp=qbsca, &
                 mpi_operation=mstm_mpi_sum, &
                 mpi_comm=mpicomm)
         endif
         end subroutine boundary_scattering

         subroutine boundary_extinction(amn,alpha,sinc,dir,qbext,common_origin)
         implicit none
         logical :: comorg
         logical, optional :: common_origin
         integer :: k,dir
         real(8) :: qt(2),targetz,qbext(2,0:number_plane_boundaries+1),alpha,sinc
         complex(8) :: amn(*)
         if(present(common_origin)) then
            comorg=common_origin
         else
            comorg=.false.
         endif
         qbext=0.d0
         do k=0,number_plane_boundaries+1
            if(k.eq.0) then
               targetz=bot_boundary
            elseif(k.eq.number_plane_boundaries+1) then
               targetz=top_boundary
            elseif(k.eq.1) then
               targetz=plane_boundary_position(k)+1.d-8
            else
               targetz=plane_boundary_position(k)-1.d-8
            endif
            if(k.eq.0.and.dimag(layer_ref_index(0)).gt.1.d-6.and.dir.eq.2) then
               qbext(:,k)=0.d0
            elseif(k.eq.number_plane_boundaries+1.and. &
             dimag(layer_ref_index(number_plane_boundaries)).gt.1.d-6.and.dir.eq.1) then
               qbext(:,k)=0.d0
            else
               call extinction_theorem(amn,sinc,dir,alpha,targetz,qt,common_origin=comorg)
               qbext(:,k)=qbext(:,k)+qt
            endif
         enddo
         end subroutine boundary_extinction

         subroutine common_origin_hemispherical_scattering(amn,qbsca)
         implicit none
         integer :: k
         real(8) :: qbsca(2,0:1),targetz
         complex(8) :: amn(2*t_matrix_order*(t_matrix_order+2),2)

         qbsca=0.d0
         do k=0,1
            if(k.eq.0) then
               targetz=bot_boundary
            elseif(k.eq.1) then
               targetz=top_boundary
            endif
            call sphere_boundary_scattering(t_matrix_order,(/0.d0,0.d0,0.d0/), &
               amn,t_matrix_order,(/0.d0,0.d0,0.d0/),amn,targetz,qbsca(:,k),lr_to_mode=.false.)
         enddo
         qbsca=qbsca/cross_section_radius**2
         end subroutine common_origin_hemispherical_scattering

         subroutine extinction_theorem(amnp,sinc,sdir,alpha,targetz,qe,common_origin)
         implicit none
         logical :: comorg
         logical, optional :: common_origin
         integer :: sdir,p,tlay
         real(8) :: sinc,alpha,qe(2),targetz,sourcez,const3(3,2)
         complex(8) :: amnp(*),s,sourceri,targetri,gfs(2,2,2),skz,tkz,sa(4),amp(2,2)
         if(present(common_origin)) then
            comorg=common_origin
         else
            comorg=.false.
         endif
         if(sdir.eq.1) then
            sourceri=layer_ref_index(0)
         else
            sourceri=layer_ref_index(number_plane_boundaries)
         endif
         sourcez=incident_field_boundary
         tlay=layer_id(targetz)
         targetri=layer_ref_index(tlay)
         s=sinc
         call layer_gf(s,sourcez,targetz,gfs,skz,tkz,include_direct=.true.)
         do p=1,2
            if(comorg) then
               call common_origin_amplitude_matrix(amnp,s,alpha,targetz,p,sa)
            else
               call multiple_origin_amplitude_matrix(amnp,s,alpha,targetz,p,sa)
            endif
            amp(p,1)=sa(2)
            amp(p,2)=sa(1)
         enddo
         const3(1,1)=16.d0*(cdabs(tkz)**2+s**2/cdabs(targetri)**2)*dble(targetri*tkz)
         const3(2,1)=-const3(1,1)
         const3(3,1)=16.d0*(cdabs(tkz)**2-s**2/cdabs(targetri)**2)*dimag(targetri*tkz)
         const3(1,2)=16.d0*dble(targetri*tkz)
         const3(2,2)=-const3(1,2)
         const3(3,2)=-16.d0*dimag(targetri*tkz)
         do p=1,2
            qe(p)=(const3(1,p)*dble(dconjg(gfs(1,sdir,p))*amp(1,p))+const3(2,p)*dble(dconjg(gfs(2,sdir,p))*amp(2,p)) &
               +const3(3,p)*dimag(dconjg(gfs(2,sdir,p))*amp(1,p)-dconjg(gfs(1,sdir,p))*amp(2,p))) &
               /incident_field_scale(p)
         enddo
         qe=qe/cross_section_radius**2
         end subroutine extinction_theorem

         subroutine periodic_lattice_scattering(amnp,qsca,scat_mat,krho_vec,num_dirs,dry_run)
         implicit none
         logical :: prop,calcsmat,dryrun
         logical, optional :: dry_run
         integer :: nmax,dir,n,i,q,p,ix,iy,nang,istrt
         integer, optional :: num_dirs(2)
         real(8) :: qsca(2,2),k0x,k0y,wx,wy,targetz,kx,ky,krho,phi,smscale
         real(8), optional :: scat_mat(32,*),krho_vec(2,*)
         complex(8) :: amnp(*),ri,s,kz,sa(4)

         if(present(dry_run)) then
            dryrun=dry_run
         else
            dryrun=.false.
         endif
         calcsmat=present(scat_mat)
         k0x=incident_lateral_vector(1)
         k0y=incident_lateral_vector(2)
         wx=cell_width(1)
         wy=cell_width(2)
         nmax=ceiling(maxval(cell_width(1:2))/pi)
         if(present(num_dirs)) num_dirs=0
         smscale=1.d0/cross_section_radius**2*64.d0*pi/wx/wy
         do dir=1,2
            if(dir.eq.1) then
               istrt=17
            else
               istrt=1
            endif
            if(dir.eq.1) then
               targetz=top_boundary
               ri=layer_ref_index(number_plane_boundaries)
            else
               targetz=bot_boundary
               ri=layer_ref_index(0)
            endif
            kx=k0x
            ky=k0y
            krho=sqrt(kx*kx+ky*ky)
            s=krho
            kz=ri*cdsqrt(1.d0-s*s/ri/ri)
            if(kx.eq.0.d0.and.ky.eq.0.d0) then
               phi=0.d0
            else
               phi=datan2(ky,kx)
            endif
            nang=1
            if(.not.dryrun) then
               call multiple_origin_amplitude_matrix(amnp,s,phi,targetz,dir,sa)
               qsca(1,dir)=kz*(cdabs(sa(2))**2+cdabs(sa(4))**2)
               qsca(2,dir)=kz*(cdabs(sa(1))**2+cdabs(sa(3))**2)
               if(calcsmat) then
                  sa=sa*cdsqrt(kz)
                  call amplitude_to_scattering_matrix(sa,scat_mat(istrt:istrt+15,nang))
                  scat_mat(istrt:istrt+15,nang)=scat_mat(istrt:istrt+15,nang)*smscale
                  if(present(krho_vec)) then
                     krho_vec(1,nang)=kx
                     krho_vec(2,nang)=ky
                  endif
               endif
            endif

            do n=1,nmax
               prop=.false.
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
                  kx=2.d0*pi*dble(ix)/wx+k0x
                  ky=2.d0*pi*dble(iy)/wy+k0y
                  krho=sqrt(kx*kx+ky*ky)
                  if(krho.le.cdabs(ri)) then
!                  if(krho.le.1.d0) then
                     prop=.true.
                     s=krho
                     kz=ri*cdsqrt(1.d0-s*s/ri/ri)
!                     kz=cdsqrt(1.d0-s*s)
                     if(kx.eq.0.d0.and.ky.eq.0.d0) then
                        phi=0.d0
                     else
                        phi=datan2(ky,kx)
                     endif
                     nang=nang+1
                     if(.not.dryrun) then
                        call multiple_origin_amplitude_matrix(amnp,s,phi,targetz,dir,sa)
                        qsca(1,dir)=qsca(1,dir)+kz*(cdabs(sa(2))**2+cdabs(sa(4))**2)
                        qsca(2,dir)=qsca(2,dir)+kz*(cdabs(sa(1))**2+cdabs(sa(3))**2)
                        if(calcsmat) then
                           sa=sa*cdsqrt(kz)
                           call amplitude_to_scattering_matrix(sa,scat_mat(istrt:istrt+15,nang))
                           scat_mat(istrt:istrt+15,nang)=scat_mat(istrt:istrt+15,nang)*smscale
                           if(present(krho_vec)) then
                              krho_vec(1,nang)=kx
                              krho_vec(2,nang)=ky
                           endif
                        endif
                     endif
                  endif
               enddo
               if(.not.prop) exit
            enddo
            if(present(num_dirs)) num_dirs(3-dir)=nang
         enddo
         qsca=qsca/cross_section_radius**2*64.d0*pi/wx/wy
         end subroutine periodic_lattice_scattering

         subroutine multiple_origin_amplitude_matrix(amnp,s,phi,targetz,dir,sa)
         implicit none
         integer :: p,i,dir
         real(8) :: phi,targetz
         complex(8) :: amnp(number_eqns,2),sa(4),sat(4),s
         complex(8), allocatable :: pmnpi(:,:),amnpi(:,:)

         sa=0.d0
         do i=1,number_spheres
            if(host_sphere(i).ne.0) cycle
            allocate(pmnpi(sphere_block(i),2),amnpi(sphere_block(i),2))
            call layervsh(s,phi,targetz,dir,sphere_position(:,i),sphere_order(i),pmnpi)
            amnpi(:,:)=amnp(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),:)
            do p=1,2
               call lr_mode_transformation(sphere_order(i),amnpi(:,p),amnpi(:,p))
            enddo
            sat(1)=sum(pmnpi(:,2)*amnpi(:,2))*0.5d0
            sat(2)=sum(pmnpi(:,1)*amnpi(:,1))*0.5d0
            sat(3)=-sum(pmnpi(:,1)*amnpi(:,2))*0.5d0
            sat(4)=-sum(pmnpi(:,2)*amnpi(:,1))*0.5d0
            sa(:)=sa(:)+sat(:)
            deallocate(pmnpi,amnpi)
         enddo
         end subroutine multiple_origin_amplitude_matrix

         subroutine common_origin_amplitude_matrix(amnp,s,phi,targetz,dir,sa)
         implicit none
         integer :: dir
         real(8) :: phi,targetz
         complex(8) :: amnp(2*t_matrix_order*(t_matrix_order+2),2),sa(4),s,pmnp(2*t_matrix_order*(t_matrix_order+2),2)
         call layervsh(s,phi,targetz,dir,cluster_origin,t_matrix_order,pmnp)
         sa(1)=0.5d0*sum(pmnp(:,2)*amnp(:,2))
         sa(2)=0.5d0*sum(pmnp(:,1)*amnp(:,1))
         sa(3)=-0.5d0*sum(pmnp(:,1)*amnp(:,2))
         sa(4)=-0.5d0*sum(pmnp(:,2)*amnp(:,1))
         end subroutine common_origin_amplitude_matrix

         subroutine amplitude_to_scattering_matrix(sa,sm)
         implicit none
         integer :: i,j
         real(8) :: sm(4,4)
         complex(8) :: sa(4),sp(4,4)
         do i=1,4
            do j=1,4
               sp(i,j)=sa(i)*dconjg(sa(j))
            enddo
         enddo
         sm(1,1)=sp(1,1)+sp(2,2)+sp(3,3)+sp(4,4)
         sm(1,2)=-sp(1,1)+sp(2,2)-sp(3,3)+sp(4,4)
         sm(2,1)=-sp(1,1)+sp(2,2)+sp(3,3)-sp(4,4)
         sm(2,2)=sp(1,1)+sp(2,2)-sp(3,3)-sp(4,4)
         sm(3,3)=2.*(sp(1,2)+sp(3,4))
         sm(3,4)=2.*dimag(sp(2,1)+sp(4,3))
         sm(4,3)=2.*dimag(sp(1,2)-sp(3,4))
         sm(4,4)=2.*(sp(1,2)-sp(3,4))
         sm(1,3)=2.*(sp(2,3)+sp(1,4))
         sm(3,1)=2.*(sp(2,4)+sp(1,3))
         sm(1,4)=-2.*dimag(sp(2,3)-sp(1,4))
         sm(4,1)=-2.*dimag(sp(4,2)+sp(1,3))
         sm(2,3)=-2.*(sp(2,3)-sp(1,4))
         sm(3,2)=-2.*(sp(2,4)-sp(1,3))
         sm(2,4)=-2.*dimag(sp(2,3)+sp(1,4))
         sm(4,2)=-2.*dimag(sp(4,2)-sp(1,3))
         end subroutine amplitude_to_scattering_matrix

         subroutine multiple_origin_scatteringmatrix(amnp,ct,phi,csca,sa,sm)
         implicit none
         integer :: dir
         real(8) :: ct,phi,sm(4,4),csca,targetz
         complex(8) :: amnp(number_eqns,2),sa(4),ri,s

         sa=0.d0
         if(ct.le.0.d0) then
            ri=layer_ref_index(0)
            dir=2
            targetz=bot_boundary
         else
            ri=layer_ref_index(number_plane_boundaries)
            dir=1
            targetz=top_boundary
         endif
         s=dble(ri)*sqrt((1.d0-ct)*(1.d0+ct))
         call multiple_origin_amplitude_matrix(amnp,s,phi,targetz,dir,sa)
         sa=sa*ri*ct/sqrt(csca/32.d0/pi)
         call amplitude_to_scattering_matrix(sa,sm)
         end subroutine multiple_origin_scatteringmatrix
!
!  scattering amplitude sa and matrix sm calculation
!
!  original: 15 January 2011
!  revised: 21 February 2011: S11 normalization changed
!  april 2013: moved things around to try to get it to work.
!
         subroutine scatteringmatrix(amn0,nodrt,ct,phi,sa,sm,rotate_plane,normalize_s11)
         implicit none
         logical, optional :: rotate_plane,normalize_s11
         logical :: rotate,norms11
         integer :: nodrt,m,n,p,m1,n1
         real(8) :: ct,phi,sm(4,4),cphi,sphi,qsca,tau(0:nodrt+1,nodrt,2)
         complex(8) :: amn0(0:nodrt+1,nodrt,2,2),sa(4),ephi,ephim(-nodrt:nodrt), &
                       ci,cin,a,b
         if(present(rotate_plane)) then
            rotate=rotate_plane
         else
            rotate=.false.
         endif
         if(present(normalize_s11)) then
            norms11=normalize_s11
         else
            norms11=.true.
         endif
         data ci/(0.d0,1.d0)/
         call taufunc(ct,nodrt,tau)
         cphi=cos(phi)
         sphi=sin(phi)
         ephi=dcmplx(cphi,sphi)
         call ephicoef(ephi,nodrt,ephim)
         sa=(0.d0,0.d0)
         qsca=0.d0
         do n=1,nodrt
            cin=(-ci)**n
            do m=-n,n
               if(m.le.-1) then
                  m1=n+1
                  n1=-m
               else
                  m1=m
                  n1=n
               endif
               do p=1,2
                  qsca=qsca+amn0(m1,n1,p,1)*dconjg(amn0(m1,n1,p,1)) &
                           + amn0(m1,n1,p,2)*dconjg(amn0(m1,n1,p,2))
                  if(rotate) then
                     a=amn0(m1,n1,p,1)*cphi+amn0(m1,n1,p,2)*sphi
                     b=amn0(m1,n1,p,1)*sphi-amn0(m1,n1,p,2)*cphi
                  else
                     a=amn0(m1,n1,p,1)
                     b=-amn0(m1,n1,p,2)
                  endif
                  sa(1)=sa(1)+cin*tau(m1,n1,3-p)*b*ephim(m)
                  sa(2)=sa(2)+ci*cin*tau(m1,n1,p)*a*ephim(m)
                  sa(3)=sa(3)+ci*cin*tau(m1,n1,p)*b*ephim(m)
                  sa(4)=sa(4)+cin*tau(m1,n1,3-p)*a*ephim(m)
               enddo
            enddo
         enddo
         qsca=qsca*2.d0
         if(.not.norms11) qsca=1.d0/pi
         sa=sa*4.d0/sqrt(qsca)
         call amplitude_to_scattering_matrix(sa,sm)
         end subroutine scatteringmatrix
!   c                                                                               c
!   c  subroutine scatexp(amn0,nodrt,nodrg,gmn) computes the expansion coefficients c
!   c  for the spherical harmonic expansion of the scattering phase function from   c
!   c  the scattering coefficients amn0.  For a complete expansion, the max. order  c
!   c  of the phase function expansion (nodrg) will be 2*nodrt, where nodrt is      c
!   c  the max. order of the scattered field expansion.   In this code nodrg is     c
!   c  typically set to 1, so that the subroutine returns the first moments         c
!   c  of the phase function; gmn(1) and gmn(2).                                    c
!   c                                                                               c
!   c  The expansion coefficients are normalized so that gmn(0)=1                   c
!   c                                                                               c
!   c  gmn(1)/3 is the asymmetry parameter.                                         c
!   c                                                                               c
         subroutine s11expansion(amn0,nodrt,mmax,nodrg,gmn)
         implicit none
         integer :: nodrt,m,n,p,ma,na,mmax,nodrg,w,w1,w2,u,uw,ww1, &
                    l1,l2,ka,la,k,l,q,ik
         real(8) :: vc1(0:nodrt*2+1),vc2(0:nodrt*2+1),g0
         complex(8) :: amn0(0:nodrt+1,nodrt,2,2),gmn(0:nodrg*(nodrg+3)/2), &
                       a(2,2),c,c2
         gmn=(0.d0,0.d0)
         do n=1,nodrt
            l1=max(1,n-nodrg)
            l2=min(nodrt,n+nodrg)
            do l=l1,l2
               c=sqrt(dble((n+n+1)*(l+l+1)))*dcmplx(0.d0,1.d0)**(l-n)
               w2=min(n+l,nodrg)
               call vcfunc(-1,l,1,n,vc2)
               do m=-n,n
                  if(m.le.-1) then
                     ma=n+1
                     na=-m
                  else
                     ma=m
                     na=n
                  endif
                  do k=-l,min(l,m)
                     if(k.le.-1) then
                        ka=l+1
                        la=-k
                     else
                        ka=k
                        la=l
                     endif
                     u=m-k
                     if(u.le.mmax) then
                        ik=(-1)**k
                        c2=ik*c
                        do p=1,2
                           do q=1,2
                              a(p,q)=c2*(amn0(ma,na,p,1)*conjg(amn0(ka,la,q,1)) &
                                    +amn0(ma,na,p,2)*conjg(amn0(ka,la,q,2)))
                           enddo
                        enddo
                        w1=max(abs(n-l),abs(u))
                        w2=min(n+l,nodrg)
                        call vcfunc(-k,l,m,n,vc1)
                        do w=w1,w2
                           uw=(w*(w+1))/2+u
                           do p=1,2
                              if(mod(n+l+w,2).eq.0) then
                                 q=p
                              else
                                 q=3-p
                              endif
                              gmn(uw)=gmn(uw)-vc1(w)*vc2(w)*a(p,q)
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo
         enddo
         g0=dble(gmn(0))
         gmn(0)=1.d0
         do w=1,nodrg
            ww1=(w*(w+1))/2
            gmn(ww1)=dcmplx(dble(gmn(ww1)),0.d0)/g0
            do u=1,min(mmax,w)
               uw=ww1+u
               gmn(uw)=(-1)**u*2.d0*gmn(uw)/g0
            enddo
         enddo
         end subroutine s11expansion
!
!  calculate azimuth--averaged scattering matrix from expansion, for cos(theta) = ct
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!  this is currently not used in v. 3.0
!
         subroutine fosmcalc(ntot,s00,s02,sp22,sm22,ct,sm,normalize_s11)
         logical :: norms11
         logical, optional :: normalize_s11
         integer :: ntot,w,ww1
         real(8) :: s00(4,4,0:ntot*2),s02(4,4,0:ntot*2),sp22(4,4,0:ntot*2),sm22(4,4,0:ntot*2), &
                    sm(4,4),dc(-2:2,0:2*ntot*(2*ntot+2)),ct,temp
         if(present(normalize_s11)) then
            norms11=normalize_s11
         else
            norms11=.true.
         endif
!if(light_up) then
!write(*,'('' foc1 '',3es13.5)') ct
!call flush(6)
!endif
         call rotcoef(ct,2,2*ntot,dc)

         sm=0.d0
         do w=0,2*ntot
            ww1=w*(w+1)
            sm(:,:)=sm(:,:)+s00(:,:,w)*dc(0,ww1)+s02(:,:,w)*dc(0,ww1+2) &
                   +sp22(:,:,w)*dc(2,ww1+2)+sm22(:,:,w)*dc(-2,ww1+2)
         enddo
!if(light_up) then
!write(*,'('' foc2 '',3es13.5)') ct,s00(1,1,0),sm(1,1)
!call flush(6)
!endif
         if(norms11) then
            if(abs(s00(1,1,0)).lt.1.d-10) then
               sm=0.d0
            else
               sm=sm/s00(1,1,0)
            endif
         endif
!
!  a patch
!
         sm(3,1)=-sm(3,1)
         sm(1,3)=-sm(1,3)
         sm(4,3)=-sm(4,3)
         temp=sm(4,1)
         sm(4,1)=sm(4,2)
         sm(4,2)=temp

         sm(1,2)=(sm(1,2)+sm(2,1))/2.d0
         sm(2,1)=sm(1,2)
         sm(3,4)=(sm(3,4)-sm(4,3))/2.d0
         sm(4,3)=-sm(3,4)

!         do i=1,4
!            do j=1,4
!               if(i.ne.1.or.j.ne.1) then
!                  sm(i,j)=sm(i,j)/sm(1,1)
!               endif
!            enddo
!         enddo
!if(light_up) then
!write(*,'('' foc3 '',i3)') mstm_global_rank
!call flush(6)
!endif
         end subroutine fosmcalc
!
!  determine the generalized spherical function expansion for the azimuth-averaged scattering matrix
!  corresponding to the target-based scattering field expansion of amnp.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: fixed flush call.
!

         subroutine fosmexpansion(ntot,amnp,s00,s02,sp22,sm22,mpi_comm)
         integer :: ntot,n,p,m,l,wmin,wmax,m1m,q,m1mq,m1mnpl,w,m1w,i,wtot,j, &
              rank,numprocs,nsend,runprintunit,mpicomm,task
         integer, optional :: mpi_comm
         real(8) :: s00(4,4,0:ntot*2),s02(4,4,0:ntot*2),sp22(4,4,0:ntot*2),sm22(4,4,0:ntot*2), &
                    cm1p1(0:ntot*2),cm1m1(0:ntot*2),cmmpm(0:ntot*2),cmmm2pm(0:ntot*2), &
                    cmmp2pm(0:ntot*2)
         complex(8) :: amnp(0:ntot+1,ntot,2,2),a1(-ntot-2:ntot+2,ntot,2),a2(-ntot-2:ntot+2,ntot,2), &
                       ci,fnl,a1122,a2112,a1p2,a1m2
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         data ci/(0.d0,1.d0)/
         call init(2*ntot)
         runprintunit=run_print_unit
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         a1=(0.d0,0.d0)
         a2=(0.d0,0.d0)
         do n=1,ntot
            do p=1,2
               do m=-n,-1
                  a1(m,n,p)=amnp(n+1,-m,p,1)
                  a2(m,n,p)=amnp(n+1,-m,p,2)
               enddo
               do m=0,n
                  a1(m,n,p)=amnp(m,n,p,1)
                  a2(m,n,p)=amnp(m,n,p,2)
               enddo
            enddo
         enddo
         s00=0.d0
         s02=0.d0
         sp22=0.d0
         sm22=0.d0
         wtot=ntot+ntot
         task=0
         do n=1,ntot
            do l=1,ntot
               wmin=abs(n-l)
               wmax=n+l
               fnl=sqrt(dble((n+n+1)*(l+l+1)))*ci**(l-n)
               call vcfunc(-1,n,1,l,cm1p1)
               call vcfunc(-1,n,-1,l,cm1m1)
               do m=-min(n,l+2),min(n,l+2)
                  m1m=(-1)**m
                  if(abs(m).le.l) then
                     call vcfunc(-m,n,m,l,cmmpm)
                  else
                     cmmpm=0.d0
                  endif
                  if(abs(-2+m).le.l) then
                     call vcfunc(-m,n,-2+m,l,cmmm2pm)
                  else
                     cmmm2pm=0.d0
                  endif
                  if(abs(2+m).le.l) then
                     call vcfunc(-m,n,2+m,l,cmmp2pm)
                  else
                     cmmp2pm=0.d0
                  endif
                  do p=1,2
                     do q=1,2
                        m1mq=(-1)**(m+q)
                        m1mnpl=(-1)**(m+n+p+l)
                        a1122=(a1(m,n,p)*conjg(a1(m,l,q)) + a2(m,n,p)*conjg(a2(m,l,q)))
                        a2112=(a2(m,n,p)*conjg(a1(m,l,q)) - a1(m,n,p)*conjg(a2(m,l,q)))
                        a1p2=(a1(m,n,p)+ci*a2(m,n,p))*conjg(a1(m-2,l,q)-ci*a2(m-2,l,q))
                        a1m2=(a1(m,n,p)-ci*a2(m,n,p))*conjg(a1(m+2,l,q)+ci*a2(m+2,l,q))
                        do w=wmin,wmax
                           task=task+1
                           if(mod(task-1,numprocs).ne.rank) cycle
                           m1w=(-1)**w
                           if(mod(n+l+w+p+q,2).eq.0) then
                              s00(1,1,w) = s00(1,1,w)-(m1m*fnl*a1122*cm1p1(w)*cmmpm(w))/2.
                           else
                              s00(4,4,w) = s00(4,4,w)+ (ci/2.*m1m*fnl*a2112*cm1p1(w)*cmmpm(w))
                           endif
                           if(w.lt.2) cycle
                           if(mod(n+l+w+p+q,2).eq.0) then
                              s02(2,1,w) = s02(2,1,w)-(m1mq*a1122*fnl*cm1m1(w)*cmmpm(w))/2.
                              s02(1,2,w) = s02(1,2,w)-(m1m*a1p2*fnl*cm1p1(w)*cmmm2pm(w))/4.
                              s02(1,2,w) = s02(1,2,w)-(m1m*a1m2*fnl*cm1p1(w)*cmmp2pm(w))/4.
                           else
                              s02(3,4,w) = s02(3,4,w)+ dimag(-ci/2.*m1mq*a2112*fnl*cm1m1(w)*cmmpm(w))
                              s02(4,3,w) = s02(4,3,w)+ dimag(m1m*a1p2*fnl*cm1p1(w)*cmmm2pm(w))/4.
                              s02(4,3,w) = s02(4,3,w)-dimag(m1m*a1m2*fnl*cm1p1(w)*cmmp2pm(w))/4.
                           endif
                           sm22(2,2,w) = sm22(2,2,w)-(m1mnpl*m1w*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                           sp22(2,2,w) = sp22(2,2,w)-(m1mq*a1p2*fnl*cm1m1(w)*cmmm2pm(w))/8.
                           sm22(2,2,w) = sm22(2,2,w)-(m1mq*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                           sp22(2,2,w) = sp22(2,2,w)-(m1mnpl*m1w*a1m2*fnl*cm1m1(w)*cmmp2pm(w))/8.
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
         call mstm_mpi(mpi_command='barrier')
         nsend=4*4*(2*ntot+1)
         call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=s00,&
              mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=s02,&
              mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=sp22,&
              mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=sm22,&
              mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_comm=mpicomm)
!
!  a patch
!
         do i=3,4
            do j=1,i
               s00(j,i,0:wtot)=-s00(j,i,0:wtot)
               s02(j,i,0:wtot)=-s02(j,i,0:wtot)
               sm22(j,i,0:wtot)=-sm22(j,i,0:wtot)
               sp22(j,i,0:wtot)=-sp22(j,i,0:wtot)
            enddo
         enddo
         sm22(3,3,:)=-sm22(2,2,:)
         sp22(3,3,:)=sp22(2,2,:)
!         deallocate(nlindex,nlnum)

         end subroutine fosmexpansion

!
!  compute the coefficients for the GSF expansion of the random orientation
!  scattering matrix.
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!  January 2012: added computation of coherent field average
!  February 2013:  added number processors option.
!
         subroutine ranorientscatmatrix(tmatrixfile,sm,smcf,beam_width,number_processors,mpi_comm)
         implicit none
         logical :: symmetrical
         integer :: nodr,nodrw,nodr2,m,n,p,k,l,q,t,v,u,w,nblk,kl,mn,nn1,tn, &
                    lmax,ll1,tvl,ku,k1,ns,ik,ik1,m1,nu,n1s,n1e,nu1,p1,n1max, &
                    in,n1,i,kt,nodrt,nodrrhs,mnm,klm,ikm, &
                    rank,numprocs, &
                    nblkrhs,numprocscalc,orig_group,new_group, &
                    new_comm,new_rank,nblkw,wv,sizedm,sizetm,nread,mpicomm
         integer, allocatable :: windex(:),vindex(:),wvindex(:),wvnum(:),group_list(:)
         integer, optional :: number_processors,mpi_comm
         real(8) :: sm(4,4,0:*),fl,xv,fl2,cbeam,gbn,wvperproc,sum, &
                    time1,time2,smcf(4,4,0:*),qextt,qscatt
         real(8), allocatable :: vc(:)
         real(8), optional :: beam_width
         complex(8) :: ci,cin,a,tct
         complex(8), allocatable :: aw(:,:,:),bw(:,:,:),cw(:), &
                       dw(:),pp(:,:,:),bm(:,:,:), &
                       am(:,:,:),fm(:,:,:,:,:),bmcf(:,:,:), &
                       amcf(:,:,:),fmcf(:,:,:,:,:),awcf(:,:,:), &
                       bwcf(:,:,:),cwcf(:),dwcf(:)
         complex(8), allocatable :: dm(:,:,:,:,:,:),dmcf(:,:,:,:,:,:)
         complex(4), allocatable :: tc(:,:,:,:)
         character*30 :: tmatrixfile
         data ci/(0.d0,1.d0)/
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         if(present(beam_width)) then
            cbeam=beam_width
         else
            cbeam=0.d0
         endif
         if(present(number_processors)) then
            numprocscalc=min(number_processors,numprocs)
         else
            numprocscalc=numprocs
         endif
         allocate(group_list(numprocscalc))
         group_list=(/(i,i=0,numprocscalc-1)/)
         call mstm_mpi(mpi_command='group',mpi_group=orig_group,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='incl',mpi_group=orig_group, &
            mpi_size=numprocscalc,mpi_new_group_list=group_list, &
            mpi_new_group=new_group,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='create',mpi_group=new_group, &
            mpi_new_comm=new_comm,mpi_comm=mpicomm)

         xv=cross_section_radius
         if(rank.le.numprocscalc-1) then

            call mstm_mpi(mpi_command='rank',mpi_comm=new_comm, &
               mpi_rank=new_rank)
            if(rank.eq.0) time1=mstm_mpi_wtime()
!
!  read the T matrix from the file
!
            do i=0,numprocscalc-1
               if(i.eq.new_rank) then
                  open(3,file=tmatrixfile)
                  read(3,*) nodr,nodrrhs
                  symmetrical=nodr.eq.nodrrhs
                  nodrt=nodr
                  nodr2=nodr+nodr
                  nodrw=nodr2
                  nblk=nodr*(nodr+2)
!                  nodrrhs=nodr
                  nblkrhs=nodrrhs*(nodrrhs+2)
                  sizetm=4*nblk*nblkrhs
                  allocate(tc(2,nblk,2,nblkrhs))
                  tc=(0.,0.)
                  do l=1,nodrrhs
                     do k=-l,l
                        kl=l*(l+1)+k
                        klm=l*(l+1)-k
                        do q=1,2
                           if(symmetrical) then
                              nread=l
                           else
                              nread=nodr
                           endif
                           do n=1,nread
                              do m=-n,n
                                 mn=n*(n+1)+m
                                 mnm=n*(n+1)-m
                                 do p=1,2
                                    read(3,*) tc(p,mn,q,kl)
                                 enddo
                                 if(n.lt.l.and.symmetrical) then
                                    ikm=(-1)**(m+k)
                                    do p=1,2
                                       tc(q,klm,p,mnm)=tc(p,mn,q,kl)*ikm
                                    enddo
                                 endif
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
                  qextt=0.d0
                  qscatt=0.d0
                  do l=1,nodrrhs
                     do k=-l,l
                        kl=l*(l+1)+k
                        do q=1,2
                           do n=1,nodr
                              do m=-n,n
                                 mn=n*(n+1)+m
                                 do p=1,2
                                    qscatt=qscatt+cabs(tc(p,mn,q,kl))**2
                                 enddo
                              enddo
                           enddo
                           if(l.le.nodr) qextt=qextt-real(tc(q,kl,q,kl))
                        enddo
                     enddo
                  enddo
                  qextt=qextt*2./cross_section_radius**2
                  qscatt=qscatt*2./cross_section_radius**2
                  close(3)
               endif
               call mstm_mpi(mpi_command='barrier',mpi_comm=new_comm)
            enddo
            if(new_rank.eq.0) then
               write(run_print_unit,'('' t matrix ext, sca:'',2e13.5)') qextt,qscatt
            endif
!!
!!  send to the other processors
!!
!            if(numprocscalc.gt.1) then
!               call mstm_mpi(mpi_command='bcast',mpi_send_buf_c=tc, &
!                    mpi_number=sizetm,mpi_rank=0,mpi_comm=new_comm)
!            endif

            allocate(vc(0:4*nodr+2))
            allocate(aw(0:2,-1:1,0:nodrw),bw(0:2,-1:1,0:nodrw),cw(0:nodrw), &
                 dw(0:nodrw),pp(nodr,2,2),bm(2,nodr*(nodr+2),2), &
                 am(2,nodrrhs+1,2),fm(3,nodr,2,nodr,2),bmcf(2,nodr*(nodr+2),2), &
                 amcf(2,nodrrhs+1,2),fmcf(3,nodr,2,nodr,2),awcf(0:2,-1:1,0:nodrw), &
                 bwcf(0:2,-1:1,0:nodrw),cwcf(0:nodrw),dwcf(0:nodrw))
            allocate(dm(-nodr-1:nodr+1,3,nodr,2,nodr,2),dmcf(-nodr-1:nodr+1,3,nodr,2,nodr,2))
            if(new_rank.eq.0) then
               time2=mstm_mpi_wtime()-time1
               call timewrite(run_print_unit,' t matrix read time:',time2)
               time1=mstm_mpi_wtime()
            endif
            dm=(0.d0,0.d0)
            dmcf=(0.d0,0.d0)
            sizedm=size(dm)
            call init(nodr2)
!
!  compute the GB modified T matrix
!
            do n=1,nodrrhs
               gbn=dexp(-((dble(n)+.5d0)*cbeam)**2.)
               cin=ci**(n+1)
               pp(n,1,1) =-.5d0*cin*fnr(n+n+1)*gbn
               pp(n,2,1) =-pp(n,1,1)
               pp(n,1,2)=-pp(n,1,1)
               pp(n,2,2)=pp(n,2,1)
            enddo
            do n=1,nodr
               nn1=n*(n+1)
               do m=-n,n
                  mn=nn1+m
                  do p=1,2
                     do l=1,nodrrhs
                        do k=-l,l
                           kl=l*(l+1)+k
                           a=tc(p,mn,1,kl)
                           tc(p,mn,1,kl)=tc(p,mn,1,kl)*pp(l,1,1)&
                              +tc(p,mn,2,kl)*pp(l,1,2)
                           tc(p,mn,2,kl)=a*pp(l,2,1)+tc(p,mn,2,kl)*pp(l,2,2)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
!
!  determine the distribution of work load among the processors
!
            nblkw=nodr2*(nodr2+2)+1
            allocate(windex(nblkw),vindex(nblkw),wvindex(0:numprocscalc-1), &
               wvnum(0:numprocscalc-1))
            w=0
            do n=0,nodr2
               do m=-n,n
                  w=w+1
                  windex(w)=n
                  vindex(w)=m
               enddo
            enddo
            wvperproc=dble(nblkw)/dble(numprocscalc)
            sum=0.
            do i=0,numprocscalc-1
               wvindex(i)=floor(sum)
               sum=sum+wvperproc
            enddo
            do i=0,numprocscalc-2
               wvnum(i)=wvindex(i+1)-wvindex(i)
            enddo
            wvnum(numprocscalc-1)=nblkw-wvindex(numprocscalc-1)
            if(new_rank.eq.0) then
               write(run_print_unit,'('' d matrix calculation, order+degree per proc.:'',f9.2)') &
                   wvperproc
               call flush(run_print_unit)
            endif
!
!  the big loop
!
            do wv=wvindex(rank)+1,wvindex(rank)+wvnum(rank)
               w=windex(wv)
               v=vindex(wv)
               bm=(0.d0,0.d0)
               bmcf=(0.d0,0.d0)
               do n=1,nodr
                  nn1=n*(n+1)
                  do l=max(1,abs(w-n)),min(nodrrhs,w+n)
                     am(1,l,1)=0.d0
                     am(1,l,2)=0.d0
                     am(2,l,1)=0.d0
                     am(2,l,2)=0.d0
                     amcf(1,l,1)=0.d0
                     amcf(1,l,2)=0.d0
                     amcf(2,l,1)=0.d0
                     amcf(2,l,2)=0.d0
                  enddo
                  do t=-n,n
                     tn=nn1+t
                     lmax=min(nodrrhs,w+n)
                     call vcfunc(v,w,-t,n,vc)
                     do l=max(1,abs(v-t),abs(n-w)),lmax
                        ll1=l*(l+1)
                        tvl=ll1+t-v
                        do k=1,2
                           do p=1,2
                              am(k,l,p)=am(k,l,p)+vc(l)*tc(p,tn,k,tvl)
                              if(l.eq.n.and.v.eq.0) then
! may 2019: orientation averaged T matrix is azimuthal independent
                                 tct=0.d0
                                 do kt=-n,n
                                    tct=tct+tc(p,nn1+kt,k,nn1+kt)
                                 enddo
                                 amcf(k,l,p)=amcf(k,l,p)+vc(l)*tct/dble(n+n+1)
!                                 amcf(k,l,p)=amcf(k,l,p)+vc(l)*tcm(p,n)
!                                 amcf(k,l,p)=amcf(k,l,p)+vc(l)*tc(p,tn,k,tvl)
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
                  do m=-n,n
                     mn=nn1+m
                     do k=1,2
                        u=m-(-3+2*k)
                        if(abs(u).le.w) then
                           lmax=min(nodrrhs,w+n)
                           call vcfunc(-u,w,m,n,vc)
                           do l=max(1,abs(w-n)),lmax
                              fl=-(-1)**l*vc(l)/dble(l+l+1)
                              do p=1,2
                                 bm(k,mn,p)=bm(k,mn,p)+am(k,l,p)*fl
                                 if(v.eq.0) then
                                    bmcf(k,mn,p)=bmcf(k,mn,p)+amcf(k,l,p)*fl
                                 endif
                              enddo
                           enddo
                        endif
                     enddo
                  enddo
               enddo

!               do k=1,2
!                  do n=1,nodr*(nodr+2)
!                     do p=1,2
!                        bmcf(k,n,p)=bm(k,n,p)-bmcf(k,n,p)
!                     enddo
!                  enddo
!               enddo

               do u=-min(w,nodr+1),min(w,nodr+1)
                  do ku=1,3
                     if(ku.eq.1) then
                        k=-1
                        k1=-1
                     elseif(ku.eq.2) then
                        k=1
                        k1=1
                     else
                        k=1
                        k1=-1
                     endif
                     m=u+k
                     ns=max(1,abs(m))
                     ik=(k+1)/2+1
                     ik1=(k1+1)/2+1
                     m1=u+k1
                     do n=ns,nodr
                        nu=n*(n+1)+m
                        n1s=max(1,abs(m1),n-nodrw)
                        n1e=min(nodr,n+nodrw)
                        do n1=n1s,n1e
                           cin=ci**(n-n1)
                           nu1=n1*(n1+1)+m1
                           fl=-fnr(n+n+1)*fnr(n1+n1+1)*dble(w+w+1)
                           do p=1,2
                              do p1=1,2
                                 a=bm(ik,nu,p)*cin*fl*conjg(bm(ik1,nu1,p1))
                                 dm(u,ku,n,p,n1,p1)=dm(u,ku,n,p,n1,p1)+a
                                 if(v.eq.0) then
                                    a=bmcf(ik,nu,p)*cin*fl*conjg(bmcf(ik1,nu1,p1))
                                    dmcf(u,ku,n,p,n1,p1)=dmcf(u,ku,n,p,n1,p1)+a
                                 endif
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            deallocate(tc)

            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=dm,&
                 mpi_number=sizedm,mpi_operation=mstm_mpi_sum,mpi_comm=new_comm)
            call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dc=dmcf,&
                 mpi_number=sizedm,mpi_operation=mstm_mpi_sum,mpi_comm=new_comm)
            if(new_rank.eq.0) then
               time2=mstm_mpi_wtime()-time1
               call timewrite(run_print_unit,' d matrix time:',time2)
               time1=mstm_mpi_wtime()
            endif
!
!  compute the expansion coefficients
!
            aw=0.d0
            bw=0.d0
            cw=0.d0
            dw=0.d0
            awcf=0.d0
            bwcf=0.d0
            cwcf=0.d0
            dwcf=0.d0
            do w=0,nodrw
               do n=1,nodr
                  n1s=max(1,abs(n-w))
                  n1e=min(nodr,n+w)
                  do n1=n1s,n1e
                     do k=1,3
                        do p=1,2
                           do p1=1,2
                              fm(k,n,p,n1,p1)=0.
                              fmcf(k,n,p,n1,p1)=0.
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               do u=-nodr-1,nodr+1
                  do k=-1,1,2
                     m=u+k
                     ik=(k+1)/2+1
                     ns=max(1,abs(m))
                     do n=ns,nodr
                        n1max=min(w+n,nodr)
                        call vcfunc(m,n,0,w,vc)
                        do n1=ns,nodr
                           if((n+n1.lt.w).or.(abs(n-n1).gt.w)) cycle
                           fl=-(-1)**n*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                           do p=1,2
                              do p1=1,2
                                 fm(ik,n,p,n1,p1)=fm(ik,n,p,n1,p1) &
                                    +dm(u,ik,n,p,n1,p1)*fl
                                 fmcf(ik,n,p,n1,p1)=fmcf(ik,n,p,n1,p1) &
                                    +dmcf(u,ik,n,p,n1,p1)*fl
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
                  if(w.lt.2) cycle
                  m=u+1
                  m1=u-1
                  ns=max(1,abs(m))
                  n1s=max(1,abs(m1))
                  do n=ns,nodr
                     n1max=min(w+n,nodr)
                     call vcfunc(m,n,-2,w,vc)
                     do n1=n1s,nodr
                        if((n+n1.lt.w).or.(abs(n-n1).gt.w)) cycle
                        fl=-(-1)**n*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                        do p=1,2
                           do p1=1,2
                              fm(3,n,p,n1,p1)=fm(3,n,p,n1,p1) &
                                 +dm(u,3,n,p,n1,p1)*fl
                              fmcf(3,n,p,n1,p1)=fmcf(3,n,p,n1,p1) &
                                 +dmcf(u,3,n,p,n1,p1)*fl
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               do n=1,nodr
                  n1s=max(1,abs(n-w))
                  n1e=min(nodr,n+w)
                  in=(-1)**n
                  n1max=min(w+n,nodr)
                  call vcfunc(1,n,0,w,vc)
                  do n1=n1s,n1e
                     fl=2.d0*in*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                     i=mod(n+n1-w,2)+1
                     do p=1,2
                        p1=(2-i)*p+(i-1)*(3-p)
                        do k=-1,1,2
                           ik=(k+1)/2+1
                           aw(0,k,w)=aw(0,k,w)+fm(ik,n,p,n1,p1)*fl
                           bw(0,k,w)=bw(0,k,w)+fm(ik,n,p,n1,3-p1)*fl
                           awcf(0,k,w)=awcf(0,k,w)+fmcf(ik,n,p,n1,p1)*fl
                           bwcf(0,k,w)=bwcf(0,k,w)+fmcf(ik,n,p,n1,3-p1)*fl
                        enddo
                        bw(2,0,w)=bw(2,0,w)+fm(3,n,p,n1,3-p1)*fl
                        aw(2,0,w)=aw(2,0,w)+fm(3,n,p,n1,p1)*fl
                        bwcf(2,0,w)=bwcf(2,0,w)+fmcf(3,n,p,n1,3-p1)*fl
                        awcf(2,0,w)=awcf(2,0,w)+fmcf(3,n,p,n1,p1)*fl
                     enddo
                  enddo
                  if(w.lt.2) cycle
                  call vcfunc(1,n,-2,w,vc)
                  do n1=n1s,n1e
                     fl=2.d0*in*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                     i=mod(n+n1-w,2)+1
                     do p=1,2
                        p1=(2-i)*p+(i-1)*(3-p)
                        do k=-1,1,2
                           ik=(k+1)/2+1
                           aw(2,k,w)=aw(2,k,w)+fm(ik,n,p,n1,p1)*fl*(-1)**p1
                           bw(2,k,w)=bw(2,k,w)+fm(ik,n,p,n1,3-p1)*fl*(-1)**(3-p1)
                           awcf(2,k,w)=awcf(2,k,w)+fmcf(ik,n,p,n1,p1)*fl*(-1)**p1
                           bwcf(2,k,w)=bwcf(2,k,w) &
                              +fmcf(ik,n,p,n1,3-p1)*fl*(-1)**(3-p1)
                        enddo
                     enddo
                     fl2=2.*(-1)**(n1+w)*vc(n1)*fnr(w+w+1)/fnr(n1+n1+1)
                     do p=1,2
                        do p1=1,2
                           cw(w)=cw(w)+fm(3,n,p,n1,p1)*fl*(-1)**p1
                           dw(w)=dw(w)+fm(3,n,p,n1,p1)*fl2*(-1)**p
                           cwcf(w)=cwcf(w)+fmcf(3,n,p,n1,p1)*fl*(-1)**p1
                           dwcf(w)=dwcf(w)+fmcf(3,n,p,n1,p1)*fl2*(-1)**p
                        enddo
                     enddo
                  enddo
               enddo
            enddo
!            do w=0,nodrw
!               do k=-1,1
!                  do i=0,2
!                     aw(i,k,w)=aw(i,k,w)*2./xv/xv
!                     bw(i,k,w)=bw(i,k,w)*2./xv/xv
!                     awcf(i,k,w)=awcf(i,k,w)*2./xv/xv
!                     bwcf(i,k,w)=bwcf(i,k,w)*2./xv/xv
!                  enddo
!               enddo
!               cw(w)=cw(w)*2./xv/xv
!               dw(w)=dw(w)*2./xv/xv
!               cwcf(w)=cwcf(w)*2./xv/xv
!               dwcf(w)=dwcf(w)*2./xv/xv
!            enddo
            do n=0,nodrw
               sm(1,1,n)=aw(0,-1,n)+aw(0,1,n)
               sm(1,2,n)=aw(2,-1,n)+aw(2,1,n)
               sm(1,3,n)=2.d0*dimag(aw(2,0,n))
               sm(1,4,n)=aw(0,1,n)-aw(0,-1,n)
               sm(2,2,n)=dw(n)
               sm(2,3,n)=dimag(dw(n))
               sm(2,4,n)=aw(2,1,n)-aw(2,-1,n)
               sm(3,1,n)=-dimag(bw(2,-1,n)+bw(2,1,n))
               sm(3,2,n)=dimag(cw(n))
               sm(3,3,n)=cw(n)
               sm(3,4,n)=dimag(bw(2,-1,n)-bw(2,1,n))
               sm(4,1,n)=bw(0,-1,n)+bw(0,1,n)
               sm(4,2,n)=2.*bw(2,0,n)
               sm(4,4,n)=bw(0,1,n)-bw(0,-1,n)

               smcf(1,1,n)=awcf(0,-1,n)+awcf(0,1,n)
               smcf(1,2,n)=awcf(2,-1,n)+awcf(2,1,n)
               smcf(1,3,n)=2.d0*dimag(awcf(2,0,n))
               smcf(1,4,n)=awcf(0,1,n)-awcf(0,-1,n)
               smcf(2,2,n)=dwcf(n)
               smcf(2,3,n)=dimag(dwcf(n))
               smcf(2,4,n)=awcf(2,1,n)-awcf(2,-1,n)
               smcf(3,1,n)=-dimag(bwcf(2,-1,n)+bwcf(2,1,n))
               smcf(3,2,n)=dimag(cwcf(n))
               smcf(3,3,n)=cwcf(n)
               smcf(3,4,n)=dimag(bwcf(2,-1,n)-bwcf(2,1,n))
               smcf(4,1,n)=bwcf(0,-1,n)+bwcf(0,1,n)
               smcf(4,2,n)=2.*bwcf(2,0,n)
               smcf(4,4,n)=bwcf(0,1,n)-bwcf(0,-1,n)

            enddo
!
!  normalization
!
!            qsca0=sm(1,1,0)
!            do n=0,nodrw
!               do i=1,4
!                  do j=1,4
!                     sm(i,j,n)=sm(i,j,n)/qsca0
!                     smcf(i,j,n)=smcf(i,j,n)/qsca0
!                  enddo
!               enddo
!            enddo
            if(new_rank.eq.0) then
               time2=mstm_mpi_wtime()-time1
               call timewrite(run_print_unit,' scat matrix coef time:',time2)
            endif
            deallocate(windex,vindex,wvindex,wvnum,dm)
            deallocate(vc)
            deallocate(aw,bw,cw,dw,pp,bm,am,fm,bmcf,amcf,fmcf,awcf, &
                 bwcf,cwcf,dwcf)
         endif
         call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
         end subroutine ranorientscatmatrix
!
!  calculation of the RO scattering matrix from the GSF expansion
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: changed normalization on S11
!
         subroutine ranorienscatmatrixcalc(ct,smc,nodrexp,sm)
         implicit none
         integer :: nodrexp,n,nn0,nnp2,nnm2
         real(8) :: smc(4,4,0:nodrexp),sm(4,4),dc(-2:2,0:nodrexp*(nodrexp+2)), &
                    ct
!
!     dc is the normalized generalized spherical function
!     dc(k,n*(n+1)+m) = ((n-k)!(n+m)!/(n+k)!/(n-m)!)^(1/2) D^k_{mn},
!     where D^k_{mn} is defined in M&M JOSA 96
!
         call rotcoef(ct,2,nodrexp,dc)
         sm=0.d0
         do n=0,nodrexp
            nn0=n*(n+1)
            nnp2=nn0+2
            nnm2=nn0-2
            sm(1,1)=sm(1,1)+dc(0,nn0)*smc(1,1,n)
            sm(1,4)=sm(1,4)+dc(0,nn0)*smc(1,4,n)
            sm(4,1)=sm(4,1)+dc(0,nn0)*smc(4,1,n)
            sm(4,4)=sm(4,4)+dc(0,nn0)*smc(4,4,n)
            if(n.ge.2) then
               sm(1,2)=sm(1,2)+dc(2,nn0)*smc(1,2,n)
               sm(2,4)=sm(2,4)+dc(2,nn0)*smc(2,4,n)
               sm(3,4)=sm(3,4)+dc(2,nn0)*smc(3,4,n)
               sm(1,3)=sm(1,3)+dc(2,nn0)*smc(1,3,n)
               sm(3,1)=sm(1,3)+dc(2,nn0)*smc(3,1,n)
               sm(4,2)=sm(4,2)+dc(2,nn0)*smc(4,2,n)
               sm(2,2)=sm(2,2)+dc(2,nnm2)*smc(2,2,n)+dc(2,nnp2)*smc(3,3,n)
               sm(2,3)=sm(2,3)+dc(2,nnp2)*smc(2,3,n)+dc(2,nnp2)*smc(3,2,n)
               sm(3,3)=sm(3,3)-dc(2,nnm2)*smc(2,2,n)+dc(2,nnp2)*smc(3,3,n)
               sm(3,2)=sm(3,2)+dc(2,nnp2)*smc(2,3,n)-dc(2,nnp2)*smc(3,2,n)
            endif
         enddo
         sm(2,1)=sm(1,2)
         sm(4,3)=-sm(3,4)
!
!  discontiued scaling option: now done in main program
!
!            if(iscale.eq.1) then
!               do j=1,4
!                  do k=j,4
!                     if(j.ne.1.or.k.ne.1) then
!                        sm(j,k,i)=sm(j,k,i)/sm(1,1,i)
!                     endif
!                  enddo
!               enddo
!            endif
!
!    here are the VV and HH differential cross sections
!
!            gvv=.25*(sm(1,1)+sm(2,2)-2.*sm(1,2))
!            ghh=.25*(sm(1,1)+sm(2,2)+2.*sm(1,2))
!
         return
         end subroutine ranorienscatmatrixcalc

      end module scatprops
!
! module nearfield contains local data and subroutines for near field calculation
!
      module nearfield
      use mpidefs
      use intrinsics
      use specialfuncs
      use spheredata
      use mie
      use numconstants
      use surface_subroutines
      use periodic_lattice_subroutines
      use scatprops
      implicit none
      type grid_info
         logical :: initialized
         integer :: cellnum
         type(cell_info), pointer :: cellinfo=>null()
      end type grid_info
      type cell_info
         logical :: outside_spheres
         integer :: ncell(3),layer,host,order,nispheres
         real(8) :: rcell(3)
         type(linked_sphere_list), pointer :: sphere_list
         complex(8), pointer :: vector(:,:,:)=>null()
         complex(8), pointer :: reg_source_vector(:,:,:)=>null()
         complex(8), pointer :: gb_vector(:,:,:)=>null()
      end type cell_info
      type linked_cell_list
         type(cell_info) :: cellinfo
         type(linked_cell_list), pointer :: next=>null()
      end type linked_cell_list
      type vector_storage
         complex(8), pointer :: vector(:,:,:)=>null()
      end type vector_storage
      type linked_sphere_data
         real(8) :: position(3), radius
         type(linked_sphere_data), pointer :: next
      end type linked_sphere_data

      type(linked_cell_list), pointer, private :: cell_info_list
      type(vector_storage), allocatable, private :: internal_field_vector(:)
      logical :: incident_gb
      logical, target :: store_surface_vector
      integer, private :: local_rank,local_numprocs,local_run_number,total_cells,number_intersecting_spheres
      integer, target :: near_field_expansion_order
      real(8), private :: grid_region(3,2),grid_spacing(3)
      real(8), target :: near_field_expansion_spacing
      type(linked_sphere_data), pointer, private :: intersecting_spheres
      complex(8), private :: vwf_0(3,3,2)
      data near_field_expansion_order,near_field_expansion_spacing,store_surface_vector/10,5.d0,.true./
      data local_run_number/1/

      contains

         subroutine near_field_calculation(amnp,alpha,sinc,dir,gridregion,griddim,incident_model,output_unit, &
            e_field_array,h_field_array,e_field_ave_array,output_header,mpi_comm)
         implicit none
         logical, optional :: output_header
         integer :: incmodel,i,p,outputunit,nblk,l,griddim(3),ix,iy,iz, &
            layer,host,ipos(3),cellnum,nsend,mpicomm,totpoints,point,dir
         integer, optional :: incident_model,output_unit,mpi_comm
         real(8) :: gridregion(3,2),rpos(3),alpha,time1,time0,rtemp,sinc
         complex(8) :: amnp(number_eqns,2),evec(3,2),hvec(3,2),evec1(3,2),hvec1(3,2)
         complex(8), allocatable :: vector(:,:,:)
         complex(8), optional :: e_field_ave_array(3,2,griddim(3))
         complex(8), target, optional :: e_field_array(3,2,griddim(1),griddim(2),griddim(3)), &
            h_field_array(3,2,griddim(1),griddim(2),griddim(3))
         complex(4) :: earray(3,2,griddim(1),griddim(2)), &
            harray(3,2,griddim(1),griddim(2))
         type(grid_info), allocatable :: gridinfo(:,:,:)
         type(cell_info), pointer :: cellinfo
         type(linked_cell_list), pointer :: clist1,clist2
!         complex(4), pointer :: earray(:,:,:,:,:),harray(:,:,:,:,:)

         if(present(incident_model)) then
            incmodel=incident_model
         else
            incmodel=1
         endif
         if(present(output_unit)) then
            outputunit=output_unit
         else
            outputunit=0
         endif
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         incident_gb=gaussian_beam_constant.ne.0.d0

         call mstm_mpi(mpi_command='size',mpi_size=local_numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=local_rank,mpi_comm=mpicomm)
         do i=1,3
            if(griddim(i).eq.1) then
               grid_spacing(i)=0.d0
               rtemp=0.5d0*(gridregion(i,1)+gridregion(i,2))
               grid_region(i,:)=rtemp
            else
               grid_spacing(i)=(gridregion(i,2)-gridregion(i,1))/dble(griddim(i))
               grid_region(i,:)=gridregion(i,:)
            endif
         enddo
         totpoints=product(griddim)

         call vwhcalc((/0.d0,0.d0,0.d0/),(/(1.d0,0.d0),(1.d0,0.d0)/),1,1,vwf_0,index_model=2)

         allocate(gridinfo(griddim(1),griddim(2),griddim(3)))
         call grid_point_initialize(griddim,gridinfo)

         if(allocated(internal_field_vector)) then
            l=ubound(internal_field_vector,1)
            do i=1,l
               if(associated(internal_field_vector(i)%vector)) deallocate(internal_field_vector(i)%vector)
            enddo
            deallocate(internal_field_vector)
         endif
         allocate(internal_field_vector(number_spheres))
         do i=1,number_spheres
            nblk=sphere_order(i)*(sphere_order(i)+2)
            allocate(vector(nblk,2,2))
            allocate(internal_field_vector(i)%vector(nblk,2,2))
            do p=1,2
               if(number_field_expansions(i).eq.1) then
                  call onemiecoeffmult(i,sphere_order(i),amnp(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),p), &
                     vector(1:nblk,1:2,p),'c')
               else
                  vector(1:nblk,1:2,p) &
                     =reshape(amnp(sphere_offset(i)+sphere_block(i)+1:sphere_offset(i)+2*sphere_block(i),p), &
                     (/nblk,2/))
               endif
            enddo
            internal_field_vector(i)%vector=vector
            deallocate(vector)
         enddo

         if(present(output_header)) then
            if(output_header.and.outputunit.ne.0) call write_output_header(griddim,outputunit)
         endif
         if(local_rank.eq.0.and.outputunit.ne.0) then
            write(run_print_unit,'('' near field calculation'')')
            write(run_print_unit,'('' grid minimum x,y,z:'',3es12.4)') grid_region(:,1)
            write(run_print_unit,'('' grid maximum x,y,z:'',3es12.4)') grid_region(:,2)
            write(run_print_unit,'('' grid x,y,z dimensions, total points:'',3i5,i12)') griddim(:),product(griddim)
            write(run_print_unit,'('' total reexpansion cells:'',i5)') total_cells
            call flush(run_print_unit)
         endif

!clist=>cell_info_list
!open(20,file='ctest.dat')
!do i=1,total_cells
!write(20,'('' cell:'',i5)') i
!write(20,'('' ncell,host,layer,order:'',6i5)') clist%ncell,clist%host,clist%layer,clist%order
!write(20,'('' rcell:'',3es12.4)') clist%rcell
!write(20,'('' vector stored:'',l2)') associated(clist%vector)
!write(20,*)
!clist=>clist%next
!if(.not.associated(clist)) exit
!enddo

         time0=mstm_mpi_wtime()
         point=0
         if(present(e_field_ave_array)) e_field_ave_array=0.d0
         do iz=1,griddim(3)
            earray=0.
            harray=0.
            do iy=1,griddim(2)
               do ix=1,griddim(1)
                  cellnum=gridinfo(ix,iy,iz)%cellnum
                  point=point+1
                  time1=mstm_mpi_wtime()
!                  if(time1-time0.gt.15.d0.and.local_rank.eq.0.and.(.not.present(e_field_ave_array))) then
                  if(time1-time0.gt.15.d0.and.mstm_global_rank.eq.0) then
                     write(run_print_unit,'('' completed field point calculation '',i8,''/'',i8)') point,totpoints
                     call flush(run_print_unit)
                     time0=time1
                  endif
                  if(mod(cellnum,local_numprocs).ne.local_rank) cycle
!
                  cellinfo=>gridinfo(ix,iy,iz)%cellinfo
                  layer=cellinfo%layer
                  host=cellinfo%host
                  ipos(:)=(/ix,iy,iz/)
                  rpos(:)=(dble(ipos(:))-(/0.5d0,0.5d0,0.5d0/))*grid_spacing(:)+grid_region(:,1)
!write(*,*) 'nf 1', layer,host
!call flush(6)
!if(local_rank.eq.0) then
!write(*,'(3i5,3es12.4)') ix,iz,host,rpos
!call flush(6)
!endif
                  if(host.eq.0.and.(.not.periodic_lattice).and.(number_plane_boundaries.eq.0)) then
                     call source_field_calculate_fast(rpos,amnp,cellinfo,evec,hvec)
                  else
                     call source_field_calculate(rpos,amnp,host,layer,evec,hvec)
                  endif
!write(*,*) 'nf 2'
!call flush(6)
                  if(host.eq.0.and.(number_plane_boundaries.gt.0.or.periodic_lattice)) then
                     call surface_field_calculate(rpos,amnp,cellinfo,evec1,hvec1)
                     evec(:,:)=evec(:,:)+evec1(:,:)
                     hvec(:,:)=hvec(:,:)+hvec1(:,:)
                  endif
!write(*,*) 'nf 3'
!call flush(6)
                  if(incmodel.ne.2.and.host.eq.0) then
                     call incident_field_calculate(rpos,layer,alpha,sinc,dir,cellinfo,evec1,hvec1)
                     evec(:,:)=evec(:,:)+evec1(:,:)
                     hvec(:,:)=hvec(:,:)+hvec1(:,:)
                  elseif(incmodel.eq.2.and.host.ne.0) then
                     call incident_field_calculate(rpos,layer,alpha,sinc,dir,cellinfo,evec1,hvec1)
                     evec(:,:)=evec(:,:)-evec1(:,:)
                     hvec(:,:)=hvec(:,:)-hvec1(:,:)
                  endif
                  earray(:,:,ix,iy)=evec
                  harray(:,:,ix,iy)=hvec
                  if(present(e_field_ave_array)) e_field_ave_array(:,:,iz)=e_field_ave_array(:,:,iz)+evec
               enddo
            enddo
            if(present(e_field_ave_array)) then
               e_field_ave_array(:,:,iz)=e_field_ave_array(:,:,iz)/dble(griddim(1)*griddim(2))
            endif

            if(local_numprocs.gt.1.and.outputunit.ne.0) then
               nsend=3*2*product(griddim(1:2))
               call mstm_mpi(mpi_command='reduce',mpi_recv_buf_c=earray, &
                    mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_rank=0,mpi_comm=mpicomm)
               call mstm_mpi(mpi_command='reduce',mpi_recv_buf_c=harray, &
                    mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_rank=0,mpi_comm=mpicomm)
            endif

            if(outputunit.ne.0.and.local_rank.eq.0) then
               do iy=1,griddim(2)
                  do ix=1,griddim(1)
                     ipos(:)=(/ix,iy,iz/)
                     rpos(:)=(dble(ipos(:))-(/0.5d0,0.5d0,0.5d0/))*grid_spacing(:)+grid_region(:,1)
                     write(outputunit,'(27es12.4)') rpos(:),earray(:,1,ix,iy),harray(:,1,ix,iy), &
                        earray(:,2,ix,iy),harray(:,2,ix,iy)
                  enddo
               enddo
            endif
            if(present(e_field_array)) then
               e_field_array(:,:,:,:,iz)=earray
            endif
            if(present(h_field_array)) then
               h_field_array(:,:,:,:,iz)=harray
            endif

         enddo

         if(present(e_field_ave_array)) then
            nsend=3*2*griddim(3)
            call mstm_mpi(mpi_command='reduce',mpi_recv_buf_dc=e_field_ave_array, &
                 mpi_number=nsend,mpi_operation=mstm_mpi_sum,mpi_rank=0,mpi_comm=mpicomm)
         endif


         clist1=>cell_info_list
         do while(associated(clist1))
            clist2=>clist1%next
            deallocate(clist1)
            nullify(clist1)
            clist1=>clist2
         enddo
         if(associated(cell_info_list)) then
            nullify(cell_info_list)
         endif
         deallocate(gridinfo)

         end subroutine near_field_calculation

         subroutine source_field_calculate_fast(rpos,sourcevec,cellinfo,evec,hvec)
         implicit none
         integer :: nodr,nblk,i,p,j
         real(8) :: rpos(3),rtran(3)
         complex(8) :: sourcevec(number_eqns,2),evec(3,2),hvec(3,2),ri2(2)
         complex(8), allocatable :: vwf(:,:,:),svec(:,:)
         type(linked_sphere_list), pointer :: slist
         type(cell_info), pointer :: cellinfo

         evec=0.d0
         hvec=0.d0
         if(.not.associated(cellinfo%reg_source_vector)) then
            call stored_source_vector_calculate(sourcevec,cellinfo)
         endif
         ri2(1:2)=layer_ref_index(cellinfo%layer)
         if(cellinfo%nispheres.gt.0) then
            slist=>cellinfo%sphere_list
            do i=1,cellinfo%nispheres
               j=slist%sphere
               if(i.lt.cellinfo%nispheres) slist=>slist%next
               nodr=sphere_order(j)
               nblk=nodr*(nodr+2)
               allocate(vwf(3,nblk,2),svec(nblk,2))
               rtran(:)=rpos(:)-sphere_position(:,j)
               call vwhcalc(rtran,ri2,nodr,3,vwf,index_model=2)
               do p=1,2
                  svec=reshape(sourcevec(sphere_offset(j)+1:sphere_offset(j)+sphere_block(j),p), &
                     (/nblk,2/))
                  evec(:,p)=evec(:,p)+(matmul(vwf(:,:,1),svec(:,1))+matmul(vwf(:,:,2),svec(:,2)))
                  hvec(:,p)=hvec(:,p)+(matmul(vwf(:,:,1),svec(:,1))-matmul(vwf(:,:,2),svec(:,2)))*ri2(1)/(0.d0,1.d0)
               enddo
               deallocate(vwf,svec)
            enddo
         endif
         if(cellinfo%outside_spheres) then
            nodr=cellinfo%order
            nblk=nodr*(nodr+2)
            allocate(vwf(3,nblk,2))
            rtran(:)=rpos(:)-cellinfo%rcell(:)
            call vwhcalc(rtran,ri2,nodr,1,vwf,index_model=2)
            do p=1,2
               evec(:,p)=evec(:,p)+(matmul(vwf(:,:,1),cellinfo%reg_source_vector(:,1,p)) &
                  +matmul(vwf(:,:,2),cellinfo%reg_source_vector(:,2,p)))
               hvec(:,p)=hvec(:,p)+(matmul(vwf(:,:,1),cellinfo%reg_source_vector(:,1,p)) &
                  -matmul(vwf(:,:,2),cellinfo%reg_source_vector(:,2,p)))*ri2(1)/(0.d0,1.d0)
            enddo
            deallocate(vwf)
         endif
         end subroutine source_field_calculate_fast

         subroutine source_field_calculate(rpos,sourcevec,host,layer,evec,hvec)
         implicit none
         integer :: layer,host,nodr,nblk,i,p,num,j,np(2),np0(2)
         real(8) :: rpos(3),rtran(3)
         complex(8) :: sourcevec(number_eqns,2),evec(3,2),hvec(3,2),ri,ri2(2),pshift,pshift0
         complex(8), allocatable :: vwf(:,:,:),svec(:,:)
         type(linked_sphere_list), pointer :: slist

         evec=0.d0
         hvec=0.d0
!return
         num=sphere_links(host,layer)%number
         if(host.eq.0) then
            ri=layer_ref_index(layer)
            ri2=ri
         else
            ri2=sphere_ref_index(:,host)
            ri=2.d0/(1.d0/ri2(1)+1.d0/ri2(2))
         endif

         np0=0
         pshift0=1.d0
         if(host.ne.0.and.periodic_lattice) then
            i=host
            do while(host_sphere(i).ne.0)
               i=host_sphere(i)
            enddo
            rtran=rpos-sphere_position(:,i)
            np0(1:2)=floor((rtran(1:2)+cell_width(1:2)/2.d0)/cell_width(1:2))
            pshift0=cdexp((0.d0,1.d0)*sum(incident_lateral_vector*dble(np0)*cell_width))
         endif

         if(num.gt.0) then
            slist=>sphere_links(host,layer)%sphere_list
            do i=1,num
               j=slist%sphere
               if(i.lt.num) slist=>slist%next
               nodr=sphere_order(j)
               nblk=nodr*(nodr+2)
               allocate(vwf(3,nblk,2),svec(nblk,2))
               rtran=rpos-sphere_position(:,j)
               if(host.eq.0.and.periodic_lattice) then
                  np(1:2)=floor((rtran(1:2)+cell_width(1:2)/2.d0)/cell_width(1:2))
                  rtran(1:2)=rtran(1:2)-cell_width(1:2)*dble(np(1:2))
                  pshift=cdexp((0.d0,1.d0)*sum(incident_lateral_vector*dble(np)*cell_width))*pshift0
               else
                  pshift=pshift0
               endif
               call vwhcalc(rtran,ri2,nodr,3,vwf,index_model=2)
               do p=1,2
                  svec=reshape(sourcevec(sphere_offset(j)+1:sphere_offset(j)+sphere_block(j),p), &
                     (/nblk,2/))
                  evec(:,p)=evec(:,p)+(matmul(vwf(:,:,1),svec(:,1))+matmul(vwf(:,:,2),svec(:,2)))*pshift
                  hvec(:,p)=hvec(:,p)+(matmul(vwf(:,:,1),svec(:,1))-matmul(vwf(:,:,2),svec(:,2)))*ri/(0.d0,1.d0)*pshift
               enddo
               deallocate(vwf,svec)
            enddo
         endif
         if(host.ne.0) then
            nodr=sphere_order(host)
            nblk=nodr*(nodr+2)
            allocate(vwf(3,nblk,2),svec(nblk,2))
            rtran=rpos-sphere_position(:,host)
            rtran(1:2)=rtran(1:2)-cell_width(1:2)*dble(np0(1:2))
            call vwhcalc(rtran,ri2,nodr,1,vwf,index_model=2)
            do p=1,2
               svec(:,:)=internal_field_vector(host)%vector(:,:,p)
               evec(:,p)=evec(:,p)+(matmul(vwf(:,:,1),svec(:,1))+matmul(vwf(:,:,2),svec(:,2)))*pshift0
               hvec(:,p)=hvec(:,p)+(matmul(vwf(:,:,1),svec(:,1))-matmul(vwf(:,:,2),svec(:,2)))*ri/(0.d0,1.d0)*pshift0
            enddo
         endif
         end subroutine source_field_calculate

         subroutine surface_field_calculate(rpos,sourcevec,cellinfo,evec,hvec)
         implicit none
         logical :: storecalc
         integer :: layer,nodr,nblk,i,p
         real(8) :: rpos(3),rtran(3),rcell(3)
         complex(8) :: sourcevec(number_eqns,2),evec(3,2),hvec(3,2),ri,rvec(3,2)
         complex(8), allocatable :: vwf(:,:,:)
         type(cell_info), pointer :: cellinfo

         layer=cellinfo%layer
         ri=layer_ref_index(layer)
         if(store_surface_vector) then
            storecalc=.true.
         else
            storecalc=.false.
         endif
!evec=0.
!hvec=0.
!return
         if(store_surface_vector) then
            nodr=cellinfo%order
            rcell=cellinfo%rcell
            nblk=nodr*(nodr+2)
            allocate(vwf(3,nblk,2))
            if(.not.associated(cellinfo%vector))  then
               call stored_surface_vector_calculate(nodr,rcell,sourcevec,cellinfo%vector)
            endif
            rtran(:)=rpos(:)-rcell(:)
            call vwhcalc(rtran,(/ri,ri/),nodr,1,vwf,index_model=2)
            do p=1,2
               evec(:,p)=matmul(vwf(:,:,1),cellinfo%vector(1:nblk,1,p)) &
                  +matmul(vwf(:,:,2),cellinfo%vector(1:nblk,2,p))
               hvec(:,p)=(matmul(vwf(:,:,1),cellinfo%vector(1:nblk,1,p)) &
                  -matmul(vwf(:,:,2),cellinfo%vector(1:nblk,2,p)))*ri/(0.d0,1.d0)
            enddo
            deallocate(vwf)
         else
            evec=0.d0
            hvec=0.d0
            do i=1,number_spheres
               if(host_sphere(i).ne.0) cycle
               rtran(:)=rpos(:)-sphere_position(:,i)
!write(*,*) rtran
!call flush(6)
               do p=1,2
                  if(periodic_lattice) then
                     call plane_boundary_lattice_interaction(1,sphere_order(i),rtran(1),rtran(2),rpos(3), &
                        sphere_position(3,i),rvec,source_vector=sourcevec(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),p), &
                        include_source=.false.,lr_transformation=.true.,index_model=2)
                  else
                     call plane_interaction(1,sphere_order(i), &
                        rtran(1),rtran(2),sphere_position(3,i),rpos(3), &
                        rvec,source_vector=sourcevec(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),p), &
                        index_model=2,lr_transformation=.true.,make_symmetric=.false.)
                  endif
                  evec(:,p)=evec(:,p)+matmul(vwf_0(:,:,1),rvec(:,1))+matmul(vwf_0(:,:,2),rvec(:,2))
                  hvec(:,p)=hvec(:,p)+(matmul(vwf_0(:,:,1),rvec(:,1))-matmul(vwf_0(:,:,2),rvec(:,2)))*ri/(0.d0,1.d0)
               enddo
            enddo
         endif
         end subroutine surface_field_calculate

         subroutine incident_field_calculate(rpos,layer,alpha,sinc,dir,cellinfo,evec,hvec)
         implicit none
         integer :: p,layer,dir,nodr,nblk
         real(8) :: alpha,sinc,rpos(3),rcell(3),rtran(3)
         complex(8) :: riinc,pmnp(3,2,2),evec(3,2),hvec(3,2)
         complex(8), allocatable :: vwf(:,:,:)
         type(cell_info), pointer :: cellinfo
         riinc=layer_ref_index(layer)
         if(incident_gb) then
            if(store_surface_vector) then
               nodr=cellinfo%order
               rcell=cellinfo%rcell
               nblk=nodr*(nodr+2)
               allocate(vwf(3,nblk,2))
!write(*,'(6f10.4)') rpos,rcell
               if(.not.associated(cellinfo%gb_vector))  then
!write(*,'('' allocating gb vec'' )')
                  allocate(cellinfo%gb_vector(nblk,2,2))
                  call layergaussbeamcoef(alpha,sinc,dir,rcell,nodr,cellinfo%gb_vector,include_direct=.true., &
                    include_indirect=.true.)
               endif
               rtran(:)=rpos(:)-rcell(:)
               call vwhcalc(rtran,(/riinc,riinc/),nodr,1,vwf,index_model=2)
               do p=1,2
                  evec(:,p)=matmul(vwf(:,:,1),cellinfo%gb_vector(1:nblk,1,p)) &
                     +matmul(vwf(:,:,2),cellinfo%gb_vector(1:nblk,2,p))
                  hvec(:,p)=(matmul(vwf(:,:,1),cellinfo%gb_vector(1:nblk,1,p)) &
                     -matmul(vwf(:,:,2),cellinfo%gb_vector(1:nblk,2,p)))*riinc/(0.d0,1.d0)
               enddo
               deallocate(vwf)
            else
               call layergaussbeamcoef(alpha,sinc,dir,rcell,1,pmnp,include_direct=.true., &
                    include_indirect=.true.)
               do p=1,2
                  evec(:,p)=matmul(vwf_0(:,:,1),pmnp(:,1,p))+matmul(vwf_0(:,:,2),pmnp(:,2,p))
                  hvec(:,p)=(matmul(vwf_0(:,:,1),pmnp(:,1,p))-matmul(vwf_0(:,:,2),pmnp(:,2,p)))*riinc/(0.d0,1.d0)
               enddo
            endif
         else
            call layerplanewavecoef(alpha,sinc,dir,rpos,1,pmnp)
            do p=1,2
               evec(:,p)=matmul(vwf_0(:,:,1),pmnp(:,1,p))+matmul(vwf_0(:,:,2),pmnp(:,2,p))
               hvec(:,p)=(matmul(vwf_0(:,:,1),pmnp(:,1,p))-matmul(vwf_0(:,:,2),pmnp(:,2,p)))*riinc/(0.d0,1.d0)
            enddo
         endif
         end subroutine incident_field_calculate

         subroutine stored_surface_vector_calculate(nodr,rc,sourcevec,storedvector)
         implicit none
         integer :: nodr,nblk,i,p
         real(8) :: rc(3),rhovec(2)
         complex(8) :: sourcevec(number_eqns,2)
         complex(8), allocatable :: vector(:,:,:)
         complex(8), pointer :: storedvector(:,:,:)

         nblk=nodr*(nodr+2)
         allocate(vector(nblk,2,2),storedvector(nblk,2,2))
         storedvector(1:nblk,1:2,1:2)=0.d0
         do i=1,number_spheres
            if(host_sphere(i).ne.0) cycle
            rhovec(1:2)=rc(1:2)-sphere_position(1:2,i)
            vector=0.d0
            do p=1,2
               if(periodic_lattice) then
                  call plane_boundary_lattice_interaction(nodr,sphere_order(i), &
                     rhovec(1),rhovec(2),rc(3),sphere_position(3,i), &
                     vector(:,:,p), &
                     source_vector=sourcevec(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),p), &
                     include_source=.false.,lr_transformation=.true.,index_model=2)
               else
                  call plane_interaction(nodr,sphere_order(i), &
                     rhovec(1),rhovec(2),sphere_position(3,i),rc(3), &
                     vector(:,:,p), &
                     source_vector=sourcevec(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),p), &
                     index_model=2,lr_transformation=.true., &
                     make_symmetric=.false.)
               endif
            enddo
            storedvector(1:nblk,1:2,1:2)=storedvector(1:nblk,1:2,1:2)+vector(1:nblk,1:2,1:2)
         enddo
         deallocate(vector)
         end subroutine stored_surface_vector_calculate

         subroutine stored_source_vector_calculate(sourcevec,cellinfo)
         implicit none
         integer :: nodr,nblk,i,p
         real(8) :: rc(3),r
         complex(8) :: sourcevec(number_eqns,2),ri2(2)
         type(cell_info), pointer :: cellinfo
         type(linked_sphere_list), pointer :: slist
         type(translation_data) :: tranmat

         ri2(1:2)=layer_ref_index(cellinfo%layer)
         nodr=cellinfo%order
         nblk=nodr*(nodr+2)
         allocate(cellinfo%reg_source_vector(nblk,2,2))
         cellinfo%reg_source_vector(1:nblk,1:2,1:2)=0.d0
         cellinfo%nispheres=0
         cellinfo%outside_spheres=.false.
         do i=1,number_spheres
            if(host_sphere(i).ne.0) cycle
            rc=cellinfo%rcell(:)-sphere_position(:,i)
            r=sqrt(sum(rc**2))
            if(r.le.2.d0*near_field_expansion_spacing) then
               if(cellinfo%nispheres.eq.0) then
                  allocate(cellinfo%sphere_list)
                  slist=>cellinfo%sphere_list
               else
                  allocate(slist%next)
                  slist=>slist%next
               endif
               slist%sphere=i
               cellinfo%nispheres=cellinfo%nispheres+1
            else
               cellinfo%outside_spheres=.true.
               tranmat%matrix_calculated=.false.
               tranmat%vswf_type=3
               tranmat%translation_vector=rc
               tranmat%refractive_index=ri2
               tranmat%rot_op=max(sphere_order(i),nodr).ge.translation_switch_order
               do p=1,2
                  call coefficient_translation(sphere_order(i),2,nodr,2, &
                     sourcevec(sphere_offset(i)+1:sphere_offset(i)+sphere_block(i),p), &
                     cellinfo%reg_source_vector(:,:,p),tranmat)
               enddo
               if(tranmat%rot_op) then
                  deallocate(tranmat%rot_mat,tranmat%phi_mat,tranmat%z_mat)
               else
                  deallocate(tranmat%gen_mat)
               endif
            endif
         enddo
         end subroutine stored_source_vector_calculate

         subroutine grid_point_initialize(griddim,gridinfo)
         implicit none
         logical :: ingrid
         integer :: depth,i,l,layer,griddim(3),ix,iy,iz,nodr,nbound,zbsign,ncell(3),jy,jx,jlim(2)
         type(grid_info) :: gridinfo(griddim(1),griddim(2),griddim(3))
         real(8) :: x,y,z,pbcellsize(1:max(1,number_plane_boundaries)),zbound,zbcell,rcell(3),cellsize,spos(3)
         type(cell_info) :: cellinfo
         type(linked_sphere_data), pointer :: slist

!         call clear_cell_info_list()
         total_cells=0
         do iz=1,griddim(3)
            do iy=1,griddim(2)
               do ix=1,griddim(1)
                  gridinfo(ix,iy,iz)%initialized=.false.
                  gridinfo(ix,iy,iz)%cellnum=0
               enddo
            enddo
         enddo

         number_intersecting_spheres=0
         allocate(intersecting_spheres)
         slist=>intersecting_spheres
         if(periodic_lattice) then
            jlim=(/-1,1/)
         else
            jlim=0
         endif
         do depth=max_sphere_depth,0,-1
            do i=1,number_spheres
               if(sphere_depth(i).ne.depth) cycle
               do jy=jlim(1),jlim(2)
                  do jx=jlim(1),jlim(2)
                     spos=sphere_position(:,i)+dble((/jx,jy,0/))*(/cell_width(1),cell_width(2),0.d0/)
                     call sphere_to_grid_points(i,spos,griddim,gridinfo,ingrid)
                     if(ingrid) then
                        slist%position=spos
                        slist%radius=sphere_radius(i)
                        number_intersecting_spheres=number_intersecting_spheres+1
                        allocate(slist%next)
                        slist=>slist%next
                     endif
                  enddo
               enddo
            enddo
         enddo
         do l=1,max(1,number_plane_boundaries)
            pbcellsize(l)=near_field_expansion_spacing
            if(plane_surface_present) then
               do i=1,number_spheres
                  if((sphere_layer(i).eq.l-1.or.sphere_layer(i).eq.l).and.host_sphere(i).eq.0) then
                     pbcellsize(l)=min(pbcellsize(l),1.d0*abs(plane_boundary_position(l)-sphere_position(3,i)))
                  endif
               enddo
            endif
         enddo

         do iz=1,griddim(3)
            z=(dble(iz)-0.5d0)*grid_spacing(3)+grid_region(3,1)
            layer=layer_id(z)
            if(plane_surface_present) then
               if(layer.eq.0) then
                  nbound=1
                  zbound=plane_boundary_position(1)-z
                  zbcell=pbcellsize(1)
               elseif(layer.eq.number_plane_boundaries) then
                  nbound=number_plane_boundaries
                  zbound=z-plane_boundary_position(number_plane_boundaries)
                  zbcell=pbcellsize(number_plane_boundaries)
               else
                  nbound=layer
                  zbound=z-plane_boundary_position(layer)
                  if(zbound.gt.plane_boundary_position(layer+1)-z) then
                     zbound=plane_boundary_position(layer+1)-z
                     nbound=layer+1
                  endif
                  zbcell=min(pbcellsize(nbound), &
                     plane_boundary_position(layer+1)-plane_boundary_position(layer))
               endif
               if(nbound.eq.layer) then
                  zbsign=1
               else
                  zbsign=-1
               endif
               if(zbound.le.pbcellsize(nbound)) then
                  cellsize=pbcellsize(nbound)
                  rcell(3)=plane_boundary_position(nbound)+0.5d0*zbsign*zbcell
               else
                  cellsize=near_field_expansion_spacing
                  zbound=plane_boundary_position(nbound)+zbsign*pbcellsize(nbound)
                  rcell(3)=zbound &
                     +(floor(abs(z-zbound)/near_field_expansion_spacing)+0.5d0) &
                     *zbsign*near_field_expansion_spacing
               endif
            else
               cellsize=near_field_expansion_spacing
               rcell(3)=(dble(floor((z-grid_region(3,1))/cellsize))+0.5d0)*cellsize+grid_region(3,1)
            endif
            if(griddim(3).eq.1) then
               ncell(3)=0
               rcell(3)=grid_region(3,1)
            else
               ncell(3)=floor(rcell(3)/grid_spacing(3))
            endif
            nodr=max(ceiling(near_field_expansion_order*0.9999d0*(cellsize/near_field_expansion_spacing)),1)
!            nodr=near_field_expansion_order

            do iy=1,griddim(2)
               y=(dble(iy)-0.5d0)*grid_spacing(2)+grid_region(2,1)
               if(griddim(2).eq.1) then
                  rcell(2)=grid_region(2,1)
                  ncell(2)=0
               else
                  rcell(2)=(dble(floor((y-grid_region(2,1))/cellsize))+0.5d0)*cellsize+grid_region(2,1)
                  ncell(2)=floor(rcell(2)/grid_spacing(2))
               endif
               do ix=1,griddim(1)
                  if(gridinfo(ix,iy,iz)%initialized) cycle
                  x=(dble(ix)-0.5d0)*grid_spacing(1)+grid_region(1,1)
                  if(griddim(1).eq.1) then
                     rcell(1)=grid_region(1,1)
                     ncell(1)=0
                  else
                     rcell(1)=(dble(floor((x-grid_region(1,1))/cellsize))+0.5d0)*cellsize+grid_region(1,1)
                     ncell(1)=floor(rcell(1)/grid_spacing(1))
                  endif
                  cellinfo%ncell=ncell
                  cellinfo%host=0
                  cellinfo%layer=layer
                  cellinfo%order=nodr
                  cellinfo%rcell(:)=rcell
                  call point_at_list_elem(cellinfo,gridinfo(ix,iy,iz)%cellinfo,cell_info_list)
                  gridinfo(ix,iy,iz)%initialized=.true.
                  gridinfo(ix,iy,iz)%cellnum=total_cells
!write(*,'(5i3,10es12.4)') ix,iz,ncell(1),ncell(3),total_cells,rcell(1),rcell(3)
!call flush(6)
               enddo
            enddo
         enddo

         end subroutine grid_point_initialize

         subroutine sphere_to_grid_points(sphere,spos,griddim,gridinfo,ingrid)
         implicit none
         logical :: ingrid,skip
         integer :: sphere,ic(3),i,limits(3,2),ix,iy,iz,griddim(3),n1,n2
         real(8) :: rpos(3),rc(3),r,r1,r2,spos(3)
         type(grid_info) :: gridinfo(griddim(1),griddim(2),griddim(3))
         type(cell_info) :: cellinfo

         ingrid=.false.
         skip=.false.
         do i=1,3
            if(grid_spacing(i).eq.0.d0) then
               if(abs(grid_region(i,1)-spos(i)).gt.sphere_radius(sphere)) then
                  skip=.true.
               else
                  limits(i,1)=1
                  limits(i,2)=1
               endif
            else
               r1=spos(i)-sphere_radius(sphere)
               if(r1.gt.grid_region(i,2)) then
                  skip=.true.
               endif
               r2=spos(i)+sphere_radius(sphere)
               if(r2.lt.grid_region(i,1)) then
                  skip=.true.
               endif
               n1=ceiling((r1-grid_region(i,1))/grid_spacing(i))
               n2=ceiling((r2-grid_region(i,1))/grid_spacing(i))
               limits(i,1)=max(1,n1)
               limits(i,2)=min(griddim(i),n2)
            endif
         enddo
         if(skip) return

         do iz=limits(3,1),limits(3,2)
            do iy=limits(2,1),limits(2,2)
               do ix=limits(1,1),limits(1,2)
                  ic=(/ix,iy,iz/)
                  rpos(:)=(dble(ic(:))-(/0.5d0,0.5d0,0.5d0/))*grid_spacing(:)+grid_region(:,1)
                  rc(:)=rpos(:)-spos(:)
                  r=sqrt(dot_product(rc(:),rc(:)))
                  if(r.gt.sphere_radius(sphere)) cycle
                  if(gridinfo(ic(1),ic(2),ic(3))%initialized) cycle
                  ingrid=.true.
                  cellinfo%ncell=(/0,0,0/)
                  cellinfo%host=sphere
                  cellinfo%layer=sphere_layer(sphere)
                  cellinfo%order=sphere_order(sphere)
                  cellinfo%rcell(:)=sphere_position(:,sphere)
                  call point_at_list_elem(cellinfo,gridinfo(ic(1),ic(2),ic(3))%cellinfo,cell_info_list)
                  gridinfo(ic(1),ic(2),ic(3))%initialized=.true.
                  gridinfo(ic(1),ic(2),ic(3))%cellnum=total_cells
               enddo
            enddo
         enddo
         end subroutine sphere_to_grid_points

         subroutine point_at_list_elem(delem,elem,list)
         implicit none
         logical :: inlist
         type(linked_cell_list), pointer :: list,tlist
         type(cell_info), pointer :: elem
         type(cell_info) :: delem

         if(associated(list)) then
            tlist=>list
            inlist=.true.
            do while(delem%ncell(1).ne.tlist%cellinfo%ncell(1) &
                    .or.delem%ncell(2).ne.tlist%cellinfo%ncell(2) &
                    .or.delem%ncell(3).ne.tlist%cellinfo%ncell(3) &
                    .or.delem%host.ne.tlist%cellinfo%host &
                    .or.delem%layer.ne.tlist%cellinfo%layer &
                    .or.delem%order.ne.tlist%cellinfo%order)
              if(.not.associated(tlist%next)) then
                  inlist=.false.
                  exit
               else
                  tlist=>tlist%next
               endif
            enddo
            if(inlist) then
               elem=>tlist%cellinfo
            else
               total_cells=total_cells+1
               allocate(tlist%next)
               tlist=>tlist%next
               tlist%cellinfo%ncell(:)=delem%ncell(:)
               tlist%cellinfo%rcell(:)=delem%rcell(:)
               tlist%cellinfo%host=delem%host
               tlist%cellinfo%layer=delem%layer
               tlist%cellinfo%order=delem%order
               elem=>tlist%cellinfo
            endif
         else
            total_cells=1
            allocate(list)
            list%cellinfo%ncell(:)=delem%ncell(:)
            list%cellinfo%rcell(:)=delem%rcell(:)
            list%cellinfo%host=delem%host
            list%cellinfo%layer=delem%layer
            list%cellinfo%order=delem%order
            elem=>list%cellinfo
         endif
         end subroutine point_at_list_elem


         subroutine write_output_header(griddim,outputunit)
         implicit none
         integer :: griddim(3),outputunit,n,j,l1,l2
         type(linked_sphere_data), pointer :: slist
         write(outputunit,'('' run number:'')')
         write(outputunit,'(i5)') local_run_number
         local_run_number=local_run_number+1
         n=number_intersecting_spheres
         write(outputunit,'(i5)') n
         slist=>intersecting_spheres
         do j=1,n
            write(outputunit,'(4es12.4)') slist%position,slist%radius
            if(j.lt.n) slist=>slist%next
         enddo
         l1=layer_id(grid_region(3,1))
         l2=layer_id(grid_region(3,2))
         write(outputunit,'(i5)') l2-l1
         do j=l1+1,l2
            write(outputunit,'(es12.4)') plane_boundary_position(j)
         enddo
         write(outputunit,'(3es12.4)') grid_region(:,1)
         write(outputunit,'(3es12.4)') grid_region(:,2)
         write(outputunit,'(3i5)') griddim(:)

         end subroutine write_output_header



      end module nearfield
!
!  module solver: subroutines for solving interaction equations for fixed orientation
!  and T matrix problems
!
!
!  last revised: 15 January 2011
!                february 2013
!
      module solver
      use mpidefs
      use intrinsics
      use numconstants
      use specialfuncs
      use spheredata
      use mie
      use translation
      use scatprops
      implicit none

      contains

         subroutine tmatrix_solution(solution_method,solution_eps,convergence_eps, &
            max_iterations,t_matrix_file,procs_per_soln,mpi_comm, &
            sphere_qeff,solution_status)
         implicit none
         logical :: firstrun,initialize,itersoln,continueloop
         integer :: iter,niter,istat,rank,maxiter,&
                    numprocs,mpicomm,i,nblkt,l,k,q,ka,la,nsolns, &
                    n,m,p,nssoln,kq,ns
         integer, save :: pcomm,prank,pgroup,ppsoln,pcomm0
         integer, optional :: mpi_comm,procs_per_soln,max_iterations,solution_status
         real(8) :: r0(3),maxerr,qeffi(3,number_spheres), &
            dqeffi(3,number_spheres),qeff(3),qeffold(3),time0,time1,timepersoln,timeleft, &
            solneps,conveps,solnerr,converr(1),qteff(3),dqeffikq(3,number_spheres), &
            ttime(0:6),dqteff(3)
         real(8), optional :: sphere_qeff(3,number_spheres), &
            solution_eps,convergence_eps
         complex(8) :: amnpkq(number_eqns),pmnpkq(number_eqns),pmnpan(number_eqns)
         complex(8), allocatable :: pmnp0(:,:,:),amnp0(:)
         character*4 :: timeunit
         character*128 :: tmatrixfile
         character*1, optional :: solution_method
         character*128, optional :: t_matrix_file
         data firstrun/.true./
         istat=0
         if(present(solution_method)) then
            itersoln=solution_method.eq.'i'
         else
            itersoln=.true.
         endif
         if(present(solution_eps)) then
            solneps=solution_eps
         else
            solneps=1.d-6
         endif
         if(present(convergence_eps)) then
            conveps=convergence_eps
         else
            conveps=1.d-6
         endif
         if(present(max_iterations)) then
            niter=max_iterations
         else
            niter=100
         endif
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(procs_per_soln)) then
            ppsoln=procs_per_soln
         else
            ppsoln=1
         endif
         if(present(t_matrix_file)) then
            tmatrixfile=t_matrix_file
         else
            tmatrixfile='tmatrix_temp.dat'
         endif

         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         if(.not.itersoln) ppsoln=1
         ppsoln=min(ppsoln,numprocs)
         nssoln=max(1,numprocs/ppsoln)


         if(firstrun) then
            firstrun=.false.
            if(numprocs.gt.1) then
               pgroup=floor(dble(rank)/dble(ppsoln))
               call mstm_mpi(mpi_command='split', &
                  mpi_color=pgroup,mpi_key=rank, &
                  mpi_new_comm=pcomm, &
                  mpi_comm=mpicomm)
               call mstm_mpi(mpi_command='rank',mpi_rank=prank,mpi_comm=pcomm)
               call mstm_mpi(mpi_command='split', &
                  mpi_color=prank,mpi_key=rank, &
                  mpi_new_comm=pcomm0, &
                  mpi_comm=mpicomm)
            else
               pgroup=0
               pcomm=mpicomm
               pcomm0=mpicomm
            endif
         endif

         nblkt=2*t_matrix_order*(t_matrix_order+2)
         r0=cluster_origin(:)
         qeffi=0.d0
         qeffold=0.d0
         qteff=0.d0
         nsolns=0
         initialize=.true.
         istat=0
         if(rank.eq.0) then
            open(20,file=tmatrixfile)
            write(20,'(2i4)') t_matrix_order,t_matrix_order
            time0=mstm_mpi_wtime()
            close(20)
         endif

         maxiter=0
         maxerr=0.d0

         if(rank.eq.0) then
            write(run_print_unit,'(''  n   # its  qext         qabs'',&
                      &''           error     sec/soln est. time rem.'')')
            call flush(run_print_unit)
         endif

         do l=1,t_matrix_order
if(rank.eq.0) then
ttime(0)=mstm_mpi_wtime()
endif
            if(rank.eq.0) then
               time1=mstm_mpi_wtime()
            endif

            allocate(pmnp0(0:l+1,l,2),amnp0(2*l*(l+2)))
            kq=0
            continueloop=.true.
            do while(continueloop)
               ns=0
               dqeffi=0.d0
               dqteff=0.d0
               do i=1,nssoln
                  ns=ns+1
                  kq=kq+1
                  k=-l+(kq-1)/2
                  q=mod(kq-1,2)+1
                  if((i-1).eq.pgroup) then
                     pmnp0=0.d0
                     if(k.le.-1) then
                        ka=l+1
                        la=-k
                     else
                        ka=k
                        la=l
                     endif
!
!  the t matrix is te-tm based; hence the following two lines
!
                     pmnp0(ka,la,1)=.5d0
                     pmnp0(ka,la,2)=-.5d0*(-1)**q
                     call distribute_from_common_origin(l,pmnp0,pmnpkq, &
                        number_rhs=1, &
                        origin_position=r0, &
                        origin_host=0, &
                        vswf_type=1, &
                        mpi_comm=pcomm)

                     call multmiecoeffmult(number_eqns,1,1,pmnpkq,pmnpan)
                     amnpkq=pmnpan
                     if(niter.ne.0) then
                        if(itersoln) then
                           call cbicg(niter,solneps,pmnpan,amnpkq,0, &
                                  iter,solnerr,initialize_solver=initialize, &
                                  mpi_comm=pcomm)
                           maxiter=max(iter,maxiter)
                           maxerr=max(solnerr,maxerr)
                           if(iter.gt.niter.or.solnerr.gt.solneps) istat=1
                        else
                           call direct_solver(pmnpan,amnpkq, &
                              initialize_solver=initialize)
                           maxerr=0.d0
                           maxiter=0
                        endif
                     endif
                     initialize=.false.
                     call qefficiencyfactors(number_spheres,1,amnpkq, &
                        pmnpkq,dqeffi,mpi_comm=pcomm)
                     dqeffi(3,:)=dqeffi(1,:)-dqeffi(2,:)
                     amnp0=0.d0
                     call merge_to_common_origin(l,amnpkq, &
                        amnp0, &
                        number_rhs=1, &
                        origin_position=r0, &
                        mpi_comm=pcomm)
                     ka=amnpaddress(k,l,q,l,2)
                     dqteff(1)=-2.d0/vol_radius**2*dble(amnp0(ka))
                     dqteff(3)=2.d0/vol_radius**2 &
                        *dble(sum(amnp0(:)*dconjg(amnp0(:))))
                     dqteff(2)=dqteff(1)-dqteff(3)

                  endif

                  if(kq.eq.2*(2*l+1)) then
                     continueloop=.false.
                     exit
                  endif
               enddo

               if(prank.eq.0) then
                  call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=dqteff, &
                     mpi_number=3,mpi_operation=mstm_mpi_sum,mpi_comm=pcomm0)
                  call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=dqeffi, &
                     mpi_number=3*number_spheres,mpi_operation=mstm_mpi_sum,mpi_comm=pcomm0)

                  qeffi=qeffi+dqeffi
                  qteff=qteff+dqteff

                  do i=1,ns
                     if(i-1.eq.pgroup) then
                        open(20,file=tmatrixfile,access='append')
                        do n=1,l
                           do m=-n,n
                              do p=1,2
                                 ka=amnpaddress(m,n,p,l,2)
                                 write(20,'(''('',e18.10,'','',e18.10,'')'')') amnp0(ka)
                              enddo
                           enddo
                        enddo
                        close(20)
                     endif
                     call mstm_mpi(mpi_command='barrier',mpi_comm=pcomm0)
                  enddo
               endif
            enddo

            if(rank.eq.0) then
               qeff=0.d0
               do i=1,number_spheres
                  qeff(:)=qeff(:)+qeffi(:,i)*sphere_radius(i)**2/vol_radius**2
               enddo

               converr=qeff(1)-qeffold(1)
               qeffold=qeff
               nsolns=nsolns+(l+l+1)*2

               timepersoln=(mstm_mpi_wtime()-time1)/dble((l+l+1)*2)
               timeleft=timepersoln*(nblkt-nsolns)
               if(timeleft.gt.3600.d0) then
                  timeleft=timeleft/3600.d0
                  timeunit=' hrs'
               elseif(timeleft.gt.60.d0) then
                  timeleft=timeleft/60.d0
                  timeunit=' min'
               else
                  timeunit=' sec'
               endif
               write(run_print_unit,'(i4,i5,3e13.5,2f8.2,a4)') l,maxiter,qeff(1:2), &
                  converr,timepersoln,timeleft,timeunit
               call flush(run_print_unit)
            endif
            deallocate(amnp0,pmnp0)

            call mstm_mpi(mpi_command='bcast',mpi_rank=0,mpi_send_buf_dp=converr,mpi_number=1)

            if(converr(1).lt.conveps) exit
         enddo

         if(l.lt.t_matrix_order) then
            t_matrix_order=l
            if(rank.eq.0) then
               open(20,file=tmatrixfile,form='formatted',access='direct',recl=8)
               write(20,'(2i4)',rec=1) l,l
               close(20)
            endif
         endif

         if(present(sphere_qeff)) sphere_qeff=qeffi
         if(present(solution_status)) solution_status=istat

         end subroutine tmatrix_solution

!
!  solution of interaction equations for a fixed orientation
!
!
!  original: 15 January 2011
!  revised: 21 February 2011: modification of efficiency calculation, to calculate
!           polarized components
!  30 March 2011: took out gbfocus argument: this is not needed since positions are defined
!  relative to the gb focus.
!  20 April 2011: used 2-group MPI formulation
!  October 2011: adapted to far field approximation.
!  December 2011: changed efficiency factor calculation, adapted to generalized sphere
!                 configuration.
! february 2013: number rhs and mpi comm options added, completely rewritten.
!
!
         subroutine fixedorsoln(alpha,sinc,dir,eps,niter,amnp,qeff, &
                    qeffdim,maxerr,maxiter,iterwrite,istat, &
                    mpi_comm,excited_spheres)
         implicit none
         logical :: firstrun,exsphere(number_spheres)
         logical, save :: inp1,inp2
         logical, optional :: excited_spheres(number_spheres)
         integer :: iter,niter,istat,rank,maxiter,iterwrite,nsend,&
                    numprocs,mpicomm,prank,oddnumproc, &
                    groupsize,pgroup,mpigroup,syncgroup,i,p,dir,qeffdim
         integer, save :: pcomm,synccomm1,synccomm2,p1,p2
         integer, allocatable :: grouplist(:)
         integer, optional :: mpi_comm
         real(8) :: alpha,sinc,eps,serr,qeff(3,qeffdim,number_spheres),maxerr
         complex(8) :: amnp(number_eqns,2)
         complex(8), allocatable :: pmnpan(:),pmnp0(:,:)
         data firstrun/.true./
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(excited_spheres)) then
            exsphere=excited_spheres
         else
            exsphere=.true.
         endif
!         if(mpicomm.eq.mpi_comm_null) return
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         if(firstrun) then
            firstrun=.false.
            if(numprocs.gt.1) then
               oddnumproc=mod(numprocs,2)
               pgroup=floor(dble(2*rank)/dble(numprocs))+1
               p1=pgroup
               p2=p1
               call mstm_mpi(mpi_command='split', &
                  mpi_color=pgroup,mpi_key=rank, &
                  mpi_new_comm=pcomm, &
                  mpi_comm=mpicomm)
               call mstm_mpi(mpi_command='rank',mpi_rank=prank,mpi_comm=pcomm)
               call mstm_mpi(mpi_command='group',mpi_group=mpigroup,mpi_comm=mpicomm)
               groupsize=numprocs/2+1
               allocate(grouplist(groupsize))
               grouplist(1)=0
               do i=1,groupsize-1
                  grouplist(i+1)=i+(numprocs/2)-1+oddnumproc
               enddo
               inp1=.false.
               do i=1,groupsize
                  if(rank.eq.grouplist(i)) then
                     inp1=.true.
                     exit
                  endif
               enddo
               call mstm_mpi(mpi_command='incl', &
                  mpi_group=mpigroup,&
                  mpi_size=groupsize,&
                  mpi_new_group_list=grouplist,&
                  mpi_new_group=syncgroup)
               call mstm_mpi(mpi_command='create',&
                  mpi_group=syncgroup,&
                  mpi_comm=mpicomm,&
                  mpi_new_comm=synccomm1)
               deallocate(grouplist)
               groupsize=groupsize+oddnumproc
               allocate(grouplist(groupsize))
               grouplist(1)=numprocs/2+oddnumproc
               do i=1,groupsize-1
                  grouplist(i+1)=i-1
               enddo
               inp2=.false.
               do i=1,groupsize
                  if(rank.eq.grouplist(i)) then
                     inp2=.true.
                     exit
                  endif
               enddo
               call mstm_mpi(mpi_command='incl', &
                  mpi_group=mpigroup,&
                  mpi_size=groupsize,&
                  mpi_new_group_list=grouplist,&
                  mpi_new_group=syncgroup)
               call mstm_mpi(mpi_command='create',&
                  mpi_group=syncgroup,&
                  mpi_comm=mpicomm,&
                  mpi_new_comm=synccomm2)
               deallocate(grouplist)
            else
               p1=1
               p2=2
               pcomm=mpicomm
            endif
         endif

!write(*,*) ' step 3'
!call flush(6)
if(light_up) then
write(*,'('' s8.2.1 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')

         allocate(pmnpan(number_eqns),pmnp0(number_eqns,2))
         call sphereplanewavecoef(alpha,sinc,dir,pmnp0, excited_spheres=exsphere)
!if(phase_shift_form) call phase_shift(pmnp0,1)

         do p=p1,p2
            istat=0
            maxiter=0
            maxerr=0.
!
!  calculate the two solutions
!
if(light_up) then
write(*,'('' s8.2.2 '',i3)') mstm_global_rank
call flush(6)
endif
            call multmiecoeffmult(number_eqns,1,1,pmnp0(:,p),pmnpan)
            amnp(:,p)=pmnpan
            if(niter.ne.0) then
!write(*,*) 'step 2, p:',p
!call flush(6)
if(light_up) then
write(*,'('' s8.2.3 '',i3)') mstm_global_rank
call flush(6)
endif
               call cbicg(niter,eps,pmnpan,amnp(:,p),iterwrite, &
                      iter,serr,mpi_comm=pcomm)
            else
               iter=0
               serr=0.d0
            endif
            maxiter=max(iter,maxiter)
            maxerr=max(serr,maxerr)
            if(iter.gt.niter.or.serr.gt.eps) istat=1
         enddo

         call mstm_mpi(mpi_command='barrier')

         if(numprocs.gt.1) then
            nsend=number_eqns
            if(inp1) then
               call mstm_mpi(mpi_command='bcast', &
                  mpi_send_buf_dc=amnp(1:nsend,1), &
                  mpi_number=nsend, &
                  mpi_rank=0, &
                  mpi_comm=synccomm1)
            endif
            call mstm_mpi(mpi_command='barrier')
            if(inp2) then
               call mstm_mpi(mpi_command='bcast', &
                  mpi_send_buf_dc=amnp(1:nsend,2), &
                  mpi_number=nsend, &
                  mpi_rank=0, &
                  mpi_comm=synccomm2)
            endif
         endif
         call mstm_mpi(mpi_command='barrier')
!
!  efficiency factor calculations
!
if(light_up) then
write(*,'('' s8.2.4 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')
         if(qeffdim.eq.1) then
            i=1
         else
            i=2
         endif
         call qefficiencyfactors(number_spheres,i,amnp,pmnp0,qeff, &
                    mpi_comm=mpicomm)

!if(phase_shift_form) call phase_shift(amnp,-1)

         deallocate(pmnp0,pmnpan)
         end subroutine fixedorsoln

         subroutine direct_solver(pnp,anp,initialize_solver)
         implicit none
         logical :: initialize
         logical, save :: firstrun
         logical, optional :: initialize_solver
         integer :: i,j,offseti,offsetj,nblki,nblkj,l,ierr,rank
         integer, allocatable, save :: indx(:)
         real(8) :: rtran(3),dsign
         complex(8)  :: pnp(number_eqns),anp(number_eqns),rimed(2)
         complex(8), allocatable :: tranmat(:,:,:)
         complex(8), allocatable, save :: amat(:,:)
         data firstrun/.true./
         if(present(initialize_solver)) then
            initialize=initialize_solver
         else
            initialize=firstrun
         endif
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         rimed=1.d0

         if(initialize) then
            if(allocated(amat)) deallocate(amat,indx)
            allocate(amat(number_eqns,number_eqns),indx(number_eqns))
            amat=0.d0
            offseti=0
            do i=1,number_spheres
               nblki=sphere_order(i)*(sphere_order(i)+2)
               offsetj=0
               do j=1,number_spheres
                  nblkj=sphere_order(j)*(sphere_order(j)+2)
                  if(j.ne.i) then
                     allocate(tranmat(nblki,nblkj,2))
                     rtran=sphere_position(:,i)-sphere_position(:,j)
                     call gentranmatrix(sphere_order(j),sphere_order(i),translation_vector=rtran, &
                          refractive_index=rimed,ac_matrix=tranmat,vswf_type=3, &
                          mode_s=2,mode_t=2)
                     amat(offseti+1:offseti+nblki,offsetj+1:offsetj+nblkj)=-tranmat(1:nblki,1:nblkj,1)
                     amat(offseti+1+nblki:offseti+2*nblki,offsetj+1+nblkj:offsetj+2*nblkj)=-tranmat(1:nblki,1:nblkj,2)
   !                  amat(offseti+1:offseti+nblki,offsetj+1+nblkj:offsetj+2*nblkj)=-tranmat(1:nblki,1:nblkj,2)
   !                  amat(offseti+1+nblki:offseti+2*nblki,offsetj+1:offsetj+nblkj)=-tranmat(1:nblki,1:nblkj,2)
                     deallocate(tranmat)
                  endif
                  offsetj=offsetj+2*nblkj
               enddo
               offseti=offseti+2*nblki
            enddo
            do l=1,number_eqns
               call multmiecoeffmult(number_eqns,1,1,amat(1:number_eqns,l),anp)
               amat(1:number_eqns,l)=anp
            enddo
            offseti=0.
            do i=1,number_spheres
               nblki=sphere_order(i)*(sphere_order(i)+2)
               do l=offseti+1,offseti+2*nblki
                  amat(l,l)=1.d0
               enddo
               offseti=offseti+2*nblki
            enddo
            call lu_decomposition(amat,number_eqns,indx,dsign,ierr)
            if(ierr.ne.0) then
               if(rank.eq.0) then
                  write(run_print_unit,'('' lu decomposition failed!!!'')')
               endif
               stop
            endif
            firstrun=.false.
         endif
         anp=pnp
         call lu_backsubstitution(amat,number_eqns,indx,anp)
         end subroutine direct_solver
!
! iteration solver
! generalized complex biconjugate gradient method
! original code: Piotr Flatau, although not much remains.
! specialized to the multiple sphere problem
!
!
!  last revised: 15 January 2011
!  october 2011: translation calls modified
!  february 2013: number rhs option added, completely rewritten.
!
         subroutine cbicg(niter,eps,pnp,anp,iterwrite,iter,errmax, &
            initialize_solver,mpi_comm)
         implicit none
         logical :: firstrun,initialize,contran2(2)
         logical, save :: inp1,inp2
         logical, optional :: initialize_solver
         integer :: neqns,niter,iter,writetime,&
                    rank,iunit,iterwrite,numprocs,i,oddnumproc, &
                    mpicomm,prank,mpigroup,syncgroup,groupsize,rank0
         integer, save :: pgroup,pcomm,synccomm1,synccomm2
         integer, allocatable :: grouplist(:)
         integer, optional :: mpi_comm
         real(8) :: eps,time1,time2,eerr,enorm,errmax,errmin,time0
         complex(8)  :: pnp(number_eqns),anp(number_eqns),cak,csk,cbk,csk2
         complex(8), allocatable :: cr(:),cp(:),cw(:),cq(:),cap(:),caw(:), &
                       capt(:),cawt(:),ctin(:,:),ctout(:,:)
         data firstrun/.true./
         data writetime/0/

         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(initialize_solver)) then
            initialize=initialize_solver
         else
            initialize=.true.
         endif
         iunit=run_print_unit
         neqns=number_eqns
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank0)

if(light_up) then
write(*,'('' s8.2.3.1 '',i3)') mstm_global_rank
call flush(6)
endif
         if(dot_product(pnp,pnp).eq.0.d0) return
         allocate(cr(neqns),cp(neqns),cw(neqns),cq(neqns), &
               cap(neqns),caw(neqns),capt(neqns), &
               cawt(neqns))

         if(firstrun) then
            firstrun=.false.
            if(numprocs.gt.1.and.niter.gt.0) then
               oddnumproc=mod(numprocs,2)
               pgroup=floor(dble(2*rank)/dble(numprocs))+1
               call mstm_mpi(mpi_command='split', &
                  mpi_color=pgroup,mpi_key=rank, &
                  mpi_new_comm=pcomm, &
                  mpi_comm=mpicomm)
               call mstm_mpi(mpi_command='rank',mpi_rank=prank,mpi_comm=pcomm)
               call mstm_mpi(mpi_command='group',mpi_group=mpigroup,mpi_comm=mpicomm)
               groupsize=numprocs/2+1
               allocate(grouplist(groupsize))
               grouplist(1)=0
               do i=1,groupsize-1
                  grouplist(i+1)=i+(numprocs/2)-1+oddnumproc
               enddo
               inp1=.false.
               do i=1,groupsize
                  if(rank.eq.grouplist(i)) then
                     inp1=.true.
                     exit
                  endif
               enddo
               call mstm_mpi(mpi_command='incl', &
                  mpi_group=mpigroup,&
                  mpi_size=groupsize,&
                  mpi_new_group_list=grouplist,&
                  mpi_new_group=syncgroup)
               call mstm_mpi(mpi_command='create',&
                  mpi_group=syncgroup,&
                  mpi_comm=mpicomm,&
                  mpi_new_comm=synccomm1)
               deallocate(grouplist)
               groupsize=groupsize+oddnumproc
               allocate(grouplist(groupsize))
               grouplist(1)=numprocs/2+oddnumproc
               do i=1,groupsize-1
                  grouplist(i+1)=i-1
               enddo
               inp2=.false.
               do i=1,groupsize
                  if(rank.eq.grouplist(i)) then
                     inp2=.true.
                     exit
                  endif
               enddo
               call mstm_mpi(mpi_command='incl', &
                  mpi_group=mpigroup,&
                  mpi_size=groupsize,&
                  mpi_new_group_list=grouplist,&
                  mpi_new_group=syncgroup)
               call mstm_mpi(mpi_command='create',&
                  mpi_group=syncgroup,&
                  mpi_comm=mpicomm,&
                  mpi_new_comm=synccomm2)
               deallocate(grouplist)
            else
               pcomm=mpicomm
            endif
         endif

         iter=0
         errmax=0.
         if(normalize_solution_error) then
            enorm=dot_product(pnp,pnp)
         else
            enorm=1.d0
         endif
         errmin=1.d10
!
!  setting niter < 0 runs the following simple order--of--scattering solution
!
         if(niter.lt.0) then
            cp=anp
            time1=mstm_mpi_wtime()
            do iter=1,-niter
               cr=0.
               if(iter.eq.1) then
                  call sphereinteraction(neqns,1,cp,cr, &
                       initial_run=initialize, &
                       mpi_comm=mpicomm)
               else
                  call sphereinteraction(neqns,1,cp,cr, &
                       mpi_comm=mpicomm)
               endif
               anp=anp+cr
               cp=cr
               eerr=dot_product(cr,cr)
               eerr=eerr/enorm
               if(eerr.lt.eps) then
                  errmax=eerr
                  exit
               endif
               time2=mstm_mpi_wtime()
               if(rank0.eq.0.and.iterwrite.eq.1.and.time2-time1.gt.5.d0) then
                  write(iunit,'('' iter,err:'',i5,e13.5)') iter,eerr
                  call flush(iunit)
                  time1=time2
               endif
            enddo
            call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
            deallocate(cr,cp,cw,cq,cap,caw,capt,cawt)
            return
         endif
!
! the following is the implementation of the complex biconjugate gradient
! iteration scheme
!
         cr=0.d0
         call sphereinteraction(neqns,1,anp,cr,initial_run=initialize, &
              mpi_comm=pcomm)
         cr=pnp-anp+cr
         cq=conjg(cr)
         cw=cq
         cp=cr
         csk=dot_product(conjg(cr),cr)
         if(cdabs(csk).eq.0.d0) then
            deallocate(cr,cp,cw,cq,cap,caw,capt,cawt)
            return
         endif
!
!  here starts the main iteration loop
!
         if(rank0.eq.0) time1=mstm_mpi_wtime()
if(light_up) then
write(*,'('' s8.2.3.2 '',i3)') mstm_global_rank
call flush(6)
endif
         do iter=1,niter
            call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
            if(rank.eq.0) time0=mstm_mpi_wtime()
            cak=0.d0
            cawt=0.d0
            capt=0.d0
            cap=0.d0
            caw=0.d0

if(light_up) then
write(*,'('' s8.2.3.3 '',3i3)') mstm_global_rank,iter
call flush(6)
endif
            if(numprocs.eq.1) then
               allocate(ctin(neqns,2),ctout(neqns,2))
               ctin(:,1)=cp(:)
               ctin(:,2)=cw(:)
               contran2=(/.false.,.true./)
               ctout=0.d0
               call sphereinteraction(neqns,2,ctin,ctout, &
                     con_tran=contran2,mpi_comm=mpicomm)
               cap(:)=ctout(:,1)
               caw(:)=ctout(:,2)
               deallocate(ctin,ctout)
            else
               if(pgroup.eq.1) then
                  call sphereinteraction(neqns,1,cp,cap, &
                        mpi_comm=pcomm)
               else
                  call sphereinteraction(neqns,1,cw,caw, &
                        con_tran=(/.true./),mpi_comm=pcomm)
               endif
               call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)
               if(inp2) then
                  call mstm_mpi(mpi_command='bcast', &
                     mpi_send_buf_dc=caw, &
                     mpi_number=neqns, &
                     mpi_rank=0, &
                     mpi_comm=synccomm2)
               endif
               if(inp1) then
                  call mstm_mpi(mpi_command='bcast', &
                     mpi_send_buf_dc=cap, &
                     mpi_number=neqns, &
                     mpi_rank=0, &
                     mpi_comm=synccomm1)
               endif
            endif

            call mstm_mpi(mpi_command='barrier',mpi_comm=mpicomm)

if(light_up) then
write(*,'('' s8.2.3.4 '',3i3)') mstm_global_rank,iter
call flush(6)
endif
            cap=cp-cap
            caw=cw-caw
            cak=dot_product(cw,cap)

            if(cdabs(cak).ne.0.d0) then
               cak=csk/cak
            else
               deallocate(cr,cp,cw,cq,cap,caw,capt,cawt)
               return
            endif

            anp=anp+cak*cp
            cr=cr-cak*cap
            cq=cq-conjg(cak)*caw
            csk2=dot_product(cq,cr)
            eerr=dot_product(cr,cr)/enorm
            errmax=eerr
            errmin=min(errmin,errmax)

            if(eerr.lt.eps.or.cdabs(csk).eq.0.d0) then
               deallocate(cr,cp,cw,cq,cap,caw,capt,cawt)
               return
            else
               cbk=csk2/csk
               csk=csk2
            endif

            cp=cr+cbk*cp
            cw=cq+conjg(cbk)*cw


if(light_up) then
write(*,'('' s8.2.3.5 '',3i3)') mstm_global_rank,iter
call flush(6)
endif
            if(rank.eq.0) time2=mstm_mpi_wtime()
            if(rank0.eq.0.and.iterwrite.eq.1.and.time2-time1.gt.5.d0) then
               write(run_print_unit,'('' iter,err,min err, tpi:'',i5,2e12.4,e12.4)') &
                     iter,errmax,errmin,time2-time0
               call flush(iunit)
               time1=time2
            endif
         enddo
         deallocate(cr,cp,cw,cq,cap,caw,capt,cawt)
         end subroutine cbicg

      end module solver
      module random_sphere_configuration
      use mpidefs
      implicit none
      type l_list
         integer :: index
         type(l_list), pointer :: next
      end type l_list
      type c_list
         integer :: number_elements
         type(l_list), pointer :: members
      end type c_list
      logical  :: target_width_specified
      logical, target :: sphere_1_fixed,periodic_bc(3)
      integer, private :: cell_dim(3)
      integer, allocatable :: sphere_cell(:,:)
      integer, target :: target_shape,wall_boundary_model,max_number_time_steps
      real(8), private :: pi,fv_crit,time_step
      real(8), private :: minimum_gap,d_cell,target_boundaries(3,2)
      real(8), target :: target_dimensions(3),psd_sigma,target_width,target_thickness,max_colls_per_sphere
      type(c_list), allocatable :: cell_list(:,:,:)
      character*1 :: c_temp
      data pi,fv_crit,time_step/3.1415926535897932385d0,0.25d0,.1d0/
      data minimum_gap,sphere_1_fixed,target_shape,psd_sigma/1.0d-3,.false.,0,0.d0/
      data periodic_bc/.true.,.true.,.true./
      data wall_boundary_model/1/
      data max_number_time_steps,max_colls_per_sphere/1000,20.d0/

      contains

         subroutine random_cluster_of_spheres(numberspheres,sphereposition,sphereradius,iunit,istatus, &
            ntsteps,skip_diffusion,use_saved_values)
         implicit none
         logical :: fitok,allin,initial0,initial1,trystage1,skipdif
         logical, optional :: skip_diffusion,use_saved_values
         logical, save :: firstrun
         integer :: i,j,maxsamp0,maxsamp1,numberspheres,ncolls,ncollstot,maxns,ntsteps,istatus,iunit,rank
         real(8) ::samppos(3),sphereposition(3,numberspheres),sphereradius(numberspheres), &
            spherevol,targetfv,u(3,numberspheres),wallboundaries(3,2),targetvol, &
            targetstretch,collspersphere,mfp,time0,time1,sum1,sum2,sdev,mean
         real(8), allocatable, save :: saved_sphereradius(:),saved_sphereposition(:,:)
         data maxsamp0,maxsamp1,firstrun/10000,100,.true./
         if(present(use_saved_values)) then
            if(use_saved_values) then
               sphereposition=saved_sphereposition
               sphereradius=saved_sphereradius
               return
            endif
         endif
         if(present(skip_diffusion)) then
            skipdif=skip_diffusion
         else
            skipdif=.false.
         endif
         if(firstrun) then
            call random_seed()
            firstrun=.false.
         endif
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         trystage1=.false.
         if(psd_sigma.eq.0.d0) then
            sphereradius=1.d0
            spherevol=dble(numberspheres)*4.d0*pi/3.d0
         else
            spherevol=0.d0
            do i=1,numberspheres
               call psdsamp(psd_sigma,2.5d0,sphereradius(i))
               spherevol=spherevol+4.d0*pi/3.d0*sphereradius(i)**3
            enddo
            if(psd_sigma.gt.0.1d0) then
               call sort_radii(numberspheres,sphereradius)
               trystage1=.true.
            endif
         endif
         call target_volume(targetvol)
         targetfv=spherevol/targetvol
         targetstretch=(1.d0/targetfv)**(1.d0/3.d0)
         targetstretch=max(targetstretch,1.02d0)
         mfp=targetvol/dble(numberspheres)/4.d0
         target_boundaries(:,1)=-target_dimensions
         target_boundaries(:,2)=target_dimensions
         wallboundaries=target_boundaries
         d_cell=2.5d0*maxval(sphereradius(1:numberspheres))
         allin=.false.
         istatus=3
         if(targetfv.le.0.25d0.or.trystage1) then
            allin=.true.
            call initialize_cells(numberspheres)
            do i=1,numberspheres
               do j=1,maxsamp0
                  call sample_position(samppos,sphereradius(i))
                  if(sphere_1_fixed.and.i.eq.1) samppos=0.d0
                  call add_sphere_to_cluster(sphereradius(i),samppos,i-1,sphereradius,sphereposition,fitok)
                  if(fitok) then
                     sphereposition(:,i)=samppos(:)
                     exit
                  endif
               enddo
               if(j.ge.maxsamp0) then
                  allin=.false.
                  exit
               endif
            enddo
            if(allin) then
               ntsteps=min(ceiling(2.d0/time_step),max_number_time_steps)
!               write(iunit,'('' target configuration computed using random sampling'')')
               istatus=0
            else
               call clear_cells()
            endif
         endif
         if(targetfv.lt.0.6d0.and.(.not.allin)) then
            allin=.true.
            sum1=0.
            sum2=0.
            do j=1,maxsamp1
               call initialize_cells(numberspheres)
               call layered_sample(numberspheres,sphereradius,sphereposition,wallboundaries,maxns)
               if(maxns.ge.numberspheres) exit
               sum1=sum1+maxns
               sum2=sum2+maxns*maxns
               mean=sum1/dble(j)
               sdev=sqrt(dble(j)*sum2-sum1*sum1)/dble(j)
!if(rank.eq.0) then
!write(*,'(3i10,es12.4)') j,maxns,numberspheres,2.d0*sdev+mean
!call flush(6)
!endif
               if(j.gt.20.and.2.d0*sdev+mean.lt.numberspheres) exit
               call clear_cells()
            enddo
            if(maxns.lt.numberspheres) then
               allin=.false.
            else
               ntsteps=min(ceiling(mfp/time_step),max_number_time_steps)
!               write(iunit,'('' target configuration computed using layered sampling + diffusion, time steps:'',i5)') ntsteps
               istatus=1
            endif
         endif
         if(.not.allin) then
            do
               call initialize_cells(numberspheres)
               call hex_position_generator(numberspheres,sphereradius,sphereposition,wallboundaries,targetstretch,allin,maxns)
!if(rank.eq.0) then
!write(*,'(i10,es12.4)') maxns,targetstretch
!call flush(6)
!endif

               if(allin) exit
               call clear_cells()
               if(targetstretch.le.1.02d0) then
                  write(iunit,'('' MC configuration sampler failed'')')
                  istatus=3
                  return
               endif
               targetstretch=targetstretch-0.001
               targetstretch=max(targetstretch,1.02)
            enddo
            istatus=2
            ntsteps=max_number_time_steps
!            write(iunit,'('' target configuration computed initial HCP + diffusion, time steps:'',i5)') ntsteps
         endif
         do i=1,numberspheres
            call check_in_target(sphereradius(i),sphereposition(:,i),wallboundaries,allin)
!            if(.not.allin) write(iunit,'('' initially outside:'',i5,3es12.4)') i,sphereposition(:,i)
         enddo
!         ntsteps=max_number_time_steps
         if(ntsteps.gt.0..and.(.not.skipdif)) then
            call samptrajectory(numberspheres,u)
            ncollstot=0
            do j=1,ntsteps

               call spheremove(numberspheres,sphereradius,sphereposition,u,time_step,wallboundaries, &
                  number_wall_hits=ncolls)
               ncollstot=ncollstot+ncolls
               collspersphere=dble(ncollstot)/dble(numberspheres)
               if(collspersphere.gt.max_colls_per_sphere) exit
            enddo
            ntsteps=min(ntsteps,j)
         endif
         do i=1,numberspheres
            call check_in_target(sphereradius(i),sphereposition(:,i),wallboundaries,allin)
            if(.not.allin) write(iunit,'('' outside:'',i5,3es12.4)') i,sphereposition(:,i)
         enddo
         if(target_shape.eq.0.or.target_shape.eq.1) then
            call sort_positions(numberspheres,sphereradius,sphereposition,3)
         else
            call sort_positions(numberspheres,sphereradius,sphereposition,0)
         endif
         if(allocated(saved_sphereposition)) then
            deallocate(saved_sphereposition,saved_sphereradius)
         endif
         allocate(saved_sphereposition(3,numberspheres),saved_sphereradius(numberspheres))
         saved_sphereposition=sphereposition
         saved_sphereradius=sphereradius
         call clear_cells()
         end subroutine random_cluster_of_spheres

         subroutine direct_overlap_test(nsphere,radius,position,overlap,distance,pair)
         implicit none
         logical :: overlap
         integer :: nsphere,i,j
         integer, optional :: pair(2)
         real(8) :: radius(nsphere),position(3,nsphere),rij
         real(8), optional :: distance
         overlap=.false.
         do i=1,nsphere-1
            do j=i+1,nsphere
               rij=sqrt(sum((position(:,i)-position(:,j))**2))
               if(rij.lt.radius(i)+radius(j)) then
                  overlap=.true.
                  if(present(distance)) distance=rij
                  if(present(pair)) pair=(/i,j/)
                  return
               endif
            enddo
         enddo
         end subroutine direct_overlap_test

         subroutine target_volume(targetvol)
         implicit none
         integer :: i,ipbc(3)
         real(8) :: targetvol
         ipbc=wall_boundary_model
         do i=1,3
            if(periodic_bc(i)) ipbc(i)=0
         enddo
         if(target_shape.eq.0) then
            targetvol=8.d0*product(target_dimensions(1:3)-dble(ipbc))
         elseif(target_shape.eq.1) then
            targetvol=2.d0*pi*((target_dimensions(1)-1)**2)*(target_dimensions(3)-ipbc(3))
         else
            targetvol=4.d0*pi*(target_dimensions(1)-1)**3/3.d0
         endif
         end subroutine target_volume

         subroutine cell_index(pos,cell)
         implicit none
         integer :: cell(3)
         real(8) :: pos(3)
         cell=floor((pos(:)-target_boundaries(:,1))/(target_boundaries(:,2)-target_boundaries(:,1))*dble(cell_dim(:)))+1
         cell=max(cell,(/1,1,1/))
         cell=min(cell,cell_dim)
         end subroutine cell_index

         subroutine sample_position(pos,rad)
         implicit none
         integer :: i
         real(8) :: pos(3),rannum(3),r,phi,ct,st,rad,wshift(3)

         call random_number(rannum)
         if(target_shape.eq.0) then
            do i=1,3
               if(periodic_bc(i)) then
                  wshift(i)=0.d0
               else
                  wshift(i)=rad*wall_boundary_model+minimum_gap
               endif
            enddo
            pos=target_boundaries(:,1)+wshift(:)+(target_boundaries(:,2)-target_boundaries(:,1)-2.d0*wshift(:))*rannum(:)
         elseif(target_shape.eq.1) then
            wshift(1)=rad*wall_boundary_model+minimum_gap
            if(periodic_bc(3)) then
               wshift(3)=0.d0
            else
               wshift(3)=rad*wall_boundary_model+minimum_gap
            endif
            r=(target_boundaries(1,2)-wshift(1))*rannum(1)**0.5d0
            phi=2.d0*pi*rannum(2)
            pos(1)=r*cos(phi)
            pos(2)=r*sin(phi)
            pos(3)=target_boundaries(3,1)+wshift(3)+(target_boundaries(3,2)-target_boundaries(3,1)-2.d0*wshift(3))*rannum(3)
         else
            wshift(1)=rad*wall_boundary_model+minimum_gap
            r=(target_boundaries(1,2)-wshift(1))*rannum(1)**0.333333d0
            phi=2.d0*pi*rannum(2)
            ct=-1.d0+2.d0*rannum(3)
            st=sqrt(1.d0-ct*ct)
            pos(1)=r*st*cos(phi)
            pos(2)=r*st*sin(phi)
            pos(3)=r*ct
         endif
         end subroutine sample_position

         subroutine clear_cells()
         implicit none
         integer :: n,ix,iy,iz,i
         type(l_list), pointer :: llist,llist2
         if(allocated(cell_list)) then
            do iz=1,cell_dim(3)
               do iy=1,cell_dim(2)
                  do ix=1,cell_dim(1)
                     n=cell_list(ix,iy,iz)%number_elements
                     if(.not.associated(cell_list(ix,iy,iz)%members)) cycle
                     llist=>cell_list(ix,iy,iz)%members
                     do i=1,n
                        llist2=>llist%next
                        deallocate(llist)
                        nullify(llist)
                        if(.not.associated(llist2)) exit
                        llist=>llist2
                     enddo
                  enddo
               enddo
            enddo
            deallocate(cell_list)
         endif
         if(allocated(sphere_cell)) deallocate(sphere_cell)
         end subroutine clear_cells

         subroutine initialize_cells(nsphere)
         implicit none
         integer :: nsphere
         if(allocated(sphere_cell)) deallocate(sphere_cell)
         allocate(sphere_cell(3,nsphere))
         sphere_cell(:,:)=0
         cell_dim(:)=floor((target_boundaries(:,2)-target_boundaries(:,1)-1.d-6)/d_cell)+1
         if(allocated(cell_list)) deallocate(cell_list)
         allocate(cell_list(cell_dim(1),cell_dim(2),cell_dim(3)))
         cell_list(:,:,:)%number_elements=0
         end subroutine initialize_cells

         subroutine swap_cell_contents(i,newcell)
         implicit none
         integer :: i,newcell(3),cell(3),n,l
         type(l_list), pointer :: llist,llist2,llistnew
         cell(:) = sphere_cell(:,i)
         n=cell_list(cell(1),cell(2),cell(3))%number_elements
         if(cell_list(cell(1),cell(2),cell(3))%members%index.eq.i) then
            llist2=>cell_list(cell(1),cell(2),cell(3))%members%next
            llistnew=>cell_list(cell(1),cell(2),cell(3))%members
            cell_list(cell(1),cell(2),cell(3))%members=>llist2
         else
            llist=>cell_list(cell(1),cell(2),cell(3))%members
            do l=1,n-1
               if(llist%next%index.eq.i) then
                  llist2=>llist%next%next
                  llistnew=>llist%next
                  llist%next=>llist2
                  exit
               endif
               llist=>llist%next
            enddo
         endif
         cell_list(cell(1),cell(2),cell(3))%number_elements=n-1
         cell=newcell
         n=cell_list(cell(1),cell(2),cell(3))%number_elements
         llist2=>cell_list(cell(1),cell(2),cell(3))%members
         cell_list(cell(1),cell(2),cell(3))%members=>llistnew
         cell_list(cell(1),cell(2),cell(3))%members%next=>llist2
         cell_list(cell(1),cell(2),cell(3))%number_elements=n+1
         sphere_cell(:,i)=cell
         end subroutine swap_cell_contents

         subroutine modify_cells(nsphere,position,start_sphere,end_sphere)
         implicit none
         integer :: isphere,nsphere,cell(3),istart,iend
         integer, optional :: start_sphere,end_sphere
         real(8) :: position(3,nsphere)
         if(present(start_sphere)) then
            istart=start_sphere
         else
            istart=1
         endif
         if(present(end_sphere)) then
            iend=end_sphere
         else
            iend=nsphere
         endif
         do isphere=istart,iend
            call cell_index(position(:,isphere),cell)
            if(any(sphere_cell(:,isphere).ne.cell)) then
               call swap_cell_contents(isphere,cell)
            endif
         enddo
         end subroutine modify_cells

         subroutine target_distribution_stats(nsphere,sdev)
         implicit none
         integer :: nsphere,n,iz,iy,ix,nt,ncell
         real(8) :: sdev,nmean
         sdev=0.d0
         ncell=product(cell_dim)
         nmean=dble(nsphere)/dble(ncell)
         nt=0
         do iz=1,cell_dim(3)
            do iy=1,cell_dim(2)
               do ix=1,cell_dim(1)
                  n=cell_list(ix,iy,iz)%number_elements
                  nt=nt+n
                  sdev=sdev+(dble(n)/nmean-1.d0)*(dble(n)/nmean-1.d0)
               enddo
            enddo
         enddo
         sdev=sqrt(sdev)
         end subroutine target_distribution_stats

         subroutine add_sphere_to_cluster(newrad,newpos,nsphere,radius,position,fitok)
         implicit none
         logical :: fitok,bndok,pbc(3)
         integer :: nsphere,i,cell(3),n,m,ccell(3),scell(3),j
         real(8) :: radius(*),position(3,*),newrad,newpos(3),rij,tpos(3)
         type(l_list), pointer :: llist
         pbc=.false.
         if(target_shape.eq.0) then
            pbc=periodic_bc
         elseif(target_shape.eq.1) then
            pbc(3)=periodic_bc(3)
         endif
         call cell_index(newpos,ccell)
         fitok=.true.
         do m=0,26
            scell(1)=mod(m,3)-1
            scell(2)=mod(m/3,3)-1
            scell(3)=mod(m/9,3)-1
            cell=ccell+scell
            bndok=.true.
            tpos=newpos
            do i=1,3
               if(cell(i).lt.1.or.cell(i).gt.cell_dim(i)) then
                  if(pbc(i)) then
                     if(cell(i).lt.1) then
                        cell(i)=cell_dim(i)
                        tpos(i)=tpos(i)+target_boundaries(i,2)-target_boundaries(i,1)
                     elseif(cell(i).gt.cell_dim(i)) then
                        cell(i)=1
                        tpos(i)=tpos(i)-target_boundaries(i,2)+target_boundaries(i,1)
                     endif
                  else
                     bndok=.false.
                     exit
                  endif
               endif
            enddo
            if(.not.bndok) cycle
            n=cell_list(cell(1),cell(2),cell(3))%number_elements
            if(n.eq.0) cycle
            llist=>cell_list(cell(1),cell(2),cell(3))%members
            do j=1,n
               i=llist%index
               rij=sqrt(sum((tpos(:)-position(:,i))**2))
               if(rij.lt.newrad+radius(i)+minimum_gap) then
                  fitok=.false.
                  return
               endif
               if(j.lt.n) llist=>llist%next
            enddo
         enddo
         sphere_cell(:,nsphere+1)=ccell(:)
         n=cell_list(ccell(1),ccell(2),ccell(3))%number_elements
         if(n.eq.0) allocate(cell_list(ccell(1),ccell(2),ccell(3))%members)
         llist=>cell_list(ccell(1),ccell(2),ccell(3))%members
         do i=1,n
            if(i.eq.n) allocate(llist%next)
            llist=>llist%next
         enddo
         llist%index=nsphere+1
         cell_list(ccell(1),ccell(2),ccell(3))%number_elements=n+1
         end subroutine add_sphere_to_cluster

         subroutine sort_positions(nsphere,radius,position,sort_elem,make_positive)
         implicit none
         logical :: makepos
         logical, optional :: make_positive
         integer :: nsphere,i,ind(nsphere),selem
         integer, optional :: sort_elem
         real(8) :: radius(nsphere),position(3,nsphere),r(nsphere), &
                    tpos(3,nsphere)
         if(present(sort_elem)) then
            selem=sort_elem
         else
            selem=0
         endif
         if(present(make_positive)) then
            makepos=make_positive
         else
            makepos=.false.
         endif
         if(selem.eq.0) then
            r(:)=sqrt(sum(position(:,:)**2,1))
         else
            if(makepos) then
               r(:)=abs(position(selem,:))
            else
               r(:)=position(selem,:)
            endif
         endif
         ind(1)=0
         call hpsort_eps_epw (nsphere, r, ind, 1.d-15)
         r=radius
         tpos=position
         do i=1,nsphere
            radius(i)=r(ind(i))
            position(:,i)=tpos(:,ind(i))
         enddo
         end subroutine sort_positions

         subroutine sort_radii(nsphere,radius)
         implicit none
         integer :: nsphere,ind(nsphere)
         real(8) :: radius(nsphere)
         radius=-radius
         ind(1)=0
         call hpsort_eps_epw (nsphere, radius, ind, 1.d-15)
         radius=-radius
         end subroutine sort_radii

         subroutine circumscribing_sphere(nsphere,radius,position,rcell)
         implicit none
         integer :: nsphere,i
         real(8) :: radius(nsphere),position(3,nsphere),ri,rcell,mpos(3)
         rcell=0.d0
         mpos=sum(position(:,:),2)/dble(nsphere)
         do i=1,nsphere
            ri=sqrt(sum((position(:,i))**2))+radius(i)
            rcell=max(rcell,ri)
         enddo
         end subroutine circumscribing_sphere

         subroutine check_in_target(rad,pos,wallbound,intarget)
         implicit none
         logical :: intarget
         integer :: i
         real(8) :: rad,pos(3),wallbound(3,2),rho,wrad
         intarget=.true.
         wrad=rad*wall_boundary_model
         if(target_shape.eq.0) then
            do i=1,3
               if(periodic_bc(i)) then
                  intarget=(pos(i).ge.wallbound(i,1).and.pos(i).le.wallbound(i,2))
               else
                  intarget=(pos(i)-wrad.ge.wallbound(i,1).and.pos(i)+wrad.le.wallbound(i,2))
               endif
               if(.not.intarget) return
            enddo
         elseif(target_shape.eq.1) then
            rho=sqrt(sum(pos(1:2)**2))
            if(rho+wrad.ge.wallbound(1,2)) then
               intarget=.false.
               return
            endif
            if(periodic_bc(3)) then
               intarget=(pos(3).ge.wallbound(3,1).and.pos(3).le.wallbound(3,2))
            else
               intarget=(pos(3)-wrad.ge.wallbound(3,1).and.pos(3)+wrad.le.wallbound(3,2))
            endif
            if(.not.intarget) return
         else
            rho=sqrt(sum(pos(1:3)**2))
            if(rho+wrad.gt.wallbound(1,2)) then
               intarget=.false.
               return
            endif
         endif
         end subroutine check_in_target

         subroutine layered_sample(nsphere,rad,pos,wallbound,nin)
         implicit none
         logical :: fitok,pbc(3)
         integer :: nsphere,nin,i,m,maxsamp
         real(8) :: rad(nsphere),pos(3,nsphere),wallbound(3,2),r2,vtot,delv,wdist(3), &
            dz,z1,z2,samp(3),rho,phi,r1,r,st,ct,rannum(3)
         data maxsamp/5000/
         if(target_shape.eq.0) then
            pbc=periodic_bc
         elseif(target_shape.eq.1) then
            pbc(1:2)=.false.
            pbc(3)=periodic_bc(3)
         else
            pbc=.false.
         endif
         if(target_shape.eq.2) then
            r2=0.d0
            vtot=0.d0
            delv=4.d0*pi*(wallbound(1,2)-dble(wall_boundary_model)-minimum_gap)**3/3.d0/dble(nsphere)
         endif
         nin=0
         do i=1,nsphere
            wdist=0.d0
            do m=1,3
               if(.not.pbc(m)) wdist(m)=rad(i)*wall_boundary_model+minimum_gap
            enddo
            if(target_shape.eq.0) then
               dz=(wallbound(3,2)-wallbound(3,1)-2.d0*wdist(3))/dble(nsphere)
               z1=wallbound(3,1)+wdist(3)+dble(i-1)*dz
               z2=z1+dz
               z2=min(z2,wallbound(3,2)-wdist(3))
               do m=1,maxsamp
                  call random_number(rannum)
                  samp=(/wallbound(1,1)+wdist(1),wallbound(2,1)+wdist(2),z1/) &
                     +((/wallbound(1,2)-wdist(1),wallbound(2,2)-wdist(2),z2/) &
                     -(/wallbound(1,1)+wdist(1),wallbound(2,1)+wdist(2),z1/))*rannum
                  call check_in_target(rad(i)*wall_boundary_model,samp,wallbound,fitok)
                  if(.not.fitok) cycle
                  call add_sphere_to_cluster(rad(i),samp,i-1,rad,pos,fitok)
                  if(fitok) then
                     pos(:,i)=samp(:)
                     exit
                  endif
               enddo
               if(.not.fitok) return
            elseif(target_shape.eq.1) then
               dz=(wallbound(3,2)-wallbound(3,1)-2.d0*wdist(3))/dble(nsphere)
               z1=wallbound(3,1)+wdist(3)+dble(i-1)*dz
               z2=z1+dz
               do m=1,maxsamp
                  call random_number(rannum)
                  rho=(wallbound(1,2)-wdist(1))*sqrt(rannum(1))
                  phi=2.d0*pi*rannum(2)
                  samp(1)=rho*cos(phi)
                  samp(2)=rho*sin(phi)
                  samp(3)=z1+dz*rannum(3)
                  call check_in_target(rad(i)*wall_boundary_model,samp,wallbound,fitok)
                  if(.not.fitok) cycle
                  call add_sphere_to_cluster(rad(i),samp,i-1,rad,pos,fitok)
                  if(fitok) then
                     pos(:,i)=samp(:)
                     exit
                  endif
               enddo
               if(.not.fitok) return
            else
               r1=r2
               vtot=vtot+delv
               r2=(3.d0*vtot/4.d0/pi)**(1.d0/3.d0)
               do m=1,maxsamp
                  if(sphere_1_fixed.and.i.eq.1) then
                     pos=0.d0
                     fitok=.true.
                  else
                     call random_number(rannum)
                     if((wallbound(1,2)-wdist(1)).le.r2) then
                        r2=wallbound(1,2)-wdist(1)
                        r2=max(r2,0.d0)
                     endif
                     r=(3.d0*delv/4.d0/pi*rannum(1)+r1**3)**(1.d0/3.d0)
                     ct=-1.d0+2.d0*rannum(2)
                     st=sqrt(1.d0-ct*ct)
                     phi=2.d0*pi*rannum(3)
                     samp(1)=r*st*cos(phi)
                     samp(2)=r*st*sin(phi)
                     samp(3)=r*ct
                     call check_in_target(rad(i)*wall_boundary_model,samp,wallbound,fitok)
                  endif
                  if(.not.fitok) cycle
                  call add_sphere_to_cluster(rad(i),samp,i-1,rad,pos,fitok)
                  if(fitok) then
                     pos(:,i)=samp(:)
                     exit
                  endif
               enddo
               if(.not.fitok) return
            endif
            nin=i
         enddo
         end subroutine layered_sample

         subroutine hex_position_generator(nsphere,rad,pos,wallbound,s,allin,ns)
         implicit none
         logical :: intarget,fitok,allin
         integer :: nsphere,i,l,m,n,ns,l0,imax,m0,i2,n2,m2,i21,l1,n0,ns0
         real(8) :: rad(nsphere),pos(3,nsphere),wallbound(3,2),s,cscale(3),tpos(3), &
                  tpos1(3),tpos2(3),trad,cscale2(3)
         data imax/200/
         cscale = (/2.d0, sqrt(3.d0), sqrt(8.d0/3.d0)/)
         cscale2=cscale*cscale
         ns=0
         allin=.true.
         do i=0,imax
            if(mod(i,2).eq.0) ns0=ns
            i2=i*i
            i21=(i+1)*(i+1)
            n0=ceiling(dble(i)/cscale(3))
            do n=-n0-1,n0+1
               n2=n*n
               m0=ceiling(sqrt(max(dble(i2)-cscale2(3)*dble(n2),0.d0))/cscale(2))
               do m=-m0-1,m0+1
                  m2=m*m
                  l0=floor(sqrt(max(dble(i2)-cscale2(3)*dble(n2)-cscale2(2)*dble(m2),0.d0))/cscale(1))
                  l0=max(l0,1)
                  l1=ceiling(sqrt(max(dble(i21)-cscale2(3)*dble(n2)-cscale2(2)*dble(m2),0.d0))/cscale(1))
                  tpos1=s*(/dble(mod(abs(m+n),2)),dble(mod(abs(n),2))/cscale(2),0.d0/)
                  do l=l0-1,l1+1
                     tpos2=s*cscale*(/dble(l),dble(m),dble(n)/)
                     tpos=tpos1+tpos2
                     trad=sqrt(sum(tpos**2))/s
                     if(trad.ge.dble(i).and.trad.lt.dble(i+1)) then
                        trad=rad(ns+1)*wall_boundary_model
                        intarget=.true.
                        call check_in_target(trad,tpos,wallbound,intarget)
                        if(intarget) then
                           call add_sphere_to_cluster(rad(ns+1),tpos,ns,rad,pos,fitok)
                           if(fitok) then
                              ns=ns+1
                              pos(:,ns)=tpos(:)
                              if(ns.eq.nsphere) return
                           endif
                        endif
                     endif
                     if(l.eq.0) cycle
                     tpos2=s*cscale*(/-dble(l),dble(m),dble(n)/)
                     tpos=tpos1+tpos2
                     trad=sqrt(sum(tpos**2))/s
                     if(trad.ge.dble(i).and.trad.lt.dble(i+1)) then
                        trad=rad(ns+1)*wall_boundary_model
                        intarget=.true.
                        call check_in_target(trad,tpos,wallbound,intarget)
                        if(intarget) then
                           call add_sphere_to_cluster(rad(ns+1),tpos,ns,rad,pos,fitok)
                           if(fitok) then
                              ns=ns+1
                              pos(:,ns)=tpos(:)
                              if(ns.eq.nsphere) return
                           endif
                        endif
                     endif
                  enddo
               enddo
            enddo
            if(ns.eq.ns0.and.mod(i,2).ne.0) then
               if(ns.lt.nsphere) allin=.false.
               return
            endif
         enddo
         end subroutine hex_position_generator

         subroutine spheremove(nsphere,radius,pos,u,maxtime,wallboundaries,number_wall_hits)
         implicit none
         logical :: collision,wallcollision,pbc(3),intarget
         integer :: nsphere,i,is,js,collisionpair(2),iwall,iswall,m,nwhits
         integer, optional :: number_wall_hits
         real(8) :: pos(3,nsphere),radius(nsphere),maxtime, &
                    tcmin,tmove,u1new(3),u2new(3),u(1:3,nsphere), &
                    twallmin,rho,cp,sp,urho,uphi,tpos(3), &
                    u1pn(3),ct,st,r,tcoll,wallboundaries(3,2),collpos(3)
         pbc=.false.
         if(target_shape.eq.0) then
            pbc=periodic_bc
         elseif(target_shape.eq.1) then
            pbc(3)=periodic_bc(3)
         endif
         tmove=maxtime
         i=1
         nwhits=0
         do while(tmove.gt.0.d0)
            call modify_cells(nsphere,pos)
            call trajectorytest(nsphere,radius,pos,u,tmove,wallboundaries, &
               collision,tcoll,collisionpair,collision_pos=collpos)
            tcmin=tcoll
            call walltest(nsphere,radius,pos,u,tmove,wallboundaries,twallmin,iswall,iwall)
            wallcollision=(twallmin.lt.tcmin)
            tcmin=min(tcmin,twallmin)
            do is=1,nsphere
               if(is.eq.1.and.sphere_1_fixed) cycle
               tpos(1:3)=pos(1:3,is)+u(1:3,is)*tcmin
               do m=1,3
                  if(pbc(m)) then
                     if(tpos(m).ge.wallboundaries(m,2)) then
                        tpos(m)=tpos(m)-(wallboundaries(m,2)-wallboundaries(m,1))
                     elseif(tpos(m).lt.wallboundaries(m,1)) then
                        tpos(m)=tpos(m)+(wallboundaries(m,2)-wallboundaries(m,1))
                     endif
                  endif
               enddo
               call check_in_target(radius(is),tpos,wallboundaries,intarget)
               if(intarget) then
                  pos(:,is)=tpos
               else
                  write(*,'('' out of target'')')
                  write(*,'(8es12.4)') tpos,tcoll,twallmin
               endif
            enddo
            if(tcmin.lt.tmove) then
               nwhits=nwhits+1
               if(wallcollision) then
                  if(target_shape.eq.0) then
                     u(iwall,iswall)=-u(iwall,iswall)
                  elseif(target_shape.eq.1) then
                     if(iwall.le.2) then
                        rho=sqrt(pos(1,iswall)*pos(1,iswall)+pos(2,iswall)*pos(2,iswall))
                        cp=pos(1,iswall)/rho
                        sp=pos(2,iswall)/rho
                        urho=cp*u(1,iswall)+sp*u(2,iswall)
                        uphi=-sp*u(1,iswall)+cp*u(2,iswall)
                        u(1,iswall)=-cp*urho-sp*uphi
                        u(2,iswall)=-sp*urho+cp*uphi
                     else
                        u(iwall,iswall)=-u(iwall,iswall)
                     endif
                  elseif(target_shape.eq.2) then
                     rho=sqrt(pos(1,iswall)*pos(1,iswall)+pos(2,iswall)*pos(2,iswall))
                     if(rho.eq.0.d0) then
                        cp=1.d0
                        sp=0.d0
                     else
                        cp=pos(1,iswall)/rho
                        sp=pos(2,iswall)/rho
                     endif
                     r=sqrt(rho*rho+pos(3,iswall)*pos(3,iswall))
                     if(r.eq.0.d0) then
                        ct=1.d0
                        st=0.d0
                     else
                        ct=pos(3,iswall)/r
                        st=rho/r
                     endif
                     u1pn(1)=(u(1,iswall)*cp+u(2,iswall)*sp)*st+u(3,iswall)*ct
                     u1pn(2)=(u(1,iswall)*cp+u(2,iswall)*sp)*ct-u(3,iswall)*st
                     u1pn(3)=u(1,iswall)*sp-u(2,iswall)*cp
                     u1pn(1)=-u1pn(1)
                     u(1,iswall)=(u1pn(1)*st+u1pn(2)*ct)*cp+u1pn(3)*sp
                     u(2,iswall)=(u1pn(1)*st+u1pn(2)*ct)*sp-u1pn(3)*cp
                     u(3,iswall)=u1pn(1)*ct-u1pn(2)*st
                  endif
               elseif(collision) then
                  is=collisionpair(1)
                  js=collisionpair(2)
                  if(is.eq.1.and.sphere_1_fixed) then
                     call collisiontrajectory(1.d20,pos(1:3,is),u(1:3,is),1.d0, &
                          pos(1:3,js),u(1:3,js),u1new,u2new)
                  else
                     call collisiontrajectory(1.d0,collpos(1:3),u(1:3,is),1.d0, &
                          pos(1:3,js),u(1:3,js),u1new,u2new)
                  endif
                  u(1:3,is)=u1new(1:3)
                  u(1:3,js)=u2new(1:3)
               endif
            endif
            tmove=tmove-abs(tcmin)
            i=i+1
            if(sphere_1_fixed) u(:,1)=0.d0
         enddo
         if(present(number_wall_hits)) number_wall_hits=nwhits
         end subroutine spheremove

         subroutine walltest(nsphere,radius,pos,u,tmove,wallboundaries,twallmin,is,iswall, &
           start_sphere,end_sphere)
         implicit none
         integer :: nsphere,is,iwall,i,iswall,i1,i2
         integer, optional :: start_sphere,end_sphere
         real(8) :: pos(3,nsphere),radius(nsphere),tmove,wallboundaries(3,2), &
                    twall,u(3,nsphere),twallmin,vel,rho,dist,cp,sp,ct,st,r,urho
         if(present(start_sphere)) then
            i1=start_sphere
         else
            i1=1
         endif
         if(present(end_sphere)) then
            i2=end_sphere
         else
            i2=nsphere
         endif
         twallmin=tmove
         if(target_shape.eq.0) then
            do iwall=1,3
               if(periodic_bc(iwall)) cycle
               do i=i1,i2
                  vel=u(iwall,i)
                  if(vel.lt.0.d0) then
                     dist=-pos(iwall,i)+wallboundaries(iwall,1)+radius(i)*wall_boundary_model+minimum_gap
                     twall=dist/vel
                  elseif(vel.gt.0.d0) then
                     dist=wallboundaries(iwall,2)-pos(iwall,i)-radius(i)*wall_boundary_model-minimum_gap
                     twall=dist/vel
                  else
                     twall=1.d6
                  endif
                  if(twall.lt.twallmin) then
                     twallmin=twall
                     is=i
                     iswall=iwall
                  endif
               enddo
            enddo
         elseif(target_shape.eq.1) then
            do iwall=2,3
               if(iwall.eq.3.and.periodic_bc(iwall)) cycle
               do i=i1,i2
                  if(iwall.lt.3) then
                     rho=sqrt(pos(1,i)*pos(1,i)+pos(2,i)*pos(2,i))
                     if(rho.ne.0.d0) then
                        vel=(pos(1,i)*u(1,i)+pos(2,i)*u(2,i))/rho
                     else
                        vel=sqrt(u(1,i)*u(1,i)+u(2,i)*u(2,i))
                     endif
                  else
                     vel=u(iwall,i)
                  endif
                  if(vel.lt.0.d0) then
                     if(iwall.lt.3) then
                        dist=-rho-wallboundaries(1,2)+radius(i)*wall_boundary_model+minimum_gap
                     else
                        dist=-pos(iwall,i)+wallboundaries(iwall,1)+radius(i)*wall_boundary_model+minimum_gap
                     endif
                     twall=dist/vel
                  elseif(vel.gt.0.d0) then
                     if(iwall.lt.3) then
                        dist=wallboundaries(1,2)-rho-radius(i)*wall_boundary_model-minimum_gap
                     else
                        dist=wallboundaries(iwall,2)-pos(iwall,i)-radius(i)*wall_boundary_model-minimum_gap
                     endif
                     twall=dist/vel
                  else
                     twall=1.d6
                  endif
                  if(twall.lt.twallmin) then
                     twallmin=twall
                     is=i
                     iswall=iwall
                  endif
               enddo
            enddo
         elseif(target_shape.eq.2) then
            do i=i1,i2
               rho=sqrt(pos(1,i)*pos(1,i)+pos(2,i)*pos(2,i))
               if(rho.eq.0.d0) then
                  cp=1.d0
                  sp=0.d0
               else
                  cp=pos(1,i)/rho
                  sp=pos(2,i)/rho
               endif
               r=sqrt(rho*rho+pos(3,i)*pos(3,i))
               if(r.eq.0.d0) then
                  vel=sqrt(dot_product(u(:,i),u(:,i)))
               else
                  ct=pos(3,i)/r
                  st=rho/r
                  urho=cp*u(1,i)+sp*u(2,i)
                  vel=urho*st+u(3,i)*ct
               endif
               if(vel.lt.0.d0) then
                  dist=-r-wallboundaries(1,2)+radius(i)*wall_boundary_model+minimum_gap
                  twall=dist/vel
               elseif(vel.gt.0.d0) then
                  dist=wallboundaries(1,2)-r-radius(i)*wall_boundary_model-minimum_gap
                  twall=dist/vel
               else
                  twall=1.d6
               endif
               if(twall.lt.twallmin) then
                  twallmin=twall
                  is=i
                  iswall=3
               endif
            enddo
         endif
         end subroutine walltest

         subroutine trajectorytest(nsphere,radius,pos,u,maxtime,wallboundaries,collision, &
            tcmin,collisionpair,start_sphere,end_sphere,minimum_distance,collision_pos)
         implicit none
         logical :: collision,bndok,loccoll,pbc(3)
         integer :: nsphere,is,j,js,cell(3),ccell(3),scell(3),collisionpair(2),m,n, &
            istart,iend,i
         integer, optional :: start_sphere,end_sphere
         real(8) :: pos(3,nsphere),radius(nsphere),maxtime,wallboundaries(3,2), &
                    tcmin,rcol,tcollision,u(1:3,nsphere),mindist,tpos(3),collisionpos(3)
         real(8), optional :: minimum_distance,collision_pos(3)
         type(l_list), pointer :: llist
         if(present(start_sphere)) then
            istart=start_sphere
         else
            istart=1
         endif
         if(present(end_sphere)) then
            iend=end_sphere
         else
            iend=nsphere
         endif
         if(present(minimum_distance)) then
            mindist=minimum_distance
         else
            mindist=minimum_gap
         endif
         if(target_shape.eq.0) then
            pbc=periodic_bc
         elseif(target_shape.eq.1) then
            pbc(3)=periodic_bc(3)
         endif
         tcmin=maxtime
         collision=.false.
         do is=istart,iend
            call cell_index(pos(:,is),ccell)
            do m=0,26
               scell(1)=mod(m,3)-1
               scell(2)=mod(m/3,3)-1
               scell(3)=mod(m/9,3)-1
               cell=ccell+scell
               bndok=.true.
               tpos=pos(:,is)
               do i=1,3
                  if(cell(i).lt.1.or.cell(i).gt.cell_dim(i)) then
                     if(pbc(i)) then
                        if(cell(i).lt.1) then
                           cell(i)=cell_dim(i)
                           tpos(i)=tpos(i)+wallboundaries(i,2)-wallboundaries(i,1)
                        elseif(cell(i).gt.cell_dim(i)) then
                           cell(i)=1
                           tpos(i)=tpos(i)-wallboundaries(i,2)+wallboundaries(i,1)
                        endif
                     else
                        bndok=.false.
                        exit
                     endif
                  endif
               enddo
               if(.not.bndok) cycle
               n=cell_list(cell(1),cell(2),cell(3))%number_elements
               if(n.eq.0) cycle
               llist=>cell_list(cell(1),cell(2),cell(3))%members
               do j=1,n
                  js=llist%index
                  if(js.ne.is) then
                     rcol=radius(is)+radius(js)+mindist
                     call paircollisiontest(tpos(1:3),u(1:3,is),pos(1:3,js),u(1:3,js), &
                          rcol,loccoll,tcollision)
                     if(loccoll) then
                        if(tcollision.lt.tcmin) then
                           collision=.true.
                           tcmin=tcollision
                           collisionpair(1:2)=(/is,js/)
                           collisionpos=tpos
                        endif
                     endif
                  endif
                  if(j.lt.n) llist=>llist%next
               enddo
            enddo
         enddo
         if(present(collision_pos)) collision_pos=collisionpos
         end subroutine trajectorytest

         subroutine paircollisiontest(pos1,u1,pos2,u2,rcol,collision,tcollision)
         implicit none
         real(8) :: pos1(3),u1(3),pos2(3),u2(3),rcol,tcollision, &
                    urel(3),posrel(3),a,b,c,d
         logical :: collision
         urel=u2-u1
         posrel=pos2-pos1
         b=2.d0*dot_product(urel,posrel)
         if(b.ge.0.d0) then
            collision=.false.
            return
         endif
         a=dot_product(urel,urel)
         c=max(dot_product(posrel,posrel)-rcol*rcol,0.d0)
         if(c.eq.0.d0) then
            collision=.true.
            tcollision=0.d0
            return
         endif
         d=b*b-4.d0*a*c
         if(d.lt.0.d0) then
            collision=.false.
            return
         endif
!         tc1=-(b+sqrt(d))/2.d0/a
!         tc2=-(b-sqrt(d))/2.d0/a
         tcollision=-(b+sqrt(b*b-4.d0*a*c))/2.d0/a
         collision=.true.
         end subroutine paircollisiontest

         subroutine collisiontrajectory(mass1,pos1,u1,mass2,pos2,u2,u1new,u2new)
         implicit none
         real(8) :: mass1,pos1(3),u1(3),mass2,pos2(3),u2(3),u1new(3),u2new(3), &
                    posrel(3),rc,cosb,sinb,alpha,cosa,sina,rotmat(3,3),u1p(3),u2p(3), &
                    u1pn(3),u2pn(3)
         posrel=pos2-pos1
         rc=sqrt(dot_product(posrel,posrel))
         cosb=posrel(3)/rc
         sinb=sqrt((1.d0-cosb)*(1.d0+cosb))
         if(posrel(1).eq.0.d0.and.posrel(2).eq.0.d0) then
            alpha=0.d0
         else
            alpha=datan2(posrel(2),posrel(1))
         endif
         cosa=cos(alpha)
         sina=sin(alpha)
         rotmat=reshape((/cosa*cosb,-sina,cosa*sinb,sina*cosb,cosa,sina*sinb,-sinb,0.d0,cosb/),(/3,3/))
         u1p=matmul(rotmat,u1)
         u2p=matmul(rotmat,u2)
         u1pn(1:2)=u1p(1:2)
         u1pn(3)=((mass1-mass2)*u1p(3)+2.d0*mass2*u2p(3))/(mass1+mass2)
         u2pn(1:2)=u2p(1:2)
         u2pn(3)=((mass2-mass1)*u2p(3)+2.d0*mass1*u1p(3))/(mass1+mass2)
         u1new=matmul(transpose(rotmat),u1pn)
         u2new=matmul(transpose(rotmat),u2pn)
         end subroutine collisiontrajectory

         subroutine samptrajectory(nsphere,u)
         implicit none
         integer :: i,nsphere
         real(8) :: u(3,nsphere),cb,sb,alpha,ca,sa,rannum(2)
         do i=1,nsphere
            call random_number(rannum)
            cb=-1.d0+2.d0*rannum(1)
            sb=sqrt((1.d0-cb)*(1.d0+cb))
            alpha=6.2831853d0*rannum(2)
            ca=cos(alpha)
            sa=sin(alpha)
            u(1:3,i)=(/ca*sb,sa*sb,cb/)
         enddo
         end subroutine samptrajectory

         subroutine psdsamp(sigma,maxradius,x)
         implicit none
         integer :: i
         real(8) :: sigma,maxradius,r2pi,f1,fd,x,fmax,s2,xmax, &
                    t1,rannum(2)
         data r2pi/2.5066282746310002d0/
         if(sigma.eq.0.d0) then
            x=1.d0
            return
         endif
         s2=sigma*sigma
         f1=1.d0
         fd=0.d0
         xmax=exp(-2.5d0*s2)
         t1=(log(xmax)+1.5d0*s2)
         fmax=exp(-t1*t1/(2.d0*s2))/r2pi/xmax/sigma
         i=0
         do while(f1.gt.fd)
            i=i+1
            call random_number(rannum)
            x=maxradius*rannum(1)
            f1=fmax*rannum(2)
            t1=(log(x)+1.5d0*s2)
            fd=exp(-t1*t1/(2.d0*s2))/r2pi/x/sigma
         enddo
         end subroutine psdsamp

!
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Adapted from flib/hpsort_eps
  !---------------------------------------------------------------------
         subroutine hpsort_eps_epw (n, ra, ind, eps)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
!           use kinds, ONLY : DP
           implicit none
           !-input/output variables
           integer, intent(in)   :: n
           real(8), intent(in)  :: eps
           integer :: ind (n)
           real(8) :: ra (n)
           !-local variables
           integer :: i, ir, j, l, iind
           real(8) :: rra
         !
           ! initialize index array
           IF (ind (1) .eq.0) then
              DO i = 1, n
                 ind (i) = i
              ENDDO
           ENDIF
           ! nothing to order
           IF (n.lt.2) return
           ! initialize indices for hiring and retirement-promotion phase
           l = n / 2 + 1

           ir = n

           sorting: do

             ! still in hiring phase
             IF ( l .gt. 1 ) then
                l    = l - 1
                rra  = ra (l)
                iind = ind (l)
                ! in retirement-promotion phase.
             ELSE
                ! clear a space at the end of the array
                rra  = ra (ir)
                !
                iind = ind (ir)
                ! retire the top of the heap into it
                ra (ir) = ra (1)
                !
                ind (ir) = ind (1)
                ! decrease the size of the corporation
                ir = ir - 1
                ! done with the last promotion
                IF ( ir .eq. 1 ) then
                   ! the least competent worker at all !
                   ra (1)  = rra
                   !
                   ind (1) = iind
                   exit sorting
                ENDIF
             ENDIF
             ! wheter in hiring or promotion phase, we
             i = l
             ! set up to place rra in its proper level
             j = l + l
             !
             DO while ( j .le. ir )
                IF ( j .lt. ir ) then
                   ! compare to better underling
                   IF ( hslt( ra (j),  ra (j + 1) ) ) then
                      j = j + 1
                   !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
                      ! this means ra(j) == ra(j+1) within tolerance
                    !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
                   ENDIF
                ENDIF
                ! demote rra
                IF ( hslt( rra, ra (j) ) ) then
                   ra (i) = ra (j)
                   ind (i) = ind (j)
                   i = j
                   j = j + j
                !else if ( .not. hslt ( ra(j) , rra ) ) then
                   !this means rra == ra(j) within tolerance
                   ! demote rra
                  ! if (iind.lt.ind (j) ) then
                  !    ra (i) = ra (j)
                  !    ind (i) = ind (j)
                  !    i = j
                  !    j = j + j
                  ! else
                      ! set j to terminate do-while loop
                  !    j = ir + 1
                  ! endif
                   ! this is the right place for rra
                ELSE
                   ! set j to terminate do-while loop
                   j = ir + 1
                ENDIF
             ENDDO
             ra (i) = rra
             ind (i) = iind

           END DO sorting
         contains

           !  internal function
           !  compare two real number and return the result

           logical function hslt( a, b )
             REAL(8) :: a, b
             IF( abs(a-b) <  eps ) then
               hslt = .false.
             ELSE
               hslt = ( a < b )
             end if
           end function hslt
         end subroutine hpsort_eps_epw

      end module random_sphere_configuration
!june 18 original
!27 july: core volume fraction added to output

      module inputinterface
      use specialfuncs
      use intrinsics
      use mpidefs
      use solver
      use spheredata
      use translation
      use mie
      use nearfield
      use scatprops
      use fft_translation
      use surface_subroutines
      use periodic_lattice_subroutines
      use random_sphere_configuration
      implicit none
      logical :: loop_job,repeat_run,first_run,data_scaled,temporary_pos_file, &
         append_near_field_output_file,incident_beta_specified,number_spheres_specified, &
         square_cell,random_configuration,use_previous_configuration,calculate_up_down_scattering
      logical, target :: append_output_file, &
                 print_scattering_matrix, &
                 copy_input_file,calculate_near_field, &
                 move_to_front,move_to_back,random_orientation, &
                 t_matrix_centered_on_1,calculate_scattering_matrix, &
                 normalize_s11,print_sphere_data,single_origin_expansion, &
                 azimuthal_average,incident_frame,configuration_average, &
                 frozen_configuration,reflection_model,input_fft_translation_option, &
                 print_random_configuration,print_timings, &
                 input_calculate_up_down_scattering,incidence_average,auto_absorption_sample_radius
      logical, allocatable :: sphere_excitation_switch(:)
      integer :: n_nest_loops,i_var_start(5),i_var_stop(5),i_var_step(5), &
                 run_number,loop_sphere_number(5),qeff_dim, &
                 scattering_map_directions,local_rank,scat_mat_ldim,scat_mat_udim,scat_mat_mdim, &
                 ran_config_stat,ran_config_time_steps,n_configuration_groups,random_configuration_number, &
                 solution_iterations,incident_direction_number,number_rl_dirs(2),max_number_rl_dirs
      integer, target :: max_iterations,t_matrix_procs_per_solution, &
         scattering_map_model,scattering_map_dimension,near_field_calculation_model, &
         incident_direction,number_configurations,min_fft_nsphere,input_node_order, &
         number_incident_directions,shifted_sphere
      real(8) :: r_var_start(5),r_var_stop(5),r_var_step(5),diffuse_scattering_ratio, &
          coherent_scattering_ratio,hemispherical_sca(2,2),evan_sca(2),prop_sca(2), &
          input_layer_thickness(max_number_plane_boundaries), &
          pl_sca(2,2),sca,scat_mat_amin,scat_mat_amax,pl_sca_ave(2,2),solution_time,solution_time_ave, &
          incident_beta,solution_error,surface_absorptance(2),surface_absorptance_ave(2), &
          position_shift(3)
      real(8), allocatable :: q_eff(:,:,:),q_vabs(:,:),q_eff_tot(:,:),scat_mat(:,:), &
         dif_scat_mat(:,:),sm_coef(:,:,:),sm_cf_coef(:,:,:),boundary_sca(:,:),boundary_ext(:,:), &
         q_eff_ave(:,:,:),q_vabs_ave(:,:),q_eff_tot_ave(:,:),scat_mat_ave(:,:), &
         boundary_sca_ave(:,:),boundary_ext_ave(:,:),sphere_position_ave(:,:),dif_boundary_sca(:,:), &
         scat_mat_exp_coef(:,:,:),scat_mat_exp_coef_ave(:,:,:),rl_vec(:,:)
      real (8), target :: incident_beta_deg,incident_alpha_deg,solution_epsilon, &
         mie_epsilon,length_scale_factor,near_field_plane_position, &
         near_field_plane_vertices(3,2),near_field_step_size, &
         translation_epsilon,t_matrix_convergence_epsilon, &
         scattering_map_increment,incident_sin_beta,input_cell_width(2), &
         sphere_volume_fraction,input_cell_width_x, &
         input_cell_volume_fraction, &
         excitation_radius,absorption_sample_radius,absorption_sample_radius_fraction,&
         x_shift,y_shift,z_shift
      complex(8) :: c_var_start(5),c_var_stop(5),c_var_step(5),effective_ref_index
      complex(8), target :: ref_index_scale_factor
      complex(8), allocatable :: amnp_s(:,:),amnp_0_ave(:,:),amnp_0(:,:),e_field(:,:,:),e_field_ave(:,:,:)
      character*1 :: loop_var_type(5)
      character*20 :: run_date_and_time
      character*256 :: loop_var_label(5),input_file
      character*256, target :: output_file,run_file,t_matrix_output_file, &
            sphere_data_input_file,near_field_output_file,solution_method
      data loop_job,repeat_run/.false.,.false./
      data append_output_file/.false./
      data copy_input_file/.false./
      data n_nest_loops/0/
      data run_number/0/
      data i_var_start,i_var_stop,i_var_step/5*0,5*0,5*0/
      data r_var_start,r_var_stop,r_var_step/5*0.d0,5*0.d0,5*0.d0/
      data c_var_start,c_var_stop,c_var_step/5*(0.d0,0.d0),5*(0.d0,0.d0),5*(0.d0,0.d0)/
      data max_iterations/10000/
      data incident_beta_deg/0.d0/
      data incident_alpha_deg/0.d0/
      data incident_sin_beta,incident_direction,incident_frame/0.d0,1,.false./
      data solution_epsilon/1.d-6/
      data mie_epsilon/1.d-6/
      data translation_epsilon/1.d-5/
      data t_matrix_convergence_epsilon/1.d-6/
      data output_file/'mstmtest.dat'/
      data run_file/'run1.dat'/
      data sphere_data_input_file/' '/
      data length_scale_factor/1.d0/
      data ref_index_scale_factor/(1.d0,0.d0)/
      data move_to_front,move_to_back/.false.,.false./
      data calculate_near_field/.false./
      data near_field_output_file/'nftest.dat'/
      data append_near_field_output_file/.false./
      data near_field_plane_vertices/-.5d0,0.d0,-.5d0,.5d0,0.d0,.5d0/
      data near_field_step_size/0.2d0/
      data data_scaled/.false./
      data temporary_pos_file/.false./
      data random_orientation/.false./
      data t_matrix_output_file/'tmattemp.dat'/
      data t_matrix_procs_per_solution/4/
      data t_matrix_centered_on_1/.false./
      data calculate_scattering_matrix/.true./
      data solution_method/'iteration'/
      data scattering_map_model/0/
      data scattering_map_dimension/15/
      data scattering_map_increment/1.d0/
      data normalize_s11,print_sphere_data,single_origin_expansion,azimuthal_average/.true.,.true.,.true.,.false./
      data number_spheres_specified,configuration_average/.true.,.false./
      data frozen_configuration,reflection_model,random_configuration_number/.false.,.false.,1/
      data min_fft_nsphere,input_fft_translation_option,input_node_order,input_cell_volume_fraction/200,.false.,-1,0.d0/
      data print_random_configuration,print_timings/.false.,.true./
      data input_calculate_up_down_scattering/.true./
      data incidence_average,number_incident_directions/.false.,16/
      data use_previous_configuration/.false./
      data absorption_sample_radius,excitation_radius/1.d10,1.d10/
      data auto_absorption_sample_radius,absorption_sample_radius_fraction/.true.,0.8d0/
      data x_shift,y_shift,z_shift,shifted_sphere/0.d0,0.d0,0.d0,0/

      contains

         subroutine variable_list_operation(varlabel, &
            var_value,var_type, &
            var_position,var_operation,var_status, &
            i_var_pointer,r_var_pointer,c_var_pointer)
         implicit none
         logical :: operate
         logical, pointer :: lvarvalue,lavarvalue(:)
         integer :: varpos,varstatus,varlen
         integer, optional :: var_position,var_status
         integer, pointer :: ivarvalue
         integer, optional, pointer :: i_var_pointer
         real(8), pointer :: rvarvalue,ravarvalue(:)
         real(8), optional, pointer :: r_var_pointer
         complex(8), pointer :: cvarvalue
         complex(8), optional, pointer :: c_var_pointer
         character*1 :: vartype
         character*1, optional :: var_type
         character*(*), optional :: var_value,var_operation
         character*256 :: varop,sentvarvalue,varlabel
         character*256, pointer :: avarvalue

         if(present(var_operation)) then
            varop=trim(var_operation)
         else
            varop=' '
         endif
         if(present(var_value)) then
            sentvarvalue=trim(var_value)
            operate=.true.
         else
            sentvarvalue=' '
            operate=.false.
         endif
         if(present(var_position)) then
            varpos=var_position
         else
            varpos=1
         endif
         varstatus=0
         vartype='n'
         varlen=1

         if(varlabel.eq.'output_file') then
            vartype='a'
            avarvalue=>output_file

         elseif(varlabel.eq.'append_output_file') then
            vartype='l'
            lvarvalue=>append_output_file

         elseif(varlabel.eq.'copy_input_file') then
            vartype='l'
            lvarvalue=>copy_input_file

         elseif(varlabel.eq.'run_file') then
            vartype='a'
            avarvalue=>run_file

         elseif(varlabel.eq.'sphere_data_input_file') then
            vartype='a'
            avarvalue=>sphere_data_input_file
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'max_iterations') then
            vartype='i'
            ivarvalue=>max_iterations

         elseif(varlabel.eq.'solution_epsilon') then
            vartype='r'
            rvarvalue=>solution_epsilon

         elseif(varlabel.eq.'normalize_solution_error') then
            vartype='l'
            lvarvalue=>normalize_solution_error

         elseif(varlabel.eq.'mie_epsilon') then
            vartype='r'
            rvarvalue=>mie_epsilon

         elseif(varlabel.eq.'translation_epsilon') then
            vartype='r'
            rvarvalue=>translation_epsilon

         elseif(varlabel.eq.'random_orientation') then
            vartype='l'
            lvarvalue=>random_orientation

         elseif(varlabel.eq.'t_matrix_centered_on_1') then
            vartype='l'
            lvarvalue=>t_matrix_centered_on_1

         elseif(varlabel.eq.'t_matrix_convergence_epsilon') then
            vartype='r'
            rvarvalue=>t_matrix_convergence_epsilon

         elseif(varlabel.eq.'solution_method') then
            vartype='a'
            avarvalue=>solution_method

         elseif(varlabel.eq.'t_matrix_procs_per_solution') then
            vartype='i'
            ivarvalue=>t_matrix_procs_per_solution

         elseif(varlabel.eq.'max_t_matrix_order') then
            vartype='i'
            ivarvalue=>max_t_matrix_order

         elseif(varlabel.eq.'fft_translation_option') then
            vartype='l'
            lvarvalue=>input_fft_translation_option

         elseif(varlabel.eq.'node_order') then
            vartype='i'
            ivarvalue=>input_node_order

         elseif(varlabel.eq.'min_fft_nsphere') then
            vartype='i'
            ivarvalue=>min_fft_nsphere

         elseif(varlabel.eq.'neighbor_node_model') then
            vartype='i'
            ivarvalue=>neighbor_node_model

         elseif(varlabel.eq.'cell_volume_fraction') then
            vartype='r'
            rvarvalue=>input_cell_volume_fraction

         elseif(varlabel.eq.'incident_beta_deg') then
            vartype='r'
            rvarvalue=>incident_beta_deg
            incident_beta_specified=.true.

         elseif(varlabel.eq.'incident_sin_beta') then
            vartype='r'
            rvarvalue=>incident_sin_beta
            incident_beta_specified=.false.

         elseif(varlabel.eq.'incident_direction') then
            vartype='i'
            ivarvalue=>incident_direction

         elseif(varlabel.eq.'incident_alpha_deg') then
            vartype='r'
            rvarvalue=>incident_alpha_deg

         elseif(varlabel.eq.'gaussian_beam_constant') then
            vartype='r'
            rvarvalue=>gaussian_beam_constant

         elseif(varlabel.eq.'excitation_radius') then
            vartype='r'
            rvarvalue=>excitation_radius

         elseif(varlabel.eq.'incidence_average') then
            vartype='l'
            lvarvalue=>incidence_average

         elseif(varlabel.eq.'number_incident_directions') then
            vartype='i'
            ivarvalue=>number_incident_directions

         elseif(varlabel.eq.'gaussian_beam_focal_point') then
            vartype='r'
            varlen=3
            ravarvalue=>gaussian_beam_focal_point(1:3)

         elseif(varlabel.eq.'calculate_scattering_matrix') then
            vartype='l'
            lvarvalue=>calculate_scattering_matrix

         elseif(varlabel.eq.'single_origin_expansion') then
            vartype='l'
            lvarvalue=>single_origin_expansion

         elseif(varlabel.eq.'scattering_map_model') then
            vartype='i'
            ivarvalue=>scattering_map_model

         elseif(varlabel.eq.'scattering_map_dimension') then
            vartype='i'
            ivarvalue=>scattering_map_dimension

         elseif(varlabel.eq.'scattering_map_increment') then
            vartype='r'
            rvarvalue=>scattering_map_increment

         elseif(varlabel.eq.'azimuthal_average') then
            vartype='l'
            lvarvalue=>azimuthal_average

         elseif(varlabel.eq.'incident_frame') then
            vartype='l'
            lvarvalue=>incident_frame

         elseif(varlabel.eq.'print_sphere_data') then
            vartype='l'
            lvarvalue=>print_sphere_data

         elseif(varlabel.eq.'number_spheres') then
            vartype='i'
            ivarvalue=>number_spheres
            recalculate_surface_matrix=.true.
            number_spheres_specified=.true.

         elseif(varlabel.eq.'length_scale_factor') then
            vartype='r'
            rvarvalue=>length_scale_factor
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'ref_index_scale_factor') then
            vartype='c'
            cvarvalue=>ref_index_scale_factor
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'number_plane_boundaries') then
            vartype='i'
            ivarvalue=>number_plane_boundaries
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'maximum_integration_subdivisions') then
            vartype='i'
            ivarvalue=>maximum_integration_subdivisions
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'integration_error_epsilon') then
            vartype='r'
            rvarvalue=>integration_error_epsilon
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'integration_limit_epsilon') then
            vartype='r'
            rvarvalue=>integration_limit_epsilon
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'gf_switch_factor') then
            vartype='r'
            rvarvalue=>gf_switch_factor
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'s_scale_constant') then
            vartype='r'
            rvarvalue=>s_scale_constant
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'real_axis_integration_limit') then
            vartype='r'
            rvarvalue=>real_axis_integration_limit
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'minimum_integration_spacing') then
            vartype='r'
            rvarvalue=>minimum_integration_spacing
            recalculate_surface_matrix=.true.

         elseif(varlabel.eq.'move_to_front') then
            vartype='l'
            lvarvalue=>move_to_front

         elseif(varlabel.eq.'move_to_back') then
            vartype='l'
            lvarvalue=>move_to_back

         elseif(varlabel.eq.'store_translation_matrix') then
            vartype='l'
            lvarvalue=>store_translation_matrix

         elseif(varlabel.eq.'store_surface_matrix') then
            vartype='l'
            lvarvalue=>store_surface_matrix

         elseif(varlabel.eq.'calculate_near_field') then
            vartype='l'
            lvarvalue=>calculate_near_field

         elseif(varlabel.eq.'store_surface_vector') then
            vartype='l'
            lvarvalue=>store_surface_vector

         elseif(varlabel.eq.'near_field_output_file') then
            vartype='a'
            avarvalue=>near_field_output_file
            append_near_field_output_file=.false.

         elseif(varlabel.eq.'near_field_calculation_model') then
            vartype='i'
            ivarvalue=>near_field_calculation_model

         elseif(varlabel.eq.'near_field_expansion_order') then
            vartype='i'
            ivarvalue=>near_field_expansion_order

         elseif(varlabel.eq.'near_field_expansion_spacing') then
            vartype='r'
            rvarvalue=>near_field_expansion_spacing

         elseif(varlabel.eq.'near_field_step_size') then
            vartype='r'
            rvarvalue=>near_field_step_size

         elseif(varlabel.eq.'near_field_minimum_border') then
            vartype='r'
            varlen=3
            ravarvalue=>near_field_plane_vertices(1:3,1)

         elseif(varlabel.eq.'near_field_maximum_border') then
            vartype='r'
            varlen=3
            ravarvalue=>near_field_plane_vertices(1:3,2)

         elseif(varlabel.eq.'normalize_s11') then
            vartype='l'
            lvarvalue=>normalize_s11

         elseif(varlabel.eq.'periodic_lattice') then
            vartype='l'
            lvarvalue=>periodic_lattice

         elseif(varlabel.eq.'phase_shift_form') then
            vartype='l'
            lvarvalue=>phase_shift_form

         elseif(varlabel.eq.'finite_lattice') then
            vartype='l'
            lvarvalue=>finite_lattice

         elseif(varlabel.eq.'cell_width') then
            vartype='r'
            varlen=2
            ravarvalue=>input_cell_width(1:2)
            square_cell=.false.

         elseif(varlabel.eq.'cell_width_x') then
            vartype='r'
            rvarvalue=>input_cell_width_x
            square_cell=.true.

         elseif(varlabel.eq.'target_dimensions') then
            vartype='r'
            varlen=3
            ravarvalue=>target_dimensions(1:3)
            target_width_specified=.false.

         elseif(varlabel.eq.'target_shape') then
            vartype='i'
            ivarvalue=>target_shape

         elseif(varlabel.eq.'target_width') then
            vartype='r'
            rvarvalue=>target_width
            target_width_specified=.true.

         elseif(varlabel.eq.'target_thickness') then
            vartype='r'
            rvarvalue=>target_thickness
            target_width_specified=.true.

         elseif(varlabel.eq.'psd_sigma') then
            vartype='r'
            rvarvalue=>psd_sigma

         elseif(varlabel.eq.'max_number_time_steps') then
            vartype='i'
            ivarvalue=>max_number_time_steps

         elseif(varlabel.eq.'number_configurations') then
            vartype='i'
            ivarvalue=>number_configurations

         elseif(varlabel.eq.'sphere_volume_fraction') then
            vartype='r'
            rvarvalue=>sphere_volume_fraction
            number_spheres_specified=.false.

         elseif(varlabel.eq.'periodic_bc') then
            vartype='l'
            varlen=3
            lavarvalue=>periodic_bc

         elseif(varlabel.eq.'sphere_1_fixed') then
            vartype='l'
            lvarvalue=>sphere_1_fixed

         elseif(varlabel.eq.'configuration_average') then
            vartype='l'
            lvarvalue=>configuration_average

         elseif(varlabel.eq.'frozen_configuration') then
            vartype='l'
            lvarvalue=>frozen_configuration

         elseif(varlabel.eq.'reflection_model') then
            vartype='l'
            lvarvalue=>reflection_model

         elseif(varlabel.eq.'absorption_sample_radius') then
            vartype='r'
            rvarvalue=>absorption_sample_radius

         elseif(varlabel.eq.'absorption_sample_radius_fraction') then
            vartype='r'
            rvarvalue=>absorption_sample_radius_fraction

         elseif(varlabel.eq.'auto_absorption_sample_radius') then
            vartype='l'
            lvarvalue=>auto_absorption_sample_radius

         elseif(varlabel.eq.'print_random_configuration') then
            vartype='l'
            lvarvalue=>print_random_configuration

         elseif(varlabel.eq.'print_timings') then
            vartype='l'
            lvarvalue=>print_timings

         elseif(varlabel.eq.'calculate_up_down_scattering') then
            vartype='l'
            lvarvalue=>input_calculate_up_down_scattering

         elseif(varlabel.eq.'x_shift') then
            vartype='r'
            rvarvalue=>x_shift

         elseif(varlabel.eq.'y_shift') then
            vartype='r'
            rvarvalue=>y_shift

         elseif(varlabel.eq.'z_shift') then
            vartype='r'
            rvarvalue=>z_shift

         elseif(varlabel.eq.'shifted_sphere') then
            vartype='i'
            ivarvalue=>shifted_sphere

         endif

         if(vartype.eq.'n') then
            varstatus=1
            if(present(var_status)) var_status=varstatus
            return
         endif
         if(present(var_type)) var_type=vartype
         if(present(i_var_pointer)) i_var_pointer=>ivarvalue
         if(present(r_var_pointer)) r_var_pointer=>rvarvalue
         if(present(c_var_pointer)) c_var_pointer=>cvarvalue

         if(operate) then
            if(vartype.eq.'i') then
               call set_string_to_int_variable(sentvarvalue, &
                  ivarvalue,var_operation=varop)
            elseif(vartype.eq.'r') then
               if(varlen.eq.1) then
                  call set_string_to_real_variable(sentvarvalue, &
                      rvarvalue,var_operation=varop)
               else
                  call set_string_to_real_array_variable(sentvarvalue, &
                      ravarvalue,var_operation=varop,var_len=varlen)
               endif
            elseif(vartype.eq.'c') then
               call set_string_to_cmplx_variable(sentvarvalue, &
                  cvarvalue,var_operation=varop)
            elseif(vartype.eq.'l') then
               if(varlen.eq.1) then
                  call set_string_to_logical_variable(sentvarvalue, &
                     lvarvalue,var_operation=varop)
               else
                  call set_string_to_logical_array_variable(sentvarvalue, &
                     lavarvalue,var_operation=varop,var_len=varlen)
               endif
            elseif(vartype.eq.'a') then
               avarvalue=sentvarvalue
            endif
         endif
         end subroutine variable_list_operation

         subroutine inputdata(inputfiledata,read_status)
         implicit none
         integer :: readok,n,spherenum,varstat,rank,stopit,istat,lines
         integer, save :: inputline
         integer, optional :: read_status
         real(8) :: rtemp(4)
         complex(8) :: ctemp(4)
         character*256 :: parmid,parmval,varop,inputfiledata(*)
         data inputline/1/

         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         readok=0
         stopit=0
         do while(readok.eq.0)
            parmid=inputfiledata(inputline)
            inputline=inputline+1
            if(trim(parmid).eq.'run_file') then
               parmval=inputfiledata(inputline)
               inputline=inputline+1
               if(trim(parmval).ne.' ') then
                  run_print_unit=3
                  if(rank.eq.0) then
                     open(3,file=trim(parmval))
                  endif
               endif
               cycle
            endif

            if(parmid(1:1).eq.'!'.or.parmid(1:1).eq.'%') then
               cycle
            endif

            if(trim(parmid).eq.'loop_variable') then
               loop_job=.true.
               n_nest_loops=n_nest_loops+1
               n=n_nest_loops
               parmid=inputfiledata(inputline)
               inputline=inputline+1
               if(trim(parmid).eq.'sphere_number') then
                  read(inputfiledata(inputline),*) spherenum
                  inputline=inputline+1
                  loop_sphere_number(n)=spherenum
                  parmid=inputfiledata(inputline)
                  inputline=inputline+1
               else
                  spherenum=1
               endif
               loop_var_label(n)=parmid
               call variable_list_operation(loop_var_label(n), &
                    var_type=loop_var_type(n),var_position=spherenum)
               if(loop_var_type(n).eq.'i') then
                  read(inputfiledata(inputline),*) i_var_start(n),i_var_stop(n),i_var_step(n)
               elseif(loop_var_type(n).eq.'r') then
                  read(inputfiledata(inputline),*) r_var_start(n),r_var_stop(n),r_var_step(n)
               elseif(loop_var_type(n).eq.'c') then
                  read(inputfiledata(inputline),*) c_var_start(n),c_var_stop(n),c_var_step(n)
               endif
               inputline=inputline+1
               cycle

            elseif(trim(parmid).eq.'sphere_data') then
               istat=0
               n=1
               if(rank.eq.0) then
                  open(20,file='temp_pos.dat')
               endif
               sphere_data_input_file='temp_pos.dat'
               do
                  parmval=inputfiledata(inputline)
                  if(trim(parmval).eq.'end_of_options') exit
                  inputline=inputline+1
                  if(trim(parmval).eq.'end_of_sphere_data') exit
                  if(parmval(1:1).eq.'!'.or.parmval(1:1).eq.'%') cycle
                  if(n.gt.number_spheres) cycle
                  read(parmval,*,iostat=istat) rtemp(1:4)
                  if(istat.ne.0) then
                     lines=3
                  else
                     read(parmval,*,iostat=istat) rtemp(1:4),ctemp(1)
                     if(istat.ne.0) then
                        lines=4
                     else
                        read(parmval,*,iostat=istat) rtemp(1:4),ctemp(1),ctemp(2)
                        if(istat.ne.0) then
                           lines=5
                        else
                           lines=6
                        endif
                     endif
                  endif
                  if(rank.eq.0) then
                     if(lines.eq.3) then
                        write(20,'(2(e20.12,'',''),e20.12)') rtemp(1:3)
                     elseif(lines.eq.4) then
                        write(20,'(3(e20.12,'',''),e20.12)') rtemp(1:4)
                     elseif(lines.eq.5) then
                        write(20,'(4(e20.12,'',''),'' ('',e20.12,'','',e20.12,'') '')') rtemp(1:4),ctemp(1)
                     else
                        write(20,'(4(e20.12,'',''),'' ('',e20.12,'','',e20.12,''), ('',e20.12,'','',e20.12,'') '')') &
                           rtemp(1:4),ctemp(1:2)
                     endif
                  endif
                  n=n+1
               enddo
               number_spheres=min(n,number_spheres)
               if(rank.eq.0) close(20)
               data_scaled=.false.
               temporary_pos_file=.true.
               recalculate_surface_matrix=.true.
               cycle

            elseif(trim(parmid).eq.'new_run') then
               repeat_run=.true.
               exit

            elseif(trim(parmid).eq.'end_of_options') then
               repeat_run=.false.
               readok=-1
               exit

            elseif(trim(parmid).eq.'layer_ref_index') then
               if(number_plane_boundaries.gt.max_number_plane_boundaries) then
                  if(rank.eq.0) write(run_print_unit,'('' max # plane boundaries exceeded:'',i3,''>'',i3)') &
                     number_plane_boundaries,max_number_plane_boundaries
                  stop
               endif
               parmval=inputfiledata(inputline)
               inputline=inputline+1
               read(parmval,*,iostat=istat) layer_ref_index(0)
               layer_ref_index(1:max(1,number_plane_boundaries))=layer_ref_index(0)
               read(parmval,*,iostat=istat) layer_ref_index(0:number_plane_boundaries)
               recalculate_surface_matrix=.true.

            elseif(trim(parmid).eq.'layer_thickness') then
               parmval=inputfiledata(inputline)
               inputline=inputline+1
               input_layer_thickness(1:max(1,number_plane_boundaries))=0.d0
               read(parmval,*,iostat=istat) input_layer_thickness(1:max(1,number_plane_boundaries))
               recalculate_surface_matrix=.true.

            else
               varstat=0
               call variable_list_operation(parmid, &
                   var_status=varstat)
               if(varstat.ne.0) then
                  if(rank.eq.0) then
                     write(run_print_unit,'('' unknown input parameter:'',a)') trim(parmid)
                     call flush(run_print_unit)
                     stopit=1
                  endif
                  cycle
               else
                  parmval=inputfiledata(inputline)
                  inputline=inputline+1
                  if(readok.ne.0) cycle
                  parmval=trim(parmval)
                  varop='assign'
                  call variable_list_operation(parmid,var_value=parmval, &
                      var_position=1,var_operation='assign', &
                      var_status=varstat)
               endif
            endif
         enddo
         if(stopit.eq.1) stop
         if(present(read_status)) read_status=varstat
         end subroutine inputdata

         subroutine main_calling_program(print_output,set_t_matrix_order,dry_run,mpi_comm)
         implicit none
         logical :: stopit,singleorigin,iframe,sett,printout,dryrun,averagerun
         logical, optional :: print_output,set_t_matrix_order,dry_run
         integer :: n,istat,niter,rank,numprocs,i,nodrw,celldim(3),itemp(6),sx,sy,maxt, &
            mpicomm
         integer, optional :: mpi_comm
         real(8) :: alpha,time1,r0(3),rtran,costheta, &
            csca,zext,targetvol,timet,tmin(3),tmax(3)
         complex(8) :: rimedium(2)
         character*256 :: timatrixfile
         if(present(dry_run)) then
            dryrun=dry_run
         else
            dryrun=.false.
         endif
         if(present(print_output)) then
            printout=print_output
         else
            printout=.true.
         endif
         if(present(set_t_matrix_order)) then
            sett=set_t_matrix_order
         else
            sett=.true.
         endif
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         averagerun=configuration_average.or.incidence_average
         calculate_up_down_scattering=input_calculate_up_down_scattering

         first_run=.false.
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
!         if(rank.ne.0) light_up=.false.
         local_rank=rank
         global_rank=rank
         if((.not.configuration_average).and.(.not.incidence_average)) n_configuration_groups=1

         random_configuration=(trim(sphere_data_input_file).eq.'random_configuration')
         if(random_configuration) then
            if(target_width_specified) then
               if(target_shape.eq.0) then
                  target_dimensions(1:2)=target_width
                  target_dimensions(3)=target_thickness
               elseif(target_shape.eq.1) then
                  target_dimensions(1:2)=target_width
                  target_dimensions(3)=target_thickness
               else
                  target_dimensions(1:3)=target_width
               endif
            endif
            call target_volume(targetvol)
            if(number_spheres_specified) then
               sphere_volume_fraction=dble(number_spheres)*4.d0*pi/3.d0/targetvol
            else
               number_spheres=ceiling(targetvol*sphere_volume_fraction)/(4.d0*pi/3.d0)
            endif
         endif
         if(allocated(sphere_radius)) then
            deallocate(sphere_radius, &
                  sphere_position, &
                  sphere_ref_index, &
                  host_sphere, &
                  number_field_expansions, &
                  sphere_excitation_switch)
         endif
         allocate(sphere_radius(number_spheres), &
                  sphere_position(3,number_spheres), &
                  sphere_ref_index(2,0:number_spheres), &
                  host_sphere(number_spheres), &
                  number_field_expansions(number_spheres), &
                  sphere_excitation_switch(number_spheres))

         if(random_configuration) then
            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' generating random configuration:'',$)')
               timet=mstm_mpi_wtime()
            endif
!            call generate_random_configuration(mpi_comm=mpicomm,skip_diffusion=dryrun)
            call generate_random_configuration(mpi_comm=mpicomm)
            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' completed, time:'',es12.5,'' s'')') mstm_mpi_wtime()-timet
            endif
            if(rank.eq.0) then
               if(ran_config_stat.ge.3) then
                  write(run_print_unit,'('' unable to generate random configuration'')')
                  stop
               endif
!               if(print_random_configuration.and.(.not.configuration_average)) then
               if(print_random_configuration.and.mstm_global_rank.eq.0) then
                  open(2,file='random_configuration.pos')
                  do i=1,number_spheres
                     write(2,'(4es13.5)') sphere_position(:,i)/length_scale_factor, &
                        sphere_radius(i)/length_scale_factor
                  enddo
                  close(2)
               endif
            endif
         else
            call read_sphere_data_input_file(mpi_comm=mpicomm)
         endif

         position_shift=(/x_shift,y_shift,z_shift/)*length_scale_factor
         if(shifted_sphere.gt.number_spheres) shifted_sphere=0
         if(any(position_shift.ne.0.d0)) then
            if(shifted_sphere.eq.0) then
               do i=1,number_spheres
                  sphere_position(:,i)=sphere_position(:,i)+position_shift(:)
               enddo
            else
               sphere_position(:,shifted_sphere)=sphere_position(:,shifted_sphere)+position_shift(:)
            endif
         endif

         sphere_ref_index(:,0)=layer_ref_index(0)
         if(periodic_lattice) then
            if(random_configuration.and.target_shape.eq.0) then
               cell_width(1:2)=target_dimensions(1:2)*2.d0*length_scale_factor
            else
               if(square_cell) then
                  cell_width=input_cell_width_x*length_scale_factor
               else
                  cell_width=input_cell_width*length_scale_factor
               endif
            endif
         endif

         plane_surface_present=number_plane_boundaries.gt.0
         layer_thickness=input_layer_thickness*length_scale_factor
         call plane_boundary_initialization()

         if(move_to_front.and.plane_surface_present) then
            zext=maxval(sphere_position(3,:)+sphere_radius(:))
            if(zext.gt.0.d0) sphere_position(3,:)=sphere_position(3,:)-zext
         endif
         if(move_to_back.and.plane_surface_present) then
            zext=minval(sphere_position(3,:)-sphere_radius(:))
            if(zext.lt.plane_boundary_position(number_plane_boundaries))  &
               sphere_position(3,:)=sphere_position(3,:)-zext+plane_boundary_position(number_plane_boundaries)
         endif

         stopit=.false.
         if(random_orientation) then
            if(number_plane_boundaries.gt.0) then
               if(rank.eq.0) write(run_print_unit,'('' random orientation requires number_plane_boundaries=0'')')
               stopit=.true.
            endif
            if(periodic_lattice) then
               if(rank.eq.0) write(run_print_unit,'('' random orientation and periodic lattice incompatible'')')
               stopit=.true.
            endif
         endif

         fft_translation_option=(input_fft_translation_option.and.number_spheres.ge.min_fft_nsphere)
         if(fft_translation_option) then
            if(number_plane_boundaries.gt.0) then
               if(rank.eq.0) write(run_print_unit,'('' fft option requires number_plane_boundaries=0'')')
               stopit=.true.
            endif
            if(periodic_lattice) then
               if(rank.eq.0) write(run_print_unit,'('' fft option and periodic lattice incompatible'')')
               stopit=.true.
            endif
         endif

!         if(fft_translation_option) single_origin_expansion=.true.

         if(stopit) return
if(light_up) then
write(*,'('' s2 '',i3)') mstm_global_rank
call flush(6)
endif
         call findhostspheres()

if(light_up) then
write(*,'('' s3 '',i3)') mstm_global_rank
call flush(6)
endif
         call sphere_layer_initialization()
         call miecoefcalc(mie_epsilon)
         call init(max_mie_order)
if(light_up) then
write(*,'('' s4 '',i3)') mstm_global_rank
call flush(6)
endif

         singleorigin=number_plane_boundaries.eq.0.and.single_origin_expansion
         iframe=singleorigin.and.incident_frame

         cluster_origin=0.d0
         if(singleorigin.or.random_orientation) then
            if(allocated(translation_order)) deallocate(translation_order)
            allocate(translation_order(number_spheres))
            translation_order(1:number_spheres)=sphere_order(1:number_spheres)
            cluster_origin=0.d0
            if((.not.configuration_average).and.(.not.incidence_average)) then
               if(gaussian_beam_constant.eq.0.d0) then
                  n=0
                  do i=1,number_spheres
                     if(host_sphere(i).eq.0) then
                        n=n+1
                        cluster_origin=cluster_origin+sphere_position(:,i)
                     endif
                  enddo
                  cluster_origin=cluster_origin/dble(n)
               else
                  cluster_origin=gaussian_beam_focal_point*length_scale_factor
               endif
            endif
            if(sett) then
               maxt=max_t_matrix_order
            else
               maxt=t_matrix_order
            endif
            t_matrix_order=max_mie_order
            do i=1,number_spheres
               if(host_sphere(i).eq.0) then
                  r0=cluster_origin
                  call exteriorrefindex(i,rimedium)
                  rtran=sqrt(sum((sphere_position(:,i)-r0(:))**2))
!                  if(rtran.gt.scattered_field_sample_length) cycle
                  call tranordertest(rtran,rimedium(1),sphere_order(i), &
                     translation_epsilon,translation_order(i))
                  translation_order(i)=min(translation_order(i),maxt)
                  t_matrix_order=max(t_matrix_order,translation_order(i))
               endif
            enddo
            if(.not.sett) t_matrix_order=maxt
         endif

if(light_up) then
write(*,'('' s5 '',i3)') mstm_global_rank
call flush(6)
endif
         one_side_only=.false.
         vol_radius=0.
         do i=1,number_spheres
            if(host_sphere(i).eq.0) then
               vol_radius=vol_radius+sphere_radius(i)**3
            endif
         enddo
         vol_radius=vol_radius**.333333

         if(periodic_lattice) then
            cross_section_radius=sqrt(product(cell_width)/pi)
         else
            if(gaussian_beam_constant.ne.0.d0) then
               cross_section_radius=1.d0/gaussian_beam_constant/sqrt(2.d0)
            elseif(reflection_model) then
               if(random_configuration) then
                  if(target_shape.eq.0) then
                     cross_section_radius=length_scale_factor*2.d0*sqrt(product(target_dimensions(1:2))/pi)
                  elseif(target_shape.ge.1) then
                     cross_section_radius=length_scale_factor*target_dimensions(1)
                  endif
               else
                  cross_section_radius=sqrt(product(sphere_max_position(1:2)-sphere_min_position(1:2))/pi)
               endif
               cross_section_radius=min(cross_section_radius,length_scale_factor*excitation_radius)
            else
               cross_section_radius=vol_radius
            endif
         endif

         if(auto_absorption_sample_radius.and.random_configuration) then
            absorption_sample_radius=absorption_sample_radius_fraction*target_dimensions(1)
         endif

         if(fft_translation_option) then
            cell_volume_fraction=input_cell_volume_fraction
            if((configuration_average.or.incidence_average).and.random_configuration) then
               tmin=-target_dimensions*length_scale_factor
               tmax=target_dimensions*length_scale_factor
               if(dryrun) call clear_fft_matrix(clear_h=.true.)
            else
               tmin=sphere_min_position
               tmax=sphere_max_position
               if(.not.averagerun) call clear_fft_matrix(clear_h=.true.)
            endif
            call node_selection(cell_volume_fraction,target_min=tmin,target_max=tmax)
            if(input_node_order.le.0) then
               node_order=-input_node_order+ceiling(d_cell)
            else
               node_order=input_node_order
            endif
            node_order=max(node_order,max_mie_order)
         endif

         sphere_excitation_switch=.true.
         do i=1,number_spheres
            if(host_sphere(i).ne.0) cycle
            if(random_configuration) then
               if(target_shape.le.1) then
                  rtran=sqrt(sum((sphere_position(1:2,i)-cluster_origin(1:2))**2))
               else
                  rtran=sqrt(sum((sphere_position(1:3,i)-cluster_origin(1:3))**2))
               endif
            else
               rtran=sqrt(sum((sphere_position(1:3,i)-cluster_origin(1:3))**2))
            endif
            sphere_excitation_switch(i)=rtran.le.excitation_radius*length_scale_factor
         enddo

if(light_up) then
write(*,'('' s6 '',i3)') mstm_global_rank
call flush(6)
endif
         if(random_orientation) then
            qeff_dim=1
            if(calculate_scattering_matrix) then
               scat_mat_udim=floor(180.00001d0/scattering_map_increment)
               scat_mat_mdim=16
               scat_mat_ldim=0
            endif
         else
            if(incident_beta_specified) then
               incident_beta=incident_beta_deg*pi/180.d0
               if(incident_beta_deg.le.90.d0) then
                  incident_direction=1
                  incident_sin_beta=dsin(incident_beta_deg*pi/180.d0)/dble(layer_ref_index(0))
               else
                  incident_direction=2
                  incident_sin_beta=dsin(incident_beta_deg*pi/180.d0) &
                     /dble(layer_ref_index(number_plane_boundaries))
               endif
            else
               incident_beta=0.d0
            endif
            if(incidence_average) then
               qeff_dim=1
            else
               qeff_dim=3
            endif
            alpha=incident_alpha_deg*pi/180.d0
            call incident_field_initialization(alpha,incident_sin_beta,incident_direction)
            if(calculate_scattering_matrix) then
               if(allocated(scat_mat)) deallocate(scat_mat)
               if(periodic_lattice) then
                  call periodic_lattice_scattering(amnp_s,pl_sca,dry_run=.true.,num_dirs=number_rl_dirs)
                  max_number_rl_dirs=maxval(number_rl_dirs)
                  if(allocated(rl_vec)) deallocate(rl_vec)
                  allocate(rl_vec(2,max_number_rl_dirs))
                  scat_mat_udim=max_number_rl_dirs
                  scat_mat_ldim=1
                  scat_mat_mdim=32
               else
                  if(scattering_map_model.eq.0) then
                     if(number_plane_boundaries.eq.0) then
                        scat_mat_udim=floor(180.00001d0/scattering_map_increment)
                        if(azimuthal_average) then
                           scat_mat_ldim=0
                        else
                           scat_mat_ldim=-scat_mat_udim
                        endif
                        scat_mat_mdim=16
                        scat_mat_amax=180.d0
                     else
                        scat_mat_udim=floor(90.00001d0/scattering_map_increment)
                        scat_mat_ldim=-scat_mat_udim
                        scat_mat_mdim=32
                        scat_mat_amax=90.d0
                     endif
                     scat_mat_amin=scat_mat_amax*(scat_mat_ldim/scat_mat_udim)
                  else
                     i=0
                     do sy=-scattering_map_dimension,scattering_map_dimension
                        do sx=-scattering_map_dimension,scattering_map_dimension
                           if(sx*sx+sy*sy.gt.scattering_map_dimension**2) cycle
                           i=i+1
                        enddo
                     enddo
                     scat_mat_udim=i
                     scat_mat_ldim=1
                     scat_mat_mdim=32
                  endif
               endif
            endif
            if(periodic_lattice.or.reflection_model.and.(target_shape.le.1)) then
               cross_section_radius=cross_section_radius*sqrt(cos(incident_beta))
            endif

         endif
if(light_up) then
write(*,'('' s7 '',i3)') mstm_global_rank
call flush(6)
endif
         if(allocated(boundary_sca)) deallocate(boundary_sca,boundary_ext)
         allocate(boundary_sca(2,0:number_plane_boundaries+1),boundary_ext(2,0:number_plane_boundaries+1))

         if(allocated(q_eff)) deallocate(q_eff,q_eff_tot,q_vabs)
         allocate(q_eff(3,qeff_dim,number_spheres),q_eff_tot(3,qeff_dim),q_vabs(qeff_dim,number_spheres))
         if(calculate_scattering_matrix) then
            if(allocated(scat_mat)) deallocate(scat_mat)
            allocate(scat_mat(scat_mat_mdim,scat_mat_ldim:scat_mat_udim))
         endif
if(light_up) then
write(*,'('' s8 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')
         if(rank.eq.0.and.printout) then
            call checkpositions()
            call print_run_variables(run_print_unit)
            open(2,file=output_file,access='append')
            call print_run_variables(2)
            close(2)
            time1=mstm_mpi_wtime()
         endif

         if(dryrun) return

         if(random_orientation) then
            niter=max_iterations
            timatrixfile='titemp.dat'
            call tmatrix_solution(solution_method=solution_method(1:1), &
               solution_eps=solution_epsilon, &
               convergence_eps=t_matrix_convergence_epsilon, &
               max_iterations=niter, &
               t_matrix_file=t_matrix_output_file, &
               procs_per_soln=t_matrix_procs_per_solution, &
               sphere_qeff=q_eff, &
               solution_status=istat, &
               mpi_comm=mpicomm)
            if(fft_translation_option) call clear_fft_matrix(clear_h=.true.)
            if(calculate_scattering_matrix) then
               if(allocated(sm_coef)) deallocate(sm_coef,sm_cf_coef)
               allocate(sm_coef(4,4,0:2*t_matrix_order),sm_cf_coef(4,4,0:2*t_matrix_order))
               nodrw=2*t_matrix_order
               call ranorientscatmatrix(t_matrix_output_file,sm_coef,sm_cf_coef, &
                  beam_width=gaussian_beam_constant, &
                  number_processors=t_matrix_procs_per_solution,mpi_comm=mpicomm)
               coherent_scattering_ratio=sm_cf_coef(1,1,0)
               do i=scat_mat_ldim,scat_mat_udim
                  costheta=cos(dble(i-scat_mat_ldim)*pi/dble(scat_mat_udim-scat_mat_ldim))
                  call ranorienscatmatrixcalc(costheta,sm_coef,nodrw,scat_mat(:,i))
               enddo
            endif
            call qtotcalc(number_spheres,qeff_dim,cross_section_radius,&
               q_eff,q_vabs,q_eff_tot)
         else
if(light_up) then
write(*,'('' s8.1 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')

            if(allocated(amnp_s)) deallocate(amnp_s)
            allocate(amnp_s(number_eqns,2))
            amnp_s=0.d0
            niter=max_iterations
            error_codes=0
            pl_error_codes=0
            recalculate_surface_matrix=.true.
if(light_up) then
write(*,'('' s8.2 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')
            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' generating solution:'',$)')
               timet=mstm_mpi_wtime()
            endif
            call fixedorsoln(alpha,incident_sin_beta,incident_direction,solution_epsilon,niter,amnp_s,q_eff, &
                    qeff_dim,solution_error,solution_iterations,1,istat, &
                    mpi_comm=mpicomm, &
                    excited_spheres=sphere_excitation_switch)
            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' completed, time:'',es12.5,'' s'')') mstm_mpi_wtime()-timet
            endif
            if(fft_translation_option) then
               call clear_fft_matrix()
            endif
            if(mstm_global_rank.eq.0..and.((.not.configuration_average).and.(.not.incidence_average))) then
               write(run_print_unit,'('' solution completed: number iterations='',i5)') solution_iterations
               call flush(run_print_unit)
            endif

if(light_up) then
write(*,'('' s8.3 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')
            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' post processing solution:'',$)')
               timet=mstm_mpi_wtime()
            endif

            call qtotcalc(number_spheres,qeff_dim,cross_section_radius,&
               q_eff,q_vabs,q_eff_tot)
            q_eff_tot(3,:)=q_eff_tot(1,:)-q_eff_tot(2,:)
            csca=(q_eff_tot(1,1)-q_eff_tot(2,1))*pi*cross_section_radius**2

            if(singleorigin) then
               if(allocated(amnp_0)) deallocate(amnp_0)
               allocate(amnp_0(2*t_matrix_order*(t_matrix_order+2),2))
               amnp_0=0.d0
if(light_up) then
write(*,'('' s8.3.1 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')
               do i=1,2
                  call merge_to_common_origin(t_matrix_order,amnp_s(:,i),amnp_0(:,i), &
                        origin_position=cluster_origin,merge_procs=.true., &
                        mpi_comm=mpicomm)
                  if(iframe) call rotvec(alpha,incident_beta,0.d0,t_matrix_order,t_matrix_order,amnp_0(:,i),1)
               enddo
if(light_up) then
write(*,'('' s8.3.2 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')
               if(calculate_up_down_scattering) call common_origin_hemispherical_scattering(amnp_0,boundary_sca)
if(light_up) then
write(*,'('' s8.3.3 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')
               if(azimuthal_average) then
                  if(allocated(scat_mat_exp_coef)) deallocate(scat_mat_exp_coef)
                  allocate(scat_mat_exp_coef(16,0:2*t_matrix_order,4))
                  call fosmexpansion(t_matrix_order,amnp_0,scat_mat_exp_coef(:,:,1),scat_mat_exp_coef(:,:,2), &
                     scat_mat_exp_coef(:,:,3),scat_mat_exp_coef(:,:,4),mpi_comm=mpicomm)
               endif
if(light_up) then
write(*,'('' s8.3.4 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')
               if(calculate_scattering_matrix) then
                  call scattering_matrix_calculation(amnp_0,scat_mat,mpi_comm=mpicomm)
               endif
if(light_up) then
write(*,'('' s8.3.5 '',i3)') mstm_global_rank
call flush(6)
endif
call mstm_mpi(mpi_command='barrier')
               if(gaussian_beam_constant.eq.0.d0) then
                  call boundary_extinction(amnp_0,alpha,incident_sin_beta,incident_direction,boundary_ext, &
                     common_origin=singleorigin)
               endif
            else
               if(calculate_up_down_scattering) call boundary_scattering(amnp_s,boundary_sca)
               if(calculate_scattering_matrix) then
                  call scattering_matrix_calculation(amnp_s,scat_mat,mpi_comm=mpicomm)
               endif
               if(gaussian_beam_constant.eq.0.d0) then
                  call boundary_extinction(amnp_s,alpha,incident_sin_beta,incident_direction,boundary_ext)
               endif
            endif
            if(gaussian_beam_constant.ne.0.d0) then
               boundary_ext=0
               boundary_ext(1:2,number_plane_boundaries+1)=-q_eff_tot(1,2:3)
            endif

            if(periodic_lattice) then
               call periodic_lattice_scattering(amnp_s,pl_sca)
            elseif(reflection_model) then
               pl_sca(:,1)=boundary_sca(:,number_plane_boundaries+1)
               pl_sca(:,2)=-boundary_sca(:,0)
            endif

            if(periodic_lattice.or.reflection_model) then
               call surface_absorptance_calculation()
            endif

            if(number_plane_boundaries.gt.0.and..not.periodic_lattice) then
               evan_sca(1:2)=q_eff_tot(1,2:3)-q_eff_tot(2,2:3)+boundary_sca(1:2,0) &
                  - boundary_sca(:,number_plane_boundaries+1)
            endif

            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' completed, time:'',es12.5,'' s'')') mstm_mpi_wtime()-timet
            endif

         endif

         if(rank.eq.0) solution_time=mstm_mpi_wtime()-time1

         if(number_plane_boundaries.gt.0) then
            itemp(1:4)=error_codes
            call mstm_mpi(mpi_command='reduce',mpi_rank=0,mpi_number=4,mpi_operation=mstm_mpi_sum, &
            mpi_recv_buf_i=error_codes,mpi_send_buf_i=itemp(1:4),mpi_comm=mpicomm)
         else
            error_codes=0
         endif
         if(periodic_lattice) then
            itemp(1:6)=pl_error_codes
            call mstm_mpi(mpi_command='reduce',mpi_rank=0,mpi_number=6,mpi_operation=mstm_mpi_sum, &
            mpi_recv_buf_i=pl_error_codes,mpi_send_buf_i=itemp(1:6),mpi_comm=mpicomm)
         else
            pl_error_codes=0
         endif
if(light_up) then
write(*,'('' s12 '',i3)') mstm_global_rank
call flush(6)
endif
         if(mstm_global_rank.eq.0.and.printout) then
            call print_calculation_results(output_file)
         endif
if(light_up) then
write(*,'('' s13 '',i3)') mstm_global_rank
call flush(6)
endif
         if(calculate_near_field) then
            celldim=ceiling((near_field_plane_vertices(:,2)-near_field_plane_vertices(:,1))/near_field_step_size)
            celldim=max(celldim,(/1,1,1/))
            if(append_near_field_output_file) then
               open(2,file=near_field_output_file,access='append')
            else
               open(2,file=near_field_output_file)
            endif
            append_near_field_output_file=.true.
            call near_field_calculation(amnp_s,alpha,incident_sin_beta,incident_direction, &
               near_field_plane_vertices,celldim, &
               incident_model=near_field_calculation_model,output_unit=2,output_header=.true.)
            close(2)
         endif
         end subroutine main_calling_program

         subroutine configuration_average_calling_program()
         implicit none
         logical :: singleorigin,iframe
         integer :: rank,numprocs, &
            numprocsperconfig,configcolor,configgroup,configcomm,configrank,config0comm,nconfigave,nsend
         real(8) :: time1,timet,diffac
         real(8), allocatable :: texpcoef(:,:,:)
         character*256 :: tmatchar1,tmatchar2
         data tmatchar1,tmatchar2/'tmat-','.tmp'/
         first_run=.false.
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
         local_rank=rank
         global_rank=rank

         if(max_iterations.lt.0) then
            numprocsperconfig=2
         else
            numprocsperconfig=4
         endif
         n_configuration_groups=numprocs/numprocsperconfig
         n_configuration_groups=max(n_configuration_groups,1)
         configcolor=floor(dble(n_configuration_groups*rank)/dble(numprocs))
         configgroup=configcolor
         call mstm_mpi(mpi_command='split', &
              mpi_color=configcolor,mpi_key=rank, &
              mpi_new_comm=configcomm)
         call mstm_mpi(mpi_command='rank', &
              mpi_rank=configrank, &
              mpi_comm=configcomm)
         configcolor=configrank
         call mstm_mpi(mpi_command='split', &
              mpi_color=configcolor,mpi_key=rank, &
              mpi_new_comm=config0comm)
         random_configuration=.true.
         singleorigin=number_plane_boundaries.eq.0.and.single_origin_expansion
!singleorigin=.true.
         iframe=singleorigin.and.incident_frame

         call main_calling_program(print_output=.false.,set_t_matrix_order=.true.,dry_run=.true.)

         if(allocated(q_eff_ave)) deallocate(q_eff_ave,q_eff_tot_ave,q_vabs_ave,sphere_position_ave,boundary_sca_ave, &
             boundary_ext_ave,dif_boundary_sca)
         allocate(q_eff_ave(3,qeff_dim,number_spheres),q_eff_tot_ave(3,qeff_dim),q_vabs_ave(qeff_dim,number_spheres), &
            sphere_position_ave(3,number_spheres),boundary_sca_ave(2,0:number_plane_boundaries+1), &
            boundary_ext_ave(2,0:number_plane_boundaries+1),dif_boundary_sca(2,0:number_plane_boundaries+1))
         q_eff_ave=0.d0
         q_eff_tot_ave=0.d0
         q_vabs_ave=0.d0
         sphere_position_ave=0.d0
         pl_sca_ave=0.d0
         boundary_sca_ave=0.d0
         boundary_ext_ave=0.d0
         solution_time_ave=0.d0
         surface_absorptance_ave=0.
         if(singleorigin) then
            if(allocated(amnp_0_ave)) deallocate(amnp_0_ave,scat_mat_exp_coef_ave)
            allocate(amnp_0_ave(2*t_matrix_order*(t_matrix_order+2),2), &
               scat_mat_exp_coef_ave(16,0:2*t_matrix_order,4))
            amnp_0_ave=0.d0
            scat_mat_exp_coef_ave=0.d0
         endif
         if(calculate_scattering_matrix) then
            if(allocated(scat_mat_ave)) deallocate(scat_mat_ave,dif_scat_mat)
            allocate(scat_mat_ave(scat_mat_mdim,scat_mat_ldim:scat_mat_udim), &
               dif_scat_mat(scat_mat_mdim,scat_mat_ldim:scat_mat_udim))
            scat_mat_ave=0.d0
         endif

         nconfigave=0
         do random_configuration_number=1,number_configurations/n_configuration_groups

            if(rank.eq.0) then
               if(random_configuration_number.eq.1) then
                  call print_run_variables(run_print_unit)
                  open(2,file=output_file,access='append')
                  call print_run_variables(2)
                  close(2)
               endif
               write(run_print_unit,'('' configuration averaging, samples:'',i5,''-'',i5)') &
                  (random_configuration_number-1)*n_configuration_groups+1, &
                  random_configuration_number*n_configuration_groups
            endif

            if(rank.eq.0) time1=mstm_mpi_wtime()

            call main_calling_program(print_output=.false.,set_t_matrix_order=.false.,mpi_comm=configcomm)

            if(singleorigin) then
               if(sphere_1_fixed) call subtract_1_from_0()
               amnp_0_ave=amnp_0_ave+amnp_0
            endif

            if(configrank.eq.0) then
               if(rank.eq.0) solution_time=mstm_mpi_wtime()-time1
               q_eff_ave=q_eff_ave+q_eff
               q_eff_tot_ave=q_eff_tot_ave+q_eff_tot
               q_vabs_ave=q_vabs_ave+q_vabs
               sphere_position_ave=sphere_position_ave+sphere_position
               pl_sca_ave=pl_sca_ave+pl_sca
               surface_absorptance_ave=surface_absorptance_ave+surface_absorptance
               if(calculate_up_down_scattering) boundary_sca_ave=boundary_sca_ave+boundary_sca
               boundary_ext_ave=boundary_ext_ave+boundary_ext
               if(singleorigin.and.azimuthal_average) scat_mat_exp_coef_ave=scat_mat_exp_coef_ave+scat_mat_exp_coef
               if(calculate_scattering_matrix) then
                  scat_mat_ave=scat_mat_ave+scat_mat
               endif
               if(rank.eq.0) solution_time_ave=solution_time_ave+solution_time
            endif

            nsend=3*number_spheres
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=sphere_position_ave, &
               mpi_recv_buf_dp=sphere_position, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=3*qeff_dim*number_spheres
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=q_eff_ave, &
               mpi_recv_buf_dp=q_eff, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=3*qeff_dim
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=q_eff_tot_ave, &
               mpi_recv_buf_dp=q_eff_tot, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=qeff_dim*number_spheres
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=q_vabs_ave, &
               mpi_recv_buf_dp=q_vabs, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=4
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=pl_sca_ave, &
               mpi_recv_buf_dp=pl_sca, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=2
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=surface_absorptance_ave, &
               mpi_recv_buf_dp=surface_absorptance, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            if(calculate_up_down_scattering) then
               nsend=2*(number_plane_boundaries+2)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=boundary_sca_ave, &
                  mpi_recv_buf_dp=boundary_sca, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=boundary_ext_ave, &
               mpi_recv_buf_dp=boundary_ext, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            if(singleorigin) then
               nsend=4*t_matrix_order*(t_matrix_order+2)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dc=amnp_0_ave, &
                  mpi_recv_buf_dc=amnp_0, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif
            if(calculate_scattering_matrix) then
               nsend=scat_mat_mdim*(scat_mat_udim-scat_mat_ldim+1)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=scat_mat_ave, &
                  mpi_recv_buf_dp=scat_mat, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif
            if(singleorigin.and.azimuthal_average) then
               nsend=16*4*(2*t_matrix_order+1)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=scat_mat_exp_coef_ave, &
                  mpi_recv_buf_dp=scat_mat_exp_coef, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif

            nconfigave=nconfigave+n_configuration_groups

            diffac=(dble(number_spheres)-1.d0)/dble(number_spheres)
            diffac=1.d0
!            diffac=(1.d0-(1.d0/dble(number_spheres))**.5d0)

            if(singleorigin) then
               if(rank.eq.0.and.print_timings) then
                  timet=mstm_mpi_wtime()
                  write(run_print_unit,'('' calculating diffuse field:'',$)')
               endif
               amnp_0=amnp_0/dble(nconfigave)
               if(azimuthal_average) then
                  allocate(texpcoef(16,0:2*t_matrix_order,4))
                  texpcoef=scat_mat_exp_coef/dble(nconfigave)
                  call fosmexpansion(t_matrix_order,amnp_0,scat_mat_exp_coef(:,:,1),scat_mat_exp_coef(:,:,2), &
                     scat_mat_exp_coef(:,:,3),scat_mat_exp_coef(:,:,4),mpi_comm=configcomm)
               endif
               if(calculate_scattering_matrix) then
                  call scattering_matrix_calculation(amnp_0,dif_scat_mat,mpi_comm=configcomm)
               endif
               if(azimuthal_average) then
                  scat_mat_exp_coef=texpcoef-scat_mat_exp_coef*diffac
                  deallocate(texpcoef)
               endif
               call common_origin_hemispherical_scattering(amnp_0,dif_boundary_sca)
               if(rank.eq.0.and.print_timings) write(run_print_unit,'('' completed, '',es12.4,'' sec'')') mstm_mpi_wtime()-timet
            endif
            if(allocated(amnp_0)) deallocate(amnp_0)

            if(rank.eq.0) then
               sphere_position=sphere_position/dble(nconfigave)
               q_eff=q_eff/dble(nconfigave)
               q_vabs=q_vabs/dble(nconfigave)
               q_eff_tot=q_eff_tot/dble(nconfigave)
               pl_sca=pl_sca/dble(nconfigave)
               surface_absorptance=surface_absorptance/dble(nconfigave)
               if(calculate_up_down_scattering) boundary_sca=boundary_sca/dble(nconfigave)
               boundary_ext=boundary_ext/dble(nconfigave)
!               if(singleorigin) dif_boundary_sca=boundary_sca-dif_boundary_sca
               if(singleorigin) dif_boundary_sca=boundary_sca &
                   -dif_boundary_sca*diffac
               if(calculate_scattering_matrix) then
                  scat_mat=scat_mat/dble(nconfigave)
! experiment

!                  if(singleorigin) dif_scat_mat=scat_mat-dif_scat_mat*(dble(number_spheres-1)/dble(number_spheres))
                  if(singleorigin) dif_scat_mat=scat_mat-dif_scat_mat*diffac

!                  dif_scat_mat=scat_mat-dif_scat_mat
               endif
               solution_time=solution_time_ave/dble(random_configuration_number)
               call print_calculation_results(output_file)
            endif
         enddo

         end subroutine configuration_average_calling_program

         subroutine incidence_average_calling_program()
         implicit none
         logical :: singleorigin,prancon,aa,soe,iframe,cuds
         integer :: rank,numprocs, &
            numprocsperconfig,configcolor,configgroup,configcomm,configrank,config0comm,nconfigave,nsend
         real(8) :: time1,timet
         real(8), allocatable :: texpcoef(:,:,:)
         character*256 :: sdatfile
         first_run=.false.
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
         local_rank=rank
         global_rank=rank
         sdatfile=sphere_data_input_file
         prancon=print_random_configuration
         print_random_configuration=.true.
         aa=azimuthal_average
         soe=single_origin_expansion
         iframe=incident_frame
         cuds=calculate_up_down_scattering
         azimuthal_average=.true.
         single_origin_expansion=.true.
         incident_frame=.true.
         calculate_up_down_scattering=.false.

         if(max_iterations.lt.0) then
            numprocsperconfig=2
         else
            numprocsperconfig=4
         endif
         n_configuration_groups=numprocs/numprocsperconfig
         n_configuration_groups=max(n_configuration_groups,1)
         configcolor=floor(dble(n_configuration_groups*rank)/dble(numprocs))
         configgroup=configcolor
         call mstm_mpi(mpi_command='split', &
              mpi_color=configcolor,mpi_key=rank, &
              mpi_new_comm=configcomm)
         call mstm_mpi(mpi_command='rank', &
              mpi_rank=configrank, &
              mpi_comm=configcomm)
         configcolor=configrank
         call mstm_mpi(mpi_command='split', &
              mpi_color=configcolor,mpi_key=rank, &
              mpi_new_comm=config0comm)
if(rank.eq.0) then
write(run_print_unit,'('' number processors, number groups:'',2i8)') numprocs,n_configuration_groups
call flush(run_print_unit)
endif
         singleorigin=number_plane_boundaries.eq.0.and.single_origin_expansion
         incident_beta_specified=.true.

         call main_calling_program(print_output=.false.,set_t_matrix_order=.true.,dry_run=.true.)

         if(trim(sphere_data_input_file).eq.'random_configuration') then
            sphere_data_input_file='random_configuration.pos'
         endif

         if(allocated(q_eff_ave)) deallocate(q_eff_ave,q_eff_tot_ave,q_vabs_ave,boundary_sca_ave, &
             boundary_ext_ave,dif_boundary_sca)
         allocate(q_eff_ave(3,qeff_dim,number_spheres),q_eff_tot_ave(3,qeff_dim),q_vabs_ave(qeff_dim,number_spheres), &
            boundary_sca_ave(2,0:number_plane_boundaries+1), &
            boundary_ext_ave(2,0:number_plane_boundaries+1),dif_boundary_sca(2,0:number_plane_boundaries+1))
         q_eff_ave=0.d0
         q_eff_tot_ave=0.d0
         q_vabs_ave=0.d0
         pl_sca_ave=0.d0
         boundary_sca_ave=0.d0
         boundary_ext_ave=0.d0
         solution_time_ave=0.d0
         if(singleorigin) then
            if(allocated(amnp_0_ave)) deallocate(amnp_0_ave,scat_mat_exp_coef_ave)
            allocate(amnp_0_ave(2*t_matrix_order*(t_matrix_order+2),2), &
               scat_mat_exp_coef_ave(16,0:2*t_matrix_order,4))
            amnp_0_ave=0.d0
            scat_mat_exp_coef_ave=0.d0
         endif
         if(calculate_scattering_matrix) then
            if(allocated(scat_mat_ave)) deallocate(scat_mat_ave,dif_scat_mat)
            allocate(scat_mat_ave(scat_mat_mdim,scat_mat_ldim:scat_mat_udim), &
               dif_scat_mat(scat_mat_mdim,scat_mat_ldim:scat_mat_udim))
            scat_mat_ave=0.d0
         endif

         nconfigave=0
         do incident_direction_number=1,ceiling(dble(number_incident_directions)/dble(n_configuration_groups))

            if(rank.eq.0) then
               if(incident_direction_number.eq.1) then
                  call print_run_variables(run_print_unit)
                  open(2,file=output_file,access='append')
                  call print_run_variables(2)
                  close(2)
               endif
               write(run_print_unit,'('' incidence averaging, samples:'',i5,''-'',i5)') &
                  (incident_direction_number-1)*n_configuration_groups+1, &
                  incident_direction_number*n_configuration_groups
            endif

            call sample_incident_direction(mpi_comm=configcomm)

            if(rank.eq.0) time1=mstm_mpi_wtime()

            call main_calling_program(print_output=.false.,set_t_matrix_order=.false.,mpi_comm=configcomm)

            if(singleorigin) then
               amnp_0_ave=amnp_0_ave+amnp_0
            endif

            if(configrank.eq.0) then
               if(rank.eq.0) solution_time=mstm_mpi_wtime()-time1
               q_eff_ave=q_eff_ave+q_eff
               q_eff_tot_ave=q_eff_tot_ave+q_eff_tot
               q_vabs_ave=q_vabs_ave+q_vabs
               pl_sca_ave=pl_sca_ave+pl_sca
               if(calculate_up_down_scattering) boundary_sca_ave=boundary_sca_ave+boundary_sca
               boundary_ext_ave=boundary_ext_ave+boundary_ext
               if(singleorigin.and.azimuthal_average) scat_mat_exp_coef_ave=scat_mat_exp_coef_ave+scat_mat_exp_coef
               if(calculate_scattering_matrix) then
                  scat_mat_ave=scat_mat_ave+scat_mat
               endif
               if(rank.eq.0) solution_time_ave=solution_time_ave+solution_time
            endif

            nsend=3*qeff_dim*number_spheres
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=q_eff_ave, &
               mpi_recv_buf_dp=q_eff, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=3*qeff_dim
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=q_eff_tot_ave, &
               mpi_recv_buf_dp=q_eff_tot, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=qeff_dim*number_spheres
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=q_vabs_ave, &
               mpi_recv_buf_dp=q_vabs, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            nsend=4
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=pl_sca_ave, &
               mpi_recv_buf_dp=pl_sca, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            if(calculate_up_down_scattering) then
               nsend=2*(number_plane_boundaries+2)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=boundary_sca_ave, &
                  mpi_recv_buf_dp=boundary_sca, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dp=boundary_ext_ave, &
               mpi_recv_buf_dp=boundary_ext, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            if(singleorigin) then
               nsend=4*t_matrix_order*(t_matrix_order+2)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dc=amnp_0_ave, &
                  mpi_recv_buf_dc=amnp_0, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif
            if(calculate_scattering_matrix) then
               nsend=scat_mat_mdim*(scat_mat_udim-scat_mat_ldim+1)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=scat_mat_ave, &
                  mpi_recv_buf_dp=scat_mat, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif
            if(singleorigin.and.azimuthal_average) then
               nsend=16*4*(2*t_matrix_order+1)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=scat_mat_exp_coef_ave, &
                  mpi_recv_buf_dp=scat_mat_exp_coef, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif

            nconfigave=nconfigave+n_configuration_groups

            if(singleorigin) then
               if(rank.eq.0.and.print_timings) then
                  timet=mstm_mpi_wtime()
                  write(run_print_unit,'('' calculating diffuse field:'',$)')
               endif
               amnp_0=amnp_0/dble(nconfigave)
               if(azimuthal_average) then
                  allocate(texpcoef(16,0:2*t_matrix_order,4))
                  texpcoef=scat_mat_exp_coef/dble(nconfigave)
                  call fosmexpansion(t_matrix_order,amnp_0,scat_mat_exp_coef(:,:,1),scat_mat_exp_coef(:,:,2), &
                     scat_mat_exp_coef(:,:,3),scat_mat_exp_coef(:,:,4),mpi_comm=configcomm)
               endif
               if(calculate_scattering_matrix) then
                  call scattering_matrix_calculation(amnp_0,dif_scat_mat,mpi_comm=configcomm)
               endif
               if(azimuthal_average) then
                  scat_mat_exp_coef=texpcoef-scat_mat_exp_coef*(dble(number_spheres-1)/dble(number_spheres))**2
                  deallocate(texpcoef)
               endif
               call common_origin_hemispherical_scattering(amnp_0,dif_boundary_sca)
               if(rank.eq.0.and.print_timings) write(run_print_unit,'('' completed, '',es12.4,'' sec'')') mstm_mpi_wtime()-timet
            endif
            if(allocated(amnp_0)) deallocate(amnp_0)

            if(rank.eq.0) then
               q_eff=q_eff/dble(nconfigave)
               q_vabs=q_vabs/dble(nconfigave)
               q_eff_tot=q_eff_tot/dble(nconfigave)
               pl_sca=pl_sca/dble(nconfigave)
               if(calculate_up_down_scattering) boundary_sca=boundary_sca/dble(nconfigave)
               boundary_ext=boundary_ext/dble(nconfigave)
!               if(singleorigin) dif_boundary_sca=boundary_sca-dif_boundary_sca
               if(singleorigin) dif_boundary_sca=boundary_sca &
                   -dif_boundary_sca*dble(number_spheres-1)/dble(number_spheres)
               if(calculate_scattering_matrix) then
                  scat_mat=scat_mat/dble(nconfigave)
! experiment
                  if(singleorigin) dif_scat_mat=scat_mat-dif_scat_mat*(dble(number_spheres-1)/dble(number_spheres))**2

!                  dif_scat_mat=scat_mat-dif_scat_mat
               endif
               solution_time=solution_time_ave/dble(random_configuration_number)
               call print_calculation_results(output_file)
            endif
         enddo
         sphere_data_input_file=sdatfile
         print_random_configuration=prancon
         azimuthal_average=aa
         single_origin_expansion=soe
         incident_frame=iframe
         calculate_up_down_scattering=cuds

         end subroutine incidence_average_calling_program

         subroutine subtract_1_from_0()
         implicit none
         integer :: i,i1,mnp0,mnp1,n,m,p
         real(8) :: fn
         complex(8) :: a(2,2),b(2,2)

         fn=dble(number_spheres-1)/dble(number_spheres)
         do i=1,number_spheres
            if(sum(sphere_position(:,i)**2).lt.1.d-7) then
               i1=i
               exit
            endif
         enddo
         do n=1,sphere_order(i1)
            do m=-n,n
               do p=1,2
                  mnp1=amnpaddress(m,n,p,sphere_order(i1),2)
                  a(p,:)=amnp_s(mnp1+sphere_offset(i1),:)
               enddo
               b(1,:)=a(1,:)+a(2,:)
               b(2,:)=a(1,:)-a(2,:)
               do p=1,2
                  mnp0=amnpaddress(m,n,p,t_matrix_order,2)
                  amnp_0(mnp0,:)=amnp_0(mnp0,:)-fn*b(p,:)
               enddo
            enddo
         enddo
         end subroutine subtract_1_from_0

         subroutine surface_absorptance_calculation()
         implicit none
         integer :: i
         real(8) :: rsamp,r,asamp
         if(periodic_lattice) then
            surface_absorptance(1:2)=q_eff_tot(2,2:3)
         else
            surface_absorptance=0.d0
            asamp=absorption_sample_radius*length_scale_factor
            rsamp=min(cross_section_radius,asamp)
            do i=1,number_spheres
               if(host_sphere(i).ne.0) cycle
               r=sqrt(sum((sphere_position(1:2,i)-cluster_origin(1:2))**2))
               if(r.le.asamp) then
                  surface_absorptance(1:2)=surface_absorptance(1:2) &
                     +q_eff(2,2:3,i)*(sphere_radius(i)/rsamp)**2
               endif
            enddo
         endif
         end subroutine surface_absorptance_calculation

         subroutine sample_incident_direction(mpi_comm)
         implicit none
         integer :: mpicomm,rank,numprocs
         integer, optional :: mpi_comm
         real(8) :: rnum(2),sbuf(2),cbeta
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         if(rank.eq.0) then
            call random_number(rnum)
            cbeta=2.d0*rnum(1)-1.d0
            sbuf(1)=180.d0/pi*dacos(cbeta)
            sbuf(2)=360.d0*rnum(2)
         endif
         if(numprocs.gt.1) then
            call mstm_mpi(mpi_command='bcast',mpi_number=2, &
            mpi_rank=0,mpi_send_buf_dp=sbuf,mpi_comm=mpicomm)
         endif
         incident_beta_deg=sbuf(1)
         incident_alpha_deg=sbuf(2)
         end subroutine sample_incident_direction

         subroutine effective_extinction_coefficient_ratio(scacoef,abscoef,srat,arat)
         implicit none
         integer :: i
         real(8) :: miesca,mieabs,scaq,absq,h,scacoef,abscoef,root0,root,func,dfunc,srat,arat,across

!         miesca=(mean_qext_mie-mean_qabs_mie)*sphere_volume_fraction*3.d0/4.d0/length_scale_factor
         mieabs=(mean_qabs_mie)*sphere_volume_fraction*3.d0/4.d0/length_scale_factor
         miesca=(mean_qext_mie)*sphere_volume_fraction*3.d0/4.d0/length_scale_factor
         if(target_shape.le.1) then
!            scaq=1.d0-0.5d0*(-sum(dif_boundary_sca(:,0))+sum(dif_boundary_sca(:,number_plane_boundaries+1)))
            scaq=1.d0-0.5d0*(-sum(dif_boundary_sca(:,0))+sum(dif_boundary_sca(:,number_plane_boundaries+1)))-q_eff_tot(2,1)
            absq=1.d0-q_eff_tot(2,1)
            scaq=max(1.d-5,scaq)
            absq=max(1.d-5,absq)
            if(target_shape.eq.0) then
               across=4.d0*product(target_dimensions(1:2))*length_scale_factor**2
            elseif(target_shape.eq.1) then
               across=pi*(target_dimensions(1)*length_scale_factor)**2
            endif
            h=4.d0*pi*dble(number_spheres)*length_scale_factor**3/(3.d0*across*sphere_volume_fraction)
            scacoef=-dlog(scaq)/h
            if(abs(mieabs).lt.1.d-7) then
               abscoef=0.d0
            else
               abscoef=-dlog(absq)/h
            endif
         else
!            scaq=0.5d0*(-sum(dif_boundary_sca(:,0))+sum(dif_boundary_sca(:,number_plane_boundaries+1)))
            scaq=0.5d0*(-sum(dif_boundary_sca(:,0))+sum(dif_boundary_sca(:,number_plane_boundaries+1)))+q_eff_tot(2,1)
            absq=q_eff_tot(2,1)
            h=length_scale_factor*(target_dimensions(1)-1.d0)**3/target_dimensions(1)**2
            root0=scaq
            do i=1,100
               root=root0
               func=1.d0-(1.d0-exp(-2.d0*root)*(1.d0+2.d0*root))/(2.d0*root*root)-scaq
               dfunc=exp(-2.d0*root)*(-1.d0+exp(2.d0*root)-2.d0*root*(1.d0+root))/root**3
               root0=root-func/dfunc
               if(abs(1.d0-root/root0).lt.1.d-6) exit
            enddo
            scacoef=root/h
            if(abs(mieabs).lt.1.d-7) then
               abscoef=0.d0
            else
               root0=absq
               do i=1,100
                  root=root0
                  func=1.d0-(1.d0-exp(-2.d0*root)*(1.d0+2.d0*root))/(2.d0*root*root)-absq
                  dfunc=exp(-2.d0*root)*(-1.d0+exp(2.d0*root)-2.d0*root*(1.d0+root))/root**3
                  root0=root-func/dfunc
                  if(abs(1.d0-root/root0).lt.1.d-6) exit
               enddo
               abscoef=root/h
            endif
         endif
         srat=scacoef/miesca
         if(abs(mieabs).lt.1.d-7) then
            arat=1.d0
         else
            arat=abscoef/mieabs
         endif
         end subroutine effective_extinction_coefficient_ratio

         subroutine read_sphere_data_input_file(mpi_comm)
         implicit none
         integer :: mpicomm,rank,istat,n
         integer, optional :: mpi_comm
         real(8) :: rtemp(4)
         complex(8) :: ctemp(2)
         character*256 :: parmval

         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)

         open(1,file=sphere_data_input_file)
         do n=1,number_spheres
            sphere_radius(n)=1.d0
            sphere_ref_index(1,n)=(1.d0,0.d0)
            sphere_ref_index(2,n)=(0.d0,0.d0)
            read(1,'(a)',iostat=istat) parmval
            if(istat.ne.0) then
               if(rank.eq.0) then
                  write(run_print_unit,'('' insufficient data in input file: '', i4,'' lines, need '',i4)') &
                     n,number_spheres
               endif
               stop
            endif
            read(parmval,*,iostat=istat) sphere_position(:,n)
            if(istat.ne.0) then
               if(rank.eq.0) then
                  write(run_print_unit,'('' read error in sphere data input file'')')
               endif
               stop
            endif
            read(parmval,*,iostat=istat) rtemp(1:4)
            if(istat.eq.0) sphere_radius(n)=rtemp(4)
            read(parmval,*,iostat=istat) rtemp(1:4),ctemp(1)
            if(istat.eq.0) sphere_ref_index(1,n)=ctemp(1)
            read(parmval,*,iostat=istat) rtemp(1:4),ctemp(1:2)
            if(istat.eq.0) then
               sphere_ref_index(2,n)=ctemp(2)
            else
               sphere_ref_index(2,n)=sphere_ref_index(1,n)
            endif
         enddo
         number_spheres=min(n,number_spheres)
         close(1)
         do n=1,number_spheres
            sphere_radius(n)=sphere_radius(n)*length_scale_factor
            sphere_position(:,n)=sphere_position(:,n)*length_scale_factor
            sphere_ref_index(:,n)=sphere_ref_index(:,n)*ref_index_scale_factor
         enddo
         end subroutine read_sphere_data_input_file

         subroutine generate_random_configuration(mpi_comm,skip_diffusion)
         implicit none
         logical :: skipdif
         logical, optional :: skip_diffusion
         integer :: mpicomm,rank,nsend
         integer, optional :: mpi_comm
         if(present(skip_diffusion)) then
            skipdif=skip_diffusion
         else
            skipdif=.false.
         endif
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         if(rank.eq.0) then
            call random_cluster_of_spheres(number_spheres,sphere_position,sphere_radius, &
               run_print_unit,ran_config_stat,ran_config_time_steps, &
               skip_diffusion=skipdif, &
               use_saved_values=use_previous_configuration)
            if(ran_config_stat.ge.3) then
               write(run_print_unit,'('' unable to generate random configuration'')')
               stop
            endif
         endif
         call mstm_mpi(mpi_command='barrier')
         nsend=number_spheres
         call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=sphere_radius, &
            mpi_number=nsend,mpi_rank=0,mpi_comm=mpicomm)
         nsend=number_spheres*3
         call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=sphere_position, &
            mpi_number=nsend,mpi_rank=0,mpi_comm=mpicomm)
         sphere_radius=sphere_radius*length_scale_factor
         sphere_position=sphere_position*length_scale_factor
         sphere_ref_index(:,1:number_spheres)=ref_index_scale_factor
         end subroutine generate_random_configuration

         subroutine scattering_matrix_calculation(amnp,scatmat,mpi_comm)
         implicit none
         logical :: singleorigin,iframe
         integer :: i,sy,sx,mpicomm
         integer, optional :: mpi_comm
         real(8) :: scatmat(scat_mat_mdim,scat_mat_ldim:scat_mat_udim),costheta,phi,csca, &
                    ky,kx,sintheta
         complex(8) :: amnp(*),ampmat(2,2)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         singleorigin=number_plane_boundaries.eq.0.and.single_origin_expansion
         iframe=singleorigin.and.incident_frame
         csca=(q_eff_tot(1,1)-q_eff_tot(2,1))*pi*cross_section_radius**2

         if(periodic_lattice) then
            call periodic_lattice_scattering(amnp,pl_sca,scat_mat=scatmat,krho_vec=rl_vec)
            return
         endif

         if(scattering_map_model.eq.0) then
            do i=scat_mat_ldim,scat_mat_udim
               costheta=cos(pi/180.d0*(scat_mat_amin+(scat_mat_amax-scat_mat_amin) &
                  *dble(i-scat_mat_ldim)/dble(scat_mat_udim-scat_mat_ldim)))
               if(costheta.eq.1.d0) costheta=0.9999999d0
               if(costheta.eq.-1.d0) costheta=-0.9999999d0
               if(i.lt.0) then
                  phi=incident_alpha_deg*pi/180.d0+pi
               else
                  phi=incident_alpha_deg*pi/180.d0
               endif
               if(number_plane_boundaries.eq.0) then
                  if(singleorigin) then
                     if(azimuthal_average) then
!                        call fosmcalc(12,s00,s02,sp22,sm22,costheta,scatmat(:,i),normalize_s11=.false.)
                        call fosmcalc(t_matrix_order,scat_mat_exp_coef(:,:,1),scat_mat_exp_coef(:,:,2), &
                           scat_mat_exp_coef(:,:,3),scat_mat_exp_coef(:,:,4), &
                           costheta,scatmat(:,i),normalize_s11=.false.)
                     else
                        call scatteringmatrix(amnp,t_matrix_order,costheta,phi,ampmat,scatmat(:,i),iframe, &
                           normalize_s11=.false.)
                     endif
                  else
                     call multiple_origin_scatteringmatrix(amnp,costheta,phi,csca,ampmat,scatmat(:,i))
                  endif
               else
                  costheta=-costheta
                  call multiple_origin_scatteringmatrix(amnp,costheta,phi,csca,ampmat,scatmat(1:16,i))
                  costheta=-costheta
                  call multiple_origin_scatteringmatrix(amnp,costheta,phi,csca,ampmat,scatmat(17:32,i))
               endif
            enddo
         else
            i=0
            do sy=-scattering_map_dimension,scattering_map_dimension
               ky=dble(sy)/dble(scattering_map_dimension)
               do sx=-scattering_map_dimension,scattering_map_dimension
                  if(sx*sx+sy*sy.gt.scattering_map_dimension**2) cycle
                  kx=dble(sx)/dble(scattering_map_dimension)
                  sintheta=kx*kx+ky*ky
                  sintheta=min(sintheta,.99999d0)
                  i=i+1
                  if(sx.eq.0.and.sy.eq.0) then
                     phi=0.d0
                  else
                     phi=datan2(ky,kx)
                  endif
                  costheta=-sqrt(1.d0-sintheta)
                  if(singleorigin) then
                     call scatteringmatrix(amnp,t_matrix_order,costheta,phi,ampmat,scatmat(1:16,i),iframe, &
                        normalize_s11=.false.)
                  else
                     call multiple_origin_scatteringmatrix(amnp,costheta,phi,csca,ampmat,scatmat(1:16,i))
                  endif
                  costheta=-costheta
                  if(singleorigin) then
                     call scatteringmatrix(amnp,t_matrix_order,costheta,phi,ampmat,scatmat(17:32,i),iframe, &
                        normalize_s11=.false.)
                  else
                     call multiple_origin_scatteringmatrix(amnp,costheta,phi,csca,ampmat, &
                        scatmat(17:32,i))
                  endif
               enddo
            enddo
         endif
         end subroutine scattering_matrix_calculation

         subroutine output_header(iunit,inputfile)
         implicit none
         integer :: iunit
         character*8 :: rundate
         character*10 :: runtime
         character*256 :: inputfile
         call date_and_time(date=rundate,time=runtime)
         run_date_and_time=trim(rundate)//' '//trim(runtime)
         write(iunit,'(''****************************************************'')')
         write(iunit,'(''****************************************************'')')
         write(iunit,'('' mstm calculation results'')')
         write(iunit,'('' date, time:'')')
         write(iunit,'(a)') run_date_and_time
         write(iunit,'('' input file:'')')
         write(iunit,'(a)') trim(inputfile)
         end subroutine output_header

         subroutine print_run_variables(iunit)
         implicit none
         integer :: iunit,i,n
         real(8) :: cb,r(2),t(2),a(2),tvol,svol

         write(iunit,'(''****************************************************'')')
         write(iunit,'('' input variables for run '',i5)') run_number
         if(trim(sphere_data_input_file).eq.'random_configuration') then
            write(iunit,'('' sphere positions randomly generated'')')
            if(target_shape.eq.0) then
               write(iunit,'('' rectangular target, unscaled half-widths in x,y,z'')')
               write(iunit,'(3es12.4)') target_dimensions
            elseif(target_shape.eq.1) then
               write(iunit,'('' cylindrical target, unscaled radius, half-thickness'')')
               write(iunit,'(3es12.4)') target_dimensions(1),target_dimensions(3)
            else
               write(iunit,'('' spherical target, unscaled radius'')')
               write(iunit,'(3es12.4)') target_dimensions(1)
            endif
            write(iunit,'('' sphere log-normal PSD sigma:'')')
            write(iunit,'(3es12.4)') psd_sigma
            write(iunit,'('' sphere volume fraction'')')
            if(.not.number_spheres_specified) then
               write(iunit,'(es12.4)') sphere_volume_fraction
            else
               call target_volume(tvol)
               tvol=tvol*length_scale_factor**3
               svol=4.d0*pi/3.d0*sum(sphere_radius(:)**3)
               write(iunit,'(es12.4)') svol/tvol
            endif
            if(ran_config_stat.eq.0) then
               write(iunit,'('' target configuration computed using random sampling + diffusion'')')
            elseif(ran_config_stat.eq.1) then
               write(iunit,'('' target configuration computed using layered sampling + diffusion'')')
            elseif(ran_config_stat.eq.2) then
               write(iunit,'('' target configuration computed initial HCP + diffusion'')')
            endif
            write(iunit,'('' number diffusion time steps:'',i5)') ran_config_time_steps
         else
            write(iunit,'('' sphere data input file:'')')
            write(iunit,'(a)') trim(sphere_data_input_file)
         endif
         write(iunit,'('' number spheres'')')
         write(iunit,'(i5)') number_spheres
         write(iunit,'('' length, ref index scale factors'')')
         write(iunit,'(3es15.7)') length_scale_factor,ref_index_scale_factor
         write(iunit,'('' volume cluster radius, area mean sphere radius, circumscribing radius, cross section radius'')')
         write(iunit,'(4es15.7)') vol_radius,area_mean_radius,circumscribing_radius,cross_section_radius
         if(print_sphere_data) then
            write(iunit,'('' sphere properties and associations'')')
            if(any_optically_active) then
               write(iunit,'(''   sphere    host   layer radius     x       y       z    '', &
                 &''     ref indx(L)             ref_indx(R)'')')
            else
               write(iunit,'(''   sphere    host   layer radius     x       y       z           ref indx'')')
            endif
            do n=1,number_spheres
               if(any_optically_active) then
                  write(iunit,'(3i8,4f8.3,4es12.4)') n,host_sphere(n),sphere_layer(n),sphere_radius(n), &
                  sphere_position(:,n),sphere_ref_index(1:2,n)
               else
                  write(iunit,'(3i8,4f8.3,4es12.4)') n,host_sphere(n),sphere_layer(n),sphere_radius(n), &
                  sphere_position(:,n),sphere_ref_index(1,n)
               endif
            enddo
         endif

         if(random_orientation) then
            write(iunit,'('' random orientation, estimated t matrix order:'')')
            write(iunit,'(i6)')  t_matrix_order
         else
            if(gaussian_beam_constant.ne.0.d0) then
               write(iunit,'('' incident Gaussian beam: 1/beam width, focal point'')')
               write(iunit,'(4es12.4)') gaussian_beam_constant,gaussian_beam_focal_point
            else
               write(iunit,'('' incident plane wave'')')
            endif
            if(incidence_average) then
               write(iunit,'('' Monte Carlo average over incident directions'')')
            else
               if(incident_beta_specified) then
                  write(iunit,'('' incident alpha, beta(deg)'')')
                  write(iunit,'(2es12.4)') incident_alpha_deg,incident_beta_deg
               else
                  write(iunit,'('' incident alpha(deg), incident sin(beta), incident direction'')')
                  write(iunit,'(2es12.4,i3)') incident_alpha_deg,incident_sin_beta,3-2*incident_direction
               endif
            endif
            if(single_origin_expansion) then
               write(iunit,'('' t matrix order:'')')
               write(iunit,'(i6)')  t_matrix_order
            endif
         endif
         if(reflection_model) then
            if(auto_absorption_sample_radius) then
               write(iunit,'('' particle layer reflectance/absorptance model, absorption sample radius fraction'')')
               write(iunit,'(es12.4)') absorption_sample_radius_fraction
            else
               write(iunit,'('' particle layer reflectance/absorptance model, absorption sample radius'')')
               write(iunit,'(es12.4)') absorption_sample_radius*length_scale_factor
            endif
         endif
         write(iunit,'('' layer 0 refractive index'')')
         write(iunit,'(4es12.4)') layer_ref_index(0)
         write(iunit,'('' number of plane boundaries '')')
         write(iunit,'(i3)') number_plane_boundaries
         if(number_plane_boundaries.gt.0) then
            write(iunit,'('' boundary, position, refractive index'')')
            do i=1,number_plane_boundaries
               write(iunit,'(i3,3es12.4)') i,plane_boundary_position(i),layer_ref_index(i)
            enddo
            if(.not.incidence_average) then
               cb=cos(incident_beta_deg*pi/180.d0)
               call boundary_energy_transfer(incident_sin_beta,incident_direction,r,t,a)
               write(iunit,'('' Fresnel boundary reflectance, transmittance, absorptance (par, perp)'')')
               write(iunit,'(3es12.4)') r(1),t(1),a(1)
               write(iunit,'(3es12.4)') r(2),t(2),a(2)
            endif
            if(number_singular_points.gt.0) then
               write(iunit,'('' GF singular points (in s)'')')
               do i=1,number_singular_points
                  write(iunit,'(2i5,es12.4,es20.10)') i,singular_point_polarization(i), &
                  singular_gf_value(i),singular_points(i)
               enddo
            endif
         endif
         if(periodic_lattice) then
            write(iunit,'('' periodic lattice cell width, incident lateral vector '')')
            write(iunit,'(4es15.7)') cell_width, incident_lateral_vector
         endif
         write(iunit,'('' max_iterations,solution_epsilon, mie_epsilon'')')
         write(iunit,'(i10,2es12.4)') max_iterations,solution_epsilon,mie_epsilon
         write(iunit,'('' maximum Mie order, number of equations:'')')
         write(iunit,'(i4,i10)') max_mie_order,number_eqns
         write(iunit,'('' mean sphere Mie extinction, absorption efficiencies'')')
         write(iunit,'(2es12.4)') mean_qext_mie,mean_qabs_mie
         if(fft_translation_option) then
            write(iunit,'('' fft translation option implemented'')')
            write(iunit,'('' cell width, cell volume fraction, cell dimension:'')')
            write(iunit,'(3i8,2es12.4)') cell_dim(1:3),cell_volume_fraction,d_cell
            write(iunit,'('' number of neighbor nodes, node order:'')')
            write(iunit,'(2i5)') number_neighbor_nodes,node_order
         endif

         write(iunit,'(''****************************************************'')')
         write(iunit,'('' calculation results for run '')')
         write(iunit,*) run_number
         call flush(iunit)
         end subroutine print_run_variables

         subroutine print_calculation_results(fout)
         implicit none
         integer :: outunit,n,i,j,smvec(6),sx,sy,s,smvec0(16),nsmat,smvecp(16)
         real(8) :: smt(16),kx,ky,s11scale,r(2),t(2),a(2),scacoef,scarat, &
            abscoef,absrat,rl(2),al(2),tl(2)
         character*2 :: smlabel(16)
         character*256 :: fout,chartemp
         smvec=(/1,5,6,11,15,16/)
         smvec0=(/1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16/)
         smlabel=(/'11','21','31','41','12','22','32','42','13','23','33','43','14','24','34','44'/)
         if(fout(1:7).eq.'console') then
            outunit=6
         else
            outunit=2
            open(2,file=trim(fout))
            do
               read(2,'(a)') chartemp
               if(trim(chartemp).eq.run_date_and_time) exit
            enddo
            do
               read(2,'(a)') chartemp
               if(chartemp(1:28).eq.' calculation results for run') then
                  read(2,*) i
                  if(i.eq.run_number) exit
               endif
            enddo
         endif
         if(configuration_average) then
            write(outunit,'('' averages collected: sample:'')')
            write(outunit,'(i5)') random_configuration_number*n_configuration_groups
         endif
         if(incidence_average) then
            write(outunit,'('' incidence averages collected: sample:'')')
            write(outunit,'(i5)') incident_direction_number*n_configuration_groups
         endif

         write(outunit,'('' number iterations, error, solution time '')')
         write(outunit,'(i6,3es12.4)') solution_iterations,solution_error,solution_time
         if(number_plane_boundaries.gt.0) then
            if(maxval(error_codes).ne.0) write(outunit,'('' warning: problems encountered with surface interaction calculations'')')
            if(error_codes(1).ne.0) write(outunit,'('' iterative surface GF algorithm did not converge'')')
            if(error_codes(2).ne.0) then
               write(outunit,'('' integration along real s axis did not converge'')')
               write(outunit,'('' increase real_axis_integration_limit from current value of '',es10.2)')  &
               real_axis_integration_limit
               write(outunit,'('' or increase integration_limit_epsilon from current value of '',es10.2, &
               &  '' and see if results change'')') integration_limit_epsilon
            endif
            if(error_codes(3).ne.0) write(outunit,'('' subdivided integration interval below 1d-12'')')
            if(error_codes(4).ne.0) then
                write(outunit,'('' maximum subdivision reached in GK integration algorithm'')')
                write(outunit,'('' increase maximum_integration_subdivisions from current value of '', i5, &
                &'' and see if results change'')') maximum_integration_subdivisions
            endif
            call boundary_energy_transfer(incident_sin_beta,incident_direction,r,t,a)
         else
            r=0.d0
            a=0.d0
            t=1.d0
         endif
         if(periodic_lattice) then
            if(maxval(pl_error_codes).ne.0) write(outunit,'('' warning: problems encountered with periodic lattice calculations'')')
            if(pl_error_codes(1).ne.0) write(outunit,'('' integration formulas for FS PL DMGF did not converge''&
            &'' in subroutine swf_lattice_sum'')')
            if(pl_error_codes(2).ne.0) write(outunit,'('' RS series for PL DMGF did not converge in subroutine''&
            &'' reciprocal_space_swf_lattice_sum'')')
            if(pl_error_codes(3).ne.0) write(outunit,'('' reciprocal space series for PL DMGF did not converge''&
                &'' in subroutine plane_boundary_lattice_interaction'')')
            if(pl_error_codes(4).ne.0) write(outunit,'('' integration did not converge in subroutine q2db'')')
            if(pl_error_codes(5).ne.0) write(outunit,'('' integration did not converge in subroutine q1dbnosource'')')
            if(pl_error_codes(6).ne.0) write(outunit,'('' series in s did not converge in subroutine swfyzlatticesum'')')
         endif

         if(random_orientation) then
            write(outunit,'('' t matrix order:'')')
            write(outunit,'(i6)')  t_matrix_order
         endif

         if(print_sphere_data) then
            write(outunit,'('' sphere extinction, absorption, volume absorption efficiencies (unpolarized incidence)'')')
            write(outunit,'(''   sphere    Qext        Qabs        Qvabs'')')
            do n=1,number_spheres
               write(outunit,'(i8,3es12.4)') n,q_eff(1:2,1,n),q_vabs(1,n)
            enddo
         endif

         if(periodic_lattice.or.reflection_model) then
            rl=pl_sca(:,2)-boundary_ext(:,0)+r(:)
            al=surface_absorptance(:)
            tl=1.d0-rl-al
            write(outunit,'('' unit cell reflectance, absorptance, transmittance (unpol, par, perp)'')')
               write(outunit,'(9es12.5)') 0.5d0*sum(rl),0.5d0*sum(al),0.5d0*sum(tl), &
                 (rl(i),al(i),tl(i),i=1,2)
!            write(outunit,'('' unit cell reflectance, absorptance, transmittance (unpol, par, perp)'')')
!               write(outunit,'(9es12.5)') 0.5d0*sum(pl_sca(:,2)-boundary_ext(:,0)+r(:)),q_eff_tot(2,1), &
!                 .5d0*sum(pl_sca(:,1)+boundary_ext(:,number_plane_boundaries+1)+t(:)), &
!                 (pl_sca(i,2)-boundary_ext(i,0)+r(i),q_eff_tot(2,i+1), &
!                  pl_sca(i,1)+boundary_ext(i,number_plane_boundaries+1)+t(i),i=1,2)

!            write(outunit,'('' down, up scattering fraction, total scattering cross section (unpol)'')')
!               write(outunit,'(9es12.5)') 0.5d0*sum(pl_sca(:,2)),.5d0*sum(pl_sca(:,1)), &
!                 pi*(0.5d0*sum(pl_sca(:,2))+.5d0*sum(pl_sca(:,1)))*cross_section_radius**2
            if(configuration_average.and.single_origin_expansion) then
               scacoef=(-0.5d0*sum(dif_boundary_sca(:,0))+0.5d0*sum(dif_boundary_sca(:,number_plane_boundaries+1))) &
                    *pi*cross_section_radius**2
               scarat=scacoef/(q_eff_tot(3,1)*pi*cross_section_radius**2)
               write(outunit,'('' down, up diffuse scattering fraction, optically thin dependent/independent ratio'')')
               write(outunit,'(9es12.5)') -0.5d0*sum(dif_boundary_sca(:,0)), &
                  0.5d0*sum(dif_boundary_sca(:,number_plane_boundaries+1)), &
                  scarat
               call effective_extinction_coefficient_ratio(scacoef,abscoef,scarat,absrat)
               write(outunit,'('' dimensionless scattering, absorption coefficients, dependent/independent ratios'')')
               write(outunit,'(4es12.5)') scacoef,abscoef,scarat,absrat
            endif
         else
            if(qeff_dim.eq.1) then
               write(outunit,'('' total extinction, absorption, scattering efficiencies (unpolarized incidence)'')')
               write(outunit,'(3es12.4)') q_eff_tot(1:3,1)
            else
               write(outunit,'('' total extinction, absorption, scattering efficiencies (unpol, par, perp incidence)'')')
               write(outunit,'(9es12.4)') q_eff_tot(1:3,1:3)
            endif
            if(.not.random_orientation) then
               if(number_plane_boundaries.gt.0) then
                  write(outunit,'(''  down and up extinction efficiencies (unpol, par, perp)'')')
                  write(outunit,'(16es12.4)')  &
                     0.5d0*sum(boundary_ext(:,0)),0.5d0*sum(-boundary_ext(:,number_plane_boundaries+1)), &
                     (boundary_ext(i,0),-boundary_ext(i,number_plane_boundaries+1),i=1,2)
               endif
               if(calculate_up_down_scattering) then
                  write(outunit,'(''  down and up hemispherical scattering efficiencies (unpol, par, perp) '')')
                  write(outunit,'(16es12.4)')  &
                     0.5d0*sum(-boundary_sca(:,0)),0.5d0*sum(boundary_sca(:,number_plane_boundaries+1)), &
                     (-boundary_sca(i,0),boundary_sca(i,number_plane_boundaries+1),i=1,2)
               endif
               if(number_plane_boundaries.gt.0.and.number_singular_points.gt.0) then
                  write(outunit,'(''  waveguide scattering efficiencies (unpol, par, perp)  '')')
                  write(outunit,'(16es12.4)')  &
                     0.5d0*sum(evan_sca(:)),(evan_sca(i),i=1,2)
               endif
            endif
         endif

         if(calculate_scattering_matrix.and..not.periodic_lattice) then
            if(normalize_s11) then
               s11scale=1.d0/((cross_section_radius**2)*2.d0*q_eff_tot(3,1))
            else
!               s11scale=(cross_section_radius**2)*q_eff_tot(3,1)
               s11scale=pi
!               if(.not.random_orientation) s11scale=pi*s11scale
            endif
            if(((.not.any_optically_active).and.random_orientation).or.azimuthal_average) then
               nsmat=6
               smvecp(1:6)=smvec(1:6)
            else
               nsmat=16
               smvecp(1:16)=smvec0(1:16)
            endif
            if(random_orientation) then
               write(outunit,'('' total scattering'')')
               call print_scat_mat_header()
               do i=scat_mat_ldim,scat_mat_udim
                  smt=scaled_scat_mat(scat_mat(1:16,i))
                  smt(1)=smt(1)*s11scale
                  call print_scat_mat_row(smt)
               enddo
            else
               if(scattering_map_model.eq.0) then
                  if(number_plane_boundaries.eq.0) then
                     if(azimuthal_average) then
                        write(outunit,'('' azimuthal averaged scattering matrix'')')
                     else
                        if(incident_frame) then
                           write(outunit,'('' scattering matrix in incident plane: 0 deg = incident direction'')')
                        else
                           write(outunit,'('' scattering matrix in incident plane: 0 deg= z axis'')')
                        endif
                     endif
                     call print_scat_mat_header()
                     do i=scat_mat_ldim,scat_mat_udim
                        smt=scaled_scat_mat(scat_mat(1:16,i))
                        smt(1)=smt(1)*s11scale
                        call print_scat_mat_row(smt)
                     enddo
                     if((configuration_average.or.incidence_average).and.single_origin_expansion) then
                        if(normalize_s11) then
                           s11scale=1.d0/(0.5d0*(-sum(dif_boundary_sca(:,0)) &
                              +sum(dif_boundary_sca(:,number_plane_boundaries+1)))*2.d0*cross_section_radius**2)
                        else
                           s11scale=pi
                        endif
                        write(outunit,'('' diffuse scattering matrix '')')
                        call print_scat_mat_header(no_numbers=.true.)
                        do i=scat_mat_ldim,scat_mat_udim
                           smt=scaled_scat_mat(dif_scat_mat(1:16,i))
                           smt(1)=smt(1)*s11scale
                           call print_scat_mat_row(smt)
                        enddo
                     endif
                  else
                     write(outunit,'('' scattering matrix in incident plane'')')
                     write(outunit,'('' reflection'')')
                     call print_scat_mat_header()
                     do i=scat_mat_ldim,scat_mat_udim
                        smt=scaled_scat_mat(scat_mat(1:16,i))
                        smt(1)=smt(1)*s11scale
                        call print_scat_mat_row(smt)
                     enddo
                     write(outunit,'('' transmission'')')
                     call print_scat_mat_header(no_numbers=.true.)
                     do i=scat_mat_ldim,scat_mat_udim
                        smt=scaled_scat_mat(scat_mat(17:32,i))
                        smt(1)=smt(1)*s11scale
                        call print_scat_mat_row(smt)
                     enddo
                  endif
               else
                  write(outunit,'('' backward hemisphere scattering'')')
                  write(outunit,'(''    kx      ky '',$)')
                  do j=1,4
                     do i=1,4
                        write(outunit,'(''     '',2i1,''     '',$)') i,j
                     enddo
                  enddo
                  write(outunit,*)
                  s=0
                  do sy=-scattering_map_dimension,scattering_map_dimension
                     ky=dble(sy)/dble(scattering_map_dimension)
                     do sx=-scattering_map_dimension,scattering_map_dimension
                        kx=dble(sx)/dble(scattering_map_dimension)
                        if(sx*sx+sy*sy.gt.scattering_map_dimension**2) cycle
                        s=s+1
                        smt=scaled_scat_mat(scat_mat(1:16,s))
                        write(outunit,'(2f9.5,$)') kx,ky
                        do j=1,16
                           write(outunit,'(es12.4,$)') smt(smvec0(j))
                        enddo
                        write(outunit,*)
                     enddo
                  enddo
                  write(outunit,'('' forward hemisphere scattering'')')
                  write(outunit,'(''    kx      ky '',$)')
                  do j=1,4
                     do i=1,4
                        write(outunit,'(''     '',2i1,''     '',$)') i,j
                     enddo
                  enddo
                  write(outunit,*)
                  s=0
                  do sy=-scattering_map_dimension,scattering_map_dimension
                     ky=dble(sy)/dble(scattering_map_dimension)
                     do sx=-scattering_map_dimension,scattering_map_dimension
                        if(sx*sx+sy*sy.gt.scattering_map_dimension**2) cycle
                        kx=dble(sx)/dble(scattering_map_dimension)
                        s=s+1
                        smt=scaled_scat_mat(scat_mat(17:32,s))
                        write(outunit,'(2f9.5,$)') kx,ky
                        do j=1,16
                           write(outunit,'(es12.4,$)') smt(smvec0(j))
                        enddo
                        write(outunit,*)
                     enddo
                  enddo
               endif
            endif
            if(azimuthal_average.and.(.not.random_orientation)) then
               s11scale=1.d0/scat_mat_exp_coef(1,0,1)
               scat_mat_exp_coef=scat_mat_exp_coef*s11scale
               write(outunit,'('' azimuthal averaged scattering matrix expansion coefficients'')')
               write(outunit,'(''    n       11          44          12          34         22p         22m'')')
               do n=0,2*t_matrix_order
                  write(outunit,'(i5,6es12.4)') n,scat_mat_exp_coef(1,n,1),scat_mat_exp_coef(16,n,1), &
                     0.5d0*(scat_mat_exp_coef(2,n,2)+scat_mat_exp_coef(5,n,2)), &
                     0.5d0*(scat_mat_exp_coef(12,n,2)+scat_mat_exp_coef(15,n,2)), &
                     scat_mat_exp_coef(6,n,3),scat_mat_exp_coef(6,n,4)
               enddo
            endif
         endif

         if(calculate_scattering_matrix.and.periodic_lattice) then
            write(outunit,'('' scattering by periodic lattice at reciprocal lattice directions'')')
            write(outunit,'('' backward hemisphere scattering'')')
            write(outunit,'('' number directions, number SM elements'')')
            write(outunit,'(2i6)') number_rl_dirs(1),16
            write(outunit,'(''    kx      ky '',$)')
            do j=1,4
               do i=1,4
                  write(outunit,'(''     '',2i1,''     '',$)') i,j
               enddo
            enddo
            write(outunit,*)
            do i=1,number_rl_dirs(1)
               smt=scaled_scat_mat(scat_mat(1:16,i))
               write(outunit,'(2f9.5,$)') rl_vec(1:2,i)/dble(layer_ref_index(0))
               do j=1,16
                  write(outunit,'(es12.4,$)') smt(smvec0(j))
               enddo
               write(outunit,*)
            enddo
            write(outunit,'('' forward hemisphere scattering'')')
            write(outunit,'('' number directions, number SM elements'')')
            write(outunit,'(2i6)') number_rl_dirs(2),16
            write(outunit,'(''    kx      ky '',$)')
            do j=1,4
               do i=1,4
                  write(outunit,'(''     '',2i1,''     '',$)') i,j
               enddo
            enddo
            write(outunit,*)
            do i=1,number_rl_dirs(2)
               smt=scaled_scat_mat(scat_mat(17:32,i))
               write(outunit,'(2f9.5,$)') rl_vec(1:2,i)/dble(layer_ref_index(number_plane_boundaries))
               do j=1,16
                  write(outunit,'(es12.4,$)') smt(smvec0(j))
               enddo
               write(outunit,*)
            enddo
         endif

         if(outunit.ne.6) close(outunit)

         contains
            subroutine print_scat_mat_header(no_numbers)
            implicit none
            logical, optional :: no_numbers
            if(present(no_numbers)) then
               if(.not.no_numbers) then
                  write(outunit,'('' number directions, number SM elements:'')')
                  write(outunit,'(2i6)') scat_mat_udim-scat_mat_ldim+1,nsmat
               endif
            else
               write(outunit,'('' number directions, number SM elements:'')')
               write(outunit,'(2i6)') scat_mat_udim-scat_mat_ldim+1,nsmat
            endif
            write(outunit,'(''   theta'',$)')
            do i=1,nsmat
               write(outunit,'(''     '',a2,''     '',$)') smlabel(smvecp(i))
            enddo
            write(outunit,*)
            end subroutine print_scat_mat_header

            subroutine print_scat_mat_row(smt)
            implicit none
            real(8) :: smt(16)
!            write(outunit,'(f8.2,$)') dble(i)*180.d0/dble(scat_mat_udim-scat_mat_ldim)
            write(outunit,'(f8.2,$)') scat_mat_amin &
               + dble(i-scat_mat_ldim)/dble(scat_mat_udim-scat_mat_ldim)*(scat_mat_amax-scat_mat_amin)
            do j=1,nsmat
               write(outunit,'(es12.4,$)') smt(smvecp(j))
            enddo
            write(outunit,*)
            end subroutine print_scat_mat_row

         end subroutine print_calculation_results

         function scaled_scat_mat(s)
         real(8) :: scaled_scat_mat(16),s(16)
         if(s(1).eq.0.d0) then
            scaled_scat_mat=0.d0
         else
            scaled_scat_mat(1)=s(1)
            scaled_scat_mat(2:16)=s(2:16)/s(1)
         endif
         end function scaled_scat_mat

         subroutine scat_mat_to_phase_mat(smat,u,up,phi,smatrot)
         implicit none
         real(8) smat(4,4),u,up,phi,smatrot(4,4),s,sp,us,ss,csig,&
            ssig,c2sig,s2sig,mat1(4,4),mat2(4,4)

         s=sqrt(1.d0-u*u)
         sp=sqrt(1.d0-up*up)
         us = u*up + s*sp*cos(phi)
         ss=sqrt(1.d0-us*us)
         if(up.eq.1.d0) then
            csig=cos(phi)
            ssig=-sin(phi)
         elseif(up.eq.-1.d0) then
            csig=-cos(phi)
            ssig=-sin(phi)
         else
            csig=(-u + up*us)/(sp*ss)
            ssig=-s*sin(phi)/ss
         endif
         c2sig = 2.d0*csig*csig - 1.d0
         s2sig = 2.d0*ssig*csig
         mat1(:,1)=(/1.d0,0.d0,0.d0,0.d0/)
         mat1(:,2)=(/0.d0, c2sig, -s2sig, 0.d0/)
         mat1(:,3)=(/0.d0, s2sig, c2sig, 0.d0/)
         mat1(:,4)=(/0.d0,0.d0,0.d0,1.d0/)
         if(u.eq.1.d0) then
            csig=cos(phi)
            ssig=-sin(phi)
         elseif(u.eq.-1.d0) then
            csig=-cos(phi)
            ssig=-sin(phi)
         else
            csig=(-up + u*us)/(s*ss)
            ssig=-sp*sin(phi)/ss
         endif
         c2sig = 2.d0*csig*csig - 1.d0
         s2sig = 2.d0*ssig*csig
         mat2(:,1)=(/1.d0,0.d0,0.d0,0.d0/)
         mat2(:,2)=(/0.d0, c2sig, -s2sig, 0.d0/)
         mat2(:,3)=(/0.d0, s2sig, c2sig, 0.d0/)
         mat2(:,4)=(/0.d0,0.d0,0.d0,1.d0/)
         mat1=matmul(smat,mat1)
         smatrot=matmul(mat2,mat1)
         end subroutine scat_mat_to_phase_mat

         subroutine checkpositions()
         implicit none
         logical :: check
         integer :: i,j,imin,jmin
         real(8) :: r,amax,amin,rmingap

         check=.true.
         rmingap=-1.d10
         do i=1,number_spheres-1
            do j=i+1,number_spheres
               amax=max(sphere_radius(i),sphere_radius(j))
               amin=min(sphere_radius(i),sphere_radius(j))
               r=sqrt(sum((sphere_position(:,i)-sphere_position(:,j))**2))
               if(r.ge.amax+amin) then
                  cycle
               else
                  if(amin+amax-r.gt.rmingap) then
                     rmingap=amin+amax-r
                     imin=i
                     jmin=j
                  endif
               endif
               if(r.le.amax-amin) cycle
               check=.false.
!               write(run_print_unit,'('' spheres '',i4,'' and '',i4,'' intersect'')') i,j
!               write(2,'('' spheres '',i4,'' and '',i4,'' intersect'')') i,j
            enddo
         enddo
         if(.not.check) then
            write(run_print_unit,'('' warning: sphere-sphere intersections detected, max overlap:'',es12.4, &
               &''  Results might be garbage!'')') rmingap
            write(run_print_unit,'('' positions:'',i5,3es12.4,i5,3es12.4)') imin,sphere_position(:,imin), &
               jmin,sphere_position(:,jmin)
            call flush(run_print_unit)
         endif
         check=.true.
         rmingap=-1.d10
         do i=1,number_spheres
            do j=1,number_plane_boundaries
               if(abs(sphere_position(3,i)-plane_boundary_position(j)).ge.sphere_radius(i)) cycle
               check=.false.
               rmingap=max(rmingap,sphere_radius(i)-abs(sphere_position(3,i)-plane_boundary_position(j)))
!               write(run_print_unit,'('' sphere '',i4,'' and plane boundary'',i4,'' intersect'')') i,j
!               write(2,'('' sphere '',i4,'' and plane boundary'',i4,'' intersect'')') i,j
            enddo
         enddo
         if(.not.check) then
            write(run_print_unit,'('' warning: sphere-plane boundary intersections detected, max overlap:'',es12.4, &
               &''  Results might be garbage!'')') rmingap
            write(run_print_unit,'('' positions:'',i5,3es12.4,i5,3es12.4)') imin,sphere_position(:,imin), &
               jmin,sphere_position(:,jmin)
            call flush(run_print_unit)
         endif
         end subroutine checkpositions

         subroutine set_string_to_int_variable(sentvarvalue, &
               ivarvalue,var_operation)
         implicit none
         integer :: itemp
         integer, pointer :: ivarvalue
         character*256 :: sentvarvalue,varop,intfile
         character*(*), optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         write(intfile,'(a)') sentvarvalue
         read(intfile,*) itemp
         if(varop(1:6).eq.'assign') then
            ivarvalue=itemp
         elseif(varop(1:3).eq.'add') then
            ivarvalue=ivarvalue+itemp
         endif
         end subroutine set_string_to_int_variable

         subroutine set_string_to_real_variable(sentvarvalue, &
               rvarvalue,var_operation)
         implicit none
         real(8) :: rtemp
         real(8), pointer :: rvarvalue
         character*256 :: sentvarvalue,varop,intfile
         character*256, optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         write(intfile,'(a)') sentvarvalue
         read(intfile,*) rtemp
         if(varop(1:6).eq.'assign') then
            rvarvalue=rtemp
         elseif(varop(1:3).eq.'add') then
            rvarvalue=rvarvalue+rtemp
         endif
         end subroutine set_string_to_real_variable

         subroutine set_string_to_real_array_variable(sentvarvalue, &
               rvarvalue,var_operation,var_len)
         implicit none
         integer :: varlen,i,ierr
         integer, optional :: var_len
         real(8) :: rtemp(4)
         real(8), pointer :: rvarvalue(:)
         character*256 :: sentvarvalue,varop,intfile
         character*256, optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         if(present(var_len)) then
            varlen=var_len
         else
            varlen=1
         endif
         write(intfile,'(a)') sentvarvalue
         do i=1,varlen
            read(intfile,*,iostat=ierr) rtemp(1:i)
            if(ierr.ne.0) then
               rtemp(i:varlen)=rtemp(i-1)
               exit
            endif
         enddo
!         read(intfile,*) rtemp(1:varlen)
         if(varop(1:6).eq.'assign') then
            rvarvalue=rtemp(1:varlen)
         elseif(varop(1:3).eq.'add') then
            rvarvalue=rvarvalue+rtemp
         endif
         end subroutine set_string_to_real_array_variable

         subroutine set_string_to_cmplx_variable(sentvarvalue, &
               cvarvalue,var_operation)
         implicit none
         complex(8) :: ctemp
         complex(8), pointer :: cvarvalue
         character*256 :: sentvarvalue,varop,intfile
         character*256, optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         write(intfile,'(a)') sentvarvalue
         read(intfile,*) ctemp
         if(varop(1:6).eq.'assign') then
            cvarvalue=ctemp
         elseif(varop(1:3).eq.'add') then
            cvarvalue=cvarvalue+ctemp
         endif
         end subroutine set_string_to_cmplx_variable

         subroutine set_string_to_logical_variable(sentvarvalue, &
               lvarvalue,var_operation)
         implicit none
         logical :: ltemp
         logical, pointer :: lvarvalue
         character*256 :: sentvarvalue,varop,intfile
         character*256, optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         write(intfile,'(a)') sentvarvalue
         read(intfile,*) ltemp
         if(varop(1:6).eq.'assign') then
            lvarvalue=ltemp
         endif
         end subroutine set_string_to_logical_variable

         subroutine set_string_to_logical_array_variable(sentvarvalue, &
            lvarvalue,var_operation,var_len)
         implicit none
         logical :: ltemp(5)
         logical, pointer :: lvarvalue(:)
         integer :: i,varlen,ierr
         integer, optional :: var_len
         character*256 :: sentvarvalue,varop,intfile
         character*256, optional :: var_operation
         if(present(var_operation)) then
            varop=var_operation(:index(var_operation,' '))
         else
            varop='assign'
         endif
         if(present(var_len)) then
            varlen=var_len
         else
            varlen=1
         endif
         write(intfile,'(a)') sentvarvalue
         do i=1,varlen
            read(intfile,*,iostat=ierr) ltemp(1:i)
            if(ierr.ne.0) then
               ltemp(i:varlen)=ltemp(i-1)
               exit
            endif
         enddo
         if(varop(1:6).eq.'assign') then
            lvarvalue=ltemp(1:varlen)
         endif
         end subroutine set_string_to_logical_array_variable

      end module inputinterface
      program main
      use inputinterface
      use solver
      use mpidefs
      use intrinsics
      use specialfuncs
      use spheredata
      implicit none
      integer :: looplevel,rank,numprocs, &
                 readstat(1),i,numargs,istat,numberinputlines,n
      character*256 :: inputfile,inputline,oldoutputfile
      character*256, allocatable :: inputfiledata(:)
      data oldoutputfile/' '/

      call mstm_mpi(mpi_command='init')
      call mstm_mpi(mpi_command='rank',mpi_rank=rank)
      call mstm_mpi(mpi_command='size',mpi_size=numprocs)

      do i=0,numprocs-1
         if(rank.eq.i) then
            numargs=mstm_nargs()
            if(numargs.eq.0) then
               inputfile='mstm.inp'
            else
               call mstm_getarg(inputfile)
            endif
            input_file=inputfile
            open(2,file=inputfile)
            numberinputlines=0
            istat=0
            do while(istat.eq.0)
               read(2,'(a)',iostat=istat) inputline
               numberinputlines=numberinputlines+1
               if(trim(inputline).eq.'end_of_options') exit
            enddo
            if(istat.ne.0) numberinputlines=numberinputlines+1
            allocate(inputfiledata(numberinputlines))
            rewind(2)
            istat=0
            n=0
            do while(istat.eq.0)
               read(2,'(a)',iostat=istat) inputline
               n=n+1
               inputfiledata(n)=inputline
               if(trim(inputline).eq.'end_of_options') exit
            enddo
            if(istat.ne.0) inputfiledata(numberinputlines)='end_of_options'
            close(2)
         endif
         call mstm_mpi(mpi_command='barrier')
      enddo

      repeat_run=.true.
      first_run=.true.

      do while(repeat_run)
         readstat=0
         do i=0,numprocs-1
            if(i.eq.rank) then
               call inputdata(inputfiledata,read_status=readstat(1))
            endif
            call mstm_mpi(mpi_command='barrier')
         enddo
         if(oldoutputfile.ne.output_file) then
            first_run=.true.
            run_number=0
            oldoutputfile=output_file
         endif
         if(rank.eq.0) then
            if(first_run) then
               if(append_output_file) then
                  open(2,file=output_file,access='append')
               else
                  open(2,file=output_file)
                  close(2,status='delete')
                  open(2,file=output_file)
               endif
               call output_header(2,inputfile)
               close(2)
            endif
         endif

         if(n_nest_loops.eq.0) then
            run_number=run_number+1
            if(configuration_average) then
               call configuration_average_calling_program()
            elseif(incidence_average) then
               call incidence_average_calling_program()
            else
               call main_calling_program()
            endif
         else
            looplevel=1
            call nested_loop(looplevel,rank)
         endif
         n_nest_loops=0
      enddo
!      if(temporary_pos_file.and.rank.eq.0) then
!         open(20,file='temp_pos.dat')
!         close(20,status='delete')
!      endif
      call mstm_mpi(mpi_command='barrier')
      call mstm_mpi(mpi_command='finalize')

      contains

         recursive subroutine nested_loop(looplevel,rank)
         implicit none
         logical :: continueloop
         integer :: looplevel,varposition,loopindex,rank
         integer, pointer :: i_loop_var_pointer
         real(8) :: maxdif,loopdif
         real(8), pointer :: r_loop_var_pointer
         complex(8), pointer :: c_loop_var_pointer
         character*256 :: varlabel
         character*1 :: vartype

         varlabel=loop_var_label(looplevel)
         vartype=loop_var_type(looplevel)
         varposition=loop_sphere_number(looplevel)
         if(vartype.eq.'i') then
            maxdif=abs(i_var_stop(looplevel)-i_var_start(looplevel))
            call variable_list_operation(varlabel, &
               var_position=varposition, &
               i_var_pointer=i_loop_var_pointer)
               i_loop_var_pointer=i_var_start(looplevel)
         elseif(vartype.eq.'r') then
            maxdif=abs(r_var_stop(looplevel)-r_var_start(looplevel))
            call variable_list_operation(varlabel, &
               var_position=varposition, &
               r_var_pointer=r_loop_var_pointer)
               r_loop_var_pointer=r_var_start(looplevel)
         elseif(vartype.eq.'c') then
            maxdif=cdabs(c_var_stop(looplevel)-c_var_start(looplevel))
            call variable_list_operation(varlabel, &
               var_position=varposition, &
               c_var_pointer=c_loop_var_pointer)
               c_loop_var_pointer=c_var_start(looplevel)
         endif

         loopindex=0
         continueloop=.true.
         do while(continueloop)
            loopindex=loopindex+1
            if(looplevel.eq.n_nest_loops) then
               run_number=run_number+1
!               if(loopindex.eq.1) then
!                  run_number=run_number+1
!               else
!                  if(.not.configuration_average) run_number=run_number+1
!               endif
               if(configuration_average) then
                  call configuration_average_calling_program()
               else
                  call main_calling_program()
               endif
            else
               call nested_loop(looplevel+1,rank)
            endif
            if(vartype.eq.'i') then
               i_loop_var_pointer=i_loop_var_pointer+i_var_step(looplevel)
               loopdif=abs(i_loop_var_pointer-i_var_start(looplevel))
            elseif(vartype.eq.'r') then
               r_loop_var_pointer=r_loop_var_pointer+r_var_step(looplevel)
               loopdif=abs(r_loop_var_pointer-r_var_start(looplevel))
            elseif(vartype.eq.'c') then
               c_loop_var_pointer=c_loop_var_pointer+c_var_step(looplevel)
               loopdif=cdabs(c_loop_var_pointer-c_var_start(looplevel))
            endif
!            if(loopindex.gt.1000) exit
            if(loopdif-maxdif.gt.1.d-6) continueloop=.false.
         enddo
         end subroutine nested_loop

      end program main
