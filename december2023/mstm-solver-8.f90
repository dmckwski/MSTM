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
            sphere_qeff,solution_status,sphere_excitation_list)
         implicit none
         logical :: firstrun,initialize,itersoln,continueloop,exlist(number_spheres)
         logical, optional :: sphere_excitation_list(number_spheres)
         integer :: iter,niter,istat,rank,maxiter,rank0,&
                    numprocs,mpicomm,i,nblkt,l,k,q,ka,la,nsolns, &
                    n,m,p,nssoln,kq,ns
         integer, save :: pcomm,prank,pgroup,ppsoln,pcomm0
         integer, optional :: mpi_comm,procs_per_soln,max_iterations,solution_status
         real(8) :: r0(3),maxerr,qeffi(3,number_spheres), &
            dqeffi(3,number_spheres),qeff(3),qeffold(3),time0,time1,timepersoln,timeleft, &
            solneps,conveps,solnerr,converr(1),qteff(3), &
            ttime(0:6),dqteff(3),converri,qeffiold(3)
         real(8), allocatable :: sexp(:,:),scexp(:,:)
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
         if(present(sphere_excitation_list)) then
            exlist=sphere_excitation_list
         else
            exlist=.true.
         endif

         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank0)
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

         if(rank0.eq.0) then
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
!if(prank.eq.0) then
!write(*,'('' pgroup,l,k,q:'',4i5)') pgroup,l,k,q
!call flush(6)
!endif
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
!                        sphere_translation_list=exlist)

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
!                     dqeffi(3,:)=dqeffi(1,:)-dqeffi(2,:)
! patch 10-22.  used new formula for qeff(3) in qefficiencyfactor SR.  Assumes ri_medium is real
!                     dqeffi(3,:)=0.d0
                     amnp0=0.d0
                     call merge_to_common_origin(l,amnpkq, &
                        amnp0, &
                        number_rhs=1, &
                        origin_position=r0, &
                        mpi_comm=pcomm, &
                        sphere_translation_list=exlist)
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
if(rank.eq.0) then
ttime(1)=mstm_mpi_wtime()
endif

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

!if(rank.eq.0) then
!ttime(2)=mstm_mpi_wtime()
!write(*,'('' timings:'',2es14.5)') ttime(1)-ttime(0),ttime(2)-ttime(1)
!endif

            if(rank.eq.0) then
               qeff=0.d0
               do i=1,number_spheres
                  if(.not.exlist(i)) cycle
                  qeff(:)=qeff(:)+qeffi(:,i)*sphere_radius(i)**2/vol_radius**2
               enddo

               converr=qeff(1)-qeffold(1)
               qeffold=qeff
               nsolns=nsolns+(l+l+1)*2

               if(rank0.eq.0) then
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
                  write(run_print_unit,'(i4,i5,4e13.5,2f8.2,a4)') l,maxiter,qeff(1:2), &
                     converr,dqeffi(1,1),timepersoln,timeleft,timeunit
                  call flush(run_print_unit)
               endif
            endif
            deallocate(amnp0,pmnp0)

!            allocate(sexp(16,0:2*l),scexp(16,0:2*l))
!            call ranorientscatmatrix(tmatrixfile,sexp, &
!               scexp,override_order=l, &
!               keep_quiet=.true., &
!               mpi_comm=mpicomm)
!            if(rank.eq.0) then
!               if(l.eq.1) then
!                  open(20,file='tmatsmexp.dat')
!               else
!                  open(20,file='tmatsmexp.dat',access='append')
!               endif
!               write(20,'(i5)') l
!               do n=0,2*l
!                  write(20,'(i5,2es14.6)') n,sexp(1,n),scexp(1,n)
!               enddo
!               close(20)
!            endif
!            deallocate(sexp,scexp)

            call mstm_mpi(mpi_command='bcast',mpi_rank=0,mpi_send_buf_dp=converr,mpi_number=1,mpi_comm=mpicomm)

            if(conveps.gt.0.d0.and.converr(1).lt.conveps) exit
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
                    mpi_comm,excited_spheres,solution_method,initialize_solver)
         implicit none
         logical :: firstrun,exsphere(number_spheres),dirsoln,initialize
         logical, save :: inp1,inp2
         logical, optional :: excited_spheres(number_spheres),initialize_solver
         integer :: iter,niter,istat,rank,maxiter,iterwrite,nsend,&
                    numprocs,mpicomm,prank,oddnumproc, &
                    groupsize,pgroup,mpigroup,syncgroup,i,p,dir,qeffdim
         integer, save :: pcomm,synccomm1,synccomm2,p1,p2
         integer, allocatable :: grouplist(:)
         integer, optional :: mpi_comm
         real(8) :: alpha,sinc,eps,serr,qeff(3,qeffdim,number_spheres),maxerr
         complex(8) :: amnp(number_eqns,2)
         complex(8), allocatable :: pmnpan(:),pmnp0(:,:)
         character*1, optional :: solution_method
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
         if(present(solution_method)) then
            dirsoln=(solution_method.ne.'i')
         else
            dirsoln=.false.
         endif
         if(present(initialize_solver)) then
            initialize=initialize_solver
         else
            initialize=firstrun
         endif

!         if(mpicomm.eq.mpi_comm_null) return
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)

         if(initialize) then
            if(numprocs.gt.1.and.(.not.dirsoln)) then
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
         firstrun=.false.

!write(*,*) ' step 3'
!call flush(6)
if(light_up) then
write(*,'('' s8.2.1 '',i3)') mstm_global_rank
call flush(6)
endif
!call mstm_mpi(mpi_command='barrier')

         allocate(pmnpan(number_eqns),pmnp0(number_eqns,2))
         call sphereplanewavecoef(alpha,sinc,dir,pmnp0,excited_spheres=exsphere,mpi_comm=mpicomm)
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
               if(dirsoln) then
                  call direct_solver(pmnpan,amnp(:,p),initialize_solver=initialize, &
                     number_iterations=0,solution_error=serr,mpi_comm=mpicomm)
                  if(numprocs.gt.1) then
                     call mstm_mpi(mpi_command='bcast', &
                        mpi_send_buf_dc=amnp(1:nsend,p), &
                        mpi_number=number_eqns, &
                        mpi_rank=0, &
                        mpi_comm=mpicomm)
                  endif
                  iter=0
               else
                  call cbicg(niter,eps,pmnpan,amnp(:,p),iterwrite, &
                         iter,serr,mpi_comm=pcomm,initialize_solver=initialize)
               endif
            else
               iter=0
               serr=0.d0
               if(number_host_spheres.gt.0) then
                  pmnpan=0.d0
                  call sphereinteraction(number_eqns,1,amnp(:,p),pmnpan, &
                       initial_run=initialize, &
                       skip_external_translation=.true., &
                       mpi_comm=pcomm)
                  call sphereinteraction(number_eqns,1,pmnpan,pmnpan, &
                       skip_external_translation=.true., &
                       mpi_comm=pcomm)
                  amnp(:,p)=amnp(:,p)+pmnpan
               endif
            endif
            initialize=.false.
            maxiter=max(iter,maxiter)
            maxerr=max(serr,maxerr)
            if(iter.gt.niter.or.serr.gt.eps) istat=1
         enddo

!         call mstm_mpi(mpi_command='barrier')

         if(numprocs.gt.1.and.(.not.dirsoln)) then
            nsend=number_eqns
            if(inp1) then
               call mstm_mpi(mpi_command='bcast', &
                  mpi_send_buf_dc=amnp(1:nsend,1), &
                  mpi_number=nsend, &
                  mpi_rank=0, &
                  mpi_comm=synccomm1)
            endif
!            call mstm_mpi(mpi_command='barrier')
            if(inp2) then
               call mstm_mpi(mpi_command='bcast', &
                  mpi_send_buf_dc=amnp(1:nsend,2), &
                  mpi_number=nsend, &
                  mpi_rank=0, &
                  mpi_comm=synccomm2)
            endif
         endif
!         call mstm_mpi(mpi_command='barrier')
!
!  efficiency factor calculations
!
if(light_up) then
write(*,'('' s8.2.4 '',i3)') mstm_global_rank
call flush(6)
endif
!
! what is this doing here?
!call mstm_mpi(mpi_command='barrier')
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

         subroutine direct_solver(pnp,anp,initialize_solver,solution_error, &
            number_iterations,solution_eps,mpi_comm)
         implicit none
         logical :: initialize
         logical, save :: firstrun
         logical, optional :: initialize_solver
         integer :: ierr,rank,niter,iter,mpicomm
         integer, allocatable, save :: indx(:)
         integer, optional :: number_iterations,mpi_comm
         real(8) :: dsign,serr,seps
         real(8), optional :: solution_error,solution_eps
         complex(8)  :: pnp(number_eqns),anp(number_eqns),rnp(number_eqns)
         complex(8), allocatable, save :: amat(:,:),lumat(:,:)
         data firstrun/.true./
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         if(present(initialize_solver)) then
            initialize=initialize_solver
         else
            initialize=firstrun
         endif
         if(present(number_iterations)) then
            niter=number_iterations
         else
            niter=3
         endif
         if(present(solution_eps)) then
            seps=solution_eps
         else
            seps=1.d-12
         endif
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)

         if(initialize) then
            if(allocated(amat)) deallocate(amat,lumat,indx)
            allocate(amat(number_eqns,number_eqns),lumat(number_eqns,number_eqns),indx(number_eqns))
            pl_error_codes=0
            call general_interaction_matrix(amat,mie_mult=.true.,mpi_comm=mpicomm)
            call mstm_mpi(mpi_command='reduce',mpi_rank=0,mpi_operation=mstm_mpi_sum, &
               mpi_recv_buf_i=pl_error_codes,mpi_number=6,mpi_comm=mpicomm)
            if(rank.eq.0) then
               if(any(pl_error_codes.ne.0)) then
                  write(run_print_unit,'('' LU decomposion pl error codes:'',10i4)') pl_error_codes,pl_rs_imax
               endif
               lumat=amat
               call lu_decomposition(lumat,number_eqns,indx,dsign,ierr)
               if(ierr.ne.0) then
                  if(rank.eq.0) then
                     write(run_print_unit,'('' lu decomposition failed!!!'')')
                  endif
                  stop
               endif
            endif
            firstrun=.false.
         endif
         if(rank.eq.0) then
            anp=pnp
            call lu_backsubstitution(lumat,number_eqns,indx,anp)
            do iter=1,niter
               rnp=matmul(amat,anp)-pnp
               serr=sqrt(sum(cdabs(rnp*dconjg(rnp))))/dble(number_eqns)
               call lu_backsubstitution(lumat,number_eqns,indx,rnp)
               anp=anp-rnp
            enddo
            if(present(solution_error)) then
               solution_error=serr
            endif
         endif
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
!         call mstm_mpi(mpi_command='rank',mpi_rank=rank0)
         rank0=mstm_global_rank

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

         if(niter.eq.0) then
            cp=anp
            cr=0.
            call sphereinteraction(neqns,1,cp,cr, &
                 initial_run=initialize, &
                 skip_external_translation=.true., &
                 mpi_comm=mpicomm)
            cp=anp
            call sphereinteraction(neqns,1,cp,anp, &
                 skip_external_translation=.true., &
                 mpi_comm=mpicomm)
            deallocate(cr,cp,cw,cq,cap,caw,capt,cawt)
            return
         endif
!
!  setting niter < 0 runs the following simple order--of--scattering solution
!
         if(niter.lt.0) then
            cp=anp
            if(rank0.eq.0) time1=mstm_mpi_wtime()
            do iter=1,-niter
               if(rank0.eq.0) time0=mstm_mpi_wtime()
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
               if(rank0.eq.0) time2=mstm_mpi_wtime()
               if(rank0.eq.0.and.iterwrite.eq.1.and.time2-time1.gt.5.d0) then
                  write(iunit,'('' iter,err,tpi:'',i5,2e13.5)') iter,eerr,time2-time0
                  call flush(iunit)
                  time1=time2
               endif
!if(rank0.eq.0) then
!write(*,'(i5,6es20.12)') iter,sum(cdabs(anp)),anp(number_eqns)
!endif
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
            if(rank0.eq.0) time0=mstm_mpi_wtime()
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
!if(rank0.eq.0) then
!write(*,'(i5,6es20.12)') iter,sum(cdabs(anp)),anp(number_eqns),cak
!endif
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
            if(rank0.eq.0) time2=mstm_mpi_wtime()
            if(rank0.eq.0.and.iterwrite.eq.1.and.time2-time1.gt.5.d0) then
               write(run_print_unit,'('' iter,err,min err, tpi:'',i5,2e12.4,e12.4)') &
                     iter,errmax,errmin,time2-time0
               call flush(iunit)
               time1=time2
            endif
         enddo
         deallocate(cr,cp,cw,cq,cap,caw,capt,cawt)
         end subroutine cbicg

         subroutine lu_decomposition(a,n,indx,d,ierr)
         implicit none
         integer :: n,indx(n),i,j,k,imax,ierr
         real(8) :: tiny,aamax,d,vv(n),dum,time1,time2
         complex(8) ::  a(n,n),sum,cdum
         data tiny/1.d-20/
         ierr=0
         d=1.d0
         time1=mstm_mpi_wtime()
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
            time2=mstm_mpi_wtime()
            if(mstm_global_rank.eq.0.and.time2-time1.gt.15.d0) then
               write(run_print_unit,'('' lu decomposition, step '', i5,''/'',i5)') j,n
               call flush(run_print_unit)
               time1=time2
            endif
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

      end module solver
