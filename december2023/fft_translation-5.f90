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
      integer, private :: neighbor_node(3,0:26),fft_local_host,fft_number_spheres
      integer, target :: node_order,neighbor_node_model
      integer, allocatable, private :: sphere_node(:,:)
!      real(8) :: d_cell
      real(8), private :: cell_origin(3),cell_boundary(3)
      real(8), private :: timedat(10)
      real(8), target :: cell_volume_fraction,d_cell
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
         complex(8), allocatable, save :: anode(:,:,:,:,:),gnode(:,:,:,:,:), &
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
timedat=0.d0
         endif
!  a test to speed up 3/23

         if(firstrun) then
            if(allocated(anode)) deallocate(anode,gnode,gout_loc,ain_t,gout_t)
            allocate(anode(cell_dim(1),cell_dim(2),cell_dim(3),node_order*(node_order+2)*2,nrhs), &
               gnode(cell_dim(1),cell_dim(2),cell_dim(3),node_order*(node_order+2)*2,nrhs), &
               gout_loc(number_eqns,nrhs),ain_t(neqns,nrhs),gout_t(neqns,nrhs))
         endif

if(light_up) then
write(*,'('' fft2 '',i3)') mstm_global_rank
call flush(6)
endif

         do rhs=1,nrhs
            noff=0
            do i=1,number_spheres
               if(host_sphere(i).ne.fft_local_host) cycle
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

!timedat(1)=mstm_mpi_wtime()
         call local_sphere_to_node_translation(nrhs,ain_t,anode, &
            store_matrix_option=smopt,initial_run=firstrun, &
            mpi_comm=mpicomm,local_host=fft_local_host,sphere_to_node=.true., &
            merge_procs=.true.)
!timedat(1)=mstm_mpi_wtime()-timedat(1)

if(light_up) then
write(*,'('' fft4 '',i3)') mstm_global_rank
call flush(6)
endif

!timedat(2)=mstm_mpi_wtime()
         do p=p1,p2
            do rhs=1,nrhs
               call fft_node_to_node_translation(anode(:,:,:,:,rhs), &
                  cell_translation_matrix(:,:,:,:,:,p), &
                  gnode(:,:,:,:,rhs),p,mpi_comm=pcomm)
            enddo
         enddo
!timedat(2)=mstm_mpi_wtime()-timedat(2)

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

!timedat(3)=mstm_mpi_wtime()

         call local_sphere_to_node_translation(nrhs,gout_t,gnode, &
                 store_matrix_option=smopt, &
                 mpi_comm=mpicomm,local_host=fft_local_host,sphere_to_node=.false., &
                 merge_procs=.false.)

!timedat(3)=mstm_mpi_wtime()-timedat(3)
!timedat(4)=mstm_mpi_wtime()

         call local_sphere_to_sphere_expansion(nrhs,ain_t,gout_loc, &
                 store_matrix_option=smopt,initial_run=firstrun, &
                 mpi_comm=mpicomm,merge_procs=.false., &
                 local_host=fft_local_host)
!timedat(4)=mstm_mpi_wtime()-timedat(4)

         gout_t=gout_t+gout_loc

if(light_up) then
write(*,'('' fft6 '',i3)') mstm_global_rank
call flush(6)
endif

         do rhs=1,nrhs
            noff=0
            do i=1,number_spheres
               if(host_sphere(i).ne.fft_local_host) cycle
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

!         deallocate(anode,gnode,gout_loc,ain_t,gout_t)
         firstrun=.false.

!if(mstm_global_rank.eq.mstm_global_numprocs-1.or.mstm_global_rank.eq.0.or..true.) then
!write(*,'(i5,4es13.5)') mstm_global_rank,timedat(1:4)
!endif

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
         gout=0.
!
!  compute offsets for scattering coefficients
!
         if(firstrun) then
            if(smopt.and.store_translation_matrix) then
               task=0
               ndim=0
               do i=1,number_spheres-1
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
         do i=1,number_spheres-1
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

         subroutine node_selection(fva,target_min,target_max,d_specified,local_host)
         implicit none
         logical :: dspec
         logical, optional :: d_specified
         integer :: nsphere,m,spherenode(3),n,i,node,ix,iy,iz,ir,j,icell,ncell,cell,lochost
         integer, optional :: local_host
         real(8) :: fva,r,fv,amean,svol,tvol,dd,targetmin(3),targetmax(3)
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
         if(present(d_specified)) then
            dspec=d_specified
         else
            dspec=.false.
         endif
         if(present(local_host)) then
            lochost=local_host
         else
            lochost=0
         endif
         fft_local_host=lochost
         cell_boundary=targetmin
         if(fft_local_host.eq.0) then
            host_ref_index=layer_ref_index(0)
         else
            host_ref_index=sphere_ref_index(:,fft_local_host)
         endif

         amean=0.
         nsphere=0
         svol=0.
         do i=1,number_spheres
            if(host_sphere(i).eq.fft_local_host) then
               amean=amean+sphere_radius(i)
               svol=svol+sphere_radius(i)**3
               nsphere=nsphere+1
            endif
         enddo
         svol=svol**(1.d0/3.d0)
         amean=amean/dble(nsphere)
         fft_number_spheres=nsphere
         tvol=1.d0
         do i=1,3
            dd=targetmax(i)-targetmin(i)
            dd=max(dd,amean)
            tvol=tvol*dd
         enddo

         if(.not.dspec) then
            if(fva.le.0.d0) then
               fv=4.d0*pi/3.d0*svol**3/tvol
               fv=min(fv,1.d0)
               fv=max(fv,.02d0)
            else
               fv=fva
            endif
            fva=fv
            d_cell=(4.d0*pi/3.d0/fv/dble(nsphere))**(1.d0/3.d0)*svol
         endif

!write(*,'('' host ri:'',2e13.5)') host_ref_index
!call flush(6)

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
            if(host_sphere(i).ne.fft_local_host) cycle
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
            if(host_sphere(i).ne.fft_local_host) cycle
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
