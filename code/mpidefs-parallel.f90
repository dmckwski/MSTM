      module mpidefs
      use mpi
      implicit none

!      include 'mpif.h'
      integer :: mstm_mpi_comm_world,mstm_mpi_sum,mstm_mpi_max,mstm_mpi_min,mstm_global_rank

      contains

         real(8) function mstm_mpi_wtime()
         implicit none
         mstm_mpi_wtime=mpi_wtime()
         end function mstm_mpi_wtime


         subroutine mstm_mpi(mpi_command,mpi_recv_buf_i,mpi_recv_buf_r,mpi_recv_buf_c,mpi_recv_buf_dp, &
                           mpi_recv_buf_dc,mpi_send_buf_i,mpi_send_buf_r,mpi_send_buf_c, &
                           mpi_send_buf_dp,mpi_send_buf_dc,mpi_number,mpi_comm,mpi_group,mpi_rank,mpi_size,&
                           mpi_new_comm,mpi_new_group,mpi_new_group_list,mpi_operation, &
                           mpi_color,mpi_key,mpi_tag,mpi_flag,mpi_recv_buf_char,mpi_send_buf_char)
         integer, optional :: mpi_number,mpi_recv_buf_i(*),mpi_send_buf_i(*),mpi_comm,mpi_group,mpi_rank, &
                              mpi_size,mpi_new_comm,mpi_new_group,mpi_new_group_list(*),mpi_operation, &
                              mpi_color,mpi_key,mpi_tag
         integer :: stat(MPI_STATUS_SIZE),mpitag,mpisource
         logical, optional :: mpi_flag
         real(4), optional :: mpi_recv_buf_r(*),mpi_send_buf_r(*)
         real(8), optional :: mpi_recv_buf_dp(*),mpi_send_buf_dp(*)
         real(8), allocatable :: dptemp(:)
         complex(4), optional :: mpi_recv_buf_c(*),mpi_send_buf_c(*)
         complex(4), allocatable :: ctemp(:)
         complex(8), optional :: mpi_recv_buf_dc(*),mpi_send_buf_dc(*)
         complex(8), allocatable :: dctemp(:)
         character, optional :: mpi_recv_buf_char(*),mpi_send_buf_char(*)
         character(*) :: mpi_command
         integer :: ntype,ierr,comm,size,rank,group,newcomm

         if(mpi_command.eq.'init') then
            call mpi_init(ierr)
            mstm_mpi_comm_world=mpi_comm_world
            mstm_mpi_sum=mpi_sum
            mstm_mpi_max=mpi_max
            mstm_mpi_min=mpi_min
            call mpi_comm_rank(mpi_comm_world,mstm_global_rank,ierr)
            return
         endif
         if(mpi_command.eq.'finalize') then
            call mpi_finalize(ierr)
            return
         endif
         if(present(mpi_comm)) then
            comm=mpi_comm
         else
            comm=mpi_comm_world
         endif
         if(mpi_command.eq.'iprobe') then
            if(present(mpi_tag)) then
               mpitag=mpi_tag
            else
               mpitag=mpi_any_tag
            endif
            if(present(mpi_rank)) then
               mpisource=mpi_rank
            else
               mpisource=mpi_any_source
            endif
            call mpi_iprobe(mpisource,mpitag,comm,mpi_flag,stat,ierr)
            return
         endif
         if(present(mpi_tag)) then
            mpitag=mpi_tag
         else
            mpitag=1
         endif
         if(mpi_command.eq.'size') then
            call mpi_comm_size(comm,size,ierr)
            mpi_size=size
            return
         endif
         if(mpi_command.eq.'rank') then
            call mpi_comm_rank(comm,rank,ierr)
            mpi_rank=rank
            return
         endif
         if(mpi_command.eq.'group') then
            call mpi_comm_group(comm,group,ierr)
            mpi_group=group
            return
         endif
         if(mpi_command.eq.'incl') then
            call mpi_group_incl(mpi_group,mpi_size,mpi_new_group_list,group,ierr)
            mpi_new_group=group
            return
         endif
         if(mpi_command.eq.'create') then
            call mpi_comm_create(comm,mpi_group,newcomm,ierr)
            mpi_new_comm=newcomm
            return
         endif
         if(mpi_command.eq.'split') then
            call mpi_comm_split(comm,mpi_color,mpi_key,newcomm,ierr)
            mpi_new_comm=newcomm
            return
         endif
         if(mpi_command.eq.'barrier') then
            call mpi_barrier (comm,ierr)
            return
         endif
         if(mpi_command.eq.'free') then
            call mpi_comm_free (comm,ierr)
            return
         endif

!         if(present(mpi_recv_buf_char).or.present(mpi_send_buf_char)) then
!            ntype=mpi_character
!            if(mpi_command.eq.'bcast') then
!               call MPI_BCAST (mpi_send_buf_char,mpi_number,ntype,mpi_rank,comm,ierr)
!               return
!            endif
!         endif

         if(present(mpi_recv_buf_i).or.present(mpi_send_buf_i)) then
            ntype=mpi_integer
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_i,mpi_number,ntype,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_i,mpi_number,ntype,mpi_rank,mpitag,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_i,mpi_number,ntype,mpi_rank,mpitag,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_i)) then
                  call mpi_reduce(mpi_send_buf_i,mpi_recv_buf_i,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_reduce(mpi_in_place,mpi_recv_buf_i,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_i)) then
                  call mpi_allreduce(mpi_send_buf_i,mpi_recv_buf_i,mpi_number,ntype,mpi_operation, &
                               comm,ierr)
               else
                  call mpi_allreduce(mpi_in_place,mpi_recv_buf_i,mpi_number,ntype,mpi_operation, &
                               comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'gather') then
               if(present(mpi_send_buf_i)) then
                  call mpi_gather(mpi_send_buf_i,mpi_number,ntype,mpi_recv_buf_i,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_gather(mpi_in_place,mpi_number,ntype,mpi_recv_buf_i,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
         endif

         if(present(mpi_recv_buf_r).or.present(mpi_send_buf_r)) then
            ntype=mpi_real
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_r,mpi_number,ntype,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_r,mpi_number,ntype,mpi_rank,mpitag,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_r,mpi_number,ntype,mpi_rank,mpitag,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_r)) then
                  call mpi_reduce(mpi_send_buf_r,mpi_recv_buf_r,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_reduce(mpi_in_place,mpi_recv_buf_r,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_r)) then
                  call mpi_allreduce(mpi_send_buf_r,mpi_recv_buf_r,mpi_number,ntype,mpi_operation, &
                               comm,ierr)
               else
                  call mpi_allreduce(mpi_in_place,mpi_recv_buf_r,mpi_number,ntype,mpi_operation, &
                               comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'gather') then
               if(present(mpi_send_buf_r)) then
                  call mpi_gather(mpi_send_buf_r,mpi_number,ntype,mpi_recv_buf_r,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_gather(mpi_in_place,mpi_number,ntype,mpi_recv_buf_r,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
         endif

         if(present(mpi_recv_buf_c).or.present(mpi_send_buf_c)) then
            ntype=mpi_complex
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_c,mpi_number,ntype,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_c,mpi_number,ntype,mpi_rank,mpitag,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_c,mpi_number,ntype,mpi_rank,mpitag,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_c)) then
                  call mpi_reduce(mpi_send_buf_c,mpi_recv_buf_c,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  allocate(ctemp(mpi_number))
                  ctemp(1:mpi_number)=mpi_recv_buf_c(1:mpi_number)
                  call mpi_reduce(ctemp,mpi_recv_buf_c,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
                  deallocate(ctemp)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_c)) then
                  call mpi_allreduce(mpi_send_buf_c,mpi_recv_buf_c,mpi_number,ntype,mpi_operation, &
                               comm,ierr)
               else
                  call mpi_allreduce(mpi_in_place,mpi_recv_buf_c,mpi_number,ntype,mpi_operation, &
                               comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'gather') then
               if(present(mpi_send_buf_c)) then
                  call mpi_gather(mpi_send_buf_c,mpi_number,ntype,mpi_recv_buf_c,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_gather(mpi_in_place,mpi_number,ntype,mpi_recv_buf_c,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif

         endif

         if(present(mpi_recv_buf_dp).or.present(mpi_send_buf_dp)) then
            ntype=mpi_double_precision
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_dp,mpi_number,ntype,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_dp,mpi_number,ntype,mpi_rank,mpitag,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_dp,mpi_number,ntype,mpi_rank,mpitag,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_dp)) then
                  call mpi_reduce(mpi_send_buf_dp,mpi_recv_buf_dp,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  allocate(dptemp(mpi_number))
                  dptemp(1:mpi_number)=mpi_recv_buf_dp(1:mpi_number)
                  call mpi_reduce(dptemp,mpi_recv_buf_dp,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
                  deallocate(dptemp)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_dp).and.present(mpi_recv_buf_dp)) then
                  call mpi_reduce(mpi_send_buf_dp,mpi_recv_buf_dp,mpi_number,ntype,mpi_operation, &
                               0,comm,ierr)
                  call MPI_BCAST (mpi_recv_buf_dp,mpi_number,ntype,0,comm,ierr)
!                  call mpi_allreduce(mpi_send_buf_dc,mpi_recv_buf_dc,mpi_number,ntype,mpi_operation, &
!                               comm,ierr)
               else
                  allocate(dptemp(mpi_number))
                  dptemp(1:mpi_number)=mpi_recv_buf_dp(1:mpi_number)
                  call mpi_reduce(dptemp,mpi_recv_buf_dp,mpi_number,ntype,mpi_operation, &
                               0,comm,ierr)
                  call MPI_BCAST (mpi_recv_buf_dp,mpi_number,ntype,0,comm,ierr)
!                  call mpi_allreduce(dctemp,mpi_recv_buf_dc,mpi_number,ntype,mpi_operation, &
!                               mpi_rank,comm,ierr)
                  deallocate(dptemp)
               endif
            endif
            if(mpi_command.eq.'gather') then
               if(present(mpi_send_buf_dp)) then
                  call mpi_gather(mpi_send_buf_dp,mpi_number,ntype,mpi_recv_buf_dp,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_gather(mpi_in_place,mpi_number,ntype,mpi_recv_buf_dp,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
         endif

         if(present(mpi_recv_buf_dc).or.present(mpi_send_buf_dc)) then
            ntype=mpi_double_complex
            if(mpi_command.eq.'bcast') then
               call MPI_BCAST (mpi_send_buf_dc,mpi_number,ntype,mpi_rank,comm,ierr)
               return
            endif
            if(mpi_command.eq.'send') then
               call mpi_send(mpi_send_buf_dc,mpi_number,ntype,mpi_rank,mpitag,comm,ierr)
               return
            endif
            if(mpi_command.eq.'recv') then
               call mpi_recv(mpi_recv_buf_dc,mpi_number,ntype,mpi_rank,mpitag,comm,stat,ierr)
               return
            endif
            if(mpi_command.eq.'reduce') then
               if(present(mpi_send_buf_dc)) then
                  call mpi_reduce(mpi_send_buf_dc,mpi_recv_buf_dc,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
               else
                  allocate(dctemp(mpi_number))
                  dctemp(1:mpi_number)=mpi_recv_buf_dc(1:mpi_number)
                  call mpi_reduce(dctemp,mpi_recv_buf_dc,mpi_number,ntype,mpi_operation, &
                               mpi_rank,comm,ierr)
                  deallocate(dctemp)
               endif
               return
            endif
            if(mpi_command.eq.'allreduce') then
               if(present(mpi_send_buf_dc).and.present(mpi_recv_buf_dc)) then
                  call mpi_reduce(mpi_send_buf_dc,mpi_recv_buf_dc,mpi_number,ntype,mpi_operation, &
                               0,comm,ierr)
                  call MPI_BCAST (mpi_recv_buf_dc,mpi_number,ntype,0,comm,ierr)
!                  call mpi_allreduce(mpi_send_buf_dc,mpi_recv_buf_dc,mpi_number,ntype,mpi_operation, &
!                               comm,ierr)
               else
                  allocate(dctemp(mpi_number))
                  dctemp(1:mpi_number)=mpi_recv_buf_dc(1:mpi_number)
                  call mpi_reduce(dctemp,mpi_recv_buf_dc,mpi_number,ntype,mpi_operation, &
                               0,comm,ierr)
                  call MPI_BCAST (mpi_recv_buf_dc,mpi_number,ntype,0,comm,ierr)
!                  call mpi_allreduce(dctemp,mpi_recv_buf_dc,mpi_number,ntype,mpi_operation, &
!                               mpi_rank,comm,ierr)
                  deallocate(dctemp)
               endif
               return
            endif
            if(mpi_command.eq.'gather') then
               if(present(mpi_send_buf_dc).and.present(mpi_recv_buf_dc)) then
                  call mpi_gather(mpi_send_buf_dc,mpi_number,ntype,mpi_recv_buf_dc,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_gather(mpi_recv_buf_dc,mpi_number,ntype,mpi_recv_buf_dc,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'allgather') then
               if(present(mpi_send_buf_dc).and.present(mpi_recv_buf_dc)) then
                  call mpi_allgather(mpi_send_buf_dc,mpi_number,ntype,mpi_recv_buf_dc,mpi_number,ntype, &
                               comm,ierr)
               else
                  call mpi_allgather(mpi_in_place,mpi_number,ntype,mpi_recv_buf_dc,mpi_number,ntype, &
                               comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'alltoall') then
               if(present(mpi_send_buf_dc).and.present(mpi_recv_buf_dc)) then
                  call mpi_alltoall(mpi_send_buf_dc,mpi_number,ntype,mpi_recv_buf_dc,mpi_number,ntype, &
                               comm,ierr)
               else
                  call mpi_alltoall(mpi_in_place,mpi_number,ntype,mpi_recv_buf_dc,mpi_number,ntype, &
                               comm,ierr)
               endif
               return
            endif
            if(mpi_command.eq.'scatter') then
               if(present(mpi_send_buf_dc).and.present(mpi_recv_buf_dc)) then
                  call mpi_scatter(mpi_send_buf_dc,mpi_number,ntype,mpi_recv_buf_dc,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               else
                  call mpi_scatter(mpi_recv_buf_dc,mpi_number,ntype,mpi_recv_buf_dc,mpi_number,ntype, &
                               mpi_rank,comm,ierr)
               endif
               return
            endif
         endif

         end subroutine mstm_mpi
      end module mpidefs
