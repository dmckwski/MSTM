      module random_sphere_configuration
      use mpidefs
      use specialfuncs
      implicit none
      type coll_list
         logical :: wallcoll
         integer :: wall,sphere
         real(8) :: time,collpos(3)
      end type coll_list
      type l_list
         integer :: index
         type(l_list), pointer :: next
      end type l_list
      type c_list
         integer :: number_elements
         type(l_list), pointer :: members
      end type c_list
      logical  :: target_width_specified
      logical, target :: sphere_1_fixed,periodic_bc(3),random_lattice_configuration
      integer, private :: cell_dim(3)
      integer, allocatable :: sphere_cell(:,:)
      integer, target :: target_shape,wall_boundary_model,max_number_time_steps,number_components
      real(8), private :: pi,fv_crit,time_step
      real(8), private :: minimum_gap,d_cell,target_boundaries(3,2)
      real(8), target :: target_dimensions(3),psd_sigma(4),target_width,target_thickness,max_collisions_per_sphere, &
                         max_diffusion_cpu_time,max_diffusion_simulation_time,component_radii(4), &
                         component_number_fraction(4)
      real(8) :: sim_timings(10),time_0
      type(c_list), allocatable :: cell_list(:,:,:)
      type(coll_list), allocatable :: coll_data(:)
      character*1 :: c_temp
      data pi,fv_crit,time_step/3.1415926535897932385d0,0.25d0,.1d0/
      data minimum_gap,sphere_1_fixed,target_shape/1.0d-3,.false.,0/
      data periodic_bc/.true.,.true.,.true./
      data wall_boundary_model/1/
      data max_number_time_steps,max_collisions_per_sphere,max_diffusion_simulation_time/100,3.d0,5.d0/
      data max_diffusion_cpu_time/100.d0/
      data number_components,component_radii,component_number_fraction,psd_sigma/1,4*1.d0, &
         1.d0,0.d0,0.d0,0.d0,4*0.d0/

      contains

         subroutine random_cluster_of_spheres(numberspheres,targetdimensions,sphereposition,sphereradius, &
            sphereindex,iunit,istatus, &
            ntsteps,skip_diffusion,use_saved_values,print_progress,simulation_file)
         implicit none
         logical :: fitok,allin,initial0,initial1,trystage1,skipdif,pprog,printsim,multicomp
         logical, optional :: skip_diffusion,use_saved_values,print_progress
         logical, save :: firstrun
         integer :: i,j,maxsamp0,maxsamp1,numberspheres,ncolls,ncollstot,maxns,ntsteps,istatus,iunit,rank,opair(2), &
            nscompi(4),sphereindex(numberspheres),n
         real(8) ::samppos(3),sphereposition(3,numberspheres),sphereradius(numberspheres), &
            spherevol,targetfv,u(3,numberspheres),wallboundaries(3,2),targetvol,targetdimensions(3), &
            targetstretch,collspersphere,mfp,time0,time1,sum1,sum2,sdev,mean,time2,rmin,rnum(3), &
            radscale
         real(8), allocatable, save :: saved_sphereradius(:),saved_sphereposition(:,:)
         character*255, optional :: simulation_file
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
         if(present(print_progress)) then
            pprog=print_progress
         else
            pprog=.false.
         endif
         printsim=(present(simulation_file).and.mstm_global_rank.eq.0)
         if(firstrun) then
            call random_seed()
            firstrun=.false.
         endif
         allocate(coll_data(numberspheres))
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         trystage1=.false.

         nscompi(1:number_components)=numberspheres*component_number_fraction(1:number_components)
         if(number_components.gt.1) nscompi(number_components)=nscompi(number_components) &
            +numberspheres-sum(nscompi(1:number_components))
         radscale=(sum(component_number_fraction(1:number_components) &
            *component_radii(1:number_components)**3))**(1.d0/3.d0)

         spherevol=0.d0
         n=0
         do j=1,number_components
            do i=1,nscompi(j)
               n=n+1
               sphereindex(n)=j
               call psdsamp(psd_sigma(j),2.5d0,sphereradius(n))
               sphereradius(n)=sphereradius(n)*component_radii(j)/radscale
               spherevol=spherevol+4.d0*pi/3.d0*sphereradius(n)**3
            enddo
   !            if(psd_sigma.gt.0.1d0) then
   !               call sort_radii(numberspheres,sphereradius)
   !               trystage1=.true.
   !            endif
         enddo
         call target_volume(targetdimensions,targetvol)
         targetfv=spherevol/targetvol
         targetstretch=(1.d0/targetfv)**(1.d0/3.d0)
         targetstretch=max(targetstretch,1.02d0)
         mfp=targetvol/dble(numberspheres)/4.d0
         target_boundaries(:,1)=-targetdimensions
         target_boundaries(:,2)=targetdimensions
         wallboundaries=target_boundaries
         d_cell=2.5d0*maxval(sphereradius(1:numberspheres))
         allin=.false.
         istatus=3
         sim_timings=0.d0
         if((targetfv.le.0.25d0.or.trystage1).and.(.not.allin).and.(.not.random_lattice_configuration)) then
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
               if(mstm_global_rank.eq.0.and.pprog) then
                  write(iunit,'('' target configuration computed using random sampling'')')
                  call flush(iunit)
               endif
               istatus=0
            else
               call clear_cells()
            endif
         endif
         if(targetfv.lt.0.6d0.and.(.not.allin).and.(.not.random_lattice_configuration)) then
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
               if(mstm_global_rank.eq.0.and.pprog) then
                  write(iunit,'('' target configuration computed using layered sampling + diffusion, time steps:'',i5)') ntsteps
                  call flush(iunit)
               endif
               istatus=1
            endif
         endif
         if((.not.allin).or.random_lattice_configuration) then
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
            if(mstm_global_rank.eq.0.and.pprog) then
               write(iunit,'('' target configuration computed initial HCP + diffusion, time steps:'',i5)') ntsteps
               call flush(iunit)
            endif
         endif
         do i=1,numberspheres
            call check_in_target(sphereradius(i),sphereposition(:,i),wallboundaries,allin)
            if(.not.allin) write(iunit,'('' initially outside:'',i5,3es12.4)') i,sphereposition(:,i)
         enddo
!call direct_overlap_test(numberspheres,sphereradius,sphereposition,allin,distance=rmin,pair=opair)
!if(allin) write(iunit,'('' initial overlap:'',2i4,f8.3)') opair,rmin
         ntsteps=max_number_time_steps
         if(printsim) then
            open(31,file=trim(simulation_file))
            write(31,'(i8)') numberspheres
            write(31,'(es13.5)') 0.d0
            do i=1,numberspheres
               write(31,'(4es13.5)') sphereposition(:,i),sphereradius(i)
            enddo
         endif
         if((.not.skipdif).and.max_number_time_steps.gt.0.and.max_diffusion_simulation_time.gt.0.) then
            call samptrajectory(numberspheres,u)
            ncollstot=0
            time1=mstm_mpi_wtime()
            time0=time1
            do j=1,max_number_time_steps
               call spheremove(numberspheres,sphereradius,sphereposition,u,time_step,wallboundaries, &
                  number_wall_hits=ncolls)
               ncollstot=ncollstot+ncolls
               collspersphere=dble(ncollstot)/dble(numberspheres)
               time2=mstm_mpi_wtime()
               if(mstm_global_rank.eq.0.and.pprog.and.time2-time1.gt.15.d0) then
                  write(iunit,'('' diffusion step, collision per sphere:'',i8,es12.4)') j,collspersphere
                  call flush(iunit)
                  time1=mstm_mpi_wtime()
               endif
               if(time2-time0.gt.max_diffusion_cpu_time) exit
               if(j*time_step.gt.max_diffusion_simulation_time) exit
               if(collspersphere.gt.max_collisions_per_sphere) exit
               if(printsim) then
                  write(31,'(es13.5)') j*time_step
                  do i=1,numberspheres
                     write(31,'(4es13.5)') sphereposition(:,i),sphereradius(i)
                  enddo
               endif
            enddo
            ntsteps=min(ntsteps,j)
            if(printsim) close(31)
            do i=1,numberspheres
               call check_in_target(sphereradius(i),sphereposition(:,i),wallboundaries,allin)
               if(.not.allin) write(iunit,'('' outside:'',i5,3es12.4)') i,sphereposition(:,i)
            enddo
!call direct_overlap_test(numberspheres,sphereradius,sphereposition,allin,distance=rmin,pair=opair)
!if(allin) write(iunit,'('' overlap:'',2i4,f8.3)') opair,rmin
         endif
         if(target_shape.eq.0.or.target_shape.eq.1) then
            call sort_positions(numberspheres,sphereradius,sphereposition,sphereindex,3)
         else
            call sort_positions(numberspheres,sphereradius,sphereposition,sphereindex,0)
         endif
         if(random_lattice_configuration) then
            call random_number(rnum(1:3))
            rnum(1)=2.d0*pi*rnum(1)
            rnum(3)=2.d0*pi*rnum(3)
            rnum(2)=dacos(-1.d0+2.d0*rnum(2))
            call eulerrotation(sphereposition(:,1:numberspheres),rnum,1, &
               sphereposition(:,1:numberspheres),numberspheres)
         endif

         sphereposition(:,1:numberspheres)=sphereposition(:,1:numberspheres)*radscale
         sphereradius(1:numberspheres)=sphereradius(1:numberspheres)*radscale

         if(allocated(saved_sphereposition)) then
            deallocate(saved_sphereposition,saved_sphereradius)
         endif
         allocate(saved_sphereposition(3,numberspheres),saved_sphereradius(numberspheres))
         saved_sphereposition=sphereposition
         saved_sphereradius=sphereradius
         call clear_cells()
         deallocate(coll_data)
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

         subroutine target_volume(targetdimensions,targetvol)
         implicit none
         integer :: i,ipbc(3)
         real(8) :: targetvol,targetdimensions(3)
         ipbc=wall_boundary_model
         do i=1,3
            if(periodic_bc(i)) ipbc(i)=0
         enddo
         if(target_shape.eq.0) then
            targetvol=8.d0*product(targetdimensions(1:3)-dble(ipbc))
         elseif(target_shape.eq.1) then
            targetvol=2.d0*pi*((targetdimensions(1)-wall_boundary_model)**2)*(targetdimensions(3)-ipbc(3))
         else
            targetvol=4.d0*pi*(targetdimensions(1)-wall_boundary_model)**3/3.d0
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

         subroutine sort_positions(nsphere,radius,position,cindex,sort_elem,make_positive)
         implicit none
         logical :: makepos
         logical, optional :: make_positive
         integer :: nsphere,i,ind(nsphere),selem,cindex(nsphere),tindex(nsphere)
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
         tindex=cindex
         do i=1,nsphere
            radius(i)=r(ind(i))
            position(:,i)=tpos(:,ind(i))
            cindex(i)=tindex(ind(i))
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
         integer :: nsphere,i,is,js,collisionpair(2),iwall,iswall,m,nwhits,it
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
time_0=mstm_mpi_wtime()
         call trajectorytest(nsphere,radius,pos,u,tmove,wallboundaries, &
            collision,tcoll,collisionpair)
sim_timings(1)=sim_timings(1)+mstm_mpi_wtime()-time_0
         do while(tmove.gt.0.d0)
time_0=mstm_mpi_wtime()
            call modify_cells(nsphere,pos)
sim_timings(2)=sim_timings(2)+mstm_mpi_wtime()-time_0
time_0=mstm_mpi_wtime()
            tcoll=tmove
            collision=.false.
            do is=1,nsphere
               if(coll_data(is)%time.lt.tcoll) then
                  tcoll=coll_data(is)%time
                  collision=.true.
                  collisionpair(1)=is
                  collisionpair(2)=coll_data(is)%sphere
                  collpos(:)=coll_data(is)%collpos(:)
               endif
            enddo
sim_timings(3)=sim_timings(3)+mstm_mpi_wtime()-time_0
time_0=mstm_mpi_wtime()
            tcmin=tcoll
            call walltest(nsphere,radius,pos,u,tmove,wallboundaries,twallmin,iswall,iwall)
sim_timings(4)=sim_timings(4)+mstm_mpi_wtime()-time_0
time_0=mstm_mpi_wtime()
            wallcollision=(twallmin.lt.tcmin)
            tcmin=min(tcmin,twallmin)
            do is=1,nsphere
               if(is.eq.1.and.sphere_1_fixed) cycle
               intarget=.false.
               do while(.not.intarget)
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
                  if(.not.intarget) tcmin=.95*tcmin
               enddo
               if(intarget) then
                  pos(:,is)=tpos
               else
                  write(*,'('' out of target'')')
                  write(*,'(8es12.4)') tpos,sqrt(sum(tpos**2)),tcmin,tcoll,twallmin
               endif
            enddo
sim_timings(5)=sim_timings(5)+mstm_mpi_wtime()-time_0
time_0=mstm_mpi_wtime()
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
            coll_data(1:nsphere)%time=coll_data(1:nsphere)%time-abs(tcmin)
            if(wallcollision) then
               call trajectorytest(nsphere,radius,pos,u,tmove,wallboundaries, &
               collision,tcoll,collisionpair,start_sphere=iswall,end_sphere=iswall)
            elseif(collision) then
               is=collisionpair(1)
               js=collisionpair(2)
               call trajectorytest(nsphere,radius,pos,u,tmove,wallboundaries, &
               collision,tcoll,collisionpair,start_sphere=is,end_sphere=is)
               call trajectorytest(nsphere,radius,pos,u,tmove,wallboundaries, &
               collision,tcoll,collisionpair,start_sphere=js,end_sphere=js)
            endif
            i=i+1
            if(sphere_1_fixed) u(:,1)=0.d0
sim_timings(6)=sim_timings(6)+mstm_mpi_wtime()-time_0
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
            coll_data(is)%wallcoll=.false.
            coll_data(is)%time=maxtime
            coll_data(is)%collpos(:)=pos(:,is)
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
                        if(coll_data(is)%time.gt.tcollision) then
                           coll_data(is)%time=tcollision
                           coll_data(is)%sphere=js
                           coll_data(is)%collpos(:)=tpos(:)
                        endif
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
