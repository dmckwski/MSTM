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
         square_cell,use_previous_configuration,calculate_up_down_scattering, &
         d_cell_specified,medium_ref_index_specified,medium_reim_ref_index_specified
      logical, target :: append_output_file, &
                 print_scattering_matrix, &
                 copy_input_file,calculate_near_field, &
                 move_to_front,move_to_back,random_orientation, &
                 t_matrix_centered_on_1,calculate_scattering_matrix, &
                 normalize_s11,print_sphere_data,single_origin_expansion, &
                 azimuthal_average,incident_frame,configuration_average, &
                 frozen_configuration,reflection_model,input_fft_translation_option, &
                 print_random_configuration,print_timings, &
                 input_calculate_up_down_scattering,incidence_average,auto_absorption_sample_radius, &
                 random_configuration,check_positions,random_configuration_host,fit_for_radius, &
                 numerical_azimuthal_average,numerical_hemispherical_integration,input_effective_medium_simulation, &
                 auto_target_radius,erase_sphere_1
      logical, allocatable :: sphere_excitation_switch(:)
      integer :: n_nest_loops,i_var_start(5),i_var_stop(5),i_var_step(5), &
                 run_number,loop_sphere_number(5),qeff_dim, &
                 scattering_map_directions,local_rank,scat_mat_ldim,scat_mat_udim,scat_mat_mdim, &
                 ran_config_stat,ran_config_time_steps,n_configuration_groups,random_configuration_number, &
                 solution_iterations,incident_direction_number,number_rl_dirs(2),max_number_rl_dirs,fit_stat
      integer, target :: max_iterations,t_matrix_procs_per_solution, &
         scattering_map_model,scattering_map_dimension,near_field_calculation_model, &
         incident_direction,number_configurations,min_fft_nsphere,input_node_order, &
         number_incident_directions,shifted_sphere,number_excited_spheres,input_number_spheres, &
         random_configuration_host_model
      integer, allocatable :: sphere_index(:)
      real(8) :: r_var_start(5),r_var_stop(5),r_var_step(5),diffuse_scattering_ratio, &
          coherent_scattering_ratio,hemispherical_sca(2,2),evan_sca(2),prop_sca(2), &
          input_layer_thickness(max_number_plane_boundaries), &
          pl_sca(2,2),scat_mat_amin,scat_mat_amax,pl_sca_ave(2,2),solution_time,solution_time_ave, &
          incident_beta,solution_error,surface_absorptance(2),surface_absorptance_ave(2), &
          position_shift(3),tot_csca_ave(1),dif_csca_ratio(1),fit_radius,tot_csca,dif_csca
      real(8), allocatable :: q_eff(:,:,:),q_vabs(:,:),q_eff_tot(:,:),scat_mat(:,:), &
         dif_scat_mat(:,:),sm_coef(:,:,:),sm_cf_coef(:,:,:),boundary_sca(:,:),boundary_ext(:,:), &
         q_eff_ave(:,:,:),q_vabs_ave(:,:),q_eff_tot_ave(:,:),scat_mat_ave(:,:), &
         boundary_sca_ave(:,:),boundary_ext_ave(:,:),sphere_position_ave(:,:),dif_boundary_sca(:,:), &
         scat_mat_exp_coef(:,:,:),scat_mat_exp_coef_ave(:,:,:),rl_vec(:,:),coh_scat_mat_exp_coef(:,:,:), &
         coh_scat_mat_exp_coef_ave(:,:,:),s_field(:,:,:,:,:),s_field_ave(:,:,:,:,:)
      real (8), target :: incident_beta_deg,incident_alpha_deg,solution_epsilon, &
         mie_epsilon,length_scale_factor,near_field_plane_position, &
         near_field_plane_vertices(3,2),near_field_step_size, &
         translation_epsilon,t_matrix_convergence_epsilon, &
         scattering_map_increment,incident_sin_beta,input_cell_width(2), &
         sphere_volume_fraction,input_cell_width_x, &
         input_cell_volume_fraction,medium_re_ref_index,medium_im_ref_index, &
         excitation_radius,absorption_sample_radius,absorption_sample_radius_fraction,&
         x_shift,y_shift,z_shift,input_d_cell,target_radius_padding
      complex(8) :: c_var_start(5),c_var_stop(5),c_var_step(5),nf_eff_ref_index
      complex(8), target :: ref_index_scale_factor,host_sphere_ref_index,medium_ref_index, &
         component_ref_index(4)
      complex(8), allocatable :: amnp_s(:,:),amnp_0_ave(:,:),amnp_0(:,:),e_field(:,:,:,:,:), &
         e_field_ave(:,:,:,:,:),h_field(:,:,:,:,:),h_field_ave(:,:,:,:,:), &
         mean_t(:,:),mean_t_ave(:,:)
      character*1 :: loop_var_type(5)
      character*20 :: run_date_and_time
      character*256 :: loop_var_label(5),input_file
      character*256, target :: output_file,run_file,t_matrix_output_file, &
            sphere_data_input_file,near_field_output_file,solution_method, &
            random_configuration_output_file
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
      data random_configuration_output_file/'random_configuration.pos'/
      data run_file/'run1.dat'/
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
      data numerical_azimuthal_average/.false./
      data number_spheres_specified,configuration_average/.true.,.false./
      data frozen_configuration,reflection_model,random_configuration_number/.false.,.false.,1/
      data random_configuration,random_configuration_output_file/.false.,'random_configuration.pos'/
      data min_fft_nsphere,input_fft_translation_option,input_node_order,input_cell_volume_fraction/200,.false.,-1,0.d0/
      data d_cell_specified/.false./
      data print_random_configuration,print_timings/.false.,.true./
      data input_calculate_up_down_scattering/.true./
      data incidence_average,number_incident_directions/.false.,16/
      data use_previous_configuration/.false./
      data absorption_sample_radius,excitation_radius/1.d10,1.d10/
      data auto_absorption_sample_radius,absorption_sample_radius_fraction/.true.,0.8d0/
      data x_shift,y_shift,z_shift,shifted_sphere/0.d0,0.d0,0.d0,0/
      data erase_sphere_1/.false./
      data check_positions/.true./
      data number_excited_spheres/1000000/
      data host_sphere_ref_index,random_configuration_host,random_configuration_host_model/(1.d0,0.d0),.false.,1/
      data fit_for_radius/.true./
      data medium_ref_index,medium_ref_index_specified/(1.d0,0.d0),.false./
      data medium_re_ref_index,medium_im_ref_index,medium_reim_ref_index_specified/1.d0,0.d0,.false./
      data numerical_hemispherical_integration/.false./
      data input_effective_medium_simulation/.false./
      data auto_target_radius,target_radius_padding/.false.,5.d0/
      data component_ref_index/4*(1.d0,0.d0)/

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
            d_cell_specified=.false.

         elseif(varlabel.eq.'d_cell') then
            vartype='r'
            rvarvalue=>input_d_cell
            d_cell_specified=.true.

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

         elseif(varlabel.eq.'interaction_radius') then
            vartype='r'
            rvarvalue=>interaction_radius

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

         elseif(varlabel.eq.'numerical_azimuthal_average') then
            vartype='l'
            lvarvalue=>numerical_azimuthal_average

         elseif(varlabel.eq.'incident_frame') then
            vartype='l'
            lvarvalue=>incident_frame

         elseif(varlabel.eq.'print_sphere_data') then
            vartype='l'
            lvarvalue=>print_sphere_data

         elseif(varlabel.eq.'number_spheres') then
            vartype='i'
            ivarvalue=>input_number_spheres
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

         elseif(varlabel.eq.'minimum_initial_segment_size') then
            vartype='r'
            rvarvalue=>minimum_initial_segment_size
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

         elseif(varlabel.eq.'fast_near_field') then
            vartype='l'
            lvarvalue=>fast_near_field

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

         elseif(varlabel.eq.'random_configuration') then
            vartype='l'
            lvarvalue=>random_configuration

         elseif(varlabel.eq.'random_configuration_output_file') then
            vartype='a'
            avarvalue=>random_configuration_output_file

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
!
!         elseif(varlabel.eq.'psd_sigma') then
!            vartype='r'
!            rvarvalue=>psd_sigma

         elseif(varlabel.eq.'number_components') then
            vartype='i'
            ivarvalue=>number_components

         elseif(varlabel.eq.'max_diffusion_simulation_time') then
            vartype='r'
            rvarvalue=>max_diffusion_simulation_time

         elseif(varlabel.eq.'max_diffusion_cpu_time') then
            vartype='r'
            rvarvalue=>max_diffusion_cpu_time

         elseif(varlabel.eq.'max_collisions_per_sphere') then
            vartype='r'
            rvarvalue=>max_collisions_per_sphere

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

         elseif(varlabel.eq.'wall_boundary_model') then
            vartype='i'
            ivarvalue=>wall_boundary_model

         elseif(varlabel.eq.'auto_target_radius') then
            vartype='l'
            lvarvalue=>auto_target_radius

         elseif(varlabel.eq.'target_radius_padding') then
            vartype='r'
            rvarvalue=>target_radius_padding

         elseif(varlabel.eq.'sphere_1_fixed') then
            vartype='l'
            lvarvalue=>sphere_1_fixed

         elseif(varlabel.eq.'random_lattice_configuration') then
            vartype='l'
            lvarvalue=>random_lattice_configuration

         elseif(varlabel.eq.'erase_sphere_1') then
            vartype='l'
            lvarvalue=>erase_sphere_1

         elseif(varlabel.eq.'configuration_average') then
            vartype='l'
            lvarvalue=>configuration_average

         elseif(varlabel.eq.'frozen_configuration') then
            vartype='l'
            lvarvalue=>frozen_configuration

         elseif(varlabel.eq.'random_configuration_host') then
            vartype='l'
            lvarvalue=>random_configuration_host

         elseif(varlabel.eq.'host_sphere_ref_index') then
            vartype='c'
            cvarvalue=>host_sphere_ref_index

         elseif(varlabel.eq.'random_configuration_host_model') then
            vartype='i'
            ivarvalue=>random_configuration_host_model

         elseif(varlabel.eq.'fit_for_radius') then
            vartype='l'
            lvarvalue=>fit_for_radius

         elseif(varlabel.eq.'effective_medium_simulation') then
            vartype='l'
            lvarvalue=>input_effective_medium_simulation

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

         elseif(varlabel.eq.'numerical_hemispherical_integration') then
            vartype='l'
            lvarvalue=>numerical_hemispherical_integration

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

         elseif(varlabel.eq.'check_positions') then
            vartype='l'
            lvarvalue=>check_positions

         elseif(varlabel.eq.'medium_ref_index') then
            vartype='c'
            cvarvalue=>medium_ref_index
            medium_ref_index_specified=.true.
            medium_reim_ref_index_specified=.false.

         elseif(varlabel.eq.'medium_re_ref_index') then
            vartype='r'
            rvarvalue=>medium_re_ref_index
            medium_ref_index_specified=.true.
            medium_reim_ref_index_specified=.true.

         elseif(varlabel.eq.'medium_im_ref_index') then
            vartype='r'
            rvarvalue=>medium_im_ref_index
            medium_ref_index_specified=.true.
            medium_reim_ref_index_specified=.true.

         elseif(varlabel.eq.'light_up') then
            vartype='l'
            lvarvalue=>light_up

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
                  if(n.gt.input_number_spheres) cycle
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
               input_number_spheres=min(n,input_number_spheres)
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
               medium_ref_index_specified=.false.

            elseif(trim(parmid).eq.'layer_thickness') then
               parmval=inputfiledata(inputline)
               inputline=inputline+1
               input_layer_thickness(1:max(1,number_plane_boundaries))=0.d0
               read(parmval,*,iostat=istat) input_layer_thickness(1:max(1,number_plane_boundaries))
               recalculate_surface_matrix=.true.

            elseif(trim(parmid).eq.'component_radii') then
               call read_real_list(component_radii,number_components)

            elseif(trim(parmid).eq.'component_number_fraction') then
               call read_real_list(component_number_fraction,number_components)

            elseif(trim(parmid).eq.'psd_sigma') then
               call read_real_list(psd_sigma,number_components)

            elseif(trim(parmid).eq.'component_ref_index') then
               call read_cmplx_list(component_ref_index,number_components)

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

         contains
            subroutine read_real_list(listvar,listnum)
            implicit none
            integer :: listnum
            real(8) :: listvar(*)
            parmval=inputfiledata(inputline)
            inputline=inputline+1
            read(parmval,*,iostat=istat) listvar(1:listnum)
            end subroutine read_real_list

            subroutine read_cmplx_list(listvar,listnum)
            implicit none
            integer :: listnum
            complex(8) :: listvar(*)
            parmval=inputfiledata(inputline)
            inputline=inputline+1
            read(parmval,*,iostat=istat) listvar(1:listnum)
            end subroutine read_cmplx_list
         end subroutine inputdata

         subroutine main_calling_program(print_output,set_t_matrix_order,dry_run,mpi_comm)
         implicit none
         logical :: stopit,singleorigin,iframe,sett,printout,dryrun,averagerun
         logical, optional :: print_output,set_t_matrix_order,dry_run
         integer :: n,istat,niter,rank,numprocs,i,nodrw,celldim(3),itemp(6),sx,sy,maxt, &
            mpicomm,lochost
         integer, optional :: mpi_comm
         real(8) :: alpha,time1,r0(3),rtran,costheta, &
            csca,zext,targetvol,timet,tmin(3),tmax(3),rannum
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
         if(reflection_model) then
            calculate_up_down_scattering=.true.
            incident_frame=.false.
         endif

         first_run=.false.
         call mstm_mpi(mpi_command='rank',mpi_rank=rank,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs,mpi_comm=mpicomm)
!         if(rank.ne.0) light_up=.false.
         local_rank=rank
         global_rank=rank
         if((.not.configuration_average).and.(.not.incidence_average)) n_configuration_groups=1
!         random_configuration=(trim(sphere_data_input_file).eq.'random_configuration')
         if(random_configuration) then
            if(auto_target_radius.and.target_shape.eq.2) then
               number_spheres=input_number_spheres
               target_dimensions(1:3)=(dble(number_spheres)/sphere_volume_fraction)**(1.d0/3.d0)
            else
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
               call target_volume(target_dimensions,targetvol)
               if(number_spheres_specified) then
                  number_spheres=input_number_spheres
                  sphere_volume_fraction=dble(input_number_spheres)*4.d0*pi/3.d0/targetvol
               else
                  number_spheres=ceiling(targetvol*sphere_volume_fraction)/(4.d0*pi/3.d0)
               endif
            endif
            if(target_shape.eq.2.and.random_configuration_host) then
               number_spheres=number_spheres+1
            endif
         else
            number_spheres=input_number_spheres
         endif
         if(medium_ref_index_specified) then
            if(medium_reim_ref_index_specified) then
               layer_ref_index(0)=dcmplx(medium_re_ref_index,medium_im_ref_index)
            else
               layer_ref_index(0)=medium_ref_index
            endif
         endif
         if(allocated(sphere_radius)) then
            deallocate(sphere_radius, &
                  sphere_position, &
                  sphere_ref_index, &
                  host_sphere, &
                  number_field_expansions, &
                  sphere_excitation_switch, &
                  sphere_index)
         endif
         allocate(sphere_radius(number_spheres), &
                  sphere_position(3,number_spheres), &
                  sphere_ref_index(2,0:number_spheres), &
                  host_sphere(number_spheres), &
                  number_field_expansions(number_spheres), &
                  sphere_excitation_switch(number_spheres), &
                  sphere_index(number_spheres))
         if(random_configuration) then
            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' generating random configuration:'',$)')
               timet=mstm_mpi_wtime()
            endif
            call generate_random_configuration(mpi_comm=mpicomm,skip_diffusion=dryrun)
!            call generate_random_configuration(mpi_comm=mpicomm)
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
                  open(2,file=trim(random_configuration_output_file))
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

         if(configuration_average.and.(target_shape.eq.2) &
            .and.random_configuration_host.and.auto_target_radius) then
            host_sphere(1:number_spheres-1)=number_spheres
            number_field_expansions(1:number_spheres-1)=1
            number_host_spheres=number_spheres
            host_sphere(number_spheres)=0
            number_field_expansions(number_spheres)=2
         endif

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
         if(singleorigin.or.random_orientation.or..true.) then
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
            t_matrix_order=min(max_mie_order,maxt)
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
            if(d_cell_specified) d_cell=input_d_cell
            lochost=0
            if(random_configuration) then
               tmin=-target_dimensions*length_scale_factor
               tmax=target_dimensions*length_scale_factor
               if(dryrun) call clear_fft_matrix(clear_h=.true.)
               if(target_shape.eq.2.and.random_configuration_host) then
                  lochost=number_spheres
               endif
            else
               tmin=sphere_min_position
               tmax=sphere_max_position
               if(.not.averagerun) call clear_fft_matrix(clear_h=.true.)
            endif
            call node_selection(cell_volume_fraction,target_min=tmin,target_max=tmax,&
               d_specified=d_cell_specified,local_host=lochost)
            if(input_node_order.le.0) then
               node_order=-input_node_order+ceiling(d_cell)
            else
               node_order=input_node_order
            endif
!            node_order=max(node_order,max_mie_order)
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
            if(excitation_radius.gt.0.d0) then
               sphere_excitation_switch(i)=rtran.le.excitation_radius*length_scale_factor
            else
               sphere_excitation_switch(i)=i.le.-int(excitation_radius)
            endif
         enddo
         if(excitation_radius.eq.0.d0) then
            sphere_excitation_switch=.false.
            if(rank.eq.0) then
               call random_seed()
               do
                  call random_number(rannum)
                  itemp(1)=1+floor(number_spheres*rannum)
                  if(host_sphere(itemp(1)).eq.0) exit
               enddo
            endif
            call mstm_mpi(mpi_command='bcast',mpi_rank=0, &
               mpi_send_buf_i=itemp(1),mpi_number=1,mpi_comm=mpicomm)
            sphere_excitation_switch(itemp(1))=.true.
         endif
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
               scat_mat_amin=0.d0
               scat_mat_amax=180.d0
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
!                        scat_mat_ldim=-scat_mat_udim
! 10-22 azimuthal_average applies to multiple plane boundaries
                        if(azimuthal_average) then
                           scat_mat_ldim=0
                        else
                           scat_mat_ldim=-scat_mat_udim
                        endif
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
!         allocate(boundary_sca(2,0:number_plane_boundaries+1),boundary_ext(2,0:number_plane_boundaries+1))
         allocate(boundary_sca(2,0:1),boundary_ext(2,0:1))

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
!call mstm_mpi(mpi_command='barrier')
         if(rank.eq.0.and.printout) then
            if(check_positions) call checkpositions()
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
            if(allocated(mean_t)) deallocate(mean_t)
            allocate(mean_t(2,t_matrix_order))
            mean_t=0.d0
            call tmatrix_solution(solution_method=solution_method(1:1), &
               solution_eps=solution_epsilon, &
               convergence_eps=t_matrix_convergence_epsilon, &
               max_iterations=niter, &
               t_matrix_file=t_matrix_output_file, &
               procs_per_soln=t_matrix_procs_per_solution, &
               sphere_qeff=q_eff, &
               solution_status=istat, &
               mpi_comm=mpicomm, &
               sphere_excitation_list=sphere_excitation_switch)
            if(fft_translation_option) call clear_fft_matrix(clear_h=.true.)
            if(calculate_scattering_matrix) then
               if(allocated(scat_mat_exp_coef)) deallocate(scat_mat_exp_coef)
               if(allocated(coh_scat_mat_exp_coef)) deallocate(coh_scat_mat_exp_coef)
               allocate(scat_mat_exp_coef(4,4,0:2*t_matrix_order),coh_scat_mat_exp_coef(4,4,0:2*t_matrix_order))
               nodrw=2*t_matrix_order
               call ranorientscatmatrix(t_matrix_output_file,scat_mat_exp_coef, &
                  coh_scat_mat_exp_coef, &
                  beam_width=gaussian_beam_constant, &
                  number_processors=t_matrix_procs_per_solution, &
                  mean_t_matrix=mean_t,mpi_comm=mpicomm)
               coherent_scattering_ratio=coh_scat_mat_exp_coef(1,1,0)
               do i=scat_mat_ldim,scat_mat_udim
                  costheta=cos(dble(i-scat_mat_ldim)*pi/dble(scat_mat_udim-scat_mat_ldim))
                  call ranorienscatmatrixcalc(costheta,scat_mat_exp_coef,nodrw,scat_mat(:,i))
               enddo
            endif
            call qtotcalc(number_spheres,qeff_dim,cross_section_radius,&
               q_eff,q_vabs,q_eff_tot)
         else
if(light_up) then
write(*,'('' s8.1 '',i3)') mstm_global_rank
call flush(6)
endif
!call mstm_mpi(mpi_command='barrier')

            if(allocated(amnp_s)) deallocate(amnp_s)
            allocate(amnp_s(number_eqns,2))
            amnp_s=0.d0
            niter=max_iterations
            error_codes=0
            pl_error_codes=0
!            recalculate_surface_matrix=.true.
if(light_up) then
write(*,'('' s8.2 '',i3)') mstm_global_rank
call flush(6)
endif
!call mstm_mpi(mpi_command='barrier')
            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' generating solution:'',$)')
               timet=mstm_mpi_wtime()
            endif
            call fixedorsoln(alpha,incident_sin_beta,incident_direction,solution_epsilon,niter,amnp_s,q_eff, &
                    qeff_dim,solution_error,solution_iterations,1,istat, &
                    mpi_comm=mpicomm, &
                    excited_spheres=sphere_excitation_switch, &
                    solution_method=solution_method(1:1), &
                    initialize_solver=.true.)
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
!call mstm_mpi(mpi_command='barrier')
            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' post processing solution:'',$)')
               timet=mstm_mpi_wtime()
            endif

            call qtotcalc(number_spheres,qeff_dim,cross_section_radius,&
               q_eff,q_vabs,q_eff_tot)
!            q_eff_tot(3,:)=q_eff_tot(1,:)-q_eff_tot(2,:)
            csca=q_eff_tot(3,1)*pi*cross_section_radius**2
            if(singleorigin) then
               if(allocated(amnp_0)) deallocate(amnp_0)
               allocate(amnp_0(2*t_matrix_order*(t_matrix_order+2),2))
               amnp_0=0.d0
if(light_up) then
write(*,'('' s8.3.1 '',i3)') mstm_global_rank
call flush(6)
endif
!call mstm_mpi(mpi_command='barrier')
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
!call mstm_mpi(mpi_command='barrier')
               if(singleorigin.and.azimuthal_average.and.(.not.numerical_azimuthal_average)) then
                  if(allocated(scat_mat_exp_coef)) deallocate(scat_mat_exp_coef)
                  allocate(scat_mat_exp_coef(16,0:2*t_matrix_order,4))
                  call fosmexpansion(t_matrix_order,amnp_0,scat_mat_exp_coef(:,:,1),scat_mat_exp_coef(:,:,2), &
                     scat_mat_exp_coef(:,:,3),scat_mat_exp_coef(:,:,4),mpi_comm=mpicomm)
               endif
if(light_up) then
write(*,'('' s8.3.4 '',i3)') mstm_global_rank
call flush(6)
endif
!call mstm_mpi(mpi_command='barrier')
               if(calculate_scattering_matrix) then
                  call scattering_matrix_calculation(amnp_0,scat_mat,mpi_comm=mpicomm)
               endif
if(light_up) then
write(*,'('' s8.3.5 '',i3)') mstm_global_rank
call flush(6)
endif
!call mstm_mpi(mpi_command='barrier')
               if(gaussian_beam_constant.eq.0.d0) then
                  call boundary_extinction(amnp_0,alpha,incident_sin_beta,incident_direction,boundary_ext, &
                     common_origin=singleorigin)
               endif
            else
if(light_up) then
write(*,'('' s8.3.5 '',i3)') mstm_global_rank
call flush(6)
endif
!call mstm_mpi(mpi_command='barrier')
if(light_up) then
write(*,'('' s8.3.5 '',i3)') mstm_global_rank
call flush(6)
endif
!call mstm_mpi(mpi_command='barrier')
               if(calculate_scattering_matrix) then
                  call scattering_matrix_calculation(amnp_s,scat_mat,mpi_comm=mpicomm)
               endif
               if(gaussian_beam_constant.eq.0.d0) then
                  call boundary_extinction(amnp_s,alpha,incident_sin_beta,incident_direction,boundary_ext)
               endif
            endif

            if(calculate_up_down_scattering) then
               if(singleorigin) then
                  call hemispherical_scattering(amnp_0,.true.,numerical_hemispherical_integration, &
                     boundary_sca,mpi_comm=mpicomm)
               else
                  call hemispherical_scattering(amnp_s,.false.,numerical_hemispherical_integration, &
                     boundary_sca,mpi_comm=mpicomm)
               endif
            endif

            if(gaussian_beam_constant.ne.0.d0) then
               boundary_ext=0
               boundary_ext(1:2,1)=-q_eff_tot(1,2:3)
            endif

            if(periodic_lattice) then
               call periodic_lattice_scattering(amnp_s,pl_sca)
            elseif(reflection_model) then
!               pl_sca(:,1)=boundary_sca(:,number_plane_boundaries+1)
               pl_sca(:,1)=boundary_sca(:,1)
               pl_sca(:,2)=-boundary_sca(:,0)
            endif

            if(periodic_lattice.or.reflection_model) then
               call surface_absorptance_calculation()
            endif

            if(number_plane_boundaries.gt.0.and..not.periodic_lattice) then
!               evan_sca(1:2)=q_eff_tot(1,2:3)-q_eff_tot(2,2:3)+boundary_sca(1:2,0) &
!                  - boundary_sca(:,number_plane_boundaries+1)
               evan_sca(1:2)=q_eff_tot(1,2:3)-q_eff_tot(2,2:3)+boundary_sca(1:2,0) &
                  - boundary_sca(:,1)
            endif

            if(print_timings.and.mstm_global_rank.eq.0) then
               write(run_print_unit,'('' completed, time:'',es12.5,'' s'')') mstm_mpi_wtime()-timet
            endif

         endif

         if(rank.eq.0) solution_time=mstm_mpi_wtime()-time1

         call gather_error_codes(mpicomm)
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
         error_codes=0
         pl_error_codes=0

         if(calculate_near_field.and.(.not.random_orientation)) then
            celldim=ceiling((near_field_plane_vertices(:,2)-near_field_plane_vertices(:,1))/near_field_step_size)
            celldim=max(celldim,(/1,1,1/))
            if(configuration_average) then
               if(allocated(e_field)) deallocate(e_field,h_field)
               if(rank.eq.0) then
                  allocate(e_field(3,2,celldim(1),celldim(2),celldim(3)), &
                     h_field(3,2,celldim(1),celldim(2),celldim(3)))
                  e_field=0.d0
                  h_field=0.d0
               endif
               call near_field_calculation(amnp_s,alpha,incident_sin_beta,incident_direction, &
                  near_field_plane_vertices,celldim, &
                  incident_model=near_field_calculation_model,output_unit=0, &
                  e_field_array=e_field,h_field_array=h_field,mpi_comm=mpicomm)
            else
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
            call gather_error_codes(mpicomm)
            if(rank.eq.0) call print_error_codes(run_print_unit)
         endif

         end subroutine main_calling_program

         subroutine gather_error_codes(mpicomm)
         implicit none
         integer :: mpicomm,itemp(6)
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
         end subroutine gather_error_codes

         subroutine configuration_average_calling_program()
         implicit none
         logical :: singleorigin,iframe
         integer :: rank,numprocs,m,n,p,mnp,griddim(3),ipos(3),ix,iy,iz, &
            numprocsperconfig,configcolor,configgroup,configcomm,configrank,config0comm,nconfigave,nsend
         real(8) :: time1,timet,diffac,csca(1),xspfit,rpos(3),rtemp(1)
         real(8), allocatable :: texpcoef(:,:,:)
         complex(8) :: ritemp(2),aneff,ctemp(1),rieff,e0
         complex(8), allocatable :: pmnp0(:,:),anp0(:,:),edat(:)
         character*256 :: tmatchar1,tmatchar2
         data tmatchar1,tmatchar2/'tmat-','.tmp'/
         first_run=.false.
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
         local_rank=rank
         global_rank=rank

         if(max_iterations.le.1) then
            numprocsperconfig=2
         else
            numprocsperconfig=4
         endif
!numprocsperconfig=2

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
            sphere_position_ave(3,number_spheres),boundary_sca_ave(2,0:1), &
            boundary_ext_ave(2,0:1),dif_boundary_sca(2,0:1))
         q_eff_ave=0.d0
         q_eff_tot_ave=0.d0
         q_vabs_ave=0.d0
         sphere_position_ave=0.d0
         pl_sca_ave=0.d0
         boundary_sca_ave=0.d0
         boundary_ext_ave=0.d0
         solution_time_ave=0.d0
         surface_absorptance_ave=0.
         effective_medium_simulation=.false.
         if(singleorigin) then
            if(allocated(amnp_0_ave)) deallocate(amnp_0_ave,scat_mat_exp_coef_ave,coh_scat_mat_exp_coef)
            allocate(amnp_0_ave(2*t_matrix_order*(t_matrix_order+2),2), &
               scat_mat_exp_coef_ave(16,0:2*t_matrix_order,4), &
               coh_scat_mat_exp_coef(16,0:2*t_matrix_order,4))
            amnp_0_ave=0.d0
            scat_mat_exp_coef_ave=0.d0
            tot_csca_ave=0.d0
            dif_csca_ratio=0.d0
            if(input_effective_medium_simulation.and.target_shape.eq.2.and. &
             (.not.random_configuration_host)) then
               effective_medium_simulation=.true.
               effective_ref_index=layer_ref_index(0)
               if(random_configuration_host_model.eq.1) then
                  effective_cluster_radius=target_dimensions(1)*length_scale_factor
               elseif(random_configuration_host_model.eq.2) then
                  effective_cluster_radius=vol_radius/(sphere_volume_fraction)**0.33333
               endif
            endif
         endif
         if(calculate_scattering_matrix) then
            if(allocated(scat_mat_ave)) deallocate(scat_mat_ave,dif_scat_mat)
            allocate(scat_mat_ave(scat_mat_mdim,scat_mat_ldim:scat_mat_udim), &
               dif_scat_mat(scat_mat_mdim,scat_mat_ldim:scat_mat_udim))
            scat_mat_ave=0.d0
         endif
         if(calculate_near_field.and.configrank.eq.0) then
            griddim=ceiling((near_field_plane_vertices(:,2)-near_field_plane_vertices(:,1))/near_field_step_size)
            griddim=max(griddim,(/1,1,1/))
            if(allocated(e_field_ave)) deallocate(e_field_ave,h_field_ave,s_field,s_field_ave)
            allocate(e_field_ave(3,2,griddim(1),griddim(2),griddim(3)), &
               h_field_ave(3,2,griddim(1),griddim(2),griddim(3)), &
               s_field(3,2,griddim(1),griddim(2),griddim(3)), &
               s_field_ave(3,2,griddim(1),griddim(2),griddim(3)))
            e_field_ave=0.d0
            h_field_ave=0.d0
            s_field_ave=0.d0
         endif

         nconfigave=0
         do random_configuration_number=1,ceiling(dble(number_configurations)/dble(n_configuration_groups))

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

            if(singleorigin.and.configrank.eq.0) then
               call common_origin_csca(t_matrix_order,amnp_0,csca)
            endif

            if(configrank.eq.0) then
               if(rank.eq.0) solution_time=mstm_mpi_wtime()-time1
               q_eff_ave=q_eff_ave+q_eff
               q_eff_tot_ave=q_eff_tot_ave+q_eff_tot
               q_vabs_ave=q_vabs_ave+q_vabs
               if(target_shape.eq.2)  &
                  call cartospherevec(number_spheres,sphere_position,sphere_position)
               sphere_position_ave=sphere_position_ave+sphere_position
               pl_sca_ave=pl_sca_ave+pl_sca
               surface_absorptance_ave=surface_absorptance_ave+surface_absorptance
               if(calculate_up_down_scattering) boundary_sca_ave=boundary_sca_ave+boundary_sca
               boundary_ext_ave=boundary_ext_ave+boundary_ext
               if(singleorigin.and.azimuthal_average.and.(.not.numerical_azimuthal_average)) &
                  scat_mat_exp_coef_ave=scat_mat_exp_coef_ave+scat_mat_exp_coef
               if(calculate_scattering_matrix) then
                  scat_mat_ave=scat_mat_ave+scat_mat
               endif
               if(calculate_near_field) then
                  e_field_ave=e_field_ave+e_field
                  h_field_ave=h_field_ave+h_field
                  s_field(1,:,:,:,:)=0.5*(e_field(2,:,:,:,:)*conjg(h_field(3,:,:,:,:)) &
                     -e_field(3,:,:,:,:)*conjg(h_field(2,:,:,:,:)))
                  s_field(2,:,:,:,:)=0.5*(-e_field(1,:,:,:,:)*conjg(h_field(3,:,:,:,:)) &
                     +e_field(3,:,:,:,:)*conjg(h_field(1,:,:,:,:)))
                  s_field(3,:,:,:,:)=0.5*(e_field(1,:,:,:,:)*conjg(h_field(2,:,:,:,:)) &
                     -e_field(2,:,:,:,:)*conjg(h_field(1,:,:,:,:)))
                  s_field_ave=s_field_ave+s_field
               endif
               if(rank.eq.0) solution_time_ave=solution_time_ave+solution_time
            endif

            if(singleorigin) then
!               if(sphere_1_fixed) call subtract_1_from_0()
!               if(sphere_1_fixed) then
!                  erase_sphere_1=.true.
!                  use_previous_configuration=.true.
!                  sphere_1_fixed=.false.
!                  call main_calling_program(print_output=.false.,set_t_matrix_order=.false.,mpi_comm=configcomm)
!                  erase_sphere_1=.false.
!                  use_previous_configuration=.false.
!                  sphere_1_fixed=.true.
!               endif
               amnp_0_ave=amnp_0_ave+amnp_0
               tot_csca_ave=tot_csca_ave+csca
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
               nsend=4
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
               nsend=1
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=tot_csca_ave, &
                  mpi_recv_buf_dp=csca, &
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
            if(singleorigin.and.azimuthal_average.and.(.not.numerical_azimuthal_average)) then
               nsend=16*4*(2*t_matrix_order+1)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=scat_mat_exp_coef_ave, &
                  mpi_recv_buf_dp=scat_mat_exp_coef, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif
            if(calculate_near_field) then
               nsend=6*product(griddim)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dc=e_field_ave, &
                  mpi_recv_buf_dc=e_field, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dc=h_field_ave, &
                  mpi_recv_buf_dc=h_field, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=s_field_ave, &
                  mpi_recv_buf_dp=s_field, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif

            nconfigave=nconfigave+n_configuration_groups
            diffac=(dble(number_spheres)-1.d0)/dble(number_spheres)
            diffac=1.d0
!            diffac=(dble(number_spheres))/dble(number_spheres-1)
!            diffac=(1.d0-(1.d0/dble(number_spheres))**.5d0)

            if(singleorigin) then
               if(rank.eq.0.and.print_timings) then
                  timet=mstm_mpi_wtime()
                  write(run_print_unit,'('' calculating diffuse field:'',$)')
               endif
               amnp_0=amnp_0/dble(nconfigave)
!
!  zero out azimuth order .ne. pm 1 for sphere targets 2/23
!
!               if(target_shape.eq.2) then
!                  do n=1,t_matrix_order
!                     do m=-n,n
!                        do p=1,2
!                           if(abs(m).ne.1) then
!                              amnp_0(amnpaddress(m,n,p,t_matrix_order,2),:)=0.d0
!                           endif
!                        enddo
!                     enddo
!                  enddo
!               endif
               csca=csca/dble(nconfigave)
               tot_csca=csca(1)
               call common_origin_csca(t_matrix_order,amnp_0,dif_csca_ratio)
               dif_csca=tot_csca-dif_csca_ratio(1)
               dif_csca_ratio=1.d0-dif_csca_ratio/csca
               if(azimuthal_average.and.(.not.numerical_azimuthal_average)) then
                  allocate(texpcoef(16,0:2*t_matrix_order,4))
                  texpcoef=scat_mat_exp_coef/dble(nconfigave)
!                  scat_mat_exp_coef=scat_mat_exp_coef/dble(nconfigave)
                  call fosmexpansion(t_matrix_order,amnp_0,coh_scat_mat_exp_coef(:,:,1),coh_scat_mat_exp_coef(:,:,2), &
                     coh_scat_mat_exp_coef(:,:,3),coh_scat_mat_exp_coef(:,:,4),mpi_comm=configcomm)
                  scat_mat_exp_coef=coh_scat_mat_exp_coef
               endif
               if(calculate_scattering_matrix) then
                  call scattering_matrix_calculation(amnp_0,dif_scat_mat,mpi_comm=configcomm)
               endif
               if(azimuthal_average.and.(.not.numerical_azimuthal_average)) then
!                  scat_mat_exp_coef=texpcoef-scat_mat_exp_coef*diffac
                  scat_mat_exp_coef=texpcoef
                  deallocate(texpcoef)
               endif
               call hemispherical_scattering(amnp_0,.true.,numerical_hemispherical_integration, &
                     dif_boundary_sca,mpi_comm=configcomm)
!               call common_origin_hemispherical_scattering(amnp_0,dif_boundary_sca)
               if(rank.eq.0.and.print_timings) write(run_print_unit,'('' completed, '',es12.4,'' sec'')') mstm_mpi_wtime()-timet
               if(rank.eq.0.and.target_shape.eq.2.and.(.not.random_configuration_host)) then
                  allocate(pmnp0(2*t_matrix_order*(t_matrix_order+2),2),anp0(2,t_matrix_order))
                  call genplanewavecoef(0.d0,(1.d0,0.d0),t_matrix_order,pmnp0,lr_tran=.false.)
                  open(30,file='anpeff.dat')
                  write(30,'(i5)') t_matrix_order
                  do n=1,t_matrix_order
                     do p=1,2
                        aneff=0.d0
                        do m=-1,1,2
                           mnp=amnpaddress(m,n,p,t_matrix_order,2)
                           aneff=aneff+0.5d0*sum(amnp_0(mnp,:)/pmnp0(mnp,:))
                        enddo
                        aneff=aneff/2.d0
                        write(30,'(2i4,2es13.5)') n,p,aneff
                        anp0(p,n)=aneff
                     enddo
                  enddo
                  close(30)
                  deallocate(pmnp0)
!                  effective_ref_index=(1.1d0,0.01d0)
!                  fit_radius=target_dimensions(1)*length_scale_factor
                  call effective_ref_index_fit(anp0,effective_ref_index,fit_radius,fit_stat)
                  deallocate(anp0)
               endif
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
!                  if(singleorigin) dif_scat_mat=scat_mat-dif_scat_mat*diffac

!                  dif_scat_mat=scat_mat-dif_scat_mat
               endif
               solution_time=solution_time_ave/dble(random_configuration_number)
               call print_calculation_results(output_file)
               if(calculate_near_field) then
                  e_field=e_field/dble(nconfigave)
                  h_field=h_field/dble(nconfigave)
                  s_field=s_field/dble(nconfigave)
                  open(2,file=near_field_output_file)
                  call write_output_header(griddim,2,print_intersecting_spheres=.false.)
                  allocate(edat(griddim(3)))
                  edat=0.d0
                  do iz=1,griddim(3)
                     do iy=1,griddim(2)
                        do ix=1,griddim(1)
                           edat(iz)=edat(iz)+e_field(1,1,ix,iy,iz)+e_field(2,2,ix,iy,iz)
                           ipos(:)=(/ix,iy,iz/)
                           rpos(:)=(dble(ipos(:))-(/0.5d0,0.5d0,0.5d0/))*grid_spacing(:)+grid_region(:,1)
                           write(2,'(33es12.4)') rpos(:),e_field(:,1,ix,iy,iz),h_field(:,1,ix,iy,iz), &
                              e_field(:,2,ix,iy,iz),h_field(:,2,ix,iy,iz), &
                              s_field(:,1,ix,iy,iz),s_field(:,2,ix,iy,iz)
                        enddo
                     enddo
                  enddo
                  edat=edat/dble(griddim(1)*griddim(2)*2.d0)
                  call effectiverefractiveindex(griddim(3),edat,grid_spacing(3),rieff,e0)
                  nf_eff_ref_index=rieff
!                  write(*,'('' field fit ri:'',2es12.5)') rieff
                  deallocate(edat)
               endif
            endif
            if(effective_medium_simulation) then
               ctemp=effective_ref_index
               call mstm_mpi(mpi_command='bcast',mpi_rank=0,mpi_number=1, &
                 mpi_send_buf_dc=ctemp)
                 effective_ref_index=ctemp(1)
            endif
         enddo

         end subroutine configuration_average_calling_program

         subroutine random_orientation_configuration_average_calling_program()
         implicit none
         integer :: rank,numprocs,itemp(1),n, &
            numprocsperconfig,configcolor,configgroup,configcomm,configrank,config0comm,nconfigave,nsend
         real(8) :: time1,timet
         real(8), allocatable :: tpos(:,:)
         complex(8) :: ctemp(1)
         character*256 :: tmatchar1,tmatchar2
         data tmatchar1,tmatchar2/'tmat-','.tmp'/
         first_run=.false.
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)
         call mstm_mpi(mpi_command='size',mpi_size=numprocs)
!         if(rank.ne.0) light_up=.false.
         local_rank=rank
         global_rank=rank

         if(max_iterations.le.0) then
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
!singleorigin=.true.

         call main_calling_program(print_output=.false.,set_t_matrix_order=.true.,dry_run=.true.)
         call groupfilename(tmatchar1,configgroup,tmatchar2,t_matrix_output_file)

         if(allocated(q_eff_ave)) deallocate(q_eff_ave,q_eff_tot_ave,q_vabs_ave,sphere_position_ave)
         if(allocated(mean_t_ave)) deallocate(mean_t_ave)
         if(allocated(scat_mat_exp_coef_ave)) deallocate(scat_mat_exp_coef_ave)
         if(allocated(coh_scat_mat_exp_coef_ave)) deallocate(coh_scat_mat_exp_coef_ave)
         allocate(q_eff_ave(3,qeff_dim,number_spheres),q_eff_tot_ave(3,qeff_dim),q_vabs_ave(qeff_dim,number_spheres), &
            sphere_position_ave(3,number_spheres),mean_t_ave(2,t_matrix_order), &
            scat_mat_exp_coef_ave(4,4,0:2*t_matrix_order),coh_scat_mat_exp_coef_ave(4,4,0:2*t_matrix_order))
         q_eff_ave=0.d0
         q_eff_tot_ave=0.d0
         q_vabs_ave=0.d0
         sphere_position_ave=0.d0
         mean_t_ave=0.d0
         scat_mat_exp_coef_ave=0.d0
         coh_scat_mat_exp_coef_ave=0.d0
         if(calculate_scattering_matrix) then
            if(allocated(scat_mat_ave)) deallocate(scat_mat_ave)
            allocate(scat_mat_ave(scat_mat_mdim,scat_mat_ldim:scat_mat_udim))
            scat_mat_ave=0.d0
         endif

         nconfigave=0
         do random_configuration_number=1,ceiling(dble(number_configurations)/dble(n_configuration_groups))

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

            if(configrank.eq.0) then
               if(rank.eq.0) solution_time=mstm_mpi_wtime()-time1
               allocate(tpos(3,number_spheres))
               call cartospherevec(number_spheres,sphere_position(:,1:number_spheres),tpos(:,1:number_spheres))
               q_eff_ave=q_eff_ave+q_eff
               q_eff_tot_ave=q_eff_tot_ave+q_eff_tot
               q_vabs_ave=q_vabs_ave+q_vabs
               sphere_position_ave=sphere_position_ave+tpos
               mean_t_ave(:,1:t_matrix_order)=mean_t_ave(:,1:t_matrix_order)+mean_t(:,1:t_matrix_order)
               if(calculate_scattering_matrix) then
                  scat_mat_ave=scat_mat_ave+scat_mat
                  scat_mat_exp_coef_ave=scat_mat_exp_coef_ave+scat_mat_exp_coef
                  coh_scat_mat_exp_coef_ave=coh_scat_mat_exp_coef_ave+coh_scat_mat_exp_coef
               endif
               if(rank.eq.0) solution_time_ave=solution_time_ave+solution_time
               deallocate(tpos)
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
            nsend=2*t_matrix_order
            call mstm_mpi(mpi_command='reduce', &
               mpi_operation=mstm_mpi_sum, &
               mpi_send_buf_dc=mean_t_ave, &
               mpi_recv_buf_dc=mean_t, &
               mpi_rank=0, &
               mpi_number=nsend, &
               mpi_comm=config0comm)
            if(calculate_scattering_matrix) then
               nsend=scat_mat_mdim*(scat_mat_udim-scat_mat_ldim+1)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=scat_mat_ave, &
                  mpi_recv_buf_dp=scat_mat, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
               nsend=16*(2*t_matrix_order+1)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=scat_mat_exp_coef_ave, &
                  mpi_recv_buf_dp=scat_mat_exp_coef, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
               call mstm_mpi(mpi_command='reduce', &
                  mpi_operation=mstm_mpi_sum, &
                  mpi_send_buf_dp=coh_scat_mat_exp_coef_ave, &
                  mpi_recv_buf_dp=coh_scat_mat_exp_coef, &
                  mpi_rank=0, &
                  mpi_number=nsend, &
                  mpi_comm=config0comm)
            endif

            nconfigave=nconfigave+n_configuration_groups

            if(rank.eq.0) then
               sphere_position=sphere_position/dble(nconfigave)
               q_eff=q_eff/dble(nconfigave)
               q_vabs=q_vabs/dble(nconfigave)
               q_eff_tot=q_eff_tot/dble(nconfigave)
               mean_t(:,1:t_matrix_order)=mean_t(:,1:t_matrix_order)/dble(nconfigave)
               if(calculate_scattering_matrix) then
                  scat_mat=scat_mat/dble(nconfigave)
                  scat_mat_exp_coef=scat_mat_exp_coef/dble(nconfigave)
                  coh_scat_mat_exp_coef=coh_scat_mat_exp_coef/dble(nconfigave)
               endif
               solution_time=solution_time_ave/dble(random_configuration_number)
               call effective_ref_index_fit(mean_t,effective_ref_index,fit_radius,fit_stat)
               call print_calculation_results(output_file)
            endif
!            if(random_configuration_host) then
!               ctemp=effective_ref_index
!               if(rank.eq.0) write(*,'('' new ri:'',2es12.4)') effective_ref_index
!               call mstm_mpi(mpi_command='bcast',mpi_rank=0,mpi_number=1, &
!                 mpi_send_buf_dc=ctemp)
!               layer_ref_index(0)=ctemp(1)
!            endif
!if(rank.eq.0) then
!write(*,'(4es12.4)') layer_ref_index(0),sphere_ref_index(1,0)
!endif
         enddo
         end subroutine random_orientation_configuration_average_calling_program

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
         singleorigin=number_plane_boundaries.eq.0.and.single_origin_expansion
         incident_beta_specified=.true.

         call main_calling_program(print_output=.false.,set_t_matrix_order=.true.,dry_run=.true.)

         if(trim(sphere_data_input_file).eq.'random_configuration') then
            sphere_data_input_file='random_configuration.pos'
         endif

         if(allocated(q_eff_ave)) deallocate(q_eff_ave,q_eff_tot_ave,q_vabs_ave,boundary_sca_ave, &
             boundary_ext_ave,dif_boundary_sca)
         allocate(q_eff_ave(3,qeff_dim,number_spheres),q_eff_tot_ave(3,qeff_dim),q_vabs_ave(qeff_dim,number_spheres), &
            boundary_sca_ave(2,0:1),boundary_ext_ave(2,0:1),dif_boundary_sca(2,0:1))
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
               if(singleorigin.and.azimuthal_average.and.(.not.numerical_azimuthal_average)) &
                   scat_mat_exp_coef_ave=scat_mat_exp_coef_ave+scat_mat_exp_coef
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
               nsend=4
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
            if(singleorigin.and.azimuthal_average.and.(.not.numerical_azimuthal_average)) then
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
               if(azimuthal_average.and.(.not.numerical_azimuthal_average)) then
                  allocate(texpcoef(16,0:2*t_matrix_order,4))
                  texpcoef=scat_mat_exp_coef/dble(nconfigave)
                  call fosmexpansion(t_matrix_order,amnp_0,scat_mat_exp_coef(:,:,1),scat_mat_exp_coef(:,:,2), &
                     scat_mat_exp_coef(:,:,3),scat_mat_exp_coef(:,:,4),mpi_comm=configcomm)
               endif
               if(calculate_scattering_matrix) then
                  call scattering_matrix_calculation(amnp_0,dif_scat_mat,mpi_comm=configcomm)
               endif
               if(azimuthal_average.and.(.not.numerical_azimuthal_average)) then
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

         subroutine common_origin_csca(n,a,c)
         implicit none
         integer :: n
         real(8) :: c(1)
         complex(8) :: a(4*n*(n+2))
         c(1)=sum(a(:)*dconjg(a(:)))
         end subroutine common_origin_csca

         subroutine subtract_1_from_0()
         implicit none
         integer :: i,i1,mnp0,mnp1,n,m,p
         real(8) :: fn
         complex(8) :: a(2,2),b(2,2)

         fn=dble(number_spheres-1)/dble(number_spheres)
         fn=1.d0
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
!                  amnp_0=amnp_0*dble(number_spheres)/dble(number_spheres-1)
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
            scaq=1.d0-0.5d0*(-sum(dif_boundary_sca(:,0))+sum(dif_boundary_sca(:,1)))-q_eff_tot(2,1)
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
            scaq=0.5d0*(-sum(dif_boundary_sca(:,0))+sum(dif_boundary_sca(:,1)))+q_eff_tot(2,1)
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

         subroutine effectiverefractiveindex(ndat,edat,d,rieff, &
                    e0)
         use mpidefs
         implicit none
         integer :: ndat,i,rank
         real(8) :: d,xdat(ndat),phase(ndat),amplitude(ndat), &
                    phaseslope,phaseintercept, &
                    ampslope,ampintercept,oldphase, &
                    newphase
         complex(8) :: edat(ndat),rieff,e0
         call mstm_mpi(mpi_command='rank',mpi_rank=rank)

         oldphase=-1.d10
         do i=1,ndat
            xdat(i)=d*dble(i-1)
            phase(i)=datan2(dimag(edat(i)),dble(edat(i)))
            amplitude(i)=dlog(cdabs(edat(i)))
         enddo
         oldphase=phase(1)
         do i=2,ndat
            newphase=phase(i)
            do while(abs(newphase-oldphase).gt.pi)
               if(newphase.gt.oldphase) then
                  newphase=newphase-2*pi
               else
                  newphase=newphase+2*pi
               endif
            enddo
            phase(i)=newphase
            oldphase=newphase
         enddo
         call linearregression(ndat,phase,xdat, &
              phaseslope,phaseintercept)
         call linearregression(ndat,amplitude,xdat, &
              ampslope,ampintercept)
         rieff=dcmplx(phaseslope,-ampslope)
         e0=dexp(ampintercept)*cdexp((0.d0,1.d0)*phaseintercept)
         end subroutine effectiverefractiveindex

         subroutine linearregression(ndat,fdat,xdat,a,b)
         implicit none
         integer :: ndat
         real(8) :: fdat(ndat),xdat(ndat),a,b,xbar,x2bar, &
                    fbar,xfbar
         fbar=sum(fdat(1:ndat))/dble(ndat)
         xbar=sum(xdat(1:ndat))/dble(ndat)
         xfbar=sum(xdat(1:ndat)*fdat(1:ndat))/dble(ndat)
         x2bar=sum(xdat(1:ndat)*xdat(1:ndat))/dble(ndat)
         a=(xfbar-xbar*fbar)/(x2bar-xbar*xbar)
         b=(fbar*x2bar-xbar*xfbar)/(x2bar-xbar*xbar)
         end subroutine linearregression

         subroutine effective_ref_index_fit(anp,rifit,xfit,info)
         use levenberg_marquardt
         implicit none
         integer :: info,n0,m0,fitorder,k
         real(8) :: parm(3),vec(4*t_matrix_order),xfit,rii,rir,fv
         complex(8) :: anp(2,t_matrix_order),rifit,ri1(2),rim(2),fdif(2)
         data k/0/

         ri1=1.d0
         rim=layer_ref_index(0)
         fitorder=min(80,t_matrix_order)
!         fitorder=2*max_mie_order
         m0=4*fitorder
         if(fit_for_radius) then
           n0=3
         else
           n0=2
         endif
         if(random_configuration_host_model.eq.1) then
            xfit=target_dimensions(1)*length_scale_factor
         elseif(random_configuration_host_model.eq.2) then
            xfit=vol_radius/(sphere_volume_fraction)**0.33333
         endif
         fv=(vol_radius/(target_dimensions(1)*length_scale_factor))**3.
         if(random_configuration_host) then
            rifit=rim(1)
         else
            rii=mean_qext_mie*area_mean_radius**2*3.d0*dble(number_spheres)/(4.d0*xfit**3)/2.d0
            rir=(1.d0-fv)+fv*dble(ref_index_scale_factor)
            rifit=dcmplx(rir,rii/2.d0)
         endif
         parm=(/dble(rifit),dimag(rifit),xfit/)
         call lmdif1 (fcn,m0,n0,parm,vec,1.d-6,info)
         rifit=dcmplx(parm(1),parm(2))
         xfit=parm(3)
         contains
            subroutine fcn(ndat,nparm,xparm,fdat,iflag)
            implicit none
            integer :: nparm,ndat,iflag,n,i
            real(8) :: xparm(nparm),fdat(ndat),xsp,qext,qabs,qsca
            complex(8) :: ri(2),anpmie(2,2,fitorder),a(2),anpold(2,2,fitorder),ao(2),ri0(2)
            ri0(:)=layer_ref_index(0)
            ri(:)=dcmplx(xparm(1),xparm(2))
            if(nparm.eq.3) then
               xsp=xparm(3)
            else
              xsp=xfit
            endif
            if(random_orientation) then
               call mieoa(xsp,ri,fitorder,0.d0,qext,qsca,qabs,anp_mie=anpmie,ri_medium=ri0)
               i=1
               do n=1,fitorder
                  a(1)=(anpmie(1,1,n)+anpmie(2,1,n))
                  a(2)=(anpmie(1,1,n)-anpmie(2,1,n))
                  fdat(i)=dble(anp(1,n)-a(1))
                  fdat(i+1)=dimag(anp(1,n)-a(1))
                  fdat(i+2)=dble(anp(2,n)-a(2))
                  fdat(i+3)=dimag(anp(2,n)-a(2))
                  i=i+4
!k=k+1
!write(*,'(i5,4es12.4)') k,fdat(1:4)
               enddo
            else
               if(effective_medium_simulation) then
                  call mieoa(xsp,ri0,fitorder,0.d0,qext,qsca,qabs,ri_medium=ri,anp_eff_mie=anpmie)
               else
                  call mieoa(xsp,ri,fitorder,0.d0,qext,qsca,qabs,anp_mie=anpmie)
               endif
               i=1
               do n=1,fitorder
                  a(1)=(anpmie(1,1,n)+anpmie(2,1,n))
                  a(2)=(anpmie(1,1,n)-anpmie(2,1,n))
                  fdat(i)=dble(anp(1,n)-a(1))
                  fdat(i+1)=dimag(anp(1,n)-a(1))
                  fdat(i+2)=dble(anp(2,n)-a(2))
                  fdat(i+3)=dimag(anp(2,n)-a(2))
                  i=i+4
               enddo
            endif
            end subroutine fcn
         end subroutine effective_ref_index_fit

         subroutine diffuse_scattering_effective_ref_index(a,qe,extrat)
         implicit none
         real(8) :: a,extrat,tvol,rc,ds,dela,qe,delextrat,extrat0
         call target_volume(target_dimensions,tvol)
         tvol=tvol*length_scale_factor**3
         rc=(tvol*3.d0/4.d0/pi)**(1.d0/3.d0)
         qe=(dif_csca_ratio(1)*q_eff_tot(3,1)+q_eff_tot(2,1))*cross_section_radius**2/rc**2
         ds=(dif_csca_ratio(1)*q_eff_tot(3,1)+q_eff_tot(2,1))*cross_section_radius**2*pi/tvol
!         a=ds/2.d0
!         dela=1.d0
!         do while(abs(dela).gt.1.d-10)
!            dela=(3.*a*(1. + 2.*a*rc)-a*exp(2.*a*rc)*(3. + 2.*a**2*rc**2*(-3. + 2.*ds*rc))) &
!               /(6. - 6.*exp(2*a*rc) + 12.*a*rc*(1. + a*rc))
!            a=a+dela
!         enddo
!         extrat=2.d0*a/((mean_qext_mie)*pi*area_mean_radius**2*dble(number_spheres)/tvol)
! patch 12/2023
!
         a=ds*rc
         dela=1.d0
         delextrat=1.d0
         extrat0=0.d0
         do while(abs(delextrat).gt.1.d-5)
            dela=-(a*(-1.d0 - 2.d0*a + exp(2.d0*a)*(1.d0 + 2.d0*a*a*(-1.d0 + qe)))) &
               /(2.d0 + 4.d0*a*(1.d0 + a) - 2.d0*exp(2.d0*a))
            a=a+dela
            extrat=a/((mean_qext_mie)*pi*area_mean_radius**2*dble(number_spheres)/tvol)/rc
            delextrat=extrat-extrat0
            extrat0=extrat
         enddo
         a=a/2.d0/rc
         end subroutine diffuse_scattering_effective_ref_index

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
         logical :: skipdif,firstrun,frozen
         logical, optional :: skip_diffusion
         integer :: mpicomm,rank,nsend,nsphere,nspheresamp,i
         integer, optional :: mpi_comm
         integer, allocatable,save :: sphereindex(:)
         real(8) :: targetdimensions(3),crad
         real(8), allocatable,save :: sphereradius(:),sphereposition(:,:)
         data firstrun/.true./
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
         nsphere=number_spheres
         if(random_configuration_host) nsphere=nsphere-1
         if(rank.eq.0) then
            frozen=((.not.firstrun).and.(frozen_configuration.or.use_previous_configuration))
            if(frozen) then
               sphere_position(:,1:nsphere)=sphereposition(:,1:nsphere)
               sphere_radius(1:nsphere)=sphereradius(1:nsphere)
            else
               if(auto_target_radius.and.target_shape.eq.2) then
                  targetdimensions=target_dimensions+target_radius_padding
!                  nspheresamp=ceiling(sphere_volume_fraction*(targetdimensions(1)-1.d0)**3.d0)
                  nspheresamp=ceiling(sphere_volume_fraction*(targetdimensions(1))**3.d0)
                  nspheresamp=max(nspheresamp,nsphere)
                  if(mstm_global_rank.eq.0) then
                     write(run_print_unit,'('' set, sampled number spheres:'',2i6)') nsphere,nspheresamp
                  endif
               else
                  targetdimensions=target_dimensions
                  nspheresamp=nsphere
               endif
               if(allocated(sphereradius)) deallocate(sphereradius,sphereposition,sphereindex)
               allocate(sphereradius(nspheresamp),sphereposition(3,nspheresamp),sphereindex(nspheresamp))
               call random_cluster_of_spheres(nspheresamp,targetdimensions,sphereposition,sphereradius, &
                  sphereindex,run_print_unit,ran_config_stat,ran_config_time_steps, &
                  skip_diffusion=skipdif,print_progress=.true.)
               firstrun=.false.
               if(ran_config_stat.ge.3) then
                  write(run_print_unit,'('' unable to generate random configuration'')')
                  stop
               endif
               sphere_position(:,1:nsphere)=sphereposition(:,1:nsphere)
               sphere_radius(1:nsphere)=sphereradius(1:nsphere)
               sphere_index(1:nsphere)=sphereindex(1:nsphere)
            endif
         endif
!         call mstm_mpi(mpi_command='barrier')
         nsend=nsphere
         call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=sphere_radius, &
            mpi_number=nsend,mpi_rank=0,mpi_comm=mpicomm)
         call mstm_mpi(mpi_command='bcast',mpi_send_buf_i=sphere_index, &
            mpi_number=nsend,mpi_rank=0,mpi_comm=mpicomm)
         nsend=nsphere*3
         call mstm_mpi(mpi_command='bcast',mpi_send_buf_dp=sphere_position, &
            mpi_number=nsend,mpi_rank=0,mpi_comm=mpicomm)
         sphere_radius(1:nsphere)=sphere_radius(1:nsphere)*length_scale_factor
         sphere_position(:,1:nsphere)=sphere_position(:,1:nsphere)*length_scale_factor
         do i=1,nsphere
            sphere_ref_index(:,i)=ref_index_scale_factor &
               *component_ref_index(sphere_index(i))
         enddo
         if(erase_sphere_1) sphere_ref_index(:,1)=(1.000001d0,0.d0)
         if(random_configuration_host) then
            crad=maxval(sqrt(sum(sphere_position(:,1:nsphere)**2,1)))+0.01d0
            if(auto_target_radius.and.target_shape.eq.2) then
               sphere_radius(number_spheres) &
                  =(sum(sphere_radius(1:number_spheres-1)**3.)/sphere_volume_fraction)**0.33333
            else
               if(random_configuration_host_model.eq.1) then
                  sphere_radius(number_spheres)=target_dimensions(1)*length_scale_factor
               elseif(random_configuration_host_model.eq.2) then
                  sphere_radius(number_spheres)=(sum(sphere_radius(1:number_spheres-1)**3.)/sphere_volume_fraction)**0.33333
               endif
            endif
            sphere_position(:,number_spheres)=0.d0
            sphere_ref_index(:,number_spheres)=host_sphere_ref_index
         endif
         end subroutine generate_random_configuration

         subroutine scattering_matrix_calculation(amnp,scatmat,mpi_comm)
         implicit none
         logical :: singleorigin,iframe
         integer :: i,sy,sx,mpicomm
         integer, optional :: mpi_comm
         real(8) :: scatmat(scat_mat_mdim,scat_mat_ldim:scat_mat_udim),costheta,phi,csca, &
                    ky,kx,sintheta,ctm
         complex(8) :: amnp(*),ampmat(2,2)
         if(present(mpi_comm)) then
            mpicomm=mpi_comm
         else
            mpicomm=mpi_comm_world
         endif
         singleorigin=number_plane_boundaries.eq.0.and.single_origin_expansion
         iframe=singleorigin.and.incident_frame
         csca=(q_eff_tot(1,1)-q_eff_tot(2,1))*pi*cross_section_radius**2
         csca=pi*2.d0

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
                        if(.not.numerical_azimuthal_average) then
!                        call fosmcalc(12,s00,s02,sp22,sm22,costheta,scatmat(:,i),normalize_s11=.false.)
                           call fosmcalc(t_matrix_order,scat_mat_exp_coef(:,:,1),scat_mat_exp_coef(:,:,2), &
                              scat_mat_exp_coef(:,:,3),scat_mat_exp_coef(:,:,4), &
                              costheta,scatmat(:,i),normalize_s11=.false.)
                        else
                           call numerical_sm_azimuthal_average_so(amnp,t_matrix_order,costheta,scatmat(:,i), &
                              rotate_plane=.true.,normalize_s11=.false.)
                        endif
                     else
                        call scatteringmatrix(amnp,t_matrix_order,costheta,phi,ampmat,scatmat(:,i), &
                           rotate_plane=iframe,normalize_s11=.false.)
                     endif
                  else
                     if(azimuthal_average) then
                        call numerical_sm_azimuthal_average_mo(amnp,costheta,scatmat(:,i),rotate_plane=.true.)
                     else
                        call multiple_origin_scatteringmatrix(amnp,costheta,phi,csca,ampmat,scatmat(:,i), &
                            rotate_plane=.true.)
                     endif
                  endif
               else
                  ctm=-costheta
                  if(azimuthal_average) then
                     call numerical_sm_azimuthal_average_mo(amnp,ctm,scatmat(1:16,i))
                     call numerical_sm_azimuthal_average_mo(amnp,costheta,scatmat(17:32,i))
                  else
                     call multiple_origin_scatteringmatrix(amnp,ctm,phi,csca,ampmat,scatmat(1:16,i))
                     call multiple_origin_scatteringmatrix(amnp,costheta,phi,csca,ampmat,scatmat(17:32,i))
                  endif
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
                     call scatteringmatrix(amnp,t_matrix_order,costheta,phi,ampmat,scatmat(1:16,i), &
                        rotate_plane=iframe,normalize_s11=.false.)
                  else
                     call multiple_origin_scatteringmatrix(amnp,costheta,phi,csca,ampmat,scatmat(1:16,i), &
                        rotate_plane=incident_frame)
                  endif
                  costheta=-costheta
                  if(singleorigin) then
                     call scatteringmatrix(amnp,t_matrix_order,costheta,phi,ampmat,scatmat(17:32,i), &
                        rotate_plane=iframe,normalize_s11=.false.)
                  else
                     call multiple_origin_scatteringmatrix(amnp,costheta,phi,csca,ampmat, &
                        scatmat(17:32,i),rotate_plane=incident_frame)
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
         if(random_configuration) then
            write(iunit,'('' sphere positions randomly generated'')')
            if(target_shape.eq.0) then
               write(iunit,'('' rectangular target, half-widths in x,y,z'')')
               write(iunit,'(3es12.4)') target_dimensions*length_scale_factor
            elseif(target_shape.eq.1) then
               write(iunit,'('' cylindrical target, radius, half-thickness'')')
               write(iunit,'(3es12.4)') target_dimensions(1)*length_scale_factor, &
                  target_dimensions(3)*length_scale_factor
            else
               write(iunit,'('' spherical target, radius'')')
               write(iunit,'(3es12.4)') target_dimensions(1)*length_scale_factor
               if(random_configuration_host) then
                  write(iunit,'('' target enclosed in host sphere w/ radius, ref index'')')
                  write(iunit,'(3es12.4)') sphere_radius(number_spheres), &
                     sphere_ref_index(1,number_spheres)
               endif
            endif
            write(iunit,'('' number of components:'')')
            write(iunit,'(i3)') number_components

            write(iunit,'('' sphere log-normal PSD sigma:'')')
            write(iunit,'(4es12.4)') psd_sigma(1:number_components)

            if(number_components.gt.1) then
               write(iunit,'('' sphere component radii:'')')
               write(iunit,'(4es12.4)') component_radii(1:number_components)
               write(iunit,'('' sphere component number_fraction:'')')
               write(iunit,'(4es12.4)') component_number_fraction(1:number_components)
               write(iunit,'('' sphere component refractive index:'')')
               write(iunit,'(8es12.4)') component_ref_index(1:number_components)
            endif

            write(iunit,'('' sphere volume fraction'')')
            if((.not.number_spheres_specified).or.auto_target_radius) then
               write(iunit,'(es12.4)') sphere_volume_fraction
            else
               call target_volume(target_dimensions,tvol)
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
         write(iunit,'(i7)') number_spheres
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
         write(iunit,'('' max_iterations,solution_epsilon, mie_epsilon, interaction radius'')')
         write(iunit,'(i10,3es12.4)') max_iterations,solution_epsilon,mie_epsilon,interaction_radius
         write(iunit,'('' maximum Mie order, number of equations:'')')
         write(iunit,'(i4,i10)') max_mie_order,number_eqns
         write(iunit,'('' mean sphere Mie extinction, absorption efficiencies, albedo'')')
         write(iunit,'(3es12.4)') mean_qext_mie,mean_qabs_mie,1.d0-mean_qabs_mie/mean_qext_mie
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

         subroutine print_error_codes(outunit)
         implicit none
         integer :: outunit
         if(number_plane_boundaries.eq.0) then
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
         end subroutine print_error_codes

         subroutine print_calculation_results(fout)
         implicit none
         integer :: outunit,n,i,j,smvec(6),sx,sy,s,smvec0(16),nsmat,smvecp(16)
         real(8) :: smt(16),kx,ky,s11scale,r(2),t(2),a(2),scacoef,scarat, &
            abscoef,absrat,rl(2),al(2),tl(2),tvol,imrieff,qeeff
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
         call print_error_codes(outunit)
         if(number_plane_boundaries.gt.0) then
            call boundary_energy_transfer(incident_sin_beta,incident_direction,r,t,a)
         else
            r=0.d0
            a=0.d0
            t=1.d0
         endif
         if(random_orientation) then
            write(outunit,'('' calculated t matrix order:'')')
            write(outunit,'(i6)')  t_matrix_order
         endif

         if(print_sphere_data) then
            write(outunit,'('' sphere extinction, absorption, volume absorption efficiencies (unpolarized incidence)'')')
            if(configuration_average) then
               if(target_shape.le.1) then
                  write(outunit,'(''   sphere z_ave       Qext        Qabs        Qvabs'')')
               else
                  write(outunit,'(''   sphere r_ave       Qext        Qabs        Qvabs'')')
               endif
            else
               write(outunit,'(''   sphere   Qext        Qabs        Qvabs'')')
            endif
            do n=1,number_spheres
               if(configuration_average) then
                  write(outunit,'(i8,4es12.4)') n,sphere_position(3,n),q_eff(1:2,1,n),q_vabs(1,n)
               else
                  write(outunit,'(i8,3es12.4)') n,q_eff(1:2,1,n),q_vabs(1,n)
               endif
            enddo
         endif

         if(periodic_lattice.or.reflection_model) then
            rl=pl_sca(:,2)-boundary_ext(:,0)+r(:)
!            tl=pl_sca(:,1)+boundary_ext(:,number_plane_boundaries+1)+t(:)
!            al=1.d0-rl-tl
            al=surface_absorptance(:)
            tl=1.d0-rl-al
            write(outunit,'('' unit cell reflectance, absorptance, transmittance (unpol, par, perp)'')')
               write(outunit,'(9es12.4)') 0.5d0*sum(rl),0.5d0*sum(al),0.5d0*sum(tl), &
                 (rl(i),al(i),tl(i),i=1,2)
!            write(outunit,'('' unit cell reflectance, absorptance, transmittance (unpol, par, perp)'')')
!               write(outunit,'(9es12.4)') 0.5d0*sum(pl_sca(:,2)-boundary_ext(:,0)+r(:)),q_eff_tot(2,1), &
!                 .5d0*sum(pl_sca(:,1)+boundary_ext(:,number_plane_boundaries+1)+t(:)), &
!                 (pl_sca(i,2)-boundary_ext(i,0)+r(i),q_eff_tot(2,i+1), &
!                  pl_sca(i,1)+boundary_ext(i,number_plane_boundaries+1)+t(i),i=1,2)

!            write(outunit,'('' down, up scattering fraction, total scattering cross section (unpol)'')')
!               write(outunit,'(9es12.4)') 0.5d0*sum(pl_sca(:,2)),.5d0*sum(pl_sca(:,1)), &
!                 pi*(0.5d0*sum(pl_sca(:,2))+.5d0*sum(pl_sca(:,1)))*cross_section_radius**2
            if(configuration_average.and.single_origin_expansion) then
               scacoef=(-0.5d0*sum(dif_boundary_sca(:,0))+0.5d0*sum(dif_boundary_sca(:,1))) &
                    *pi*cross_section_radius**2
               scarat=scacoef/(q_eff_tot(3,1)*pi*cross_section_radius**2)
               write(outunit,'('' down, up diffuse scattering fraction, optically thin dependent/independent ratio'')')
               write(outunit,'(9es12.4)') -0.5d0*sum(dif_boundary_sca(:,0)), &
                  0.5d0*sum(dif_boundary_sca(:,number_plane_boundaries+1)), &
                  scarat
               call effective_extinction_coefficient_ratio(scacoef,abscoef,scarat,absrat)
               write(outunit,'('' dimensionless extinction, absorption coefficients, dependent/independent ratios'')')
               write(outunit,'(4es12.4)') scacoef,abscoef,scarat,absrat
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
                     0.5d0*sum(boundary_ext(:,0)),0.5d0*sum(-boundary_ext(:,1)), &
                     (boundary_ext(i,0),-boundary_ext(i,1),i=1,2)
               endif
               if(calculate_up_down_scattering) then
                  write(outunit,'(''  down and up hemispherical scattering efficiencies (unpol, par, perp) '')')
                  write(outunit,'(16es12.4)')  &
                     0.5d0*sum(-boundary_sca(:,0)),0.5d0*sum(boundary_sca(:,1)), &
                     (-boundary_sca(i,0),boundary_sca(i,1),i=1,2)
               endif
               if(number_plane_boundaries.gt.0.and.number_singular_points.gt.0) then
                  write(outunit,'(''  waveguide scattering efficiencies (unpol, par, perp)  '')')
                  write(outunit,'(16es12.4)')  &
                     0.5d0*sum(evan_sca(:)),(evan_sca(i),i=1,2)
               endif
            endif
         endif
         if(configuration_average.and.single_origin_expansion.and..not.random_orientation) then
            call diffuse_scattering_effective_ref_index(imrieff,qeeff,scarat)
            write(outunit,'('' diffuse/total scattering ratio, extinction coefficient, extinction prob., albedo, RT ratio'')')
            write(outunit,'(6es12.4)') dif_csca_ratio,2.d0*imrieff, qeeff,&
               dif_csca_ratio*q_eff_tot(3,1)/(dif_csca_ratio*q_eff_tot(3,1)+q_eff_tot(2,1)),scarat
            write(outunit,'('' NI formulation ratio'')')
            write(outunit,'(5es12.4)') tot_csca,dif_csca, &
               pi*(dif_csca)/(pi*cross_section_radius**2*q_eff_tot(1,1))
         endif
         if(configuration_average.and.target_shape.eq.2.and.(.not.random_configuration_host)) then
            scarat=mean_qext_mie*area_mean_radius**2*3.d0*dble(number_spheres)/(4.d0*fit_radius**3)
            scarat=2*dimag(effective_ref_index)/scarat
            write(outunit,'('' mie fit effective ri, radius, RT ratio, status'')')
            write(outunit,'(4es12.4,i5)') effective_ref_index,fit_radius,scarat,fit_stat
         endif
         if(configuration_average.and.target_shape.eq.2.and.calculate_near_field) then
            write(outunit,'('' field fit effective ri'')')
            write(outunit,'(2es12.4)') nf_eff_ref_index
         endif

         if(calculate_scattering_matrix.and..not.periodic_lattice) then
            if(normalize_s11) then
               s11scale=1.d0/((cross_section_radius**2)*pi*q_eff_tot(3,1))
            else
!               s11scale=(cross_section_radius**2)*q_eff_tot(3,1)
               s11scale=2.d0*pi
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
                  call print_scat_mat_row(i,smt)
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
                        call print_scat_mat_row(i,smt)
                     enddo
                     if((configuration_average.or.incidence_average).and.single_origin_expansion) then
!write(*,'(es12.4)') s11scale
                        write(outunit,'('' diffuse scattering matrix '')')
                        call print_scat_mat_header(no_numbers=.true.)
                        do i=scat_mat_ldim,scat_mat_udim
                           smt=scaled_scat_mat(dif_scat_mat(1:16,i))
                           smt(1)=smt(1)*s11scale
                           call print_scat_mat_row(i,smt)
                        enddo
                     endif
                  else
                     write(outunit,'('' scattering matrix in incident plane'')')
                     write(outunit,'('' reflection'')')
                     call print_scat_mat_header()
                     do i=scat_mat_ldim,scat_mat_udim
                        smt=scaled_scat_mat(scat_mat(1:16,i))
                        smt(1)=smt(1)*s11scale
                        call print_scat_mat_row(i,smt)
                     enddo
                     write(outunit,'('' transmission'')')
                     call print_scat_mat_header(no_numbers=.true.)
                     do i=scat_mat_ldim,scat_mat_udim
                        smt=scaled_scat_mat(scat_mat(17:32,i))
                        smt(1)=smt(1)*s11scale
                        call print_scat_mat_row(i,smt)
                     enddo
                  endif
               else
                  write(outunit,'('' 2-D scattering in backward and forward hemispheres.   Number points:'')')
                  write(outunit,'(i10)') scat_mat_udim
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
                        smt(1)=smt(1)*s11scale
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
                        smt(1)=smt(1)*s11scale
                        write(outunit,'(2f9.5,$)') kx,ky
                        do j=1,16
                           write(outunit,'(es12.4,$)') smt(smvec0(j))
                        enddo
                        write(outunit,*)
                     enddo
                  enddo
               endif
            endif
            if(azimuthal_average.and.(.not.random_orientation).and.(.not.numerical_azimuthal_average)) then
!               s11scale=1.d0/scat_mat_exp_coef(1,0,1)
!               scat_mat_exp_coef=scat_mat_exp_coef*s11scale
               write(outunit,'('' azimuthal averaged scattering matrix expansion coefficients, total field'')')
               write(outunit,'(''    n         11            44            12            34           22p           22m'')')
               do n=0,2*t_matrix_order
                  write(outunit,'(i5,6es14.6)') n,scat_mat_exp_coef(1,n,1),scat_mat_exp_coef(16,n,1), &
                     0.5d0*(scat_mat_exp_coef(2,n,2)+scat_mat_exp_coef(5,n,2)), &
                     0.5d0*(scat_mat_exp_coef(12,n,2)+scat_mat_exp_coef(15,n,2)), &
                     scat_mat_exp_coef(6,n,3),scat_mat_exp_coef(6,n,4)
               enddo
               if(configuration_average) then
                  write(outunit,'('' azimuthal averaged scattering matrix expansion coefficients, coherent field'')')
                  write(outunit,'(''    n         11            44            12            34           22p           22m'')')
                  do n=0,2*t_matrix_order
                     write(outunit,'(i5,6es14.6)') n,coh_scat_mat_exp_coef(1,n,1),coh_scat_mat_exp_coef(16,n,1), &
                        0.5d0*(coh_scat_mat_exp_coef(2,n,2)+coh_scat_mat_exp_coef(5,n,2)), &
                        0.5d0*(coh_scat_mat_exp_coef(12,n,2)+coh_scat_mat_exp_coef(15,n,2)), &
                        coh_scat_mat_exp_coef(6,n,3),coh_scat_mat_exp_coef(6,n,4)
                  enddo
               endif
            endif
            if(random_orientation) then
               write(outunit,'('' orientation averaged scattering matrix expansion coefficients'')')
!               write(outunit,'(''    w  a11           a22           a33           '',&
!                &''a23           a32           a44           a12         '',&
!                &''a34           a13           a24           a14'')')
!               do n=0,2*t_matrix_order
!                  write(outunit,'(i5,11es14.6)') w,scat_mat_exp_coef(1,1,n),scat_mat_exp_coef(2,2,n),&
!                    scat_mat_exp_coef(3,3,n),scat_mat_exp_coef(2,3,n),scat_mat_exp_coef(3,2,n),scat_mat_exp_coef(4,4,n),&
!                    scat_mat_exp_coef(1,2,n),scat_mat_exp_coef(3,4,n),scat_mat_exp_coef(1,3,n),scat_mat_exp_coef(2,4,n),&
!                    scat_mat_exp_coef(1,4,n)
!               enddo
               write(outunit,'(''    n         11            44            12            34           22            33'')')
               do n=0,2*t_matrix_order
                  write(outunit,'(i5,11es14.6)') n,scat_mat_exp_coef(1,1,n),scat_mat_exp_coef(4,4,n),&
                    scat_mat_exp_coef(1,2,n),scat_mat_exp_coef(3,4,n), &
                    scat_mat_exp_coef(2,2,n),scat_mat_exp_coef(3,3,n)
               enddo
               write(outunit,'('' orientation averaged coherent scattering matrix expansion coefficients'')')
               write(outunit,'(''    n         11            44            12            34           22            33'')')
               do n=0,2*t_matrix_order
                  write(outunit,'(i5,11es14.6)') n,coh_scat_mat_exp_coef(1,1,n),coh_scat_mat_exp_coef(4,4,n),&
                    coh_scat_mat_exp_coef(1,2,n),coh_scat_mat_exp_coef(3,4,n), &
                    coh_scat_mat_exp_coef(2,2,n),coh_scat_mat_exp_coef(3,3,n)
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

         if(random_orientation.and.configuration_average) then
            write(outunit,'('' mean t matrix elements (p=1,2)'')')
            do i=1,t_matrix_order
               write(outunit,'(i5,4es12.4)') i,mean_t(1,i),mean_t(2,i)
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

            subroutine print_scat_mat_row(i,smt)
            implicit none
            integer :: i
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
!            write(run_print_unit,'('' positions:'',i5,3es12.4,i5,3es12.4)') imin,sphere_position(:,imin), &
!               jmin,sphere_position(:,jmin)
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
