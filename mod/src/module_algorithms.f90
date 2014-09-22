module module_algorithms

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_box, only: num_dimensions
use module_types_micro, only: Box_Parameters, Particle_Index, Particle_Energy
use module_geometry, only: geometry
use module_physics_micro, only: random_surface, markov_surface
use class_external_field, only: External_Field
use class_hard_spheres, only: Hard_Spheres, Dipolar_Hard_Spheres
use class_neighbour_cells, only: Neighbour_Cells
use class_hard_spheres_potential, only: Between_Hard_Spheres_Potential_Energy
use module_types_macro, only: Hard_Spheres_Macro, Dipolar_Hard_Spheres_Macro
use class_discrete_observable, only: Discrete_Observables
use class_hard_spheres_observables, only: Hard_Spheres_Monte_Carlo_Observables, &
                                          Dipolar_Hard_Spheres_Monte_Carlo_Observables, &
                                          Hard_Spheres_Post_Processing_Observables

implicit none
private
public move, widom, switch, rotate

contains

    !> Particle move
    
    subroutine move(Box, &
                    this_spheres, this_macro, this_observables, &
                    other_spheres, other_between_cells, &
                    between_spheres_potential, mix_potential_energy)
    
        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this_spheres, other_spheres
        class(Hard_Spheres_Macro), intent(inout) :: this_macro
        class(Neighbour_Cells), intent(inout) :: other_between_cells
        class(Hard_Spheres_Monte_Carlo_Observables), intent(inout) :: this_observables
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: between_spheres_potential
        real(DP), intent(inout) :: mix_potential_energy
        
        real(DP) :: random
        real(DP), dimension(num_dimensions) :: random_vector
        type(Particle_Index) :: old, new
        logical :: overlap
        real(DP) :: energy_delta
        real(DP) :: this_energy_delta, mix_energy_delta
        real(DP) :: this_energy_new, this_energy_old
        real(DP) :: mix_energy_new, mix_energy_old
        
        real(DP) :: this_energy_real_new, this_energy_real_old
        
        this_observables%move%num_hits = this_observables%move%num_hits + 1
        
        call random_number(random)
        old%number = int(random*this_spheres%get_num_particles()) + 1
        old%position(:) = this_spheres%get_position(old%number)
        
        new%number = old%number
        call random_number(random_vector)
        new%position(:) = old%position(:) + (random_vector(:)-0.5_DP) * this_macro%move%get_delta()

        if (geometry%bulk) then
            new%position(:) = modulo(new%position(:), Box%size(:))
        else if (geometry%slab) then
            if (new%position(3) < this_spheres%get_diameter()/2._DP .or. &
                new%position(3) > Box%height-this_spheres%get_diameter()/2._DP) then
                this_observables%move%num_rejections = this_observables%move%num_rejections + 1
                return
            end if
            new%position(1:2) = modulo(new%position(1:2), Box%size(1:2))
        end if
        
        if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
            new%same_i_cell = this_macro%same_cells%index_from_position(new%position)
            call this_macro%hard_potential%neighbours(Box%size, this_spheres, this_macro%same_cells, &
                                                      new, overlap, this_energy_new)
        else
            new%between_i_cell = other_between_cells%index_from_position(new%position)
            call between_spheres_potential%neighbours(Box%size, other_spheres, this_macro%between_cells, &
                                                      new, overlap, mix_energy_new)
        end if
        
        if (.not. overlap) then
        
            if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                new%between_i_cell = other_between_cells%index_from_position(new%position)
                call between_spheres_potential%neighbours(Box%size, other_spheres, &
                                                          this_macro%between_cells, &
                                                          new, overlap, mix_energy_new)
            else
                new%same_i_cell = this_macro%same_cells%index_from_position(new%position)
                call this_macro%hard_potential%neighbours(Box%size, this_spheres, &
                                                          this_macro%same_cells, new, overlap, &
                                                          this_energy_new)
            end if
                        
            if (.not. overlap) then
    
                old%same_i_cell = this_macro%same_cells%index_from_position(old%position)
                select type (this_spheres)
                    type is (Dipolar_Hard_Spheres)
                        old%orientation(:) = this_spheres%get_orientation(old%number)
                        new%orientation(:) = old%orientation(:)
                        select type (this_macro)
                            type is (Dipolar_Hard_Spheres_Macro)
                                this_energy_real_new = this_macro%ewald_real%solo_energy(Box%size, &
                                                           this_spheres, new)
                                this_energy_real_old = this_macro%ewald_real%solo_energy(Box%size, &
                                                           this_spheres, old)
                                this_energy_delta = (this_energy_real_new - this_energy_real_old) + &
                                                    this_macro%ewald_reci%move_energy(Box, old, new)
                                if (geometry%slab) then
                                    this_energy_delta = this_energy_delta - &
                                                        this_macro%elc%move_energy(Box, old, new)
                                end if
                        end select
                    class default
                        call this_macro%hard_potential%neighbours(Box%size, this_spheres, &
                                                                  this_macro%same_cells, &
                                                                  old, overlap, this_energy_old)
                        this_energy_delta = this_energy_new - this_energy_old
                end select
                    
                old%between_i_cell = other_between_cells%index_from_position(old%position)
                call between_spheres_potential%neighbours(Box%size, other_spheres, &
                                                          this_macro%between_cells, &
                                                          old, overlap, mix_energy_old)
                
                mix_energy_delta = mix_energy_new - mix_energy_old

                energy_delta = this_energy_delta + mix_energy_delta
                
                call random_number(random)
                if (random < exp(-energy_delta/Box%temperature)) then
                
                    select type (this_macro)
                        type is (Dipolar_Hard_Spheres_Macro)
                            call this_macro%ewald_reci%update_structure_move(Box, old, new)
                            if (geometry%slab) then
                                call this_macro%elc%update_structure_move(Box, old, new)
                            end if
                    end select
                
                    call this_spheres%set_position(old%number, new%position)
                    this_observables%potential_energy = this_observables%potential_energy + &
                                                        this_energy_delta
                    mix_potential_energy = mix_potential_energy + mix_energy_delta
                    
                    if (old%same_i_cell /= new%same_i_cell) then
                        call this_macro%same_cells%remove_particle_from_cell(old%number, old%same_i_cell)
                        call this_macro%same_cells%add_particle_to_cell(new%number, new%same_i_cell)
                    end if
                    if (old%between_i_cell /= new%between_i_cell) then
                        call other_between_cells%remove_particle_from_cell(old%number, old%between_i_cell)
                        call other_between_cells%add_particle_to_cell(new%number, new%between_i_cell)
                    end if
                    
                else
                    this_observables%move%num_rejections = this_observables%move%num_rejections + 1
                end if
         
            else
                this_observables%move%num_rejections = this_observables%move%num_rejections + 1
            end if
            
        else
            this_observables%move%num_rejections = this_observables%move%num_rejections + 1
        end if
    
    end subroutine move
    
    !> Widom's method

    subroutine widom(Box, ext_field, &
                     this_spheres, this_macro, this_observables, &
                     other_spheres, other_between_cells, &
                     between_spheres_potential)
        
        type(Box_Parameters), intent(in) :: Box
        class(External_Field), intent(in) :: ext_field
        class(Hard_Spheres), intent(in) :: this_spheres
        class(Hard_Spheres_Macro), intent(in) :: this_macro
        class(Neighbour_Cells), intent(in) ::  other_between_cells
        class(Hard_Spheres_Post_Processing_Observables), intent(inout) :: this_observables
        class(Hard_Spheres), intent(in) :: other_spheres
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: between_spheres_potential
        
        integer :: i_widom_particule
        real(DP) :: inv_activity_sum
        real(DP), dimension(num_dimensions) :: random_vector
        type(Particle_Index) :: test
        logical :: overlap
        real(DP) :: energy_test
        real(DP) :: this_energy_test, mix_energy_test
        
        inv_activity_sum = 0._DP
        test%number = 0
        
        do i_widom_particule = 1, this_spheres%get_widom_num_particles()
            
            call random_number(random_vector)
            if (geometry%bulk) then
                test%position(:) = Box%size(:) * random_vector(:)
            else if (geometry%slab) then
                test%position(1:2) = Box%size(1:2) * random_vector(1:2)
                test%position(3) = (Box%height - this_spheres%get_diameter()) * random_vector(3) + &
                                   this_spheres%get_diameter()/2._DP
            end if

            if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                test%same_i_cell = this_macro%same_cells%index_from_position(test%position)
                call this_macro%hard_potential%neighbours(Box%size, this_spheres, &
                                                          this_macro%same_cells, test, &
                                                          overlap, this_energy_test)
            else
                test%between_i_cell = other_between_cells%index_from_position(test%position)
                call between_spheres_potential%neighbours(Box%size, other_spheres, &
                                                          this_macro%between_cells, &
                                                          test, overlap, mix_energy_test)
            end if
            
            if (.not. overlap) then
            
                if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                    test%between_i_cell = other_between_cells%index_from_position(test%position)
                    call between_spheres_potential%neighbours(Box%size, other_spheres, &
                                                              this_macro%between_cells, test, overlap, &
                                                              mix_energy_test)
                else
                    test%same_i_cell = this_macro%same_cells%index_from_position(test%position)
                    call this_macro%hard_potential%neighbours(Box%size, this_spheres, &
                                                              this_macro%same_cells, test, &
                                                              overlap, this_energy_test)
                end if
                
                if (.not. overlap) then
                
                    select type (this_spheres)
                        type is (Dipolar_Hard_Spheres)
                            test%add = .true.
                            test%orientation(:) = random_surface()
                            select type (this_macro)
                                type is (Dipolar_Hard_Spheres_Macro)
                                    this_energy_test = &
                                        this_macro%ewald_real%solo_energy(Box%size, &
                                                                          this_spheres, test) + &
                                        this_macro%ewald_reci%exchange_energy(Box, test) - &
                                        this_macro%ewald_self%solo_energy(test%orientation) + &
                                        this_macro%ewald_bound%exchange_energy(Box%size, test) + &
                                        ext_field%exchange_energy(test)
                                       
                                    if (geometry%slab) then
                                        this_energy_test = this_energy_test - &
                                                           this_macro%elc%exchange_energy(Box, test)
                                    end if
                                                    
                            end select
                    end select
                
                    energy_test = this_energy_test + mix_energy_test
                    inv_activity_sum = inv_activity_sum + exp(-energy_test/Box%temperature)
                    
                end if
                
            end if
            
        end do
        
        this_observables%inv_activity = inv_activity_sum/real(this_spheres%get_widom_num_particles(), DP)
        
    end subroutine widom
    
    !> Particle switch
    
    subroutine switch(Box, &
                      type1_spheres, type1_macro, type1_observables, &
                      type2_spheres, type2_macro, type2_observables, &
                      between_spheres_potential, mix_potential_energy, &
                      switch_observable)
    
        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: type1_spheres, type2_spheres
        class(Hard_Spheres_Macro), intent(inout) :: type1_macro, type2_macro
        class(Hard_Spheres_Monte_Carlo_Observables), intent(inout) :: type1_observables, type2_observables
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: between_spheres_potential
        real(DP), intent(inout) :: mix_potential_energy
        type(Discrete_Observables), intent(inout) :: switch_observable
        
        real(DP) :: random
        type(Particle_Index) :: old1, old2
        type(Particle_Index) :: new1, new2
        logical :: overlap
        real(DP) :: energy_delta, type1_energy_delta, type2_energy_delta
        real(DP) :: type1_mix_energy_delta, type2_mix_energy_delta
        type(Particle_Energy) :: type1_energy_old, type1_energy_new
        type(Particle_Energy) :: type2_energy_old, type2_energy_new
        
        switch_observable%num_hits = switch_observable%num_hits + 1
        
        if (type1_spheres%get_num_particles()==0 .or. type2_spheres%get_num_particles()==0) then
            switch_observable%num_rejections = switch_observable%num_rejections + 1
            return
        end if
        
        call random_number(random)
        old1%number = int(random*type1_spheres%get_num_particles()) + 1
        new1%number = old1%number
        call random_number(random)
        old2%number = int(random*type2_spheres%get_num_particles()) + 1
        new2%number = old2%number
        
        old1%other_number = old2%number; new1%other_number = new2%number
        old2%other_number = old1%number; new2%other_number = new1%number
                
        call before_switch_energy(Box%size, &
                                  type1_spheres, type1_macro, old1, &
                                  type2_spheres, type2_macro%between_cells, &
                                  between_spheres_potential, &
                                  type1_energy_old)
        call before_switch_energy(Box%size, &
                                  type2_spheres, type2_macro, old2, &
                                  type1_spheres, type1_macro%between_cells, &
                                  between_spheres_potential, &
                                  type2_energy_old)
             
        call after_switch_energy(Box, &
                                 type1_spheres, type1_macro, old1, new1, &
                                 type2_spheres, type2_macro%between_cells, &
                                 between_spheres_potential, &
                                 overlap, &
                                 type1_energy_new)
        
        if (.not. overlap) then
        
            call after_switch_energy(Box, &
                                     type2_spheres, type2_macro, old2, new2, &
                                     type1_spheres, type1_macro%between_cells, &
                                     between_spheres_potential, &
                                     overlap, &
                                     type2_energy_new)
            
            if (.not. overlap) then

                type1_energy_delta = type1_energy_new%same - type1_energy_old%same
                type1_mix_energy_delta = type1_energy_new%mix - type1_energy_old%mix
                type2_energy_delta = type2_energy_new%same - type2_energy_old%same
                type2_mix_energy_delta = type2_energy_new%mix - type2_energy_old%mix
                energy_delta = type1_energy_delta + type1_mix_energy_delta + type2_energy_delta + &
                               type2_mix_energy_delta
                
                call random_number(random)
                if (random < exp(-energy_delta/Box%temperature)) then
                
                    call after_switch_update(Box, &
                                             type1_spheres, type1_macro, old1, new1, &
                                             type2_macro%between_cells)
                    call after_switch_update(Box, &
                                             type2_spheres, type2_macro, old2, new2, &
                                             type1_macro%between_cells)
                                             
                    type1_observables%potential_energy = type1_observables%potential_energy + &
                                                         type1_energy_delta
                    type2_observables%potential_energy = type2_observables%potential_energy + &
                                                         type2_energy_delta
                    mix_potential_energy = mix_potential_energy + type1_mix_energy_delta + &
                                                                  type2_mix_energy_delta
                    
                else
                    switch_observable%num_rejections = switch_observable%num_rejections + 1
                end if
                
            else
                switch_observable%num_rejections = switch_observable%num_rejections + 1
            end if
            
        else
            switch_observable%num_rejections = switch_observable%num_rejections + 1
        end if
        
    end subroutine switch
    
    subroutine before_switch_energy(Box_size, &
                                    this_spheres, this_macro, old, &
                                    other_spheres, other_between_cells, &
                                    between_spheres_potential, &
                                    energy_old)
        
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Hard_Spheres_Macro), intent(in) :: this_macro
        class(Neighbour_Cells), intent(in) :: other_between_cells
        type(Particle_Index), intent(inout) :: old
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: between_spheres_potential
        type(Particle_Energy), intent(out) :: energy_old
        logical :: overlap
        
        old%position(:) = this_spheres%get_position(old%number)
        
        old%same_i_cell = this_macro%same_cells%index_from_position(old%position)
        select type (this_spheres)
            type is (Dipolar_Hard_Spheres)
                old%orientation(:) = this_spheres%get_orientation(old%number)
                select type (this_macro)
                    type is (Dipolar_Hard_Spheres_Macro)
                        energy_old%same = this_macro%ewald_real%solo_energy(Box_size, this_spheres, old)
                               ! potential_energy_reci: cf. after_switch_energy
                end select
            type is (Hard_Spheres)
                energy_old%same = 0._DP
        end select
        
        old%between_i_cell = other_between_cells%index_from_position(old%position)
        call between_spheres_potential%neighbours(Box_size, other_spheres, this_macro%between_cells, &
                                                  old, overlap, energy_old%mix)
        
    end subroutine before_switch_energy
    
    subroutine after_switch_energy(Box, &
                                   this_spheres, this_macro, old, new, &
                                   other_spheres, other_between_cells, &
                                   between_spheres_potential, &
                                   overlap, &
                                   energy_new)

        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Hard_Spheres_Macro), intent(in) :: this_macro
        class(Neighbour_Cells), intent(in) :: other_between_cells
        type(Particle_Index), intent(in) :: old
        type(Particle_Index), intent(inout) :: new
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: between_spheres_potential
        logical, intent(out) :: overlap
        type(Particle_Energy), intent(out) :: energy_new
        
        new%position(:) = other_spheres%get_position(new%other_number)
        
        if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
            new%same_i_cell = this_macro%same_cells%index_from_position(new%position)
            call this_macro%hard_potential%neighbours(Box%size, this_spheres, this_macro%same_cells, &
                                                      new, overlap, energy_new%same)
        else
            new%between_i_cell = other_between_cells%index_from_position(new%position)
            call between_spheres_potential%neighbours(Box%size, other_spheres, this_macro%between_cells, &
                                                      new, overlap, energy_new%mix)
        end if
        
        if (.not. overlap) then
        
            if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                new%between_i_cell = other_between_cells%index_from_position(new%position)
                call between_spheres_potential%neighbours(Box%size, other_spheres, &
                                                          this_macro%between_cells, &
                                                          new, overlap, energy_new%mix)
            else
                new%same_i_cell = this_macro%same_cells%index_from_position(new%position)
                call this_macro%hard_potential%neighbours(Box%size, this_spheres, &
                                                          this_macro%same_cells, new, &
                                                          overlap, energy_new%same)
            end if
            
            if (.not. overlap) then
            
                select type (this_spheres)
                    type is (Dipolar_Hard_Spheres)
                        new%orientation(:) = this_spheres%get_orientation(new%number)
                        select type (this_macro)
                            type is (Dipolar_Hard_Spheres_Macro)
                                energy_new%same = this_macro%ewald_real%solo_energy(Box%size, &
                                                      this_spheres, new) + &
                                                  this_macro%ewald_reci%move_energy(Box, old, new)
                                if (geometry%slab) then
                                    energy_new%same = energy_new%same - &
                                                      this_macro%elc%move_energy(Box, old, new)
                                end if
                        end select
                                       
                end select
            
            end if
                
        end if
    
    end subroutine after_switch_energy
    
    subroutine after_switch_update(Box, &
                                   this_spheres, this_macro, old, new, &
                                   other_between_cells)

        type(Box_Parameters), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this_spheres
        class(Hard_Spheres_Macro), intent(inout) :: this_macro
        class(Neighbour_Cells), intent(inout) :: other_between_cells
        type(Particle_Index), intent(in) :: old, new
        
        call this_spheres%set_position(old%number, new%position)
        
        select type (this_macro)
            type is (Dipolar_Hard_Spheres_Macro)
                call this_macro%ewald_reci%update_structure_move(Box, old, new)
                if (geometry%slab) then
                    call this_macro%elc%update_structure_move(Box, old, new)
                end if
        end select
        
        if (old%same_i_cell /= new%same_i_cell) then
            call this_macro%same_cells%remove_particle_from_cell(old%number, old%same_i_cell)
            call this_macro%same_cells%add_particle_to_cell(new%number, new%same_i_cell)
        end if
        if (old%between_i_cell /= new%between_i_cell) then
            call other_between_cells%remove_particle_from_cell(old%number, old%between_i_cell)
            call other_between_cells%add_particle_to_cell(new%number, new%between_i_cell)
        end if
        
    end subroutine after_switch_update
    
    !> Dipole rotation
    
    subroutine rotate(Box, ext_field, &
                      this_spheres, this_macro, this_observables)
    
        type(Box_Parameters), intent(in) :: Box
        class(External_Field), intent(in) :: ext_field
        class(Dipolar_Hard_Spheres), intent(inout) :: this_spheres
        class(Dipolar_Hard_Spheres_Macro), intent(inout) :: this_macro
        class(Dipolar_Hard_Spheres_Monte_Carlo_Observables), intent(inout) :: this_observables
        
        real(DP) :: random
        type(Particle_Index) :: old, new
        real(DP) :: energy_delta
        real(DP) :: energy_real_delta, energy_self_delta
        
        this_observables%rotation%num_hits = this_observables%rotation%num_hits + 1

        call random_number(random)
        old%number = int(random*this_spheres%get_num_particles()) + 1
        old%position(:) = this_spheres%get_position(old%number)
        old%orientation(:) = this_spheres%get_orientation(old%number)
        
        new%number = old%number
        new%position(:) = old%position(:)
        new%orientation(:) = old%orientation(:)
        call markov_surface(new%orientation, this_macro%rotation%get_delta())
        
        energy_real_delta = this_macro%ewald_real%solo_energy(Box%size, this_spheres, new) - &
                            this_macro%ewald_real%solo_energy(Box%size, this_spheres, old)
        
        energy_self_delta = this_macro%ewald_self%solo_energy(new%orientation) - &
                            this_macro%ewald_self%solo_energy(old%orientation)
        
        energy_delta = energy_real_delta + this_macro%ewald_reci%rotation_energy(Box, old, new) - &
                       energy_self_delta + this_macro%ewald_bound%rotation_energy(Box%size, old, new) + &
                       ext_field%rotation_energy(old, new)
        if (geometry%slab) then
            energy_delta = energy_delta - this_macro%elc%rotation_energy(Box, old, new)
        end if
        
        call random_number(random)
        if (random < exp(-energy_delta/Box%temperature)) then
        
            call this_macro%ewald_reci%update_structure_rotation(Box, old, new)
            call this_macro%ewald_bound%update_total_moment_rotation(old, new)
            call this_spheres%set_orientation(old%number, new%orientation)
            if (geometry%slab) then
                call this_macro%elc%update_structure_rotation(Box, old, new)
            end if
            
            this_observables%potential_energy = this_observables%potential_energy + energy_delta
            
        else
            this_observables%rotation%num_rejections = this_observables%rotation%num_rejections + 1
        end if
    
    end subroutine rotate

end module module_algorithms
