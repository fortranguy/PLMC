module module_algorithms

use data_precisions, only: DP
use data_box, only: Ndim
use module_types_micro, only: Box_Dimensions, Particle_Index, Particle_Energy
use module_physics_micro, only: random_surface, markov_surface
use class_neighbour_cells, only: Neighbour_Cells
use class_hard_spheres, only: Hard_Spheres, Dipolar_Hard_Spheres
use class_small_move, only: Small_Move
use class_small_rotation, only: Small_Rotation
use module_types_macro, only: Hard_Spheres_Macro, Dipolar_Hard_Spheres_Macro
use class_between_spheres_potential, only: Between_Hard_Spheres_Potential_Energy
use class_hard_spheres_observables, only: Hard_Spheres_Observables, Dipolar_Hard_Spheres_Observables

implicit none
private
public move, widom, switch, rotate

contains

    !> Particle move
    
    subroutine move(Box, &
                    this_spheres, this_macro, this_observables, &
                    other_spheres, other_mix_cells, &
                    mix, mix_potential_energy)
    
        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this_spheres, other_spheres
        class(Hard_Spheres_Macro), intent(inout) :: this_macro
        class(Neighbour_Cells), intent(inout) :: other_mix_cells
        class(Hard_Spheres_Observables), intent(inout) :: this_observables
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: mix
        real(DP), intent(inout) :: mix_potential_energy
        
        real(DP) :: random
        real(DP), dimension(Ndim) :: xRand
        type(Particle_Index) :: old, new
        logical :: overlap
        real(DP) :: deltaEpot
        real(DP) :: this_deltaEpot, mix_deltaEpot
        real(DP) :: this_EpotNew, this_EpotOld
        real(DP) :: mix_EpotNew, mix_EpotOld
        
        real(DP) :: this_EpotNew_real, this_EpotOld_real
        
        this_observables%move_num_hits = this_observables%move_num_hits + 1
        
        call random_number(random)
        old%number = int(random*this_spheres%get_num_particles()) + 1
        old%position(:) = this_spheres%get_position(old%number)
        
        new%number = old%number
        call random_number(xRand)
        new%position(:) = old%position(:) + (xRand(:)-0.5_DP) * this_macro%move%get_delta()
        new%position(:) = modulo(new%position(:), Box%size(:))
        
        if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
            new%same_iCell = this_macro%same_cells%index_from_position(new%position)
            call this_macro%hard_potential%neighCells(Box%size, this_spheres, this_macro%same_cells, &
                                                      new, overlap, this_EpotNew)
        else
            new%mix_iCell = other_mix_cells%index_from_position(new%position)
            call mix%neighCells(Box%size, new, this_macro%mix_cells, other_spheres, overlap, &
                                     mix_EpotNew)
        end if
        
        if (.not. overlap) then
        
            if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                new%mix_iCell = other_mix_cells%index_from_position(new%position)
                call mix%neighCells(Box%size, new, this_macro%mix_cells, other_spheres, overlap, &
                                         mix_EpotNew)
            else
                new%same_iCell = this_macro%same_cells%index_from_position(new%position)
                call this_macro%hard_potential%neighCells(Box%size, this_spheres, &
                                                          this_macro%same_cells, new, overlap, &
                                                          this_EpotNew)
            end if
                        
            if (.not. overlap) then
    
                old%same_iCell = this_macro%same_cells%index_from_position(old%position)
                select type (this_spheres)
                    type is (Dipolar_Hard_Spheres)
                        old%orientation(:) = this_spheres%get_orientation(old%number)
                        new%orientation(:) = old%orientation(:)
                        select type (this_macro)
                            type is (Dipolar_Hard_Spheres_Macro)
                                this_EpotNew_real = this_macro%ewald_real%solo(Box%size, this_spheres, new)
                                this_EpotOld_real = this_macro%ewald_real%solo(Box%size, this_spheres, old)
                                this_deltaEpot = (this_EpotNew_real - this_EpotOld_real) + &
                                                  this_macro%ewald_reci%move(Box, old, new)
                        end select
                    class default
                        call this_macro%hard_potential%neighCells(Box%size, this_spheres, this_macro%same_cells, &
                                                                  old, overlap, this_EpotOld)
                        this_deltaEpot = this_EpotNew - this_EpotOld
                end select
                    
                old%mix_iCell = other_mix_cells%index_from_position(old%position)
                call mix%neighCells(Box%size, old, this_macro%mix_cells, other_spheres, overlap, &
                                         mix_EpotOld)
                
                mix_deltaEpot = mix_EpotNew - mix_EpotOld

                deltaEpot = this_deltaEpot + mix_deltaEpot
                
                call random_number(random)
                if (random < exp(-deltaEpot/Box%temperature)) then
                
                    select type (this_macro)
                        type is (Dipolar_Hard_Spheres_Macro)
                            call this_macro%ewald_reci%update_structure_move(Box, old, new)
                    end select
                
                    call this_spheres%set_position(old%number, new%position)
                    this_observables%potential_energy = this_observables%potential_energy + this_deltaEpot
                    mix_potential_energy = mix_potential_energy + mix_deltaEpot
                    
                    if (old%same_iCell /= new%same_iCell) then
                        call this_macro%same_cells%remove_col_from_cell(old%number, old%same_iCell)
                        call this_macro%same_cells%add_col_to_cell(new%number, new%same_iCell)
                    end if
                    if (old%mix_iCell /= new%mix_iCell) then
                        call other_mix_cells%remove_col_from_cell(old%number, old%mix_iCell)
                        call other_mix_cells%add_col_to_cell(new%number, new%mix_iCell)
                    end if
                    
                else
                    this_observables%move_num_rejections = this_observables%move_num_rejections + 1
                end if
         
            else
                this_observables%move_num_rejections = this_observables%move_num_rejections + 1
            end if
            
        else
            this_observables%move_num_rejections = this_observables%move_num_rejections + 1
        end if
    
    end subroutine move
    
    !> Widom's method

    subroutine widom(Box, &
                     this_spheres, this_macro, this_observables, &
                     other_spheres, other_mix_cells, &
                     mix)
        
        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(in) :: this_spheres
        class(Hard_Spheres_Macro), intent(in) :: this_macro
        class(Neighbour_Cells), intent(in) ::  other_mix_cells
        class(Hard_Spheres_Observables), intent(inout) :: this_observables
        class(Hard_Spheres), intent(in) :: other_spheres
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: mix
        
        integer :: iWidom
        real(DP) :: widTestSum
        real(DP), dimension(Ndim) :: xRand
        type(Particle_Index) :: test
        logical :: overlap
        real(DP) :: EpotTest
        real(DP) :: this_EpotTest, mix_EpotTest
        
        widTestSum = 0._DP
        test%number = 0
        
        do iWidom = 1, this_spheres%get_widom_num_particles()
            
            call random_number(xRand)
            test%position(:) = Box%size(:) * xRand(:)

            if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                test%same_iCell = this_macro%same_cells%index_from_position(test%position)
                call this_macro%hard_potential%neighCells(Box%size, this_spheres, this_macro%same_cells, test, &
                                                          overlap, this_EpotTest)
            else
                test%mix_iCell = other_mix_cells%index_from_position(test%position)
                call mix%neighCells(Box%size, test, this_macro%mix_cells, other_spheres, overlap, &
                                         mix_EpotTest)
            end if
            
            if (.not. overlap) then
            
                if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                    test%mix_iCell = other_mix_cells%index_from_position(test%position)
                    call mix%neighCells(Box%size, test, this_macro%mix_cells, other_spheres, overlap, &
                                             mix_EpotTest)
                else
                    test%same_iCell = this_macro%same_cells%index_from_position(test%position)
                    call this_macro%hard_potential%neighCells(Box%size, this_spheres, this_macro%same_cells, test, &
                                                              overlap, this_EpotTest)
                end if
                
                if (.not. overlap) then
                
                    select type (this_spheres)
                        type is (Dipolar_Hard_Spheres)
                            test%add = .true.
                            test%orientation(:) = random_surface()
                            select type (this_macro)
                                type is (Dipolar_Hard_Spheres_Macro)                                
                                    this_EpotTest = this_macro%ewald_real%solo(Box%size, this_spheres, &
                                                                               test) + &
                                                    this_macro%ewald_reci%exchange(Box, test) - &
                                                    this_macro%ewald_self%solo(test%orientation) + &
                                                    this_macro%ewald_bound%exchange(Box%size, &
                                                                                  test%orientation)
                                                    
                            end select                                            
                    end select
                
                    EpotTest = this_EpotTest + mix_EpotTest
                    widTestSum = widTestSum + exp(-EpotTest/Box%temperature)
                    
                end if
                
            end if
            
        end do
        
        this_observables%inv_activity = widTestSum/real(this_spheres%get_widom_num_particles(), DP)
        
    end subroutine widom
    
    !> Particle switch
    
    subroutine switch(Box, &
                      type1_spheres, type1_macro, type1_observables, &
                      type2_spheres, type2_macro, type2_observables, &
                      mix, mix_potential_energy, &
                      switch_num_rejections)
    
        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: type1_spheres, type2_spheres
        class(Hard_Spheres_Macro), intent(inout) :: type1_macro, type2_macro
        class(Hard_Spheres_Observables), intent(inout) :: type1_observables, type2_observables
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: mix
        real(DP), intent(inout) :: mix_potential_energy
        integer, intent(inout) :: switch_num_rejections
        
        real(DP) :: random
        type(Particle_Index) :: old1, old2
        type(Particle_Index) :: new1, new2
        logical :: overlap
        real(DP) :: deltaEpot, type1_deltaEpot, type2_deltaEpot
        real(DP) :: type1_mix_deltaEpot, type2_mix_deltaEpot
        type(Particle_Energy) :: type1_EpotOld, type1_EpotNew
        type(Particle_Energy) :: type2_EpotOld, type2_EpotNew
        
        if (type1_spheres%get_num_particles()==0 .or. type2_spheres%get_num_particles()==0) then
            switch_num_rejections = switch_num_rejections + 1
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
                                  type2_spheres, type2_macro%mix_cells, &
                                  mix, &
                                  type1_EpotOld)
        call before_switch_energy(Box%size, &
                                  type2_spheres, type2_macro, old2, &
                                  type1_spheres, type1_macro%mix_cells, &
                                  mix, &
                                  type2_EpotOld)
             
        call after_switch_energy(Box, &
                                 type1_spheres, type1_macro, old1, new1, &
                                 type2_spheres, type2_macro%mix_cells, &
                                 mix, &
                                 overlap, &
                                 type1_EpotNew)
        
        if (.not. overlap) then
        
            call after_switch_energy(Box, &
                                     type2_spheres, type2_macro, old2, new2, &
                                     type1_spheres, type1_macro%mix_cells, &
                                     mix, &
                                     overlap, &
                                     type2_EpotNew)
            
            if (.not. overlap) then

                type1_deltaEpot = type1_EpotNew%same - type1_EpotOld%same
                type1_mix_deltaEpot = type1_EpotNew%mix - type1_EpotOld%mix
                type2_deltaEpot = type2_EpotNew%same - type2_EpotOld%same
                type2_mix_deltaEpot = type2_EpotNew%mix - type2_EpotOld%mix
                deltaEpot = type1_deltaEpot + type1_mix_deltaEpot + type2_deltaEpot + &
                            type2_mix_deltaEpot
                
                call random_number(random)
                if (random < exp(-deltaEpot/Box%temperature)) then
                
                    call after_switch_update(Box, &
                                             type1_spheres, type1_macro, old1, new1, &
                                             type2_macro%mix_cells)
                    call after_switch_update(Box, &
                                             type2_spheres, type2_macro, old2, new2, &
                                             type1_macro%mix_cells)
                                             
                    type1_observables%potential_energy = type1_observables%potential_energy + type1_deltaEpot
                    type2_observables%potential_energy = type2_observables%potential_energy + type2_deltaEpot
                    mix_potential_energy = mix_potential_energy + type1_mix_deltaEpot + type2_mix_deltaEpot
                    
                else
                    switch_num_rejections = switch_num_rejections + 1
                end if
                
            else
                switch_num_rejections = switch_num_rejections + 1
            end if
            
        else
            switch_num_rejections = switch_num_rejections + 1
        end if
        
    end subroutine switch
    
    subroutine before_switch_energy(Box_size, &
                                    this_spheres, this_macro, old, &
                                    other_spheres, other_mix_cells, &
                                    mix, &
                                    EpotOld)
        
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Hard_Spheres_Macro), intent(in) :: this_macro
        class(Neighbour_Cells), intent(in) :: other_mix_cells
        type(Particle_Index), intent(inout) :: old
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: mix
        type(Particle_Energy), intent(out) :: EpotOld
        logical :: overlap
        
        old%position(:) = this_spheres%get_position(old%number)
        
        old%same_iCell = this_macro%same_cells%index_from_position(old%position)
        select type (this_spheres)
            type is (Dipolar_Hard_Spheres)
                old%orientation(:) = this_spheres%get_orientation(old%number)
                select type (this_macro)
                    type is (Dipolar_Hard_Spheres_Macro)
                        EpotOld%same = this_macro%ewald_real%solo(Box_size, this_spheres, old)
                               ! potential_energy_reci: cf. after_switch_energy
                end select
            type is (Hard_Spheres)
                EpotOld%same = 0._DP
        end select
        
        old%mix_iCell = other_mix_cells%index_from_position(old%position)
        call mix%neighCells(Box_size, old, this_macro%mix_cells, other_spheres, overlap, EpotOld%mix)
        
    end subroutine before_switch_energy
    
    subroutine after_switch_energy(Box, &
                                   this_spheres, this_macro, old, new, &
                                   other_spheres, other_mix_cells, &
                                   mix, &
                                   overlap, &
                                   EpotNew)

        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Hard_Spheres_Macro), intent(in) :: this_macro
        class(Neighbour_Cells), intent(in) :: other_mix_cells
        type(Particle_Index), intent(in) :: old
        type(Particle_Index), intent(inout) :: new
        class(Between_Hard_Spheres_Potential_Energy), intent(in) :: mix
        logical, intent(out) :: overlap
        type(Particle_Energy), intent(out) :: EpotNew
        
        new%position(:) = other_spheres%get_position(new%other_number)
        
        if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
            new%same_iCell = this_macro%same_cells%index_from_position(new%position)
            call this_macro%hard_potential%neighCells(Box%size, this_spheres, this_macro%same_cells, new, overlap, &
                                                      EpotNew%same)
        else
            new%mix_iCell = other_mix_cells%index_from_position(new%position)
            call mix%neighCells(Box%size, new, this_macro%mix_cells, other_spheres, overlap, EpotNew%mix)
        end if
        
        if (.not. overlap) then
        
            if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                new%mix_iCell = other_mix_cells%index_from_position(new%position)
                call mix%neighCells(Box%size, new, this_macro%mix_cells, other_spheres, overlap, &
                                         EpotNew%mix)
            else
                new%same_iCell = this_macro%same_cells%index_from_position(new%position)
                call this_macro%hard_potential%neighCells(Box%size, this_spheres, this_macro%same_cells, new, &
                                                    overlap, EpotNew%same)
            end if
            
            if (.not. overlap) then
            
                select type (this_spheres)
                    type is (Dipolar_Hard_Spheres)
                        new%orientation(:) = this_spheres%get_orientation(new%number)
                        select type (this_macro)
                            type is (Dipolar_Hard_Spheres_Macro)
                                EpotNew%same = this_macro%ewald_real%solo(Box%size, this_spheres, new) + &
                                               this_macro%ewald_reci%move(Box, old, new)
                        end select
                                       
                end select
            
            end if
                
        end if
    
    end subroutine after_switch_energy
    
    subroutine after_switch_update(Box, &
                                   this_spheres, this_macro, old, new, &
                                   other_mix_cells)

        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this_spheres
        class(Hard_Spheres_Macro), intent(inout) :: this_macro
        class(Neighbour_Cells), intent(inout) :: other_mix_cells
        type(Particle_Index), intent(in) :: old, new
        
        call this_spheres%set_position(old%number, new%position)
        
        select type (this_macro)
            type is (Dipolar_Hard_Spheres_Macro)
                call this_macro%ewald_reci%update_structure_move(Box, old, new)
        end select
        
        if (old%same_iCell /= new%same_iCell) then
            call this_macro%same_cells%remove_col_from_cell(old%number, old%same_iCell)
            call this_macro%same_cells%add_col_to_cell(new%number, new%same_iCell)
        end if
        if (old%mix_iCell /= new%mix_iCell) then
            call other_mix_cells%remove_col_from_cell(old%number, old%mix_iCell)
            call other_mix_cells%add_col_to_cell(new%number, new%mix_iCell)
        end if
        
    end subroutine after_switch_update
    
    !> Dipole rotation
    
    subroutine rotate(Box, &
                      this_spheres, this_macro, this_observables)
    
        type(Box_Dimensions), intent(in) :: Box
        class(Dipolar_Hard_spheres), intent(inout) :: this_spheres
        class(Dipolar_Hard_spheres_Macro), intent(inout) :: this_macro
        class(Dipolar_Hard_spheres_Observables), intent(inout) :: this_observables
        
        real(DP) :: random
        type(Particle_Index) :: old, new
        real(DP) :: deltaEpot
        real(DP) :: deltaEpot_real, deltaEpot_self
        real(DP) :: real_EpotNew, real_EpotOld
        
        this_observables%rotate_num_hits = this_observables%rotate_num_hits + 1

        call random_number(random)
        old%number = int(random*this_spheres%get_num_particles()) + 1
        old%position(:) = this_spheres%get_position(old%number)
        old%orientation(:) = this_spheres%get_orientation(old%number)
        
        new%number = old%number
        new%position(:) = old%position(:)
        new%orientation(:) = old%orientation(:)
        call markov_surface(new%orientation, this_macro%rotation%get_delta())
        
        real_EpotOld = this_macro%ewald_real%solo(Box%size, this_spheres, old)
        real_EpotNew = this_macro%ewald_real%solo(Box%size, this_spheres, new)
        deltaEpot_real = real_EpotNew - real_EpotOld
        
        deltaEpot_self = this_macro%ewald_self%solo(new%orientation) - &
                         this_macro%ewald_self%solo(old%orientation)
        
        deltaEpot = deltaEpot_real + this_macro%ewald_reci%rotation(Box, old, new) - &
                    deltaEpot_self + this_macro%ewald_bound%rotation(Box%size, old%orientation, &
                                                                               new%orientation)
        
        call random_number(random)
        if (random < exp(-deltaEpot/Box%temperature)) then
        
            call this_macro%ewald_reci%update_structure_rotation(Box, old, new)
            call this_macro%ewald_bound%update_total_moment_rotation(old%orientation, new%orientation)
            call this_spheres%set_orientation(old%number, new%orientation)
            
            this_observables%potential_energy = this_observables%potential_energy + deltaEpot
            
        else
            this_observables%rotate_num_rejections = this_observables%rotate_num_rejections + 1
        end if
    
    end subroutine rotate

end module module_algorithms
