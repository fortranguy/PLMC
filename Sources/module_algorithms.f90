module module_algorithms

use data_precisions, only: DP
use data_box, only: Ndim
use data_monte_carlo, only: Temperature
use module_types_micro, only: Box_Dimensions, Particle_Index, Particle_Energy
use module_physics_micro, only: random_surface, markov_surface
use class_neighbour_cells
use class_hard_spheres
use class_small_move
use class_dipolar_spheres
use class_hard_spheres_potential
use class_small_rotation
use class_mixing_potential
use class_observables

implicit none
private
public move, widom, switch, rotate

contains

    !> Particle move
    
    subroutine move(Box, this_spheres, this_same_cells, this_mix_cells, this_move, this_hard_potential, &
                    this_obs, other_spheres, other_mix_cells, mix, mix_Epot)
    
        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this_spheres, other_spheres
        class(Neighbour_Cells), intent(inout) :: this_same_cells, this_mix_cells, other_mix_cells
        class(Small_Move), intent(in) :: this_move
        class(Hard_Spheres_Potential), intent(in) :: this_hard_potential
        class(Observables), intent(inout) :: this_obs
        class(Mixing_Potential), intent(in) :: mix
        real(DP), intent(inout) :: mix_Epot
        
        real(DP) :: random
        real(DP), dimension(Ndim) :: xRand
        type(Particle_Index) :: old, new
        logical :: overlap
        real(DP) :: deltaEpot
        real(DP) :: this_deltaEpot, mix_deltaEpot
        real(DP) :: this_EpotNew, this_EpotOld
        real(DP) :: mix_EpotNew, mix_EpotOld
        
        real(DP) :: this_EpotNew_real, this_EpotOld_real
        
        this_obs%move_Nhit = this_obs%move_Nhit + 1
        
        call random_number(random)
        old%number = int(random*this_spheres%get_num_particles()) + 1
        old%position(:) = this_spheres%get_position(old%number)
        
        new%number = old%number
        call random_number(xRand)
        new%position(:) = old%position(:) + (xRand(:)-0.5_DP)*this_move%get_delta()
        new%position(:) = modulo(new%position(:), Box%size(:))
        
        if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
            new%same_iCell = this_same_cells%index_from_position(new%position)
            call this_hard_potential%neighCells(Box%size, this_spheres, this_same_cells, new, &
                                                overlap, this_EpotNew)
        else
            new%mix_iCell = other_mix_cells%index_from_position(new%position)
            call mix%Epot_neighCells(Box%size, new, this_mix_cells, other_spheres, overlap, &
                                     mix_EpotNew)
        end if
        
        if (.not. overlap) then
        
            if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                new%mix_iCell = other_mix_cells%index_from_position(new%position)
                call mix%Epot_neighCells(Box%size, new, this_mix_cells, other_spheres, overlap, &
                                         mix_EpotNew)
            else
                new%same_iCell = this_same_cells%index_from_position(new%position)
                call this_hard_potential%neighCells(Box%size, this_spheres, this_same_cells, new, &
                                                    overlap, this_EpotNew)
            end if
                        
            if (.not. overlap) then
    
                old%same_iCell = this_same_cells%index_from_position(old%position)
                select type (this_spheres)
                    type is (Dipolar_Spheres)
                        old%orientation(:) = this_spheres%get_orientation(old%number)
                        new%orientation(:) = old%orientation(:)
                        this_EpotNew_real = this_spheres%Epot_real_solo(Box%size, new)
                        this_EpotOld_real = this_spheres%Epot_real_solo(Box%size, old)
                        this_deltaEpot = (this_EpotNew_real - this_EpotOld_real) + &
                                         this_spheres%deltaEpot_reci_move(Box, old, new)
                    class default
                        call this_hard_potential%neighCells(Box%size, this_spheres, this_same_cells, &
                                                            old, overlap, this_EpotOld)
                        this_deltaEpot = this_EpotNew - this_EpotOld
                end select
                    
                old%mix_iCell = other_mix_cells%index_from_position(old%position)
                call mix%Epot_neighCells(Box%size, old, this_mix_cells, other_spheres, overlap, &
                                         mix_EpotOld)
                
                mix_deltaEpot = mix_EpotNew - mix_EpotOld

                deltaEpot = this_deltaEpot + mix_deltaEpot
                
                call random_number(random)
                if (random < exp(-deltaEpot/Temperature)) then
                
                    select type (this_spheres)
                        type is (Dipolar_Spheres)
                            call this_spheres%reci_update_structure_move(Box, old, new)
                    end select
                
                    call this_spheres%set_position(old%number, new%position)
                    this_obs%Epot = this_obs%Epot + this_deltaEpot
                    mix_Epot = mix_Epot + mix_deltaEpot
                    
                    if (old%same_iCell /= new%same_iCell) then
                        call this_same_cells%remove_col_from_cell(old%number, old%same_iCell)
                        call this_same_cells%add_col_to_cell(new%number, new%same_iCell)
                    end if
                    if (old%mix_iCell /= new%mix_iCell) then
                        call other_mix_cells%remove_col_from_cell(old%number, old%mix_iCell)
                        call other_mix_cells%add_col_to_cell(new%number, new%mix_iCell)
                    end if
                    
                else
                    this_obs%move_Nreject = this_obs%move_Nreject + 1
                end if
         
            else
                this_obs%move_Nreject = this_obs%move_Nreject + 1
            end if
            
        else
            this_obs%move_Nreject = this_obs%move_Nreject + 1
        end if
    
    end subroutine move
    
    !> Widom's method

    subroutine widom(Box, this_spheres, this_same_cells, this_mix_cells, this_hard_potential, this_obs, &
                     other_spheres, other_mix_cells, mix)
        
        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(in) :: this_spheres
        class(Neighbour_Cells), intent(inout) :: this_same_cells, this_mix_cells, other_mix_cells
        class(Hard_Spheres_Potential), intent(in) :: this_hard_potential
        class(Observables), intent(inout) :: this_obs
        class(Hard_Spheres), intent(in) :: other_spheres
        class(Mixing_Potential), intent(in) :: mix
        
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
                test%same_iCell = this_same_cells%index_from_position(test%position)
                call this_hard_potential%neighCells(Box%size, this_spheres, this_same_cells, test, &
                                                    overlap, this_EpotTest)
            else
                test%mix_iCell = other_mix_cells%index_from_position(test%position)
                call mix%Epot_neighCells(Box%size, test, this_mix_cells, other_spheres, overlap, &
                                         mix_EpotTest)
            end if
            
            if (.not. overlap) then
            
                if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                    test%mix_iCell = other_mix_cells%index_from_position(test%position)
                    call mix%Epot_neighCells(Box%size, test, this_mix_cells, other_spheres, overlap, &
                                             mix_EpotTest)
                else
                    test%same_iCell = this_same_cells%index_from_position(test%position)
                    call this_hard_potential%neighCells(Box%size, this_spheres, this_same_cells, test, &
                                                        overlap, this_EpotTest)
                end if
                
                if (.not. overlap) then
                
                    select type (this_spheres)
                        type is (Dipolar_Spheres)
                            test%add = .true.
                            test%orientation(:) = random_surface()
                            this_EpotTest = this_spheres%Epot_real_solo(Box%size, test) + &
                                            this_spheres%deltaEpot_reci_exchange(Box, test) - &
                                            this_spheres%Epot_self_solo(test%orientation) + &
                                            this_spheres%deltaEpot_bound_exchange(Box%size, &
                                                                                  test%orientation)
                    end select
                
                    EpotTest = this_EpotTest + mix_EpotTest
                    widTestSum = widTestSum + exp(-EpotTest/Temperature)
                    
                end if
                
            end if
            
        end do
        
        this_obs%activ = widTestSum/real(this_spheres%get_widom_num_particles(), DP)
        
    end subroutine widom
    
    !> Particle switch
    
    subroutine switch(Box, type1_spheres, type1_same_cells, type1_mix_cells, type1_hard_potential, &
                      type1_obs, type2_spheres, type2_same_cells, type2_mix_cells, type2_hard_potential, &
                      type2_obs, mix, mix_Epot, switch_Nreject)
    
        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: type1_spheres, type2_spheres
        class(Neighbour_Cells), intent(inout) :: type1_same_cells, type1_mix_cells
        class(Neighbour_Cells), intent(inout) :: type2_same_cells, type2_mix_cells
        class(Hard_Spheres_Potential), intent(in) :: type1_hard_potential, type2_hard_potential
        class(Observables), intent(inout) :: type1_obs, type2_obs
        class(Mixing_Potential), intent(in) :: mix
        real(DP), intent(inout) :: mix_Epot
        integer, intent(inout) :: switch_Nreject
        
        real(DP) :: random
        type(Particle_Index) :: old1, old2
        type(Particle_Index) :: new1, new2
        logical :: overlap
        real(DP) :: deltaEpot, type1_deltaEpot, type2_deltaEpot
        real(DP) :: type1_mix_deltaEpot, type2_mix_deltaEpot
        type(Particle_Energy) :: type1_EpotOld, type1_EpotNew
        type(Particle_Energy) :: type2_EpotOld, type2_EpotNew
        
        if (type1_spheres%get_num_particles()==0 .or. type2_spheres%get_num_particles()==0) then
            switch_Nreject = switch_Nreject + 1
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
                
        call before_switch_energy(Box%size, type1_spheres, type1_same_cells, type1_mix_cells, old1, &
                                  type2_spheres, type2_mix_cells, mix, type1_EpotOld)
        call before_switch_energy(Box%size, type2_spheres, type2_same_cells, type2_mix_cells, old2, &
                                  type1_spheres, type1_mix_cells, mix, type2_EpotOld)
             
        call after_switch_energy(Box, type1_spheres, type1_same_cells, type1_mix_cells, &
                                 type1_hard_potential, old1, new1, type2_spheres, type2_mix_cells, &
                                 mix, overlap, type1_EpotNew)
        
        if (.not. overlap) then
        
            call after_switch_energy(Box, type2_spheres, type2_same_cells, type2_mix_cells, &
                                     type2_hard_potential, old2, new2, type1_spheres, type1_mix_cells, & 
                                     mix, overlap, type2_EpotNew)
            
            if (.not. overlap) then

                type1_deltaEpot = type1_EpotNew%same - type1_EpotOld%same
                type1_mix_deltaEpot = type1_EpotNew%mix - type1_EpotOld%mix
                type2_deltaEpot = type2_EpotNew%same - type2_EpotOld%same
                type2_mix_deltaEpot = type2_EpotNew%mix - type2_EpotOld%mix
                deltaEpot = type1_deltaEpot + type1_mix_deltaEpot + type2_deltaEpot + &
                            type2_mix_deltaEpot
                
                call random_number(random)
                if (random < exp(-deltaEpot/Temperature)) then
                
                    call after_switch_update(Box, type1_spheres, type1_same_cells, old1, new1, &
                                             type2_mix_cells)
                    call after_switch_update(Box, type2_spheres, type2_same_cells, old2, new2, &
                                             type1_mix_cells)
                                             
                    type1_obs%Epot = type1_obs%Epot + type1_deltaEpot
                    type2_obs%Epot = type2_obs%Epot + type2_deltaEpot
                    mix_Epot = mix_Epot + type1_mix_deltaEpot + type2_mix_deltaEpot
                    
                else
                    switch_Nreject = switch_Nreject + 1
                end if
                
            else
                switch_Nreject = switch_Nreject + 1
            end if
            
        else
            switch_Nreject = switch_Nreject + 1
        end if
        
    end subroutine switch
    
    subroutine before_switch_energy(Box_size, this_spheres, this_same_cells, this_mix_cells, old, &
                                    other_spheres, other_mix_cells, mix, EpotOld)
        
        real(DP), dimension(:), intent(in) :: Box_size
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Neighbour_Cells), intent(inout) :: this_same_cells, this_mix_cells, other_mix_cells
        type(Particle_Index), intent(inout) :: old
        class(Mixing_Potential), intent(in) :: mix
        type(Particle_Energy), intent(out) :: EpotOld
        logical :: overlap
        
        old%position(:) = this_spheres%get_position(old%number)
        
        old%same_iCell = this_same_cells%index_from_position(old%position)
        select type (this_spheres)
            type is (Dipolar_Spheres)
                old%orientation(:) = this_spheres%get_orientation(old%number)
                EpotOld%same = this_spheres%Epot_real_solo(Box_size, old)
                               ! Epot_reci: cf. after_switch_energy
            type is (Hard_Spheres)
                EpotOld%same = 0._DP
        end select
        
        old%mix_iCell = other_mix_cells%index_from_position(old%position)
        call mix%Epot_neighCells(Box_size, old, this_mix_cells, other_spheres, overlap, EpotOld%mix)
        
    end subroutine before_switch_energy
    
    subroutine after_switch_energy(Box, this_spheres, this_same_cells, this_mix_cells, &
                                   this_hard_potential, old, new, other_spheres, other_mix_cells, &
                                   mix, overlap, EpotNew)

        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(in) :: this_spheres, other_spheres
        class(Neighbour_Cells), intent(inout) :: this_same_cells, this_mix_cells, other_mix_cells
        class(Hard_Spheres_Potential), intent(in) :: this_hard_potential
        type(Particle_Index), intent(in) :: old
        type(Particle_Index), intent(inout) :: new
        class(Mixing_Potential), intent(in) :: mix
        logical, intent(out) :: overlap
        type(Particle_Energy), intent(out) :: EpotNew
        
        new%position(:) = other_spheres%get_position(new%other_number)
        
        if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
            new%same_iCell = this_same_cells%index_from_position(new%position)
            call this_hard_potential%neighCells(Box%size, this_spheres, this_same_cells, new, overlap, &
                                                EpotNew%same)
        else
            new%mix_iCell = other_mix_cells%index_from_position(new%position)
            call mix%Epot_neighCells(Box%size, new, this_mix_cells, other_spheres, overlap, EpotNew%mix)
        end if
        
        if (.not. overlap) then
        
            if (this_spheres%get_num_particles() >= other_spheres%get_num_particles()) then
                new%mix_iCell = other_mix_cells%index_from_position(new%position)
                call mix%Epot_neighCells(Box%size, new, this_mix_cells, other_spheres, overlap, &
                                         EpotNew%mix)
            else
                new%same_iCell = this_same_cells%index_from_position(new%position)
                call this_hard_potential%neighCells(Box%size, this_spheres, this_same_cells, new, &
                                                    overlap, EpotNew%same)
            end if
            
            if (.not. overlap) then
            
                select type (this_spheres)
                    type is (Dipolar_Spheres)
                        new%orientation(:) = this_spheres%get_orientation(new%number)
                        EpotNew%same = this_spheres%Epot_real_solo(Box%size, new) + &
                                       this_spheres%deltaEpot_reci_move(Box, old, new)
                end select
            
            end if
                
        end if
    
    end subroutine after_switch_energy
    
    subroutine after_switch_update(Box, this_spheres, this_same_cells, old, new, &
                                   other_mix_cells)

        type(Box_Dimensions), intent(in) :: Box
        class(Hard_Spheres), intent(inout) :: this_spheres
        class(Neighbour_Cells), intent(inout) :: this_same_cells, other_mix_cells
        type(Particle_Index), intent(in) :: old, new
        
        call this_spheres%set_position(old%number, new%position)
        
        select type (this_spheres)
            type is (Dipolar_Spheres)
                call this_spheres%reci_update_structure_move(Box, old, new)
        end select
        
        if (old%same_iCell /= new%same_iCell) then
            call this_same_cells%remove_col_from_cell(old%number, old%same_iCell)
            call this_same_cells%add_col_to_cell(new%number, new%same_iCell)
        end if
        if (old%mix_iCell /= new%mix_iCell) then
            call other_mix_cells%remove_col_from_cell(old%number, old%mix_iCell)
            call other_mix_cells%add_col_to_cell(new%number, new%mix_iCell)
        end if
        
    end subroutine after_switch_update
    
    !> Dipole rotation
    
    subroutine rotate(Box, spheres, rotation, obs)
    
        type(Box_Dimensions), intent(in) :: Box
        class(Dipolar_Spheres), intent(inout) :: spheres
        class(Small_rotation), intent(in) :: rotation
        class(MoreObservables), intent(inout) :: obs
        
        real(DP) :: random
        type(Particle_Index) :: old, new
        real(DP) :: deltaEpot
        real(DP) :: deltaEpot_real, deltaEpot_self
        real(DP) :: real_EpotNew, real_EpotOld
        
        obs%rotate_Nhit = obs%rotate_Nhit + 1

        call random_number(random)
        old%number = int(random*spheres%get_num_particles()) + 1
        old%position(:) = spheres%get_position(old%number)
        old%orientation(:) = spheres%get_orientation(old%number)
        
        new%number = old%number
        new%position(:) = old%position(:)
        new%orientation(:) = old%orientation(:)
        call markov_surface(new%orientation, rotation%get_delta())
        
        real_EpotOld = spheres%Epot_real_solo(Box%size, old)
        real_EpotNew = spheres%Epot_real_solo(Box%size, new)
        deltaEpot_real = real_EpotNew - real_EpotOld
        
        deltaEpot_self = spheres%Epot_self_solo(new%orientation) - spheres%Epot_self_solo(old%orientation)
        
        deltaEpot = deltaEpot_real + spheres%deltaEpot_reci_rotate(Box, old, new) - &
                    deltaEpot_self + spheres%deltaEpot_bound_rotate(Box%size, old%orientation, &
                                                                    new%orientation)
        
        call random_number(random)
        if (random < exp(-deltaEpot/Temperature)) then
        
            call spheres%reci_update_structure_rotate(Box, old, new)
            call spheres%update_totalMoment_rotate(old%orientation, new%orientation)
            call spheres%set_orientation(old%number, new%orientation)
            
            obs%Epot = obs%Epot + deltaEpot
            
        else
            obs%rotate_Nreject = obs%rotate_Nreject + 1
        end if
    
    end subroutine rotate

end module module_algorithms
