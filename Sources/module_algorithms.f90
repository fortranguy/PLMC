module module_algorithms

use data_precisions, only: DP
use data_box, only: Ndim
use data_monteCarlo, only: Temperature
use module_types, only: Box_dimensions, particle_index, particle_energy
use module_physics_micro, only: random_surface, markov_surface
use class_hardSpheres
use class_dipolarSpheres
use class_mixingPotential
use class_observables

implicit none
private
public move, widom, switch, rotate

contains

    !> Particle move
    
    subroutine move(Box, this, this_obs, other, mix, mix_Epot)
    
        type(Box_dimensions), intent(in) :: Box
        class(HardSpheres), intent(inout) :: this
        class(Observables), intent(inout) :: this_obs
        class(HardSpheres), intent(inout) :: other
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: mix_Epot
        
        real(DP) :: random
        real(DP), dimension(Ndim) :: xRand
        type(particle_index) :: old, new
        logical :: overlap
        real(DP) :: deltaEpot
        real(DP) :: this_deltaEpot, mix_deltaEpot
        real(DP) :: this_EpotNew, this_EpotOld
        real(DP) :: mix_EpotNew, mix_EpotOld
        
        real(DP) :: this_EpotNew_real, this_EpotOld_real
        
        this_obs%move_Nhit = this_obs%move_Nhit + 1
        
        call random_number(random)
        old%number = int(random*this%get_num_particles()) + 1
        old%xCol(:) = this%all_positions(:, old%number)
        
        new%number = old%number
        call random_number(xRand)
        new%xCol(:) = old%xCol(:) + (xRand(:)-0.5_DP)*this%get_move_delta()
        new%xCol(:) = modulo(new%xCol(:), Box%size(:))
        
        if (this%get_num_particles() >= other%get_num_particles()) then
            new%same_iCell = this%sameCells%index_from_position(new%xCol)
            call this%Epot_neighCells(Box%size, new, overlap, this_EpotNew)
        else
            new%mix_iCell = other%mixCells%index_from_position(new%xCol)
            call mix%Epot_neighCells(Box%size, new, this%mixCells, other%all_positions, overlap, &
                                     mix_EpotNew)
        end if
        
        if (.not. overlap) then
        
            if (this%get_num_particles() >= other%get_num_particles()) then
                new%mix_iCell = other%mixCells%index_from_position(new%xCol)
                call mix%Epot_neighCells(Box%size, new, this%mixCells, other%all_positions, overlap, &
                                         mix_EpotNew)
            else
                new%same_iCell = this%sameCells%index_from_position(new%xCol)
                call this%Epot_neighCells(Box%size, new, overlap, this_EpotNew)
            end if
                        
            if (.not. overlap) then
    
                old%same_iCell = this%sameCells%index_from_position(old%xCol)
                select type (this)
                    type is (DipolarSpheres)
                        old%mCol(:) = this%orientations(:, old%number)
                        new%mCol(:) = old%mCol(:)
                        this_EpotNew_real = this%Epot_real_solo(Box%size, new)
                        this_EpotOld_real = this%Epot_real_solo(Box%size, old)
                        this_deltaEpot = (this_EpotNew_real - this_EpotOld_real) + &
                                         this%deltaEpot_reci_move(Box, old, new)
                    class default
                        call this%Epot_neighCells(Box%size, old, overlap, this_EpotOld)
                        this_deltaEpot = this_EpotNew - this_EpotOld
                end select
                    
                old%mix_iCell = other%mixCells%index_from_position(old%xCol)
                call mix%Epot_neighCells(Box%size, old, this%mixCells, other%all_positions, overlap, &
                                         mix_EpotOld)
                
                mix_deltaEpot = mix_EpotNew - mix_EpotOld

                deltaEpot = this_deltaEpot + mix_deltaEpot
                
                call random_number(random)
                if (random < exp(-deltaEpot/Temperature)) then
                
                    select type (this)
                        type is (DipolarSpheres)
                            call this%reci_update_structure_move(Box, old, new)
                    end select
                
                    this%all_positions(:, old%number) = new%xCol(:)
                    this_obs%Epot = this_obs%Epot + this_deltaEpot
                    mix_Epot = mix_Epot + mix_deltaEpot
                    
                    if (old%same_iCell /= new%same_iCell) then
                        call this%sameCells%remove_col_from_cell(old%number, old%same_iCell)
                        call this%sameCells%add_col_to_cell(new%number, new%same_iCell)
                    end if
                    if (old%mix_iCell /= new%mix_iCell) then
                        call other%mixCells%remove_col_from_cell(old%number, old%mix_iCell)
                        call other%mixCells%add_col_to_cell(new%number, new%mix_iCell)
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

    subroutine widom(Box, this, this_obs, other, mix)
        
        type(Box_dimensions), intent(in) :: Box
        class(HardSpheres), intent(in) :: this
        class(Observables), intent(inout) :: this_obs
        class(HardSpheres), intent(in) :: other
        class(MixingPotential), intent(in) :: mix
        
        integer :: iWidom
        real(DP) :: widTestSum
        real(DP), dimension(Ndim) :: xRand
        type(particle_index) :: test
        logical :: overlap
        real(DP) :: EpotTest
        real(DP) :: this_EpotTest, mix_EpotTest
        
        widTestSum = 0._DP
        test%number = 0
        
        do iWidom = 1, this%get_Nwidom()
            
            call random_number(xRand)
            test%xCol(:) = Box%size(:) * xRand(:)

            if (this%get_num_particles() >= other%get_num_particles()) then
                test%same_iCell = this%sameCells%index_from_position(test%xCol)
                call this%Epot_neighCells(Box%size, test, overlap, this_EpotTest)
            else
                test%mix_iCell = other%mixCells%index_from_position(test%xCol)
                call mix%Epot_neighCells(Box%size, test, this%mixCells, other%all_positions, overlap, &
                                         mix_EpotTest)
            end if
            
            if (.not. overlap) then
            
                if (this%get_num_particles() >= other%get_num_particles()) then
                    test%mix_iCell = other%mixCells%index_from_position(test%xCol)
                    call mix%Epot_neighCells(Box%size, test, this%mixCells, other%all_positions, overlap, &
                                             mix_EpotTest)
                else
                    test%same_iCell = this%sameCells%index_from_position(test%xCol)
                    call this%Epot_neighCells(Box%size, test, overlap, this_EpotTest)
                end if
                
                if (.not. overlap) then
                
                    select type (this)
                        type is (DipolarSpheres)
                            test%add = .true.
                            test%mCol(:) = random_surface()
                            this_EpotTest = this%Epot_real_solo(Box%size, test) + &
                                            this%deltaEpot_reci_exchange(Box, test) - &
                                            this%Epot_self_solo(test%mCol) + &
                                            this%deltaEpot_bound_exchange(Box%size, test%mCol)
                    end select
                
                    EpotTest = this_EpotTest + mix_EpotTest
                    widTestSum = widTestSum + exp(-EpotTest/Temperature)
                    
                end if
                
            end if
            
        end do
        
        this_obs%activ = widTestSum/real(this%get_Nwidom(), DP)
        
    end subroutine widom
    
    !> Particle switch
    
    subroutine before_switch_energy(Box_size, this, old, other, mix, EpotOld)
        
        real(DP), dimension(:), intent(in) :: Box_size
        class(HardSpheres), intent(in) :: this, other
        type(particle_index), intent(inout) :: old
        class(MixingPotential), intent(in) :: mix
        type(particle_energy), intent(out) :: EpotOld
        logical :: overlap
        
        old%xCol(:) = this%all_positions(:, old%number)
        
        old%same_iCell = this%sameCells%index_from_position(old%xCol)
        select type (this)
            type is (DipolarSpheres)
                old%mCol(:) = this%orientations(:, old%number)
                EpotOld%same = this%Epot_real_solo(Box_size, old) ! Epot_reci: cf. after_switch_energy
            type is (HardSpheres)
                EpotOld%same = 0._DP
        end select
        
        old%mix_iCell = other%mixCells%index_from_position(old%xCol)
        call mix%Epot_neighCells(Box_size, old, this%mixCells, other%all_positions, overlap, EpotOld%mix)
        
    end subroutine before_switch_energy
    
    subroutine after_switch_energy(Box, this, old, new, other, mix, overlap, EpotNew)

        type(Box_dimensions), intent(in) :: Box
        class(HardSpheres), intent(in) :: this, other
        type(particle_index), intent(in) :: old
        type(particle_index), intent(inout) :: new
        class(MixingPotential), intent(in) :: mix
        logical, intent(out) :: overlap
        type(particle_energy), intent(out) :: EpotNew
        
        new%xCol(:) = other%all_positions(:, new%other_number)
        
        if (this%get_num_particles() >= other%get_num_particles()) then ! optimisation: more chance to overlap
            new%same_iCell = this%sameCells%index_from_position(new%xCol)
            call this%Epot_neighCells(Box%size, new, overlap, EpotNew%same)
        else
            new%mix_iCell = other%mixCells%index_from_position(new%xCol)
            call mix%Epot_neighCells(Box%size, new, this%mixCells, other%all_positions, overlap, &
                                     EpotNew%mix)
        end if
        
        if (.not. overlap) then
        
            if (this%get_num_particles() >= other%get_num_particles()) then
                new%mix_iCell = other%mixCells%index_from_position(new%xCol)
                call mix%Epot_neighCells(Box%size, new, this%mixCells, other%all_positions, overlap, &
                                         EpotNew%mix)
            else
                new%same_iCell = this%sameCells%index_from_position(new%xCol)
                call this%Epot_neighCells(Box%size, new, overlap, EpotNew%same)
            end if
            
            if (.not. overlap) then
            
                select type (this)
                    type is (DipolarSpheres)
                        new%mCol(:) = this%orientations(:, new%number)
                        EpotNew%same = this%Epot_real_solo(Box%size, new) + &
                                      this%deltaEpot_reci_move(Box, old, new)
                end select
            
            end if
                
        end if
    
    end subroutine after_switch_energy
    
    subroutine after_switch_update(Box, this, old, new, other)

        type(Box_dimensions), intent(in) :: Box
        class(HardSpheres), intent(inout) :: this, other
        type(particle_index), intent(in) :: old, new
        
        this%all_positions(:, old%number) = new%xCol(:)
        
        select type (this)
            type is (DipolarSpheres)
                call this%reci_update_structure_move(Box, old, new)
        end select
        
        if (old%same_iCell /= new%same_iCell) then
            call this%sameCells%remove_col_from_cell(old%number, old%same_iCell)
            call this%sameCells%add_col_to_cell(new%number, new%same_iCell)
        end if
        if (old%mix_iCell /= new%mix_iCell) then
            call other%mixCells%remove_col_from_cell(old%number, old%mix_iCell)
            call other%mixCells%add_col_to_cell(new%number, new%mix_iCell)
        end if
        
    end subroutine after_switch_update
    
    subroutine switch(Box, type1, type1_obs, type2, type2_obs, mix, mix_Epot, switch_Nreject)
    
        type(Box_dimensions), intent(in) :: Box
        class(HardSpheres), intent(inout) :: type1, type2
        class(Observables), intent(inout) :: type1_obs, type2_obs
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inout) :: mix_Epot
        integer, intent(inout) :: switch_Nreject
        
        real(DP) :: random
        type(particle_index) :: old1, old2
        type(particle_index) :: new1, new2
        logical :: overlap
        real(DP) :: deltaEpot, type1_deltaEpot, type2_deltaEpot
        real(DP) :: type1_mix_deltaEpot, type2_mix_deltaEpot
        type(particle_energy) :: type1_EpotOld, type1_EpotNew
        type(particle_energy) :: type2_EpotOld, type2_EpotNew
        
        if (type1%get_num_particles()==0 .or. type2%get_num_particles()==0) then
            switch_Nreject = switch_Nreject + 1
            return
        end if
        
        call random_number(random)
        old1%number = int(random*type1%get_num_particles()) + 1
        new1%number = old1%number
        call random_number(random)
        old2%number = int(random*type2%get_num_particles()) + 1
        new2%number = old2%number
        
        old1%other_number = old2%number; new1%other_number = new2%number
        old2%other_number = old1%number; new2%other_number = new1%number
                
        call before_switch_energy(Box%size, type1, old1, type2, mix, type1_EpotOld)
        call before_switch_energy(Box%size, type2, old2, type1, mix, type2_EpotOld)        
             
        call after_switch_energy(Box, type1, old1, new1, type2, mix, overlap, type1_EpotNew)
        
        if (.not. overlap) then
        
            call after_switch_energy(Box, type2, old2, new2, type1, mix, overlap, type2_EpotNew)
            
            if (.not. overlap) then

                type1_deltaEpot = type1_EpotNew%same - type1_EpotOld%same
                type1_mix_deltaEpot = type1_EpotNew%mix - type1_EpotOld%mix
                type2_deltaEpot = type2_EpotNew%same - type2_EpotOld%same
                type2_mix_deltaEpot = type2_EpotNew%mix - type2_EpotOld%mix
                deltaEpot = type1_deltaEpot + type1_mix_deltaEpot + type2_deltaEpot + &
                            type2_mix_deltaEpot
                
                call random_number(random)
                if (random < exp(-deltaEpot/Temperature)) then
                
                    call after_switch_update(Box, type1, old1, new1, type2)
                    call after_switch_update(Box, type2, old2, new2, type1)
                                             
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
    
    !> Dipole rotation
    
    subroutine rotate(Box, this, obs)
    
        type(Box_dimensions), intent(in) :: Box
        class(DipolarSpheres), intent(inout) :: this
        class(MoreObservables), intent(inout) :: obs
        
        real(DP) :: random
        type(particle_index) :: old, new
        real(DP) :: deltaEpot
        real(DP) :: deltaEpot_real, deltaEpot_self
        real(DP) :: real_EpotNew, real_EpotOld
        
        obs%rotate_Nhit = obs%rotate_Nhit + 1

        call random_number(random)
        old%number = int(random*this%get_num_particles()) + 1
        old%xCol(:) = this%all_positions(:, old%number)
        old%mCol(:) = this%orientations(:, old%number)
        
        new%number = old%number
        new%xCol(:) = old%xCol(:)
        new%mCol(:) = old%mCol(:)
        call markov_surface(new%mCol, this%get_rotation_delta())
        
        real_EpotOld = this%Epot_real_solo(Box%size, old)
        real_EpotNew = this%Epot_real_solo(Box%size, new)
        deltaEpot_real = real_EpotNew - real_EpotOld
        
        deltaEpot_self = this%Epot_self_solo(new%mCol) - this%Epot_self_solo(old%mCol)
        
        deltaEpot = deltaEpot_real + this%deltaEpot_reci_rotate(Box, old, new) - &
                    deltaEpot_self + this%deltaEpot_bound_rotate(Box%size, old%mCol, new%mCol)
        
        call random_number(random)
        if (random < exp(-deltaEpot/Temperature)) then
        
            call this%reci_update_structure_rotate(Box, old, new)
            call this%update_totalMoment_rotate(old%mCol, new%mCol)
            this%orientations(:, old%number) = new%mCol(:)
            
            obs%Epot = obs%Epot + deltaEpot
            
        else
            obs%rotate_Nreject = obs%rotate_Nreject + 1
        end if
    
    end subroutine rotate

end module module_algorithms
