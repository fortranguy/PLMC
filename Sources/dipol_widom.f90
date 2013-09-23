    subroutine DipolarSpheres_widom(this, other, mix, activ)
        
        class(DipolarSpheres), intent(in) :: this
        class(Spheres), intent(in) :: other
        class(MixingPotential), intent(in) :: mix
        real(DP), intent(inOut) :: activ 
        
        integer :: iWidom
        real(DP) :: widTestSum
        real(DP), dimension(Ndim) :: xTest, xRand
        real(DP), dimension(Ndim) :: mTest
        integer :: this_iCellTest, mix_iCellTest
        logical :: overlap
        real(DP) :: EpotTest
        real(DP) :: this_EpotTest
        real(DP) :: mix_EpotTest
        
        widTestSum = 0._DP
        
        do iWidom = 1, this%Nwidom
            
            call random_number(xRand)
            xTest(:) = Lsize(:) * xRand(:)

            if (this%Ncol >= other%Ncol) then
                this_iCellTest = this%sameCells%index_from_position(xTest)
                call this%Epot_real_test_overlap(0, xTest, this_iCellTest, overlap)
            else            
                mix_iCellTest = this%mixCells%index_from_position(xTest)
                call mix%Epot_neighCells(xTest, mix_iCellTest, this%mixCells, other%positions, &
                                         overlap, mix_EpotTest)
            end if
            
            if (.not. overlap) then
                               
                if (this%Ncol >= other%Ncol) then
                    mix_iCellTest = this%mixCells%index_from_position(xTest)
                    call mix%Epot_neighCells(xTest, mix_iCellTest, this%mixCells, other%positions, &
                                            overlap, mix_EpotTest)
                else
                    this_iCellTest = this%sameCells%index_from_position(xTest)
                    call this%Epot_real_test_overlap(0, xTest, this_iCellTest, overlap)
                end if
                
                if (.not. overlap) then
                
                    mTest(:) = random_surface()
                                        
                    this_EpotTest = this%Epot_real_solo(0, xTest, mTest) + &
                                    this%deltaEpot_reci_exchange(xTest, +mTest) - &
                                    this%Epot_self_solo(mTest) + this%deltaEpot_bound_exchange(+mTest)
                
                    EpotTest = this_EpotTest + mix_EpotTest
                    widTestSum = widTestSum + exp(-EpotTest/Temperature)
                    
                end if
            
            end if
            
        end do
        
        activ = widTestSum/real(this%Nwidom, DP)
        
    end subroutine DipolarSpheres_widom
