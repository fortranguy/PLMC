    subroutine DipolarSpheres_rotate(this, obs)
    
        class(DipolarSpheres), intent(inout) :: this
        class(MoreObservables), intent(inout) :: obs
        
        real(DP) :: random
        integer :: iOld
        real(DP), dimension(Ndim) :: xCol
        real(DP), dimension(Ndim) :: mOld, mNew
        real(DP) :: deltaEpot, deltaEpot_real, deltaEpot_self
        real(DP) :: real_EpotNew, real_EpotOld
        integer :: iTotalCell
        
        call random_number(random)
        iOld = int(random*this%Ncol) + 1
        
        xCol(:) = this%positions(:, iOld)
        mOld(:) = this%orientations(:, iOld)
        mNew(:) = mOld(:)
        call markov_surface(mNew, this%rotate_delta)
        
        iTotalCell = this%sameCells%index_from_position(xCol)
        real_EpotNew = this%Epot_real_solo(iOld, xCol, mNew)
        real_EpotOld = this%Epot_real_solo(iOld, xCol, mOld)
        deltaEpot_real = real_EpotNew - real_EpotOld        
        
        deltaEpot_self = this%Epot_self_solo(mNew) - this%Epot_self_solo(mOld)
        
        deltaEpot = deltaEpot_real + this%deltaEpot_reci_rotate(xCol, mOld, mNew) - deltaEpot_self + &
                    this%deltaEpot_bound_rotate(mOld, mNew)
        
        call random_number(random)
        if (random < exp(-deltaEpot/Temperature)) then
        
            call this%deltaEpot_reci_rotate_update_structure(xCol, mOld, mNew)
            call this%deltaEpot_bound_rotate_update_totalMoment(mOld, mNew)
            this%orientations(:, iOld) = mNew(:)
            
            obs%Epot = obs%Epot + deltaEpot
            
        else
            obs%rotate_Nreject = obs%rotate_Nreject + 1
        end if
    
    end subroutine DipolarSpheres_rotate
