    ! External field ------------------------------------------------------------------------------

    !> \f[ U_\text{field} = -(\vec{M} \cdot \vec{E}) \f]
    
    pure function DipolarSpheres_Epot_field(this) result(Epot_field)

        class(DipolarSpheres), intent(in) :: this
        real(DP) :: Epot_field

        Epot_field = -dot_product(this%totalMoment, Field)

    end function DipolarSpheres_Epot_field
    + this%Epot_Field()
        
        
    real(DP), dimension(Ndim), parameter :: Field = [0._DP, 0._DP, 0._DP]
    
    use module_physics_micro, only: random_surface, markov_surface, deltaEpot_field_exchange, &


deltaEpot_field_exchange(+mTest)

+ deltaEpot_field_rotate(mOld, mNew)

    !> \Delta U_{N\rightarrow\N+1} = -(\vec{\mu}_{N+1} \cdot \vec{E})

    pure function deltaEpot_field_exchange(mCol)

        real(DP), dimension(:), intent(in) :: mCol
        real(DP) :: deltaEpot_field_exchange

        deltaEpot_field_exchange = -dot_product(mCol, Field)

    end function deltaEpot_field_exchange

    pure function deltaEpot_field_rotate(mOld, mNew)
    
        real(DP), dimension(:), intent(in) :: mOld, mNew
        real(DP) :: deltaEpot_field_rotate

        deltaEpot_field_rotate = -dot_product(mNew - mOld, Field)

    end function deltaEpot_field_rotate

GC
xNew(1:2) = Lsize(1:2) * xRand(1:2)
        xNew(3) = (Height - this%get_sigma()) * xRand(3) + this%get_sigma()/2._DP
        this%Volume (GC)
