    !> Field
    !> \f[
    !>      \vec{E}(\vec{r}_l) = -\frac{4\pi}{V} \sum_{\vec{k} \neq 0} w(\alpha, \vec{k})
    !>                                      \Re(S(\vec{k}) e^{-i(\vec{k}\cdot \vec{r}_l)}) \vec{k}
    !> \f]

    pure function Abstract_Weighted_Structure_solo_field(this, Box, particle) result(solo_field)

        class(Abstract_Weighted_Structure), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        real(DP), dimension(num_dimensions) :: solo_field

        real(DP), dimension(num_dimensions) :: position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx

        real(DP), dimension(num_dimensions) :: wave_vector

        integer :: kx, ky, kz

        position_div_box(:) = 2._DP*PI * particle%position(:) / Box%size(:)
        call set_fourier(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call set_fourier(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_fourier(Box%wave(3), position_div_box(3), exp_Ikx_3)

        solo_field(:) = 0._DP

        do kz = -Box%wave(3), Box%wave(3)
            wave_vector(3) = 2._DP*PI * real(kz, DP) / Box%size(3)

        do ky = -Box%wave(2), Box%wave(2)
            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

        do kx = -Box%wave(1), Box%wave(1)
            wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

            exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)

            solo_field(:) = solo_field(:) + &
                this%weight(kx, ky, kz) * real(this%structure(kx, ky, kz) * conjg(exp_Ikx), DP) * &
                wave_vector(:)

        end do

        end do

        end do

        solo_field(:) =-4._DP*PI / product(Box%size) * solo_field(:)

    end function Abstract_Weighted_Structure_solo_field

    !> Field
    !> \f[
    !>      \vec{E}(\vec{r}_{N+1}) = -\frac{2\pi}{V} \sum_{\vec{k} \neq 0} w(\alpha, \vec{k})
    !>                               [
    !>                                  2\Re(S(\vec{k}) e^{-i(\vec{k}\cdot \vec{r}_{N+1})}) +
    !>                                  (\vec{\mu}_{N+1} \cdot \vec{k})
    !>                               ] \vec{k}
    !> \f]

    pure function Abstract_Weighted_Structure_test_field(this, Box, particle) result(test_field)

        class(Abstract_Weighted_Structure), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        real(DP), dimension(num_dimensions) :: test_field

        real(DP), dimension(num_dimensions) :: position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx

        real(DP), dimension(num_dimensions) :: wave_vector

        integer :: kx, ky, kz

        position_div_box(:) = 2._DP*PI * particle%position(:) / Box%size(:)
        call set_fourier(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call set_fourier(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_fourier(Box%wave(3), position_div_box(3), exp_Ikx_3)

        test_field(:) = 0._DP

        do kz = -Box%wave(3), Box%wave(3)
            wave_vector(3) = 2._DP*PI * real(kz, DP) / Box%size(3)

        do ky = -Box%wave(2), Box%wave(2)
            wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

        do kx = -Box%wave(1), Box%wave(1)
            wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

            exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)

            test_field(:) = test_field(:) + &
                this%weight(kx, ky, kz) * &
                (2._DP*real(this%structure(kx, ky, kz) * conjg(exp_Ikx), DP) + &
                dot_product(particle%orientation, wave_vector)) * wave_vector(:)

        end do

        end do

        end do

        test_field(:) =-2._DP*PI / product(Box%size) * test_field(:)

    end function Abstract_Weighted_Structure_test_field


    !> Move

    !> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta S^2
    !> w(\alpha, \vec{k}) \f]
    !> \f[
    !>  \Delta S^2 = 2\Re[
    !>                  (\vec{\mu}_l\cdot\vec{k})
    !>                  (e^{-i\vec{k}\cdot\vec{x}^\prime_l} - e^{-i\vec{k}\cdot\vec{x}_l})
    !>                  S_\underline{l}(\vec{k})
    !>               ]
    !> \f]

    pure function Abstract_Weighted_Structure_move_energy(this, Box, old, new) result(move_energy)

        class(Abstract_Weighted_Structure), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: move_energy

        real(DP), dimension(num_dimensions) :: new_position_div_box, old_position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP) :: real_part

        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz

        new_position_div_box(:) = 2._DP*PI * new%position(:) / Box%size(:)
        call set_fourier(Box%wave(1), new_position_div_box(1), exp_IkxNew_1)
        call set_fourier(Box%wave(2), new_position_div_box(2), exp_IkxNew_2)
        call set_fourier(Box%wave(3), new_position_div_box(3), exp_IkxNew_3)

        old_position_div_box(:) = 2._DP*PI * old%position(:) / Box%size(:)
        call set_fourier(Box%wave(1), old_position_div_box(1), exp_IkxOld_1)
        call set_fourier(Box%wave(2), old_position_div_box(2), exp_IkxOld_2)
        call set_fourier(Box%wave(3), old_position_div_box(3), exp_IkxOld_3)

        move_energy = 0._DP

        do kz = 0, Box%wave(3) ! symmetry: half wave vectors -> double Energy
            wave_vector(3) = 2._DP*PI * real(kz, DP) / Box%size(3)

            do ky = -reciprocal_size_2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

                do kx = -reciprocal_size_1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                    wave_dot_orientation = dot_product(wave_vector, new%orientation)

                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)

                    real_part = wave_dot_orientation * real((conjg(exp_IkxNew) - conjg(exp_IkxOld)) * &
                                (this%structure(kx, ky, kz) - cmplx(wave_dot_orientation, 0._DP, DP) * &
                                exp_IkxOld), DP)

                    move_energy = move_energy + 2._DP * this%weight(kx, ky, kz) * real_part

                end do

            end do

        end do

        move_energy = 4._DP*PI / product(Box%size) * move_energy

    end function Abstract_Weighted_Structure_move_energy

    !> Update position -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot\vec{\mu}_l)
    !>                          (e^{+i\vec{k}\cdot\vec{x}^\prime_l} - e^{+i\vec{k}\cdot\vec{x}_l})
    !>  \f]
    !>

    pure subroutine Abstract_Weighted_Structure_update_structure_move(this, Box, old, new)

        class(Abstract_Weighted_Structure), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new

        real(DP), dimension(num_dimensions) :: new_position_div_box, old_position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxNew_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxNew_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxNew_3
        complex(DP) :: exp_IkxNew

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_IkxOld_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_IkxOld_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_IkxOld_3
        complex(DP) :: exp_IkxOld

        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz

        new_position_div_box(:) = 2._DP*PI * new%position(:) / Box%size(:)
        call set_fourier(Box%wave(1), new_position_div_box(1), exp_IkxNew_1)
        call set_fourier(Box%wave(2), new_position_div_box(2), exp_IkxNew_2)
        call set_fourier(Box%wave(3), new_position_div_box(3), exp_IkxNew_3)

        old_position_div_box(:) = 2._DP*PI * old%position(:) / Box%size(:)
        call set_fourier(Box%wave(1), old_position_div_box(1), exp_IkxOld_1)
        call set_fourier(Box%wave(2), old_position_div_box(2), exp_IkxOld_2)
        call set_fourier(Box%wave(3), old_position_div_box(3), exp_IkxOld_3)

        do kz = 0, Box%wave(3)
            wave_vector(3) = 2._DP*PI * real(kz, DP) / Box%size(3)

            do ky = -reciprocal_size_2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

                do kx = -reciprocal_size_1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                    wave_dot_orientation = dot_product(wave_vector, new%orientation)
                    exp_IkxNew = exp_IkxNew_1(kx) * exp_IkxNew_2(ky) * exp_IkxNew_3(kz)
                    exp_IkxOld = exp_IkxOld_1(kx) * exp_IkxOld_2(ky) * exp_IkxOld_3(kz)

                    this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                                 cmplx(wave_dot_orientation, 0._DP, DP) * &
                                                 (exp_IkxNew - exp_IkxOld)

                end do

            end do

        end do

    end subroutine Abstract_Weighted_Structure_update_structure_move

    !> Rotation

    !> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta S^2
    !>                                       w(\alpha, \vec{k}) \f]
    !> \f[
    !>  \Delta S^2 = (\vec{k} \cdot \vec{\mu}_l^\prime)^2 - (\vec{k} \cdot \vec{\mu}_l)^2 +
    !>               2\Re\{
    !>                  [(\vec{k} \cdot \vec{\mu}_l^\prime) - (\vec{k} \cdot \vec{\mu}_l)]
    !>                  e^{-i \vec{k} \cdot \vec{x}_l} S_\underline{l}(\vec{k})
    !>               \}
    !> \f]

    pure function Abstract_Weighted_Structure_rotation_energy(this, Box, old, new) result(rotation_energy)

        class(Abstract_Weighted_Structure), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new
        real(DP) :: rotation_energy

        real(DP), dimension(num_dimensions) :: position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx

        real(DP) :: real_part

        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_new_orientation, wave_dot_old_orientation
        integer :: kx, ky, kz

        position_div_box(:) = 2._DP*PI * new%position(:) / Box%size(:)
        call set_fourier(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call set_fourier(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_fourier(Box%wave(3), position_div_box(3), exp_Ikx_3)

        rotation_energy = 0._DP

        do kz = 0, Box%wave(3)
            wave_vector(3) = 2._DP*PI * real(kz, DP) / Box%size(3)

            do ky = -reciprocal_size_2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

                do kx = -reciprocal_size_1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                    wave_dot_new_orientation = dot_product(wave_vector, new%orientation)
                    wave_dot_old_orientation = dot_product(wave_vector, old%orientation)
                    exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)

                    real_part = wave_dot_new_orientation**2 - wave_dot_old_orientation**2
                    real_part = real_part + &
                                2._DP * (wave_dot_new_orientation - wave_dot_old_orientation) * &
                                real(conjg(exp_Ikx) * &
                                (this%structure(kx, ky, kz) - &
                                wave_dot_old_orientation * exp_Ikx), DP)

                    rotation_energy = rotation_energy + this%weight(kx, ky, kz) * real_part

                end do

            end do

        end do

        rotation_energy = 4._DP*PI / product(Box%size) * rotation_energy

    end function Abstract_Weighted_Structure_rotation_energy

    !> Update moment -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = [\vec{k}\cdot(\vec{\mu}_l^\prime - \vec{\mu}_l)]
    !>                          e^{+i\vec{k}\cdot\vec{x}_l}
    !>  \f]
    !>

    pure subroutine Abstract_Weighted_Structure_update_structure_rotation(this, Box, old, new)

        class(Abstract_Weighted_Structure), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: old, new

        real(DP), dimension(num_dimensions) :: position_div_box

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx

        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: k_dot_deltaMcol
        integer :: kx, ky, kz

        position_div_box(:) = 2._DP*PI * new%position(:)/Box%size(:)
        call set_fourier(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call set_fourier(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_fourier(Box%wave(3), position_div_box(3), exp_Ikx_3)

        do kz = 0, Box%wave(3)
            wave_vector(3) = 2._DP*PI * real(kz, DP) / Box%size(3)

            do ky = -reciprocal_size_2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

                do kx = -reciprocal_size_1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                    k_dot_deltaMcol = dot_product(wave_vector, &
                                                  new%orientation - old%orientation)
                    exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)

                    this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                                 cmplx(k_dot_deltaMcol, 0._DP, DP) * exp_Ikx

                end do

            end do

        end do

    end subroutine Abstract_Weighted_Structure_update_structure_rotation

    !> Energy of 1 dipole with others

    !> Addition :

    !> Difference of Energy
    !> \f[ \Delta U^{N+1} = \frac{2\pi}{V} \sum_{\vec{k} \neq \vec{0}}
    !>                          (\vec{k} \cdot +\vec{\mu}_{N+1}) w(\alpha, \vec{k})
    !>                          \{
    !>                              (\vec{k} \cdot +\vec{\mu}_{N+1}) +
    !>                              2\Re[S(\vec{k}) e^{-i \vec{k} \cdot \vec{x}_{N+1}}]
    !>                          \}
    !> \f]

    !> Summary: only the sign of \f$\vec{\mu}\f$ changes.

    pure function Abstract_Weighted_Structure_exchange_energy(this, Box, particle) result(exchange_energy)

        class(Abstract_Weighted_Structure), intent(in) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle
        real(DP) :: exchange_energy

        real(DP), dimension(num_dimensions) :: position_div_box
        real(DP), dimension(num_dimensions) :: orientation

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx

        real(DP) :: real_part

        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz

        position_div_box(:) = 2._DP*PI * particle%position(:) / Box%size(:)
        call set_fourier(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call set_fourier(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_fourier(Box%wave(3), position_div_box(3), exp_Ikx_3)

        orientation(:) = exchange_sign(particle%add) * particle%orientation(:)

        exchange_energy = 0._DP

        do kz = 0, Box%wave(3)
            wave_vector(3) = 2._DP*PI * real(kz, DP) / Box%size(3)

            do ky = -reciprocal_size_2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

                do kx = -reciprocal_size_1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                    wave_dot_orientation = dot_product(wave_vector, orientation)
                    exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)

                    real_part = wave_dot_orientation * (wave_dot_orientation + 2._DP * &
                                real(this%structure(kx, ky, kz) * conjg(exp_Ikx), DP))

                    exchange_energy = exchange_energy + this%weight(kx, ky, kz) * real_part

                end do

            end do

        end do

        exchange_energy = 4._DP*PI / product(Box%size) * exchange_energy

    end function Abstract_Weighted_Structure_exchange_energy

    !> Exchange a particle -> update the ``structure factor''

    !> Add particle
    !>  \f[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot+\vec{\mu}_{N+1}) e^{+i\vec{k}\cdot\vec{x}_{N+1}}
    !>  \f]
    !>

    !> Remove particle -> update the ``structure factor''
    !>  \f[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot-\vec{\mu}_{N}) e^{+i\vec{k}\cdot\vec{x}_{N}}
    !>  \f]
    !>

    pure subroutine Abstract_Weighted_Structure_update_structure_exchange(this, Box, particle)

        class(Abstract_Weighted_Structure), intent(inout) :: this
        type(Box_Parameters), intent(in) :: Box
        type(Particle_Index), intent(in) :: particle

        real(DP), dimension(num_dimensions) :: position_div_box
        real(DP), dimension(num_dimensions) :: orientation

        complex(DP), dimension(-Box%wave(1):Box%wave(1)) :: exp_Ikx_1
        complex(DP), dimension(-Box%wave(2):Box%wave(2)) :: exp_Ikx_2
        complex(DP), dimension(-Box%wave(3):Box%wave(3)) :: exp_Ikx_3
        complex(DP) :: exp_Ikx

        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_orientation
        integer :: kx, ky, kz

        position_div_box(:) = 2._DP*PI * particle%position(:) / Box%size(:)
        call set_fourier(Box%wave(1), position_div_box(1), exp_Ikx_1)
        call set_fourier(Box%wave(2), position_div_box(2), exp_Ikx_2)
        call set_fourier(Box%wave(3), position_div_box(3), exp_Ikx_3)

        orientation(:) = exchange_sign(particle%add) * particle%orientation(:)

        do kz = 0, Box%wave(3)
            wave_vector(3) = 2._DP*PI * real(kz, DP) / Box%size(3)

            do ky = -reciprocal_size_2_sym(Box%wave, kz), Box%wave(2)
                wave_vector(2) = 2._DP*PI * real(ky, DP) / Box%size(2)

                do kx = -reciprocal_size_1_sym(Box%wave, ky, kz), Box%wave(1)
                    wave_vector(1) = 2._DP*PI * real(kx, DP) / Box%size(1)

                    wave_dot_orientation = dot_product(wave_vector, orientation)
                    exp_Ikx = exp_Ikx_1(kx) * exp_Ikx_2(ky) * exp_Ikx_3(kz)

                    this%structure(kx, ky, kz) = this%structure(kx, ky, kz) + &
                                                 cmplx(wave_dot_orientation, 0._DP, DP) * exp_Ikx

                end do

            end do

        end do

    end subroutine Abstract_Weighted_Structure_update_structure_exchange
