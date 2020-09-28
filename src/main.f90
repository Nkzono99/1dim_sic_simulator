!> 1次元静電SICシミュレータ.
program main
    use parameters
    use commons
    use particles
    use efield
    use status
    use argparse
    use progressbar
    use sic
    use utils

    implicit none

    integer istep

    ! コマンドライン引数の処理
    call argparse_init
    call argparse_add('inputfilename', description='parameter namelist filename (*.in)')
    call argparse_add('--output', subflagname='-o', default='.', description='data output directory')
    call argparse_parse

    call parameters_init(argparse_get('inputfilename'))
    call commons_init
    call utils_set_random_seed(random_seed)

    call sic_init

    call sic_scatter_on_grid

    call efield_update

    call show_simulation_settings

    call status_start(argparse_get('output'))

    do istep = 1, nsteps
        if (mod(istep, 100) == 0) then
            call progressbar_show(istep, nsteps, n=50)
        end if

        call sic_update

        ! if (mod(istep, 50) == 0) then
        !     call sic_correct_temprature
        ! end if
!
!         if (mod(istep, 1000) == 0) then
!             call particles_sort
!         end if

        call sic_scatter_on_grid

        call efield_update

        if (refinement_interval /= 0 .and. mod(istep, refinement_interval) == 0) then
            call sic_refinement(refinement_threshold*dx)
        end if

        ! call show_temprature
        ! call check_rho_and_phi
        ! print *, nsimp(1:nspec)

        if (mod(istep - 1, output_skips) == 0) then
            call status_write
        end if
    end do

    call status_close

contains

    !> @brief シミュレーション設定を出力する.
    subroutine show_simulation_settings
        print *, '---- Simulation Parameter --------'
        print *, 'total steps =', nsteps
        print *, 'initial number of super particles', npcl_init(1:nspec)
        print *, 'max number of super particles', max_npcl
        print *, 'number of grids =', ngrid
        print *, ''
        print *, 'dt =', dt
        print *, 'dx =', dx
        print *, 'npcl_per_super =', npcl_per_super
        print *, '----------------------------------'
        print *, '---- Plasma Parameter ------------'
        print *, 'debye length =', lambda
        print *, 'wpe =', 1/lambda*sqrt(kB*Ts(1)/me)
        print *, '----------------------------------'
        print *, '---- Parallel Parameter ----------'
        print *, 'number of thread =', nthreads
        print *, '----------------------------------'
        print *, '---- Refinement Parameter --------'
        print *, 'refinement interval =', refinement_interval
        print *, 'refinement threthold =', refinement_threshold
        print *, '----------------------------------'
    end subroutine

    !> @brief 電位が正しく計算されているかチェックする.
    subroutine check_rho_and_phi
        integer i
        real(8) error, tmp_rho
        error = 0

        do i = 1, ngrid
            tmp_rho = -eps0*(phi(i + 1) - 2*phi(i) + phi(i - 1))/(dx*dx)
            error = error + abs(rho(i) - tmp_rho)
        end do

        print *, istep, ': Error =', error/ngrid, 'max =', maxval(abs(rho))
    end subroutine

    !> 現在の粒子温度を出力する.
    subroutine show_temprature
        integer :: ispec, isimp, ipcl1, ipcl2
        real(8) :: mean_pvx
        real(8) :: Ts_current
        real(8) :: energy

        do ispec = 1, nspec
            if (Ts(ispec) == 0.0) then
                pvx(:, ispec) = 0.0
                cycle
            end if

            ! 各粒子の運動エネルギーの和を計算する.
            energy = 0.0
            do isimp = 1, nsimp(ispec)
                ipcl1 = simplices(isimp, ispec)%ipcl1
                ipcl2 = simplices(isimp, ispec)%ipcl2
                mean_pvx = 0.5*(pvx(ipcl1, ispec) + pvx(ipcl2, ispec))

                energy = energy + 0.5d0*ms(ispec)*mean_pvx*mean_pvx
            end do

            ! 3/2kBT = 1/2mv^2より、Tを求める.
            Ts_current = 2.0d0/kB*energy/(npcl(ispec)*npcl_per_super)

            print *, ispec, Ts_current
        end do
    end subroutine

end program
