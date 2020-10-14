!> @brief 現在の状態の出力処理モジュール.
module status
    use parameters
    use commons
    implicit none

    private

    public status_start, status_write, status_close

contains

    !> @brief 出力ファイルを作成する.
    !>
    !> @param[in] output_dir 出力ディレクトリ
    subroutine status_start(output_dir)
        character(:), allocatable, intent(in) :: output_dir
        character(6) str
        integer ispec

        do ispec = 1, nspec
            write (str, '(i0)') ispec
            open (200 + ispec, file=output_dir//'/rho'//trim(str)//'.csv', status='replace')
        end do
        open (103, file=output_dir//'/phi.csv', status='replace')
        open (104, file=output_dir//'/ex.csv', status='replace')

        if (output_electrostatic_energy) open (105, file=output_dir//'/es_energy.csv', status='replace')
        if (output_kinetic_energy) open (106, file=output_dir//'/kinetic_energy.csv', status='replace')
        if (output_distance_between_tracers) open (107, file=output_dir//'distance_between_tracers.csv', status='replace')
        if (output_npcl) open (108, file=output_dir//'npcl.csv', status='replace')
    end subroutine

    !> @brief 現在の状況をファイルに書き込む.
    subroutine status_write(istep)
        integer, intent(in) :: istep
        integer ispec
        if (mod(istep - 1, output_skips) /= 0) then
            return
        end if

        do ispec = 1, nspec
            write (200 + ispec, '(*(F20.6, :, ","))') rhospec(1:ngrid, ispec)
        end do
        write (103, '(*(F20.6, :, ","))') phi(1:ngrid)
        write (104, '(*(F20.6, :, ","))') ex(1:ngrid)

        if (output_electrostatic_energy) call write_electrostatic_energy(105)
        if (output_kinetic_energy) call write_kinetic_energy(106)
        if (output_distance_between_tracers) call write_distance_between_tracers(107)
        if (output_npcl) call write_npcl(108)
    end subroutine

    !> @brief 出力を終了する.
    subroutine status_close
        integer ispec

        do ispec = 1, nspec
            close (200 + ispec)
        end do
        close (103)
        close (104)

        if (output_electrostatic_energy) close (105)
        if (output_kinetic_energy) close (106)
        if (output_distance_between_tracers) close (107)
        if (output_npcl) close (108)
    end subroutine

    !> @brief 静電エネルギーの総和をファイルに書き込む.
    !>
    !> @param[in] output_number 出力番号
    subroutine write_electrostatic_energy(output_number)
        integer, intent(in) :: output_number
        integer i
        real(8) energy

        energy = 0
        do i = 1, ngrid
            energy = energy + ex(i)*ex(i)
        end do
        energy = 0.5*eps0*energy*dx

        write (output_number, '(E20.6)') energy
    end subroutine

    !> @brief 粒子の運動エネルギーの総和をファイルに書き込む.
    !>
    !> @param[in] output_number 出力番号
    subroutine write_kinetic_energy(output_number)
        integer, intent(in) :: output_number
        integer ispec, isimp, ipcl1, ipcl2
        real(8) energy, mean_pvx, rate

        energy = 0
        do ispec = 1, nspec
            do isimp = 1, nsimp(ispec)
                ipcl1 = simplices(isimp, ispec)%ipcl1
                ipcl2 = simplices(isimp, ispec)%ipcl2
                rate = simplices(isimp, ispec)%rate
                mean_pvx = 0.5*(pvx(ipcl1, ispec) + pvx(ipcl2, ispec))
                energy = energy + 0.5*rate*ms(ispec)*mean_pvx*mean_pvx
            end do
        end do

        write (output_number, '(E20.6)') energy
    end subroutine

    !> @brief トレーサー間の平均距離をファイルに書き込む
    !>
    !> @param[in] output_number 出力番号
    subroutine write_distance_between_tracers(output_number)
        integer, intent(in) :: output_number
        integer :: ispec, isimp, ipcl1, ipcl2
        real(8) :: offset1, offset2, px1, px2
        real(8) :: mean_dist(nspec)

        do ispec = 1, nspec
            mean_dist(ispec) = 0

            do isimp = 1, nsimp(ispec)
                ipcl1 = simplices(isimp, ispec)%ipcl1
                ipcl2 = simplices(isimp, ispec)%ipcl2
                offset1 = simplices(isimp, ispec)%offset1
                offset2 = simplices(isimp, ispec)%offset2
                px1 = px(ipcl1, ispec) + ngrid*dx*ncycles(ipcl1, ispec) + offset1
                px2 = px(ipcl2, ispec) + ngrid*dx*ncycles(ipcl2, ispec) + offset2
                mean_dist(ispec) = mean_dist(ispec) + abs(px1 - px2)
            end do

            mean_dist(ispec) = mean_dist(ispec)/npcl(ispec)
        end do

        write (output_number, '(*(F20.6, :, ","))') mean_dist
    end subroutine

    !> @brief 粒子の個数をファイルに書き込む.
    !>
    !> @param[in] output_number 出力番号
    subroutine write_npcl(output_number)
        integer, intent(in) :: output_number

        write (output_number, '(*(I16, :, ","))') npcl(1:nspec)
    end subroutine

end module
