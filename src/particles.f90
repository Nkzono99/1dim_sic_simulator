!> @brief 粒子に対する処理モジュール.
module particles
    use parameters
    use commons
    use boundary
    use utils
    implicit none

    private

    public particles_add_particle
    public particles_add_uniformly, particles_add_randomly
    public particles_sort

contains

    !> @brief 粒子を空間に一様に配置する.
    !>
    !> @param[in] npcl_add 追加する粒子の個数
    !> @param[in] ispec 粒子の種類
    subroutine particles_add_uniformly(npcl_add, ispec)
        integer, intent(in) :: npcl_add, ispec
        integer :: i, j
        real(8) :: x, vx
        integer :: npcl_prev

        ! 粒子を追加する前に存在していた粒子の個数を覚えておく
        npcl_prev = npcl(ispec)

        ! 各グリッドに同じ数だけ粒子を配置する.
        do i = 1, ngrid
            do j = 1, int((npcl_add + i - 1)/ngrid)
                call random_number(x)
                x = x*dx + (i - 1.5)*dx

                ! 粒子の初期速度を正規分布に従うように設定
                call utils_rand_bm(vx)

                call particles_add_particle(ispec, x, vx, 0)
            end do
        end do

        call sort_for_species(ispec, npcl_prev+1, npcl(ispec))
        call boundary_correct_pcl
    end subroutine

    !> @brief 粒子をランダムに配置する.
    !>
    !> @param[in] npcl_add 追加する粒子の個数
    !> @param[in] ispec 粒子の種類
    !> @param[in] start_x 粒子を追加する範囲の初めの位置
    !> @param[in] end_x 粒子を追加する範囲の最後の位置
    subroutine particles_add_randomly(npcl_add, ispec, start_x, end_x)
        integer, intent(in) :: npcl_add, ispec
        real(8), intent(in) :: start_x, end_x
        integer :: ipcl
        real(8) :: x, vx
        integer :: npcl_prev

        ! 粒子を追加する前に存在していた粒子の個数を覚えておく
        npcl_prev = npcl(ispec)

        ! ランダムに粒子を配置する
        do ipcl = 1, npcl_add
            call random_number(x)
            x = start_x + (end_x - start_x) * x
            
            ! 粒子の初期速度を正規分布に従うように設定
            call utils_rand_bm(vx)
            call particles_add_particle(ispec, x, vx, 0)
        end do

        call sort_for_species(ispec, npcl_prev+1, npcl(ispec))
        call boundary_correct_pcl
    end subroutine

    !> @brief 粒子を一つ追加する.
    !>
    !> @param[in] ispec 粒子の種類
    !> @param[in] x 追加する位置
    !> @param[in] vx 設定する速度
    !> @param[in] ncycle 設定する周期境界をまたいだ数
    subroutine particles_add_particle(ispec, x, vx, ncycle)
        integer, intent(in) :: ispec
        real(8), intent(in) :: x, vx
        integer, intent(in) :: ncycle

        if (npcl(ispec) + 1 > max_npcl) then
            return
        end if

        npcl(ispec) = npcl(ispec) + 1
        px(npcl(ispec), ispec) = x
        pvx(npcl(ispec), ispec) = vx
        ncycles(npcl(ispec), ispec) = ncycle
    end subroutine

    !> @brief !> @brief 粒子位置で粒子のインデックスをソートする.
    subroutine particles_sort
        integer ispec

        do ispec = 1, nspec
            call sort_for_species(ispec, 1, npcl(ispec))
        end do
    end subroutine

    !> @brief 粒子位置で粒子のインデックスをソートするためのヘルパ関数.
    !>
    !> @param[in] ispec 粒子の種類
    !> @param[in] left 左インデックス
    !> @param[in] right 右インデックス
    recursive subroutine sort_for_species(ispec, left, right)
        integer, intent(in) :: ispec, left, right
        integer :: i, j
        real(8) :: pivot

        if (left >= right) then
            return
        end if

        pivot = 0.5*(px(left, ispec) + px(right, ispec))

        i = left
        j = right
        do while (.true.)
            do while (px(i, ispec) < pivot)
                i = i + 1
            end do
            do while (pivot < px(j, ispec))
                j = j - 1
            end do

            if (i >= j) then
                exit
            end if

            call swap_particles(ispec, i, j)
            i = i + 1
            j = j - 1
        end do

        call sort_for_species(ispec, left, i - 1)
        call sort_for_species(ispec, j + 1, right)
    end subroutine

    !> 粒子インデックスを入れ替える.
    !>
    !> @param[in] ispec 粒子の種類
    !> @param[in] ipcl1 1つ目の粒子インデックス
    !> @param[in] ipcl2 2つ目の粒子インデックス
    subroutine swap_particles(ispec, ipcl1, ipcl2)
        integer ispec, ipcl1, ipcl2
        real(8) px_buf, pvx_buf

        px_buf = px(ipcl1, ispec)
        pvx_buf = pvx(ipcl1, ispec)

        px(ipcl1, ispec) = px(ipcl2, ispec)
        pvx(ipcl1, ispec) = pvx(ipcl2, ispec)

        px(ipcl2, ispec) = px_buf
        pvx(ipcl2, ispec) = pvx_buf
    end subroutine

end module
