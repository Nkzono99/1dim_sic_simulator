!> @brief 粒子に対する処理モジュール.
module particles
    use parameters
    use commons
    use boundary
    use utils
    implicit none

    private

    public particles_add_particle
    public particles_distribute
    public particles_sort

contains

    !> @brief 粒子を分配する.
    !>
    !> @details
    !!> 分配方法
    !!>   - 'uniform' : startからendまでの範囲に一様に分配する
    !!>   - 'random' : startからendまでの範囲にランダムに分配する
    !!>   - 'normal' : 平均mean, 分散stdの正規分布に従うように分配する
    !>
    !> @param[in] npcl_add 追加する粒子数
    !> @param[in] ispec 追加する粒子の種類
    !> @param[in] x_mode 位置の分配方法('uniform', 'random', 'normal')
    !> @param[in] vx_mode 速度の分配方法('uniform', 'random', 'normal')
    !> @param[in] x_start 初めの位置 (default: 0)
    !> @param[in] x_end 最後の位置 (default: ngrid*dx)
    !> @param[in] x_mean 平均位置 (default: 0.5*ngrid*dx)
    !> @param[in] x_std 位置の分散 (default: dx)
    !> @param[in] vx_start 最小速度 (default: -1)
    !> @param[in] vx_end 最大速度 (default: 1)
    !> @param[in] vx_mean 平均速度 (default: 0)
    !> @param[in] vx_std 速度の分散 (default: 1)
    subroutine particles_distribute(npcl_add, ispec, x_mode, vx_mode, &
                                    x_start, x_end, x_mean, x_std, &
                                    vx_start, vx_end, vx_mean, vx_std)
        integer, intent(in) :: npcl_add
        integer, intent(in) :: ispec
        character(*), intent(in) :: x_mode, vx_mode
        real(8), intent(in), optional :: x_start, x_end
        real(8), intent(in), optional :: x_mean, x_std
        real(8), intent(in), optional :: vx_start, vx_end
        real(8), intent(in), optional :: vx_mean, vx_std

        integer ipcl
        integer :: npcl_prev
        real(8) :: x, vx
        real(8) :: x_start_, x_end_
        real(8) :: x_mean_, x_std_
        real(8) :: vx_start_, vx_end_
        real(8) :: vx_mean_, vx_std_

        if (present(x_start)) then
            x_start_ = x_start
        else
            x_start_ = 0
        end if

        if (present(x_end)) then
            x_end_ = x_end
        else
            x_end_ = ngrid * dx
        end if

        if (present(x_mean)) then
            x_mean_ = x_mean
        else
            x_mean_ = 0.5 * ngrid * dx
        end if

        if (present(x_std)) then
            x_std_ = x_std
        else
            x_std_ = dx
        end if

        if (present(vx_start)) then
            vx_start_ = vx_start
        else
            vx_start_ = -1
        end if

        if (present(vx_end)) then
            vx_end_ = vx_end
        else
            vx_end_ = 1
        end if

        if (present(vx_mean)) then
            vx_mean_ = vx_mean
        else
            vx_mean_ = 0
        end if

        if (present(vx_std)) then
            vx_std_ = vx_std
        else
            vx_std_ = 1
        end if

        ! 粒子を追加する前に存在していた粒子の個数を覚えておく
        npcl_prev = npcl(ispec)

        do ipcl = 1, npcl_add
            if (x_mode == 'uniform') then
                x = uniform(ipcl, npcl_add, x_start_, x_end_)
            else if (x_mode == 'random') then
                x = random(x_start_, x_end_)
            else if (x_mode == 'normal') then
                x = normal(x_mean_, x_std_)
            end if

            if (vx_mode == 'uniform') then
                vx = uniform(ipcl, npcl_add, vx_start_, vx_end_)
            else if (vx_mode == 'random') then
                vx = random(vx_start_, vx_end_)
            else if (vx_mode == 'normal') then
                vx = normal(vx_mean_, vx_std_)
            end if

            call particles_add_particle(ispec, x, vx, 0)
        end do
        
        ! 追加した粒子について位置でソートする
        call sort_for_species(ispec, npcl_prev + 1, npcl(ispec))
        call boundary_correct_pcl
    end subroutine

    !> @brief 一様な分布に従う値を返す.
    !> 
    !> @param[in] i 番号
    !> @param[in] n 最大番号
    !> @param[in] min_val 最小値
    !> @param[in] max_val 最大値
    function uniform(i, n, min_val, max_val)
        integer, intent(in) :: i, n
        real(8), intent(in) :: min_val, max_val
        real(8) :: uniform

        uniform = min_val + (max_val - min_val)*i/n
    end function

    !> @brief ランダムな分布に従う値を返す.
    !> 
    !> @param[in] min_val 最小値
    !> @param[in] max_val 最大値
    function random(min_val, max_val)
        real(8), intent(in), optional :: min_val, max_val
        real(8) :: random

        call random_number(random)
        random = min_val + (max_val - min_val)*random
    end function

    !> @brief 正規分布に従う値を返す
    !>
    !> @param[in] val_mean 平均値
    !> @param[in] val_std 分散
    function normal(val_mean, val_std)
        real(8), intent(in), optional :: val_mean, val_std
        real(8) :: normal

        call utils_rand_bm(normal)
        normal = val_mean + normal*val_std
    end function

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

    !> @brief 粒子位置で粒子のインデックスをソートする.
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
