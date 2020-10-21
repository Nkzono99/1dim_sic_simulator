!> @brief 粒子に対する処理モジュール.
module particles
    use parameters
    use commons
    use boundary
    use utils
    implicit none

    private

    public particles_add_particle
    public particles_delete_particle
    public particles_distribute
    public particles_sort

contains

    !> @brief 粒子を一つ追加する.
    !>
    !> @param[in] ispec 粒子の種類
    !> @param[in] x 追加する位置
    !> @param[in] vx 設定する速度
    !> @param[in] ncycle 設定する周期境界をまたいだ数
    !> @param[out] status 追加に成功したら.true.
    subroutine particles_add_particle(ispec, x, vx, ncycle, status)
        integer, intent(in) :: ispec
        real(8), intent(in) :: x, vx
        integer, intent(in) :: ncycle
        logical, intent(out), optional :: status

        if (npcl(ispec) + 1 > max_npcl) then
            if (present(status)) status = .false.
            return
        end if

        npcl(ispec) = npcl(ispec) + 1
        px(npcl(ispec), ispec) = x
        pvx(npcl(ispec), ispec) = vx
        ncycles(npcl(ispec), ispec) = ncycle
        pcl2simp(npcl(ispec), ispec, 1) = -1
        pcl2simp(npcl(ispec), ispec, 2) = -1

        if (present(status)) status = .true.
    end subroutine

    !> @brief 粒子を削除し、削除した位置に最後の粒子を移動する.
    !>
    !> @param[in] ispec 粒子の種類
    !> @param[in] ipcl 粒子番号
    subroutine particles_delete_particle(ispec, ipcl)
        integer, intent(in) :: ispec
        integer, intent(in) :: ipcl
        integer :: isimp1, isimp2

        px(ipcl, ispec) = px(npcl(ispec), ispec)
        pvx(ipcl, ispec) = pvx(npcl(ispec), ispec)
        ncycles(ipcl, ispec) = ncycles(npcl(ispec), ispec)
        pcl2simp(ipcl, ispec, 1) = pcl2simp(npcl(ispec), ispec, 1)
        pcl2simp(ipcl, ispec, 2) = pcl2simp(npcl(ispec), ispec, 2)

        isimp1 = pcl2simp(npcl(ispec), ispec, 1)
        isimp2 = pcl2simp(npcl(ispec), ispec, 2)
        if (isimp1 /= -1) then
            if (simplices(isimp1, ispec)%ipcl1 == npcl(ispec)) then
                simplices(isimp1, ispec)%ipcl1 = ipcl
            else
                simplices(isimp1, ispec)%ipcl2 = ipcl
            end if
        end if

        if (isimp2 /= -1) then
            if (simplices(isimp2, ispec)%ipcl1 == npcl(ispec)) then
                simplices(isimp2, ispec)%ipcl1 = ipcl
            else
                simplices(isimp2, ispec)%ipcl2 = ipcl
            end if
        end if

        npcl(ispec) = npcl(ispec) - 1
    end subroutine

    !> @brief 粒子を分配する.
    !>
    !> @details
    !!> 分配方法(x)
    !!>   - 'uniform' : startからendまでの範囲に一様に分配する
    !!>         x1 : 初めの位置
    !!>         x2 : 最後の位置
    !!>   - 'random' : startからendまでの範囲にランダムに分配する
    !!>         x1 : 初めの位置
    !!>         x2 : 最後の位置
    !!>   - 'normal' : 平均mean, 分散stdの正規分布に従うように分配する
    !!>         x1 : 平均位置
    !!>         x2 : 位置の分散
    !!>
    !!> 分配方法(vx)
    !!>   - 'uniform' : startからendまでの範囲に一様に分配する
    !!>         vx1 : 初めの位置
    !!>         vx2 : 最後の位置
    !!>   - 'random' : startからendまでの範囲にランダムに分配する
    !!>         vx1 : 初めの位置
    !!>         vx2 : 最後の位置
    !!>   - 'normal' : 平均mean, 分散stdの正規分布に従うように分配する
    !!>         vx1 : 平均位置
    !!>         vx2 : 位置の分散
    !>
    !> @param[in] npcl_add 追加する粒子数
    !> @param[in] ispec 追加する粒子の種類
    !> @param[in] x_mode 位置の分配方法('uniform', 'random', 'normal')
    !> @param[in] vx_mode 速度の分配方法('uniform', 'random', 'normal')
    !> @param[in] x1 位置パラメータ1(詳細は分配方法(x)を参照)
    !> @param[in] x2 位置パラメータ2(詳細は分配方法(x)を参照)
    !> @param[in] vx1 速度パラメータ1(詳細は分配方法(vx)を参照)
    !> @param[in] vx2 速度パラメータ2(詳細は分配方法(vx)を参照)
    subroutine particles_distribute(npcl_add, ispec, x_mode, vx_mode, x1, x2, vx1, vx2)
        integer, intent(in) :: npcl_add
        integer, intent(in) :: ispec
        character(*), intent(in) :: x_mode, vx_mode
        real(8), intent(in), optional :: x1, x2
        real(8), intent(in), optional :: vx1, vx2

        integer ipcl
        integer :: npcl_prev
        real(8) :: x, vx

        ! 粒子を追加する前に存在していた粒子の個数を覚えておく
        npcl_prev = npcl(ispec)

        do ipcl = 1, npcl_add
            if (x_mode == 'uniform') then
                x = uniform(ipcl, npcl_add, x1, x2)
            else if (x_mode == 'random') then
                x = random(x1, x2)
            else if (x_mode == 'normal') then
                x = normal(x1, x2)
            end if

            if (vx_mode == 'uniform') then
                vx = uniform(ipcl, npcl_add, vx1, vx2)
            else if (vx_mode == 'random') then
                vx = random(vx1, vx2)
            else if (vx_mode == 'normal') then
                vx = normal(vx1, vx2)
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
