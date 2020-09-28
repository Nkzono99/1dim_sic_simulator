!> simplexに対する処理モジュール.
module sic
    use parameters
    use commons
    use utils
    use particles
    use boundary

    implicit none

    private

    public sic_init
    public sic_correct_temprature
    public sic_refinement
    public sic_scatter_on_grid
    public sic_update

contains

    !> simplexを初期化する.
    subroutine sic_init
        integer ispec

        do ispec = 1, nspec
            call particles_distribute(npcl_init(ispec), ispec, x_mode='uniform', vx_mode='normal')
        end do

        do ispec = 1, nspec
            call create_mesh(ispec, 1, npcl(ispec), 1.0d0, periodic=.true.)
        end do

        call sic_correct_temprature
    end subroutine

    !> 粒子メッシュを作成する.
    !>
    !> @param[in] ispec 粒子の種類
    !> @param[in] start_ipcl 初めの粒子インデックス
    !> @param[in] end_ipcl 最後の粒子インデックス
    !> @param[in] rate 粒子の割合
    !> @param[in] periodic 周期的なメッシュを作成する場合.true.
    subroutine create_mesh(ispec, start_ipcl, end_ipcl, rate, periodic)
        integer, intent(in) :: ispec, start_ipcl, end_ipcl
        real(8), intent(in) :: rate
        logical, intent(in) :: periodic
        integer :: ipcl

        do ipcl = start_ipcl, end_ipcl - 1
            call sic_add_simplex(ispec, ipcl, ipcl + 1, rate)
        end do

        if (periodic) then
            call sic_add_simplex(ispec, start_ipcl, end_ipcl, rate, offset1=dx*ngrid)
        end if
    end subroutine

    !> @brief 新しいsimplexを追加する.
    !>
    !> @param[in] ispec 粒子の種類
    !> @param[in] ipcl1 1つ目の粒子のインデックス
    !> @param[in] ipcl2 2つ目の粒子のインデックス
    !> @param[in] offset1 1つ目の粒子のオフセット(default: 0)
    !> @param[in] offset2 2つ目の粒子のオフセット(default: 0)
    !> @param[in] isimp simplexを上書きする場合はそのインデックスを指定する(default:追加)
    subroutine sic_add_simplex(ispec, ipcl1, ipcl2, rate, offset1, offset2, isimp)
        integer, intent(in) :: ispec
        integer, intent(in) :: ipcl1, ipcl2
        real(8), intent(in) :: rate
        real(8), intent(in), optional :: offset1, offset2
        integer, intent(in), optional :: isimp
        integer :: isimp_

        if (present(isimp)) then
            isimp_ = isimp
        else
            ! 最大数に達していた場合終了
            if (nsimp(ispec) >= max_npcl) then
                return
            end if

            ! simplexの個数をインクリメント
            nsimp(ispec) = nsimp(ispec) + 1

            isimp_ = nsimp(ispec)
        end if

        ! パラメータの設定
        simplices(isimp_, ispec)%ipcl1 = ipcl1
        simplices(isimp_, ispec)%ipcl2 = ipcl2
        simplices(isimp_, ispec)%rate = rate

        if (present(offset1)) then
            simplices(isimp_, ispec)%offset1 = offset1
        else
            simplices(isimp_, ispec)%offset1 = 0
        end if

        if (present(offset2)) then
            simplices(isimp_, ispec)%offset2 = offset2
        else
            simplices(isimp_, ispec)%offset2 = 0
        end if

        ! particleからsimplexへのリンクを作る
        if (pcl2simp(ipcl1, ispec, 1) == -1) then
            pcl2simp(ipcl1, ispec, 1) = isimp_
        else
            pcl2simp(ipcl1, ispec, 2) = isimp_
        end if

        if (pcl2simp(ipcl2, ispec, 1) == -1) then
            pcl2simp(ipcl2, ispec, 1) = isimp_
        else
            pcl2simp(ipcl2, ispec, 2) = isimp_
        end if
    end subroutine

    !> @brief Refinementを適用する.
    !>
    !> @details トレーサー間がしきい値より離れた場合そのsimplexをしきい値未満になるまで分割する.
    !>
    !> @param[in] threshold しきい値
    subroutine sic_refinement(threshold)
        real(8), intent(in) :: threshold
        integer :: isimp, ispec
        integer :: ipcl1, ipcl2
        real(8) :: offset1, offset2
        real(8) :: px1, px2

        do ispec = 1, nspec
            isimp = 1
            do while (isimp <= nsimp(ispec))
                if (nsimp(ispec) >= max_npcl) then
                    exit
                end if
                ipcl1 = simplices(isimp, ispec)%ipcl1
                ipcl2 = simplices(isimp, ispec)%ipcl2
                offset1 = simplices(isimp, ispec)%offset1
                offset2 = simplices(isimp, ispec)%offset2
                px1 = px(ipcl1, ispec) + ngrid*dx*ncycles(ipcl1, ispec) + offset1
                px2 = px(ipcl2, ispec) + ngrid*dx*ncycles(ipcl2, ispec) + offset2

                if (abs(px2 - px1) > threshold) then
                    call sic_split_simplex(isimp, ispec)

                    ! まだsimplexを追加できるなら分割後のsimplexをチェックするためcontinue
                    continue
                end if

                isimp = isimp + 1
            end do
        end do
    end subroutine

    !> @brief simplexを分割する.
    !>
    !> @param[in] isimp simplexのインデックス
    !> @param[in] ispec 粒子の種類
    subroutine sic_split_simplex(isimp, ispec)
        integer, intent(in) :: isimp, ispec
        integer :: ipcl1, ipcl2
        real(8) :: rate, offset1, offset2
        real(8) :: px1, px2, pvx1, pvx2
        real(8) :: px_new, pvx_new
        integer :: ncycle_new

        ! simplex数が最大数に達していたら終了
        if (nsimp(ispec) >= max_npcl) then
            return
        end if

        ipcl1 = simplices(isimp, ispec)%ipcl1
        ipcl2 = simplices(isimp, ispec)%ipcl2
        rate = simplices(isimp, ispec)%rate
        offset1 = simplices(isimp, ispec)%offset1
        offset2 = simplices(isimp, ispec)%offset2
        px1 = px(ipcl1, ispec) + ngrid*dx*ncycles(ipcl1, ispec) + offset1
        px2 = px(ipcl2, ispec) + ngrid*dx*ncycles(ipcl2, ispec) + offset2
        pvx1 = pvx(ipcl1, ispec)
        pvx2 = pvx(ipcl2, ispec)

        ! 中点に新しいトレーサーを追加
        ! 位置計算
        px_new = 0.5*(px1 + px2)
        ncycle_new = int(0.5*(px1 + px2)/(dx*ngrid))
        px_new = px_new - ncycle_new*dx*ngrid
        ! 速度計算
        pvx_new = 0.5*(pvx1 + pvx2)
        call particles_add_particle(ispec, px_new, pvx_new, ncycle_new)

        ! particleからsimplexへのリンクを削除
        if (pcl2simp(ipcl1, ispec, 1) == isimp) then
            pcl2simp(ipcl1, ispec, 1) = -1
        else
            pcl2simp(ipcl1, ispec, 2) = -1
        end if

        if (pcl2simp(ipcl2, ispec, 1) == isimp) then
            pcl2simp(ipcl2, ispec, 1) = -1
        else
            pcl2simp(ipcl2, ispec, 2) = -1
        end if

        ! 左側のsimplexを作成
        call sic_add_simplex(ispec, ipcl1, npcl(ispec), 0.5*rate, offset1=offset1, isimp=isimp)

        ! 右側のsimplexを作成
        call sic_add_simplex(ispec, npcl(ispec), ipcl2, 0.5*rate, offset2=offset2)
    end subroutine

    !> @brief 粒子電荷をグリッドに分配する.
    subroutine sic_scatter_on_grid
        integer :: ispec, i

        rho(0:ngrid + 1) = 0
        rhospec(0:ngrid + 1, 1:nspec) = 0

        call sic_scatter_b0

        ! apply period boundary
        if (boundary_type == 0) then
            do ispec = 1, nspec
                rhospec(1, ispec) = rhospec(1, ispec) + rhospec(ngrid + 1, ispec)
                rhospec(ngrid, ispec) = rhospec(ngrid, ispec) + rhospec(0, ispec)
            end do
        end if

        do ispec = 1, nspec
            do i = 1, ngrid
                rho(i) = rho(i) + qs(ispec)*rhospec(i, ispec)
            end do
        end do

        call boundary_correct_rho
    end subroutine

    !> @brief 粒子の電荷を分配する.
    subroutine sic_scatter_b0
        integer :: ispec, isimp
        integer :: i

        ! simplex構造体の値を格納する一時変数
        integer :: ipcl1, ipcl2
        real(8) :: rate, offset1, offset2

        integer :: il, ir  ! 粒子の所属するグリッド番号
        real(8) :: dl, dr  ! 粒子の所属するグリッドでの相対位置 (0.0 ~ 1.0)

        integer :: dcycle  ! シミュレーション空間を何周しているか

        real(8) :: px1, px2  ! 粒子の実際の位置 (-inf ~ +inf)
        real(8) :: lpx1, lpx2  ! 粒子のシミュレーション空間での位置 (0.0 ~ ngrid*dx)

        real(8) :: dpx  ! 粒子の実際の距離 (0.0 ~ +inf)
        real(8) :: ldpx  ! 粒子のシミュレーション空間での距離 (0.0 ~ ngrid*dx)

        ! グリッド全体に足す値を格納する変数
        real(8) :: all_sum

        do ispec = 1, nspec
            all_sum = 0
            do isimp = 1, nsimp(ispec)
                ipcl1 = simplices(isimp, ispec)%ipcl1
                ipcl2 = simplices(isimp, ispec)%ipcl2
                rate = simplices(isimp, ispec)%rate
                offset1 = simplices(isimp, ispec)%offset1
                offset2 = simplices(isimp, ispec)%offset2

                px1 = px(ipcl1, ispec) + ngrid*dx*ncycles(ipcl1, ispec) + offset1
                px2 = px(ipcl2, ispec) + ngrid*dx*ncycles(ipcl2, ispec) + offset2
                lpx1 = utils_pmod(px1, ngrid*dx)
                lpx2 = utils_pmod(px2, ngrid*dx)

                dpx = abs(px2 - px1)
                ldpx = abs(lpx2 - lpx1)

                dcycle = int(dpx/ngrid*dx)
                all_sum = all_sum + dcycle*dx/dpx*rate

                if (lpx1 < lpx2) then
                    il = int(lpx1/dx + 0.5)
                    dl = (lpx1/dx + 0.5) - il
                    ir = int(lpx2/dx + 0.5)
                    dr = (lpx2/dx + 0.5) - ir
                else
                    il = int(lpx2/dx + 0.5)
                    dl = (lpx2/dx + 0.5) - il
                    ir = int(lpx1/dx + 0.5)
                    dr = (lpx1/dx + 0.5) - ir
                end if

                if ((lpx2 - lpx1)*(px2 - px1) >= 0) then
                    ! 以下のようなsimplexの場合(ここで足し込むのは+マーク部分)
                    ! |                      |
                    ! |  case 1              |
                    ! |    o+++++++++++o     |
                    ! |  case 2              |
                    ! |    o-----------------|
                    ! |-----++++++++++++o    |
                    ! |                      |
                    if (il == ir) then
                        rhospec(il, ispec) = rhospec(il, ispec) + rate
                    else
                        rhospec(il, ispec) = rhospec(il, ispec) + (1 - dl)*dx/dpx*rate
                        do i = il + 1, ir - 1
                            rhospec(i, ispec) = rhospec(i, ispec) + dx/dpx*rate
                        end do
                        rhospec(ir, ispec) = rhospec(ir, ispec) + dr*dx/dpx*rate
                    end if
                else
                    ! 以下のようなsimplexの場合(ここで足し込むのは+マーク部分)
                    ! |                     |
                    ! |  case 1             |
                    ! |++++o           o++++|
                    ! |  case 2             |
                    ! |++++o                |
                    ! |---------------------|
                    ! |                o++++|
                    ! |                     |
                    if (il == 0) then
                        rhospec(0, ispec) = rhospec(0, ispec) + (dl - 0.5)*dx/dpx*rate
                    else
                        rhospec(0, ispec) = rhospec(0, ispec) + 0.5*dx/dpx*rate
                        do i = 1, il - 1
                            rhospec(i, ispec) = rhospec(i, ispec) + dx/dpx*rate
                        end do
                        rhospec(il, ispec) = rhospec(il, ispec) + dl*dx/dpx*rate
                    end if

                    if (ir == ngrid) then
                        rhospec(ngrid, ispec) = rhospec(ngrid, ispec) + (0.5 - dr)*dx/dpx*rate
                    else
                        rhospec(ir, ispec) = rhospec(ir, ispec) + (1 - dr)*dx/dpx*rate
                        do i = ir + 1, ngrid - 1
                            rhospec(i, ispec) = rhospec(i, ispec) + dx/dpx*rate
                        end do
                        rhospec(ngrid, ispec) = rhospec(ngrid, ispec) + 0.5*dx/dpx*rate
                    end if
                end if
            end do

            do i = 1, ngrid
                rhospec(i, ispec) = rhospec(i, ispec) + all_sum
            end do
        end do
    end subroutine

    !> @brief トレーサー位置・速度を更新する.
    subroutine sic_update
        integer :: ispec, ipcl
        integer :: isimp1, isimp2
        real(8) :: ex_p

        do ispec = 1, nspec
            do ipcl = 1, npcl(ispec)
                isimp1 = pcl2simp(ipcl, ispec, 1)
                isimp2 = pcl2simp(ipcl, ispec, 2)

                ex_p = 0
                if (isimp1 /= -1) then
                    ex_p = ex_p + 0.5*calc_ex(isimp1, ispec)
                end if
                if (isimp2 /= -1) then
                    ex_p = ex_p + 0.5*calc_ex(isimp2, ispec)
                end if

                pvx(ipcl, ispec) = pvx(ipcl, ispec) + dt*qs(ispec)/ms(ispec)*ex_p
            end do

            do ipcl = 1, npcl(ispec)
                px(ipcl, ispec) = px(ipcl, ispec) + dt*pvx(ipcl, ispec)
            end do
        end do

        call boundary_correct_pcl
    end subroutine

    !> @brief simplexに及ぼされる電場の大きさを計算する.
    !>
    !> @param[in] isimp simplexのインデックス
    !> @param[in] ispec 粒子の種類
    function calc_ex(isimp, ispec) result(ex_p)
        integer, intent(in) :: isimp, ispec
        real(8) :: ex_p
        integer :: i

        ! simplex構造体の値を格納する一時変数
        integer :: ipcl1, ipcl2
        real(8) :: rate, offset1, offset2

        integer :: il, ir  ! 粒子の所属するグリッド番号
        real(8) :: dl, dr  ! 粒子の所属するグリッドでの相対位置 (0.0 ~ 1.0)

        real(8) :: px1, px2  ! 粒子の実際の位置 (-inf ~ +inf)
        real(8) :: lpx1, lpx2  ! 粒子のシミュレーション空間での位置 (0.0 ~ ngrid*dx)

        real(8) :: dpx  ! 粒子の実際の距離 (0.0 ~ +inf)
        real(8) :: ldpx  ! 粒子のシミュレーション空間での距離 (0.0 ~ ngrid*dx)

        ipcl1 = simplices(isimp, ispec)%ipcl1
        ipcl2 = simplices(isimp, ispec)%ipcl2
        offset1 = simplices(isimp, ispec)%offset1
        offset2 = simplices(isimp, ispec)%offset2

        px1 = px(ipcl1, ispec) + ngrid*dx*ncycles(ipcl1, ispec) + offset1
        px2 = px(ipcl2, ispec) + ngrid*dx*ncycles(ipcl2, ispec) + offset2
        lpx1 = utils_pmod(px1, ngrid*dx)
        lpx2 = utils_pmod(px2, ngrid*dx)

        dpx = abs(px2 - px1)
        ldpx = abs(lpx2 - lpx1)

        if (lpx1 < lpx2) then
            il = int(lpx1/dx + 0.5)
            dl = (lpx1/dx + 0.5) - il
            ir = int(lpx2/dx + 0.5)
            dr = (lpx2/dx + 0.5) - ir
        else
            il = int(lpx2/dx + 0.5)
            dl = (lpx2/dx + 0.5) - il
            ir = int(lpx1/dx + 0.5)
            dr = (lpx1/dx + 0.5) - ir
        end if

        ex_p = 0.0
        if ((lpx2 - lpx1)*(px2 - px1) >= 0) then
            ! 以下のようなsimplexの場合(ここで足し込むのは+マーク部分)
            ! |                      |
            ! |  case 1              |
            ! |    o+++++++++++o     |
            ! |  case 2              |
            ! |    o-----------------|
            ! |-----++++++++++++o    |
            ! |                      |
            if (il == ir) then
                ex_p = ex_p + ex(il)
            else
                ex_p = ex_p + (1 - dl)*dx/dpx*ex(il)
                do i = il + 1, ir - 1
                    ex_p = ex_p + dx/dpx*ex(i)
                end do
                ex_p = ex_p + dr*dx/dpx*ex(ir)
            end if
        else
            ! 以下のようなsimplexの場合(ここで足し込むのは+マーク部分)
            ! |                     |
            ! |  case 1             |
            ! |++++o           o++++|
            ! |  case 2             |
            ! |++++o                |
            ! |---------------------|
            ! |                o++++|
            ! |                     |
            if (il == 0) then
                ex_p = ex_p + (dl - 0.5)*dx/dpx*ex(0)
            else
                ex_p = ex_p + 0.5*dx/dpx*ex(0)
                do i = 1, il - 1
                    ex_p = ex_p + dx/dpx*ex(i)
                end do
                ex_p = ex_p + dl*dx/dpx*ex(il)
            end if

            if (ir == ngrid) then
                ex_p = ex_p + (0.5 - dr)*dx/dpx*rate
            else
                ex_p = ex_p + (1 - dr)*dx/dpx*ex(ir)
                do i = ir + 1, ngrid - 1
                    ex_p = ex_p + dx/dpx*ex(i)
                end do
                ex_p = ex_p + 0.5*dx/dpx*ex(ngrid)
            end if
        end if
    end function

    !> 設定した粒子温度になるように速度を補正する.
    subroutine sic_correct_temprature
        integer :: ispec, isimp, ipcl, ipcl1, ipcl2
        real(8) :: Ts_current, mean_pvx
        real(8) :: energy
        real(8) :: xs

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

            ! 1/2kBT = 1/2mv^2より、Tを求める.
            Ts_current = 2.0/kB*energy/(npcl(ispec)*npcl_per_super)

            ! 温度補正係数を計算する.
            xs = sqrt(Ts(ispec)/Ts_current)

            do ipcl = 1, npcl(ispec)
                pvx(ipcl, ispec) = pvx(ipcl, ispec)*xs
            end do
        end do

        call boundary_correct_pcl
    end subroutine

end module
