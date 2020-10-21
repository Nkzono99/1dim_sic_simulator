module simplex
    use parameters
    use commons
    use particles
    implicit none

contains
    !> @brief 新しいsimplexを追加する.
    !>
    !> @param[in] ispec 粒子の種類
    !> @param[in] ipcl1 1つ目の粒子のインデックス
    !> @param[in] ipcl2 2つ目の粒子のインデックス
    !> @param[in] offset1 1つ目の粒子のオフセット(default: 0)
    !> @param[in] offset2 2つ目の粒子のオフセット(default: 0)
    !> @param[in] isimp simplexを上書きする場合はそのインデックスを指定する(default:追加)
    subroutine simplex_add_simplex(ispec, ipcl1, ipcl2, rate, offset1, offset2, isimp)
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

    !> @brief simplexを削除し、削除したインデックスに最後のsimplexを移動する.
    !>
    !> @param[in] ispec 粒子の種類
    !> @param[in] isimp 削除するsimplex番号
    subroutine simplex_delete_simplex(ispec, isimp)
        integer, intent(in) :: ispec
        integer, intent(in) :: isimp
        integer ipcl1, ipcl2

        ipcl1 = simplices(isimp, ispec)%ipcl1
        ipcl2 = simplices(isimp, ispec)%ipcl2

        ! particleからsimplexへのリンクがまだ残っている場合、そのリンクを削除する
        if (pcl2simp(ipcl1, ispec, 1) == isimp) then
            pcl2simp(ipcl1, ispec, 1) = -1
        else if (pcl2simp(ipcl1, ispec, 2) == isimp) then
            pcl2simp(ipcl1, ispec, 2) = -1
        end if

        if (pcl2simp(ipcl2, ispec, 1) == isimp) then
            pcl2simp(ipcl2, ispec, 1) = -1
        else if (pcl2simp(ipcl2, ispec, 2) == isimp) then
            pcl2simp(ipcl2, ispec, 2) = -1
        end if

        ! インデックスを移動するためparticleからsimiplexへのリンクを書き換える
        ipcl1 = simplices(nsimp(ispec), ispec)%ipcl1
        ipcl2 = simplices(nsimp(ispec), ispec)%ipcl2
        if (pcl2simp(ipcl1, ispec, 1) == nsimp(ispec)) then
            pcl2simp(ipcl1, ispec, 1) = isimp
        else
            pcl2simp(ipcl1, ispec, 2) = isimp
        end if

        if (pcl2simp(ipcl2, ispec, 1) == nsimp(ispec)) then
            pcl2simp(ipcl2, ispec, 1) = isimp
        else
            pcl2simp(ipcl2, ispec, 2) = isimp
        end if

        ! 最後のsimplexを移動する
        simplices(isimp, ispec)%ipcl1 = simplices(nsimp(ispec), ispec)%ipcl1
        simplices(isimp, ispec)%ipcl2 = simplices(nsimp(ispec), ispec)%ipcl2
        simplices(isimp, ispec)%offset1 = simplices(nsimp(ispec), ispec)%offset1
        simplices(isimp, ispec)%offset2 = simplices(nsimp(ispec), ispec)%offset2
        simplices(isimp, ispec)%rate = simplices(nsimp(ispec), ispec)%rate

        ! simplex数をデクリメントする
        nsimp(ispec) = nsimp(ispec) - 1
    end subroutine

    !> @brief simplexを分割する.
    !>
    !> @param[in] isimp simplexのインデックス
    !> @param[in] ispec 粒子の種類
    subroutine simplex_split_simplex(isimp, ispec)
        integer, intent(in) :: isimp, ispec
        integer :: ipcl1, ipcl2
        real(8) :: rate, offset1, offset2
        real(8) :: px1, px2, pvx1, pvx2
        real(8) :: px_new, pvx_new
        integer :: ncycle_new
        logical :: status

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
        call particles_add_particle(ispec, px_new, pvx_new, ncycle_new, status=status)

        ! 追加できなかった場合終了
        if (.not. status) then
            return
        end if

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
        call simplex_add_simplex(ispec, ipcl1, npcl(ispec), 0.5*rate, offset1=offset1, isimp=isimp)

        ! 右側のsimplexを作成
        call simplex_add_simplex(ispec, npcl(ispec), ipcl2, 0.5*rate, offset2=offset2)
    end subroutine

    !> @brief simplexを縮約する.
    !>
    !> @param[in] isimp simplexのインデックス
    !> @param[in] ispec 粒子の種類
    subroutine simplex_contract_simplex(isimp, ispec)
        integer, intent(in) :: isimp, ispec
        integer :: ipcl1, ipcl2, isimp1, isimp2
        real(8) :: add_rate, offset1, offset2
        real(8) :: px1, px2, pvx1, pvx2
        real(8) :: px_new, pvx_new
        integer :: ncycle_new
        logical :: status

        ipcl1 = simplices(isimp, ispec)%ipcl1
        ipcl2 = simplices(isimp, ispec)%ipcl2
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
        call particles_add_particle(ispec, px_new, pvx_new, ncycle_new, status=status)

        ! 追加できなかった場合終了
        if (.not. status) then
            return
        end if

        ! 左側のsimplex番号を取得(存在しない場合-1)
        if (pcl2simp(ipcl1, ispec, 1) /= isimp) then
            isimp1 = pcl2simp(ipcl1, ispec, 1)
        else
            isimp1 = pcl2simp(ipcl1, ispec, 2)
        end if

        ! 右側のsimplex番号を取得(存在しない場合-1)
        if (pcl2simp(ipcl2, ispec, 1) /= isimp) then
            isimp2 = pcl2simp(ipcl2, ispec, 1)
        else
            isimp2 = pcl2simp(ipcl2, ispec, 2)
        end if

        ! 両側に粒子密度を分配する
        if (isimp1 /= -1 .and. isimp2 /= -1) then
            add_rate = simplices(isimp, ispec)%rate*0.5d0
        else
            add_rate = simplices(isimp, ispec)%rate
        end if

        ! 左側のsimplexに分配する
        if (isimp1 /= -1) then
            simplices(isimp1, ispec)%rate = simplices(isimp1, ispec)%rate + add_rate
        end if

        ! 右側のsimplexに分配する
        if (isimp2 /= -1) then
            simplices(isimp2, ispec)%rate = simplices(isimp2, ispec)%rate + add_rate
        end if

        ! 両側のsimplexのトレーサーを付け替える
        if (isimp1 /= -1) then
            ! 左側のsimplexのトレーサーを新しく追加したトレーサーに変更する
            if (simplices(isimp1, ispec)%ipcl1 == ipcl1) then
                simplices(isimp1, ispec)%ipcl1 = npcl(ispec)
            else
                simplices(isimp1, ispec)%ipcl2 = npcl(ispec)
            end if
        end if

        if (isimp2 /= -1) then
            ! 右側のsimplexのトレーサーを新しく追加したトレーサーに変更する
            if (simplices(isimp2, ispec)%ipcl1 == ipcl2) then
                simplices(isimp2, ispec)%ipcl1 = npcl(ispec)
            else
                simplices(isimp2, ispec)%ipcl2 = npcl(ispec)
            end if
        end if

        ! トレーサーからsimplexへのリンクを作る
        pcl2simp(npcl(ispec), ispec, 1) = isimp1
        pcl2simp(npcl(ispec), ispec, 2) = isimp2

        call simplex_delete_simplex(ispec, isimp)
        call particles_delete_particle(ispec, ipcl1)
        call particles_delete_particle(ispec, ipcl2)
    end subroutine
end module
