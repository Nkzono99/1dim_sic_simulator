!> @brief 境界処理モジュール.
module boundary
    use commons
    use parameters

    implicit none

    private

    public boundary_correct_pcl, boundary_correct_rho
    public boundary_correct_phi, boundary_correct_ex

contains

    !> @brief 粒子に境界条件を適用する.
    subroutine boundary_correct_pcl
        if (boundary_type == 0) call correct_pcl_periodic
        if (boundary_type == 1) call correct_pcl_reflective
    end subroutine

    !> @brief 粒子密度に境界条件を適用する.
    subroutine boundary_correct_rho
        if (boundary_type == 0) call correct_rho_periodic
        if (boundary_type == 1) call correct_rho_reflective
    end subroutine

    !> @brief 電位に境界条件を適用する.
    subroutine boundary_correct_phi
        if (boundary_type == 0) call correct_phi_periodic
        if (boundary_type == 1) call correct_phi_reflective
    end subroutine

    !> @brief 電場に境界条件を適用する.
    subroutine boundary_correct_ex
        if (boundary_type == 0) call correct_ex_periodic
        if (boundary_type == 1) call correct_ex_reflective
    end subroutine

    !> @brief 粒子に周期境界条件を適用する.
    subroutine correct_pcl_periodic
        integer :: ispec, ipcl

        do ispec = 1, nspec
            do ipcl = 1, npcl(ispec)
                if (px(ipcl, ispec) < 0) then
                    px(ipcl, ispec) = px(ipcl, ispec) + ngrid*dx
                    ncycles(ipcl, ispec) = ncycles(ipcl, ispec) - 1
                end if
                if (px(ipcl, ispec) >= ngrid*dx) then
                    px(ipcl, ispec) = px(ipcl, ispec) - ngrid*dx
                    ncycles(ipcl, ispec) = ncycles(ipcl, ispec) + 1
                end if
            end do
        end do
    end subroutine

    !> @brief 粒子密度に周期境界条件を適用する.
    subroutine correct_rho_periodic
        integer :: ispec

        do ispec = 1, nspec
            rhospec(0, ispec) = rhospec(ngrid, ispec)
            rhospec(ngrid + 1, ispec) = rhospec(1, ispec)
        end do

        rho(0) = rho(ngrid)
        rho(ngrid + 1) = rho(1)
    end subroutine

    !> @brief 電位に周期境界条件を適用する.
    subroutine correct_phi_periodic
        phi(0) = phi(ngrid)
        phi(ngrid + 1) = phi(1)

        phibuf(0) = phibuf(ngrid)
        phibuf(ngrid + 1) = phi(1)
    end subroutine

    !> @brief 電場に周期境界条件を適用する.
    subroutine correct_ex_periodic
        ex(0) = ex(ngrid)
        ex(ngrid + 1) = ex(1)
    end subroutine

    !> @brief 粒子に完全反射境界条件を適用する.
    subroutine correct_pcl_reflective
        integer :: ispec, ipcl

        do ispec = 1, nspec
            do ipcl = 1, npcl(ispec)
                if (px(ipcl, ispec) < dx) then
                    px(ipcl, ispec) = 2*dx - px(ipcl, ispec)
                    pvx(ipcl, ispec) = -pvx(ipcl, ispec)
                end if
                if (px(ipcl, ispec) >= ngrid*dx) then
                    px(ipcl, ispec) = 2*ngrid*dx - px(ipcl, ispec)
                    pvx(ipcl, ispec) = -pvx(ipcl, ispec)
                end if
            end do
        end do
    end subroutine

    !> @brief 粒子密度に完全反射境界条件を適用する.
    subroutine correct_rho_reflective
        integer :: ispec

        do ispec = 1, nspec
            rhospec(0, ispec) = 0
            rhospec(ngrid + 1, ispec) = 0
        end do

        rho(0) = 0
        rho(ngrid + 1) = 0
    end subroutine

    !> @brief 電位に完全反射境界条件を適用する.
    subroutine correct_phi_reflective
        phi(0) = 0
        phi(ngrid + 1) = 0

        phibuf(0) = 0
        phibuf(ngrid + 1) = 0
    end subroutine

    !> @brief 電場に完全反射境界条件を適用する.
    subroutine correct_ex_reflective
        ex(0) = 0
        ex(ngrid + 1) = 0
    end subroutine
end module
