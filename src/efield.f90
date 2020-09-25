!> @brief 電場・電位に対する処理モジュール.
module efield
    use parameters
    use commons
    use boundary
    use fft
    implicit none

    private calc_phi, calc_ex

    public efield_update

contains

    !> @brief 電場と電位を更新する.
    subroutine efield_update
        call calc_phi
        call calc_ex
    end subroutine

    !> 電位を計算する.
    subroutine calc_phi
        if (solve_method == 0) then
            call calc_phi_fft
        else if (solve_method == 1) then
            call calc_phi_iterative
        end if
    end subroutine

    !> FFTを用いてポアソン方程式を解き電位を計算する.
    subroutine calc_phi_fft
        integer :: i
        real(8) :: k

        if (.not. fft_initialized) then
            call fft_init(ngrid)
        end if

        ! rho(x) to rho(k)
        !$omp parallel do
        do i = 1, ngrid
            fft_buf(i) = -rho(i)/eps0
        end do
        !$omp end parallel do
        call fft_transform(ngrid)

        ! phi(k) = dx^2 / (2 * (cos(k*dx) - 1)) * rho(k)
        fft_buf(1) = 0
        !$omp parallel do private(k)
        do i = 2, ngrid
            k = (i - 1)*2*pi/(dx*ngrid)
            fft_buf(i) = fft_buf(i)*(dx*dx)/2/(cos(k*dx) - 1)
        end do
        !$omp end parallel do

        ! phi(k) to phi(x)
        call fft_inverse_transform(ngrid)

        !$omp parallel do
        do i = 1, ngrid
            phi(i) = real(fft_buf(i))
        end do
        !$omp end parallel do
        call boundary_correct_phi
    end subroutine

    !> 反復法(ガウスザイデル法)を用いてポアソン方程式を解き電位を計算する.
    subroutine calc_phi_iterative
        integer i, loop
        real(8) b_norm
        logical is_converged

        if (find_exact_solution) then
            b_norm = 0
            do i = 1, ngrid
                b_norm = b_norm + rho(i)*rho(i)
            end do
            b_norm = sqrt(b_norm)/eps0
        end if

        call boundary_correct_phi

        ! 反復法によりポアソン方程式を解き電位分布を求める
        is_converged = .false.
        do loop = 1, max_loop_count
            phibuf(1) = 0.5*(phi(2) + phi(0) + dx*dx/eps0*rho(1))
            do i = 2, ngrid
                phibuf(i) = 0.5*(phi(i + 1) + phibuf(i - 1) + dx*dx/eps0*rho(i))
            end do

            ! 収束判定
            if (mod(loop, convergence_judge_interval) == 0) then
                if (find_exact_solution) then
                    is_converged = exactly_converged(b_norm)
                else
                    is_converged = converged()
                end if
            end if

            ! phibufからphiにコピーする.
            do i = 1, ngrid
                phi(i) = phibuf(i)
            end do
            call boundary_correct_phi

            ! 収束していたら反復を終了
            if (is_converged) exit
        end do
    end subroutine

    !> 相対誤差により収束判定を行う.
    logical function converged()
        real(8) max_error, cur_error
        integer i

        max_error = 0
        do i = 1, ngrid
            cur_error = abs((phibuf(i) - phi(i))/phi(i))
            max_error = max(max_error, cur_error)
        end do

        converged = max_error < error_limit
    end function

    !> 真の解との残差により収束判定を行う.
    !>
    !> @param[in] b_norm ノルム
    function exactly_converged(b_norm)
        real(8), intent(in) :: b_norm
        logical :: exactly_converged
        real(8) :: val, diff
        integer :: i

        call boundary_correct_phi

        val = 0
        do i = 1, ngrid
            diff = -rho(i)/eps0 - (phibuf(i + 1) - 2*phibuf(i) + phibuf(i - 1))/(dx*dx)
            val = val + diff*diff
        end do
        val = sqrt(val)

        exactly_converged = (val/b_norm) < error_limit
    end function

    !> 電場を計算する.
    subroutine calc_ex
        integer i

        ! 電位分布から電場分布を計算する.
        do i = 0, ngrid
            ex(i) = (phi(i) - phi(i + 1))/dx
        end do

        ! self-forceを回避するため電荷密度と同じグリッドに再分配する.
        do i = ngrid, 1, -1
            ex(i) = (ex(i - 1) + ex(i))*0.5
        end do

        call boundary_correct_ex
    end subroutine

end module
