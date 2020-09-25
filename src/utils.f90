!> @brief Utilityモジュール
module utils
    implicit none

    !> 円周率
    real, parameter :: pi = 4.0*atan(1.0)

    private pi
    public utils_rand_bm, utils_rand_bm2
    public utils_pmod

contains

    !> @brief 乱数のシード値を設定する.
    !>
    !> @param[in] seed シード値
    subroutine utils_set_random_seed(seed)
        integer, intent(in), optional :: seed
        integer :: i
        integer :: seedsize
        integer, allocatable :: seeds(:)
        real(8) :: value

        call random_seed(size=seedsize)
        allocate (seeds(seedsize))
        call random_seed(get=seeds)

        if (present(seed) .and. seed /= -1) then
            seeds = seed
            call random_seed(put=seeds)

            do i = 1, seedsize
                call random_number(value)
                seeds(i) = int((2*(value - 0.5))*2147483647)
            end do
        else
            call system_clock(count=seeds(1))
        end if

        call random_seed(put=seeds)
    end subroutine

    !> @brief Box-Muller法により正規分布に従う乱数を生成する.
    !>
    !> @param[out] z 生成した乱数を代入する変数
    subroutine utils_rand_bm(z)
        real(8), intent(out) :: z
        real(8) :: x, y

        call random_number(x)
        call random_number(y)

        z = sqrt(-2*log(x))*cos(2*pi*y)
    end subroutine

    !> @brief Box-Muller法により正規分布に従う2つの独立な乱数を生成する.
    !>
    !> @param[out] z1 生成した乱数を代入する変数
    !> @param[out] z2 生成した乱数を代入する変数
    subroutine utils_rand_bm2(z1, z2)
        real(8), intent(out) :: z1, z2
        real(8) :: x, y

        call random_number(x)
        call random_number(y)

        z1 = sqrt(-2*log(x))*cos(2*pi*y)
        z2 = sqrt(-2*log(x))*sin(2*pi*y)
    end subroutine

    !> @brief 正剰余を返す.
    !>
    !> @param[in] a 割られる数
    !> @param[in] b 割る数
    !> @retval utils_pmod 正剰余
    function utils_pmod(a, b)
        real(8), intent(in) :: a
        real(8), intent(in) :: b
        real(8) :: utils_pmod
        utils_pmod = a - floor(a/b)*b
    end function

end module
