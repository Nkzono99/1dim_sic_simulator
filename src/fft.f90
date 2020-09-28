!> @brief 2の乗数長の配列に対する1次元FFT(Fast Fourier Transform)モジュール. 
!> @attention fftを行う前にfft_initを呼び出す必要あり.
!>
!> @details
!> 使用例:
!> @code
!> integer, parameter :: n = 32
!> real(8) :: xs(n)  ! x(t)
!>
!> call random_number(xs(1:n))
!>
!> call fft_init(n)
!> fft_buf(1:n) = xs(1:n)
!>
!> ! x(t) to x(k)
!> call fft_transform
!> xs(1:n) = fft_buf(1:n)
!>
!> ! x(k) to x(t)
!> call fft_inverse_transform
!> xs(1:n) = fft_buf(1:n)
!> @endcode
module fft
    implicit none
    !> 虚数単位
    complex(8), parameter :: j = (0.0, 1.0)
    !> 円周率
    real(8), parameter :: pi = 4.0*atan(1.0d0)

    !> fftモジュールが初期化されているか
    logical :: fft_initialized = .false.
    !> fftを行う配列
    complex(8), allocatable :: fft_buf(:)
    !> fftを行う際の一時的なバッファ
    complex(8), allocatable :: fft_buf2(:)

    private

    public fft_buf
    public fft_initialized
    public fft_init
    public fft_transform, fft_inverse_transform

contains

    !> モジュールを初期化する.
    !> @param[in] n 配列長
    subroutine fft_init(n)
        integer, intent(in) :: n
        integer :: bufsize

        if (allocated(fft_buf)) then
            deallocate (fft_buf)
            deallocate (fft_buf2)
        end if

        bufsize = find_near_pow2(n)

        allocate (fft_buf(bufsize))
        allocate (fft_buf2(bufsize))

        fft_initialized = .true.
    end subroutine

    !> @brief FFT(Fast Fourier Transform)を行う.
    !> @param[in] n 配列長
    subroutine fft_transform(n)
        integer, intent(in) :: n
        integer :: i, istart, num, skip

        do i = n + 1, size(fft_buf)
            fft_buf(i) = 0
        end do

        !$omp parallel private(i, num, skip)
        do i = 1, ceiling(log(real(n))/log(2.0)) + 1
            num = 2**(i - 1)
            skip = int(n/num)
            !$omp do
            do istart = 1, skip
                call sub_transform(istart, num, skip)
            end do
            !$omp end do
        end do
        !$omp end parallel
    end subroutine

    !> @brief IFFT(Inverse Fast Fourier Transform)を行う.
    !> @param[in] n 配列長
    subroutine fft_inverse_transform(n)
        integer, intent(in) :: n
        integer :: i

        do i = 1, n
            fft_buf(i) = conjg(fft_buf(i))
        end do

        call fft_transform(n)

        do i = 1, n
            fft_buf(i) = conjg(fft_buf(i))/n
        end do
    end subroutine

    !> @brief指定された範囲でバタフライ演算を行う.
    !> @param[in] istert 範囲の初め
    !> @param[in] n 個数
    !> @param[in] skip スキップ数
    subroutine sub_transform(istart, n, skip)
        integer, intent(in) :: istart, n, skip
        integer :: n_half, i, k, ieven, iodd
        complex(8) :: wn

        if (n == 1) then
            return
        end if

        n_half = int(n/2)

        do i = 1, n_half
            k = istart + (i - 1)*skip
            ieven = istart + (i - 1)*skip*2
            iodd = istart + (i - 1)*skip*2 + skip

            wn = exp(-j*(2*pi*(i - 1))/n)

            fft_buf2(k) = fft_buf(ieven) + wn*fft_buf(iodd)
            fft_buf2(k + n_half*skip) = fft_buf(ieven) - wn*fft_buf(iodd)
        end do

        do i = 1, n
            k = istart + (i - 1)*skip
            fft_buf(k) = fft_buf2(k)
        end do
    end subroutine

    !> @brief n以上の最小の2の乗数を返す(nが0以下の場合は0).
    !>
    !> 例:
    !!>   2 => 2
    !!>   9 => 16
    !!>  -5 => 0
    function find_near_pow2(n)
        integer, intent(in) :: n
        integer :: find_near_pow2

        if (n <= 0) then
            find_near_pow2 = 0
            return
        end if

        find_near_pow2 = 2**ceiling(log(real(n))/log(2.0))
        return
    end function
end module
