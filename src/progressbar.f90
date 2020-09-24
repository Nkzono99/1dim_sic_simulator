!> @brief プログレスバーモジュール.
!>
!> @details
!> 次のようなプログレスバーを表示する.
!> @code
!> |############               |       5500    /     10000 (10.3 sec)
!> @endcode
module progressbar
    use omp_lib
    implicit none

    !> 初めにprogressbar_showが呼ばれたときの時間
    real(8) :: start_time = -1

    private start_time
    private time

    public progressbar_show

contains

    !> プログレスバーを表示する.
    !>
    !> @param[in] current_step 現在のステップ数
    !> @param[in] max_step 最大ステップ数
    !> @param[in] n プログレスバーの長さ (default: 20)
    !> @param[in] rewrite プログレスバーを上書き出力するなら.true. (default: .true.)
    subroutine progressbar_show(current_step, max_step, n, rewrite)
        integer, intent(in) :: current_step, max_step
        integer, intent(in), optional :: n
        logical, intent(in), optional :: rewrite

        integer :: i
        integer :: n_
        logical :: rewrite_
        real(8) :: elapsed_time_sec

        if (present(rewrite)) then
            rewrite_ = rewrite
        else
            rewrite_ = .true.
        end if

        if (present(n)) then
            n_ = n
        else
            n_ = 20
        end if

        call time(elapsed_time_sec)

        if (rewrite_) then
            write (*, fmt='(a)', advance='no') char(13)
        end if

        write (*, fmt='(a)', advance='no') '|'

        do i = 1, int(current_step*n_/max_step)
            write (*, fmt='(a)', advance='no') '#'
        end do
        do i = 1, n_ - int(current_step*n_/max_step)
            write (*, fmt='(a)', advance='no') ' '
        end do

        write (*, fmt='(a)', advance='no') '|'

        write (*, fmt='(I10)', advance='no') current_step
        write (*, fmt='(a)', advance='no') ' / '
        write (*, fmt='(I10)', advance='no') max_step

        write (*, fmt='(a)', advance='no') '('
        write (*, fmt='(F0.2)', advance='no') elapsed_time_sec
        write (*, fmt='(a)', advance='no') ' sec)'
    end subroutine

    !> 現在の時間[sec]を返す.
    !>
    !> @param[out] time_sec 現在時間を格納する変数
    subroutine time(time_sec)
        real(8), intent(out) :: time_sec

        if (start_time == -1) then
            start_time = omp_get_wtime()
        end if

        time_sec = omp_get_wtime() - start_time
    end subroutine

end module
