!> @brief コマンドライン引数パーサーモジュール.
!>
!> @attention 不正な引数を指定した場合プログラムを終了する.
!>
!> @details
!> 使用可能な引数の種類
!!>   - 位置引数 (./program.exe [arg])
!!>   - オプション引数 (./program.exe --flag [arg])
!!>   - 省略オプション引数 (./program.exe -f [arg])
!>
!> 使用例:(program.exe inputfile.txt -o outputfile.out)
!> @code
!> character(30) :: input_filename
!> call argparse_init
!> call argparse_add('inputfilename', description='input filename')
!> call argparse_add('--output', subflagname='-o', description='output filename')
!> call argparse_add('--directory', subflagname='-d', default='data', description='data directory')
!> call argparse_parse
!>
!> print *, argparse_get('inputfilename')  ! "inputfile.txt"
!> print *, argparse_get('output')  ! "outputfile.out"
!> print *, argparse_get('directory') ! "data"
!> @endcode
module argparse
    implicit none

    !> 引数情報を格納する構造体
    type argument
        !> 名前
        character(30) :: name
        !> フラグ名(--flag)
        character(30) :: flagname
        !> フラグの省略形(-f)
        character(30) :: subflagname
        !> 説明
        character(360) :: description
        !> デフォルト値
        character(120) :: default
        !> 実際の値
        character(120) :: value
        !> デフォルト値を持つかコマンドラインにより指定されたら.true.
        logical :: configured
        !> オプション引数(--flag, -f など)の場合.true.
        logical :: is_option
    end type

    !> 値が存在しないことを示す文字列
    character(*), parameter :: none_value = "__NoneValue__"

    !> 引数の個数
    integer :: nargs
    !> プログラムの説明
    character(360) :: program_description
    !> 引数配列
    type(argument), allocatable :: arguments(:)

    private

    public argparse_init
    public argparse_add
    public argparse_parse
    public argparse_get

contains

    !> @brief パーサを初期化する.
    !>
    !> @param[in] max_nargs 設定できる引数の最大個数 (default: 10個)
    !> @param[in] description プログラムの説明 (default: none) 
    subroutine argparse_init(max_nargs, description)
        integer, intent(in), optional :: max_nargs
        character(*), intent(in), optional :: description

        if (present(max_nargs)) then
            allocate(arguments(max_nargs))
        else
            allocate(arguments(10))
        end if

        if (present(description)) then
            program_description = description
        else
            program_description = none_value
        end if
        
        nargs = 0
    end subroutine

    !> @brief 引数を追加する.
    !>
    !> @param[in] name 引数の名前 (先頭に--をつけた場合オプション引数として扱う)
    !> @param[in] subflagname オプション引数の省略形 (default: none)
    !> @param[in] description 引数の説明 (default: none)
    !> @param[in] default デフォルト値 (default: none)
    subroutine argparse_add(name, subflagname, description, default)
        character(*), intent(in) :: name
        character(*), intent(in), optional :: subflagname
        character(*), intent(in), optional :: description
        character(*), intent(in), optional :: default

        if (name(1:2) == "--") then
            arguments(nargs+1)%name = name(3:len(name))
            arguments(nargs+1)%is_option = .true.
        else
            arguments(nargs+1)%name = name
            arguments(nargs+1)%is_option = .false.
        end if

        arguments(nargs+1)%flagname = name

        if (present(subflagname)) then
            arguments(nargs+1)%subflagname = subflagname
        else
            arguments(nargs+1)%subflagname = none_value
        end if

        if (present(description)) then
            arguments(nargs+1)%description = description
        else
            arguments(nargs+1)%description = none_value
        end if
        
        if (present(default)) then
            arguments(nargs+1)%value = default
            arguments(nargs+1)%default = default
            arguments(nargs+1)%configured = .false.
        else
            arguments(nargs+1)%value = none_value
            arguments(nargs+1)%default = none_value
        end if

        nargs = nargs + 1
    end subroutine

    !> @brief コマンドライン引数をパースする.
    subroutine argparse_parse
        integer i, j, length, status
        character(:), allocatable :: arg
        character(120) :: flagname = none_value

        ! ヘルプフラグが立っているか判定
        do i = 1, command_argument_count()
            ! コマンドライン引数の文字列長を取得
            call get_command_argument(i, length=length, status=status)
            if (status /= 0) then
                print *, 'error ', status, 'on arguments', i
            end if

            ! コマンドライン引数を取得
            allocate(character(length) :: arg)
            call get_command_argument(i, arg, status=status)
            if (status /= 0) then
                print *, 'error ', status, 'on arguments', i
            end if

            if (arg == '--help' .or. arg == '-h') then
                call argparse_show_help
                call exit(1)
            end if
            deallocate(arg)
        end do

        ! コマンドライン引数の処理
        do i = 1, command_argument_count()
            ! コマンドライン引数の文字列長を取得
            call get_command_argument(i, length=length, status=status)
            if (status /= 0) then
                print *, 'error ', status, 'on arguments', i
            end if

            ! コマンドライン引数を取得
            allocate(character(length) :: arg)
            call get_command_argument(i, arg, status=status)
            if (status /= 0) then
                print *, 'error ', status, 'on arguments', i
            end if

            ! オプション引数の処理
            if(flagname /= none_value) then
                call argparse_set_option_arg(flagname, arg)
                flagname = none_value
                deallocate(arg)
                cycle
            end if

            ! --flagnameの処理
            if (length > 2) then
                if (arg(1:2) == "--") then
                    flagname = arg
                end if
            end if

            ! -subflagnameの処理
            if (flagname == none_value .and. length > 1) then
                if (arg(1:1) == '-') then
                    do j = 1, nargs
                        if (arguments(j)%subflagname == arg) then
                            flagname = arguments(j)%flagname
                            exit
                        end if
                    end do
                end if
            end if

            ! 引数がフラグなら次の引数を処理
            if (flagname /= none_value) then
                deallocate(arg)
                cycle
            end if

            ! ポジショナル引数の処理
            call argparse_set_positional_arg(arg)
            deallocate(arg)
        end do

        ! まだ指定されていない引数があれば使用方法を出力し終了する.
        do i = 1, nargs
            if(arguments(i)%value /= none_value) then
                cycle
            end if

            if(arguments(i)%default /= none_value) then
                arguments(i)%configured = .true.
                cycle
            end if

            call argparse_show_usage
            call exit(1)
        end do
    end subroutine

    !> 引数の値を取得する.
    function argparse_get(name) result(res)
        character(*), intent(in) :: name  !> 引数名
        character(:), allocatable :: res
        integer i

        do i = 1, nargs
            if (name == arguments(i)%name) then
                res = trim(arguments(i)%value)
                return
            end if
        end do

        call argparse_show_usage
        call exit(1)
    end function

    !> パースされたオプション引数に値を設定する.
    subroutine argparse_set_option_arg(flagname, value)
        character(*), intent(in) :: flagname  !> フラグ名 (例: --flag)
        character(:), allocatable, intent(in) :: value  !> 値
        integer i

        do i = 1, nargs
            if (arguments(i)%flagname == flagname) then
                arguments(i)%value = value
                arguments(i)%configured = .true.
                return
            end if
        end do

        ! もし指定されたフラグが存在しない場合エラーを出力し終了する.
        print *, 'error: unrecognized arguments:', flagname
        call exit(1)
    end subroutine

    !> パースされた位置引数に値を設定する.
    subroutine argparse_set_positional_arg(value)
        character(*), intent(in) :: value  !> 値
        integer i

        ! まだ設定されていない位置引数に設定する.
        do i = 1, nargs            
            if (arguments(i)%is_option) then
                cycle
            end if

            if(arguments(i)%configured) then
                cycle
            end if

            arguments(i)%value = value
            arguments(i)%configured = .true.
            return
        end do

        ! もしすべての位置引数が設定されていた場合エラーを出力し終了する.
        print *, 'error: unrecognized arguments:', value
        call exit(1)
    end subroutine

    !> ヘルプメッセージを出力する.
    subroutine argparse_show_help
        integer i
        type(argument) :: arg

        call argparse_show_usage

        if (program_description /= none_value) then
            print *, ''
            print *, program_description
        end if

        print *, ''

        print *, 'positional arguments:'
        do i = 1, nargs
            arg = arguments(i)
            if (.not. arg%is_option) then
                if(arg%description == none_value) then
                    print *, ' ', arg%name
                else
                    print *, ' ', arg%name, arg%description
                end if
            end if
        end do

        print *, ''

        print *, 'optional arguments:'
        do i = 1, nargs
            arg = arguments(i)
            if (arg%is_option) then
                if(arg%description == none_value) then
                    print *, ' ', arg%flagname
                else
                    print *, ' ', arg%flagname, arg%description
                end if
            end if
        end do
        print *, ' ', '--help', 'show this help message and exit'
    end subroutine

    !> プログラムの利用方法を出力する.
    subroutine argparse_show_usage
        integer i
        type(argument) :: arg

        write (*, fmt='(a)', advance='no') 'usage: ./program [--help] '
        do i = 1, nargs
            arg = arguments(i)
            if (arg%is_option) then
                if (arg%subflagname == none_value) then
                    write (*, fmt='(a)', advance='no') '[' // trim(arg%flagname) // ']'
                else
                    write (*, fmt='(a)', advance='no') '[' // trim(arg%subflagname) // ']'
                end if
            else
                write (*, fmt='(a)', advance='no') trim(arg%name)
            end if
            write (*, fmt='(a)', advance='no') ' '
        end do
        print *, ''
    end subroutine
end module
