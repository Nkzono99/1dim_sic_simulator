!> @brief パラメータモジュール.
module parameters
    use omp_lib
    implicit none

    !> 真空の誘電率[F/m]
    real(8), parameter :: eps0 = 8.8541878128e-12
    !> ボルツマン定数[J/K]
    real(8), parameter :: kB = 1.380649e-23
    !> 電気素量 [C]
    real(8), parameter :: e = 1.6021766343e-19
    !> 電子質量 [kg]
    real(8), parameter :: me = 9.10938356e-31
    !> 円周率
    real(8), parameter :: pi = 3.14159265359

    !> シミュレーションできる粒子の最大種類数
    integer, parameter :: max_nspec = 5

    !> 乱数のシード値(-1の場合実行ごとに異なるシード値を与える)
    integer :: random_seed = -1

    !> シミュレーションステップ数
    integer :: nsteps

    !> 出力ステップ数
    integer :: output_steps
    !> 出力ステップ間隔
    integer :: output_skips

    !> 出力フラグ
    !> 運動エネルギーを出力するか
    logical :: output_kinetic_energy = .false.
    !> 静電エネルギーを出力するか
    logical :: output_electrostatic_energy = .false.
    !> 平均トレーサー間距離を出力するか
    logical :: output_distance_between_tracers = .false.
    !> トレーサー数を出力するか
    logical :: output_npcl = .false.
    !> Simplexの位置を出力するステップ数(0なら出力しない)
    integer :: output_simplex_steps = 0
    !> 出力するSimplex数
    integer :: output_nsimp = 100

    !> 境界条件(0: 周期境界条件, 1: 反射境界条件)
    integer :: boundary_type = 0

    ! 静電場を解くためのパラメータ
    !> 静電場の解き方(0: スペクトル法, 1: 反復法)
    integer :: solve_method = 0
    ! 反復法パラメータ
    !> 反復法の収束条件(平均相対変化率がこれ以下になれば終了)
    real(8) :: error_limit = 1.0d-6
    !> 最大反復回数
    integer :: max_loop_count = 1e5
    !> 厳密に収束判定を行うか
    logical :: find_exact_solution = .false.
    !> 何ループごとに収束判定を行うか
    integer :: convergence_judge_interval = 1

    !> 初期のトレーサー数
    integer :: npcl_init(max_nspec)
    !> 現在のトレーサー数
    integer :: npcl(max_nspec) = 0
    !> 最大トレーサー数
    integer :: max_npcl

    !> 空間グリッド数
    integer :: ngrid
    !> グリッド幅 [m]
    real(8) :: dx
    !> 時間幅 [s]
    real(8) :: dt

    !> デバイ長 [m]
    real(8) :: lambda

    !> 粒子の種類数
    integer :: nspec
    !> 各粒子の質量と電子の質量の比
    real(8) :: q_ratio(max_nspec)
    !> 各粒子の電荷と電気素量の比
    real(8) :: m_ratio(max_nspec)
    !> 各粒子の温度 [K]
    real(8) :: Ts(max_nspec)

    !> 一つのsimplexに割り当てる粒子数
    real(8) :: npcl_per_super
    !> 一つのsimplexの電荷
    real(8) :: qs(max_nspec)
    !> 一つのsimplexの質量
    real(8) :: ms(max_nspec)

    ! Refinement Parameter
    !> refinementを実行する間隔(0ならrefinementを行わない)
    integer :: refinement_interval = 0
    !> refinementを実行するしきい値 [grid]
    real(8) :: refinement_threshold = 1

    ! Simplification Parameter
    !> simplificationを実行する間隔(0ならsimplificationを行わない)
    integer :: simplification_interval = 0
    !> simplificationを実行するしきい値 [grid]
    real(8) :: simplification_threshold = 0.1

    !> スレッド数 (プログラム実行時に格納される)
    integer :: nthreads = 1

    private

    public eps0, kB, e, me, pi
    public max_nspec
    public random_seed
    public nsteps, output_steps, output_skips
    public output_kinetic_energy, output_electrostatic_energy
    public output_distance_between_tracers, output_npcl
    public output_simplex_steps, output_nsimp
    public boundary_type
    public solve_method
    public error_limit, max_loop_count
    public find_exact_solution, convergence_judge_interval
    public npcl_init, max_npcl, npcl
    public ngrid, dx, dt
    public lambda
    public nspec, q_ratio, m_ratio, Ts
    public npcl_per_super, qs, ms
    public refinement_interval, refinement_threshold
    public simplification_interval, simplification_threshold
    public nthreads

    public parameters_init

contains

    !> パラメータをファイルから読み込み初期化する.
    !>
    !> @param[in] inputfilename パラメータファイル名
    !> @param[in] output_number 出力番号(default: 11)
    subroutine parameters_init(inputfilename, output_number)
        character(*), intent(in) :: inputfilename
        integer, intent(in), optional :: output_number
        integer :: output_number_

        namelist /random/ random_seed
        namelist /simulation/ nsteps, npcl_init, max_npcl, ngrid, dx, dt, boundary_type
        namelist /solver/ solve_method, error_limit, max_loop_count, find_exact_solution, convergence_judge_interval, &
                & refinement_interval, refinement_threshold, simplification_interval, simplification_threshold
        namelist /output/ output_steps, output_kinetic_energy, output_electrostatic_energy, output_distance_between_tracers, &
                & output_npcl, output_simplex_steps, output_nsimp
        namelist /plasma/ nspec, lambda, q_ratio, m_ratio, Ts

        if (present(output_number)) then
            output_number_ = output_number
        else
            output_number_ = 11
        end if

        open (output_number_, file=inputfilename)
        read (output_number_, nml=random)
        read (output_number_, nml=simulation)
        read (output_number_, nml=solver)
        read (output_number_, nml=output)
        read (output_number_, nml=plasma)
        close (output_number_)

        output_skips = nsteps/output_steps

        if (lambda == -1.0) then
            npcl_per_super = 1
        else
            npcl_per_super = eps0*kB*Ts(1)/(lambda*lambda*e*e)*dx*ngrid/npcl_init(1)
        end if

        qs = e*q_ratio*npcl_per_super
        ms = me*m_ratio*npcl_per_super

        !$omp parallel
        !$omp single
        nthreads = omp_get_num_threads()
        !$omp end single
        !$omp end parallel
    end subroutine
end module
