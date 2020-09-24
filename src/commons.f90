!> @brief 共有変数定義モジュール.
!>
!> @attention 初めにcommon_init関数を呼び出し、メモリを確保すること.
module commons
    use parameters
    implicit none

    !> 電荷分布 [C]
    real(8), allocatable :: rho(:)
    !> 粒子分布 [super particles]
    real(8), allocatable :: rhospec(:, :)
    !> 電位分布 [V]
    real(8), allocatable :: phi(:)
    !> 電位分布バッファ [V]
    real(8), allocatable :: phibuf(:)
    !> 電場分布 [N/C]
    real(8), allocatable :: ex(:)
    !> 粒子位置 [m]
    real(8), allocatable :: px(:, :)
    !> 粒子速度 [m/s]
    real(8), allocatable :: pvx(:, :)
    !> 粒子が周期境界をまたいだ数
    integer, allocatable :: ncycles(:, :)

    !> シンプレックス構造体
    type :: simplex
        !> 頂点となる1つ目の粒子のインデックス
        integer :: ipcl1
        !> 頂点となる2つ目の粒子のインデックス
        integer :: ipcl2
        !> 粒子密度(初期が1.0, シンプレックス分割に応じて減少)
        real(8) :: rate
        !> 1つ目の頂点粒子のオフセット
        real(8) :: offset1
        !> 2つ目の頂点粒子のオフセット
        real(8) :: offset2
    end type

    !> シンプレックス配列
    type(simplex), allocatable :: simplices(:, :)
    !> シンプレックスの個数
    integer :: nsimp(max_nspec)
    !> particleからsimplexへのリンク配列
    integer, allocatable :: pcl2simp(:, :, :)

contains

    !> @brief 共有変数を初期化する.
    subroutine commons_init
        allocate(rho(0:ngrid+1))
        allocate(rhospec(0:ngrid+1, nspec))
        allocate(phi(0:ngrid+1))
        allocate(phibuf(0:ngrid+1))
        allocate(ex(0:ngrid+1))

        allocate(px(max_npcl, nspec))
        allocate(pvx(max_npcl, nspec))

        allocate(ncycles(max_npcl, nspec))
        allocate(simplices(max_npcl, nspec))
        allocate(pcl2simp(max_npcl, nspec, 2))

        rho(0:ngrid+1) = 0
        rhospec(0:ngrid+1, 1:nspec) = 0
        phi(0:ngrid+1) = 0
        phibuf(0:ngrid+1) = 0
        ex(0:ngrid+1) = 0
        px(1:max_npcl, 1:nspec) = 0
        pvx(1:max_npcl, 1:nspec) = 0
        ncycles(1:max_npcl, 1:nspec) = 0
    
        nsimp(1:max_nspec) = 0
        pcl2simp(1:max_npcl, 1:nspec, 1:2) = -1
    end subroutine

end module