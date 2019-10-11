
include 'Global.f90'                                    !全局变量模块
include 'IO.f90'                                        !输入输出模块
include 'FunctionList.f90'                              !自定义函数模块


 program MF
    use Global
    use IO
    use FunctionList
    
    !=================================声明区================================
    implicit none
    integer i,j
    real(kind=8),allocatable ::IESM(:,:,:)                !局部坐标系下的单元刚度矩阵
    real(kind=8),allocatable ::ESM(:,:,:)                 !总体坐标系下的单元刚度矩阵
    real(kind=8),allocatable ::TSM(:,:)                   !总体刚度矩阵
    real(kind=8),allocatable ::PK(:,:)                    !过渡荷载向量
    real(kind=8),allocatable ::LT(:,:)                    !移植荷载向量
    real(kind=8),allocatable ::P(:)                       !荷载向量
    real(kind=8),allocatable ::U(:)                       !位移向量
    real(kind=8),allocatable :: InternalDeformationForce(:,:)!节点内力
    real(kind=8),allocatable :: EndReaction(:,:)           !支座反力
    real(kind=8),allocatable :: SumForce(:,:)              !节点总内力
    real(kind=8) :: TOL                                    !迭代精度
    
    !=================================界面区================================
    write(*,*)
    write(*,*)'***********************************************************'
    write(*,*)'*                 结构力学计算程序                        *'
    write(*,*)'*                                                         *'
    write(*,*)'*                  Power by Matte                         *'
    write(*,*)'*                     2018.1.1                            *'
    write(*,*)'*              2017港口、海岸及近海工程                   *'
    write(*,*)'*                                                         *'
    write(*,*)'***********************************************************'
    write(*,*)'*本程序可以进行以下计算：                                 *'
    write(*,*)'*1.平面桁架  2.平面刚架                                   *'
    write(*,*)'*                                                         *'
    write(*,*)'***********************************************************'
    write(*,*)'请输入您要进行的计算种类：                                 '     
    read(*,*)CalType
    write(*,*)'请输入弹性模量(KN/m^2)：                                 ' 
    !=================================输入区================================
    read(*,*)E

    do while(CalType/='1'.AND.CalType/='2')
       write(*,*)'错误的指令，请重新输入计算种类'  
       read(*,*)CalType
    end do
    !=================================加载区================================
    call InputKP()
    call InputEM()
    !================================初始化区================================
    if(CalType=='1')then
        allocate(IESM(4,4,num2))
        allocate(ESM(4,4,num2))
        allocate(TSM(2*num1,2*num1)) 
        allocate(PK(2*num1,num1))
        allocate(P(2*num1)) 
        allocate(U(2*num1)) 
    else if(CalType=='2')then
        allocate(IESM(6,6,num2))
        allocate(ESM(6,6,num2))
        allocate(TSM(3*num1,3*num1)) 
        allocate(PK(3*num1,num1))
        allocate(P(3*num1)) 
        allocate(U(3*num1))
    end if
    !=================================计算区================================  
     do i=1,num2
        IESM(:,:,i)=IESMC(i)                                      !局部坐标系下的单元刚度矩阵计算
        ESM(:,:,i)=ESMC(IESM(:,:,i),size(IESM,1),i)               !总体坐标系下的单元刚度矩阵计算
     end do
     
    TSM=TSMC(ESM,size(ESM,1))                                     !总体刚度矩阵计算
    
    do i=1,size(P)                                                !外力荷载向量初始化
        P(i)=0
    end do
    if(CalType=='2')then 
    allocate(LT(6,num2))   
    do i=1,6
        do j=1,num2
            LT(i,j)=0                                             !单元等效荷载初始化            
        end do
    end do
    
    do i=1,num2
        LT(:,i)=LoadEQ(i)                                         !局部坐标系下单元等效荷载计算
        P=P+LoadTrans(LT(:,i),i)                                  !将单元等效荷载转换至总体坐标系并升阶
    end do
    
    end if
    
    do i=1,num1
    PK(:,i)=PC(i)                                                 !节点荷载计算并升阶
    end do
    
    P=P+sum(PK,2)                                                 !组集外力荷载向量
    
    

    do i=1,num1                                                   !根据约束类型进行消行修正
      select case(cons(i))
      case(1)
          TSM=TSMM(TSM,i,CD(i),size(TSM,1))               !总体刚度矩阵消行修正（TSM总体刚度矩阵，i节点号,CD为约束方向）
          P=PM(P,i,CD(i),size(P))                         !荷载向量消行修正  （P荷载向量，i节点号,CD为约束方向）
          
      case default
      continue
     end select 
    end do
   !=================================迭代区================================   
   
    TOL=0.00001                                                                   !迭代精度
    U=Jacobi(TSM,P,TOL,size(P))                                                   !雅克比迭代法解线性方程组
    
    !===============================后处理区===============================
 
    allocate(InternalDeformationForce(size(ESM,1),num2))                          !节点变形内力
    if(CalType=='1')then
    do i=1,num2
        InternalDeformationForce(:,i)=InFor(ESM(:,:,i),U,i,size(ESM,1),size(U))   !桁架节点变形内力
    end do
    
    allocate(EndReaction((size(U)/num1)+1,num1))
    EndReaction=EndReac(InternalDeformationForce,size(InternalDeformationForce,1),size(U))!支座反力求解
    
    else if(CalType=='2')then
    allocate(SumForce(6,num2))    
     
    do i=1,num2                                                                   !刚架节点变形内力求解
      InternalDeformationForce(:,i)=ELMLoad(U,IESM(:,:,i),i)              
    end do
    
    do i=1,num2
    SumForce(:,i)=SumForceCal(InternalDeformationForce(:,i),LT(:,i))              !刚架总内力求解
    end do
        
    end if    
    
    
    
    !===============================输出区===============================
    Call Outputdis(U,size(U))                                                     !输出位移
    Call OutputIDF(InternalDeformationForce,size(InternalDeformationForce,1))     !输出节点变形内力
    if(CalType=='1')then
    Call OutputENR(EndReaction,size(EndReaction,1))                               !若为桁架，输出支座反力
    else if(CalType=='2')then
    Call OutputSUF(SumForce,size(SumForce,1))                                     !若为钢架，输出总内力
        
    end if
    
    write(*,*)'结果已写入根目录中的Result文件夹'
    write(*,*)'按任意键结束程序...'
    read(*,*)
    
    stop
    end program MF
    
