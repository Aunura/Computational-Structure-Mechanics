module Global
    implicit none
    integer,parameter :: icol=50
    
    
    integer:: num1=0               !节点数
    integer:: num2=0               !单元数
    real(kind=8):: X(icol)=0       !节点的X坐标
    real(kind=8):: Y(icol)=0       !节点的Y坐标
    real(kind=8):: LX(icol)=0      !节点上的沿X轴正向荷载大小（负值则沿X轴反向）
    real(kind=8):: LY(icol)=0      !节点上的沿Y轴正向荷载大小（负值则沿Y轴反向）
    real(kind=8):: MZ(icol)=0       !节点上的弯矩荷载大小（正值为逆时针，负值为顺时针）
    integer:: cons(icol)=0         !节点上的约束种类  0：无约束 1：刚性支座 
    integer:: CD(icol)=0           !节点上约束方向    0: 无约束 1：X方向约束 2:Y方向约束 3:XY均约束 4:弯矩约束  5:X与弯矩约束 6:Y与弯矩约束 7:刚性固定(全约束）
    integer:: NA(icol)=0           !单元A端节点编号
    integer:: NB(icol)=0           !单元B端节点编号
    integer:: EMLT(icol)=0         !单元上的荷载种类  0:无荷载 1:集中力  2:集中力偶  3：分布力
    real:: EML(icol)=0             !单元上的荷载大小（此处取局部坐标系，局部坐标系x轴正方向永远为NA指向NB，y轴正方向为x轴正方向逆时针旋转90°,弯矩以逆时针为正）
    real:: EMP(icol)=0             !单元上的荷载位置（局部坐标系下的X坐标）
    real(kind=8):: A(icol)=0       !单元横截面积
    real(kind=8):: MOI(icol)=0     !截面惯性矩
    real(kind=8):: E=1      !弹性模量
    character(len=1) :: CalType     !计算种类
    
    common  num1,num2,X,Y,LX,LY,MZ,cons,CD,NA,NB,EMLT,EML,EMP,A,MOI,E,CalType
    end module