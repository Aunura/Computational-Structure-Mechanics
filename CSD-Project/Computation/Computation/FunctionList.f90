module FunctionList
    use Global
    implicit none
    
    
    contains
function TransMatrix(ELN)                                 !单元坐标转换矩阵计算
    use Global                                          
    implicit none
    integer ELN,i,j 
    real(kind=8),allocatable ::T(:,:)                     !坐标转换矩阵
    real(kind=8),allocatable ::TransMatrix(:,:)           !坐标转换矩阵
    real(kind=8) ::AL                                     !单元长度
    real(kind=8) ::st=0                                   !sin(theta)
    real(kind=8) ::ct=0                                   !cos(theta)
    
    
    
  
    if(CalType=='1')then
    allocate(T(4,4))
    allocate(TransMatrix(4,4))
    else if(CalType=='2')then
    allocate(T(6,6))
    allocate(TransMatrix(6,6))
    end if   
    
     do i=1,size(T,1)
       do j=1,size(T,1)
        T(i,j)=0
       end do
     end do
    
    AL=sqrt((X(NA(ELN))-X(NB(ELN)))**2+(Y(NA(ELN))-Y(NB(ELN)))**2)  !单元长度计算
    st=(Y(NB(ELN))-Y(NA(ELN)))/AL                                   !sin(theta)
    ct=(X(NB(ELN))-X(NA(ELN)))/AL                                   !cos(theta)
     
    if(CalType=='1')then
    T(1,1)=ct
    T(1,2)=-st
    T(2,1)=st
    T(2,2)=ct
    T(3,3)=ct
    T(3,4)=-st
    T(4,3)=st
    T(4,4)=ct
    else if(CalType=='2')then    
    T(1,1)=ct
    T(1,2)=-st
    T(2,1)=st
    T(2,2)=ct
    T(3,3)=1
    do i=1,3
        do j=1,3
        T(i+3,j+3)=T(i,j)    
        end do
    end do    
    end if
    TransMatrix=T
        
end function
    

function LoadEQ(ELN)                !局部坐标系下单元等效荷载计算
    use Global                                          
    implicit none
    integer ELN,i,j    !单元号
    real(kind=8) :: Vecc(6),Vecg(3*num1),T(6,6)
    real(kind=8) :: FA,FB,MA,MB,AL,C1,C2
    real(kind=8),allocatable :: LoadEQ(:)
    
    allocate(LoadEQ(6))
    do i=1,3*num1
        Vecg(i)=0
    end do
    do i=1,6
        Vecc(i)=0
    end do
    
    
    AL=sqrt((X(NA(ELN))-X(NB(ELN)))**2+(Y(NA(ELN))-Y(NB(ELN)))**2)  !长度计算
    C1=EMP(ELN)                                                     !荷载作用位移离A节点距离
    C2=AL-EMP(ELN)                                                  !荷载作用位移离B节点距离
    select case(EMLT(ELN))
    case(1)                                                         !集中力
        MA=-EML(ELN)*C1*(C2**2)/(AL**2)
        MB=EML(ELN)*C2*(C1**2)/(AL**2)
        FA=-EML(ELN)*(C2**2)*(1+2*C1/AL)/(AL**2)
        FB=-EML(ELN)*(C1**2)*(1+2*C2/AL)/(AL**2)
    case(2)                                                         !集中力偶
        MA=(2*C1-C2)*EML(ELN)*C2/(AL**2)
        MB=(2*C2-C1)*EML(ELN)*C1/(AL**2)
        FA=6*EML(ELN)*C1*C2/(AL**2)
        FB=-6*EML(ELN)*C1*C2/(AL**2)
    case(3)                                                         !分布力
        MA=-EML(ELN)*(AL**2)/12
        MB=EML(ELN)*(AL**2)/12
        FA=-EML(ELN)*AL/2
        FB=-EML(ELN)*AL/2
    case default
    MA=0
    MB=0
    FA=0
    FB=0
    end select
    Vecc(2)=-FA
    Vecc(3)=-MA
    Vecc(5)=-FB
    Vecc(6)=-MB
    

    
    LoadEQ=Vecc
end function
 
function LoadTrans(LoadEQ,ELN)                          !将单元等效荷载转换至总体坐标系并升阶
    use Global                                          
    implicit none
    integer ELN,i,j    !单元号
    real(kind=8) :: T(6,6)
    real(kind=8):: Vecg(3*num1)
    real(kind=8):: Vecc(3*num1)
    real(kind=8):: LoadTrans(3*num1)
    real(kind=8):: LoadEQ(3*num1)
    T=TransMatrix(ELN)
    
    do i=1,3*num1
        Vecg(i)=0
    end do
    
    Vecc=matmul(T,LoadEQ)                     !坐标转换至整体坐标系      
    
    Vecg(1+(NA(ELN)-1)*3)=Vecc(1)             !升阶
    Vecg(2+(NA(ELN)-1)*3)=Vecc(2)
    Vecg(3+(NA(ELN)-1)*3)=Vecc(3)
    Vecg(1+(NB(ELN)-1)*3)=Vecc(4)
    Vecg(2+(NB(ELN)-1)*3)=Vecc(5)
    Vecg(3+(NB(ELN)-1)*3)=Vecc(6)
    
    LoadTrans=Vecg
end function



function IESMC(ELN)                          !局部坐标系下的单元刚度矩阵计算，返回对应的单元刚度矩阵
    use Global                                          
    implicit none
    integer ELN,i,j    !单元号
    real(kind=8),allocatable ::IESM(:,:)                  !初始刚度矩阵
    real(kind=8),allocatable ::IESMC(:,:)                  !单元刚度矩阵
    real(kind=8),allocatable ::T(:,:)                     !坐标转换矩阵

    real(kind=8) ::AL                                     !单元长度
    real(kind=8) ::C(4)                                   !单元刚度矩阵基础元素
    
 
    
    
    AL=sqrt((X(NA(ELN))-X(NB(ELN)))**2+(Y(NA(ELN))-Y(NB(ELN)))**2)  !单元长度
    if(CalType=='1')then
    allocate(IESM(4,4))
    allocate(IESMC(4,4))
    allocate(T(4,4))
    
    else if(CalType=='2')then
    allocate(IESM(6,6))
    allocate(IESMC(6,6))
    allocate(T(6,6))
    end if   
        
        
    
    do i=1,size(IESM,1)
       do j=1,size(IESM,1)
        IESM(i,j)=0
       
       end do
    end do
    
   
    if(CalType=='1')then                              !桁架单元刚度矩阵计算
    C(1)=E*A(ELN)/AL
    IESM(1,1)=C(1)                          
    IESM(1,3)=-C(1)
    IESM(3,1)=-C(1)
    IESM(3,3)=C(1)
    
   
    
    
    
    else if(CalType=='2')then                   !刚架单元刚度矩阵计算
    C(1)=E*A(ELN)/AL    
    C(2)=12*E*MOI(ELN)/(AL**3)
    C(3)=6*E*MOI(ELN)/(AL**2)
    C(4)=4*E*MOI(ELN)/(AL)    
        
    IESM(1,1)=C(1)
    IESM(1,4)=-C(1)
    IESM(4,1)=-C(1)
    IESM(4,4)=C(1)
    IESM(2,2)=C(2)
    IESM(2,5)=-C(2)
    IESM(5,2)=-C(2)
    IESM(5,5)=C(2)
    IESM(2,3)=C(3)
    IESM(2,6)=C(3)
    IESM(3,2)=C(3)
    IESM(3,5)=-C(3)
    IESM(5,3)=-C(3)
    IESM(5,6)=-C(3)
    IESM(6,2)=C(3)
    IESM(6,5)=-C(3)
    IESM(3,3)=C(4)
    IESM(3,6)=0.5*C(4)
    IESM(6,3)=0.5*C(4)
    IESM(6,6)=C(4)
    end if
    IESMC=IESM
    
    
end function
function ESMC(IESM,order,ELN)              !总体坐标系下的单元刚度矩阵计算
    use Global                                          
    implicit none
    integer order,ELN   !单元号
    real(kind=8) :: IESM(order,order)
    real(kind=8) :: ESMC(order,order)
    real(kind=8) :: T(order,order)
    T=TransMatrix(ELN)
    ESMC=matmul((matmul(T,IESM)),(transpose(T)))  

end function


function TSMC(ESM,order)           !总体刚度矩阵计算
   use Global
   implicit none
   integer:: i,j,k,bu,order,COE
   real(kind=8) :: ESM(order,order,num2)
   real(kind=8),allocatable :: TSMC(:,:)
   real(kind=8),allocatable :: TM(:,:,:)
   
   
   
!========================初始化=============================!
if(CalType=='1')then   
COE=2                               !矩阵组集系数，根据运算种类决定后面局部坐标系→总体坐标系中，单元刚度矩阵如何升阶
bu=2*num1
else if(CalType=='2')then
COE=3
bu=3*num1
end if

allocate(TSMC(bu,bu))
allocate(TM(bu,bu,num2))
do i=1,bu
   do j=1,bu
    TSMC(i,j)=0
   end do 
end do
do k=1,num2
   do i=1,bu
    do j=1,bu
        TM(i,j,k)=0
    end do
end do
end do
     
!========================单元刚度矩阵升阶=============================!
do k=1,num2                                   !k为单元号
    do i=1,COE                                !i为单元刚度矩阵的行号
        do j=1,COE                            !j为单元刚度矩阵的列号
            TM(i+(NA(k)-1)*COE,j+(NA(k)-1)*COE,k)=ESM(i,j,k)                             !左上角
        end do
        do j=COE+1,2*COE
            TM(i+(NA(k)-1)*COE,j+(NA(k)-1)*COE+(NB(k)-NA(k)-1)*COE,k)=ESM(i,j,k)         !右上角
        end do
    end do
    do i=COE+1,2*COE                                                                     
       do j=1,COE                                                                        !左下角
            TM(i+(NA(k)-1)*COE+(NB(k)-NA(k)-1)*COE,j+(NA(k)-1)*COE,k)=ESM(i,j,k)
        end do  
       do j=COE+1,2*COE                                                                  !右下角
            TM(i+(NA(k)-1)*COE+(NB(k)-NA(k)-1)*COE,j+(NA(k)-1)*COE+(NB(k)-NA(k)-1)*COE,k)=ESM(i,j,k)       
        end do 
    end do
end do

TSMC=sum(TM,3)                                 !组集总体刚度矩阵
  
end function
function PC(KPN)                               !节点荷载计算并升阶
    use Global
    implicit none
    integer:: KPN,i,COE
    real(kind=8),allocatable :: PC(:)
    
    if(CalType=='1')then     !计算种类检查
    COE=2
    else if(CalType=='2')then
    COE=3
    end if
    
    allocate(PC(COE*num1))
    
    do i=1,COE*num1
        PC(i)=0
    end do
    
    if(CalType=='1')then        
    if(LX(KPN)/=0)then
    PC(KPN*2-1)=LX(KPN)
    end if
    
    if(LY(KPN)/=0)then
    PC(KPN*2)=LY(KPN)
    end if
    
    
    else if(CalType=='2')then
        
    if(LX(KPN)/=0)then
    PC(KPN*COE-2)=LX(KPN)
    end if
    
    if(LY(KPN)/=0)then
    PC(KPN*COE-1)=LY(KPN)
    end if
    
    if(MZ(KPN)/=0)then
    PC(KPN*COE)=MZ(KPN)
    end if
    end if

end function
   
   
   
function TSMM(TSM,KPN,consdir,order)           !消行修正法  TSM:总体刚度矩阵  KPN：节点号  consdir:约束方向
    use Global
    implicit none
    integer :: order
    real(kind=8) :: TSM(order,order)
    real(kind=8),allocatable :: TSMT(:,:)    
    real(kind=8),allocatable :: TSMM(:,:)  
    integer :: SwitchRes(5)
    integer :: KPN,consdir,i,COE,switch1,switch2,switch3,switch4,switch5
    allocate(TSMM(order,order))
    allocate(TSMT(order,order))
    TSMT=TSM
    COE=3
    SwitchRes=switchtest(consdir)           !控制开关
    switch1=SwitchRes(1)
    switch2=SwitchRes(2)
    switch3=SwitchRes(3)
    switch4=SwitchRes(4)
    switch5=SwitchRes(5)
   
 
 !==============================消行修正==================================   
if(switch1==1)then                        !消去除对角元上的元素FX相关（桁架）
     do i=1,order
      if(i/=KPN*2-1)then                 
          TSMT(KPN*2-1,i)=0
          TSMT(i,KPN*2-1)=0
          end if
      end do
end if
if(switch2==1)then                        !消去除对角元上的元素FX相关（桁架）
     do i=1,order
      if(i/=KPN*2)then                 
          TSMT(KPN*2,i)=0
          TSMT(i,KPN*2)=0
          end if
      end do
end if
if(switch3==1)then                             !消去除对角元上的元素FX相关（刚架）
do i=1,order  
      if(i/=KPN*COE-2)then                   
          TSMT(KPN*COE-2,i)=0
          TSMT(i,KPN*COE-2)=0
          
          end if
      end do    
end if
if(switch3==1)then
     do i=1,order 
         if(i/=KPN*COE-1)then                   !消去除对角元上的元素FY相关（刚架）
     
          TSMT(KPN*COE-1,i)=0 
          TSMT(i,KPN*COE-1)=0
          end if
          end do   
end if
if(switch5==1)then
     do i=1,order 
         if(i/=KPN*COE)then                        !消去除对角元上的元素M相关 （刚架）
          TSMT(KPN*COE,i)=0 
          TSMT(i,KPN*COE)=0
          end if
          end do   
end if    
TSMM=TSMT
end function
      
function PM(P,KPN,consdir,order)               !消行修正法  P:外力荷载向量    KPN：节点号  consdir:约束方向
    use Global
    implicit none
    integer :: order   
    real(kind=8) :: P(order)
    integer :: KPN,consdir,i,switch1,switch2,switch3,switch4,switch5
    integer :: SwitchRes(5)
    integer :: COE=3
    real(kind=8),allocatable :: PM(:)
    
    allocate(PM(order))
    PM=P
    
    SwitchRes=switchtest(consdir)
    switch1=SwitchRes(1)
    switch2=SwitchRes(2)
    switch3=SwitchRes(3)
    switch4=SwitchRes(4)
    switch5=SwitchRes(5)
    
    

if(switch1==1)then
do i=1,order
      if(i/=KPN*2-1)then                
          PM(KPN*2-1)=0                      !消除荷载向量上的对应X向荷载（桁架）
          end if
      end do

end if

if(switch2==1)then
    do i=1,order
 if(i/=KPN*2)then                
          PM(KPN*2)=0                       !消除荷载向量上的对应Y向荷载（桁架）
          end if
     end do        

end if
if(switch3==1)then
do i=1,order
      if(i/=KPN*COE-2)then                
          PM(KPN*COE-2)=0                       !消除荷载向量上的对应X向荷载（刚架）
          end if
    end do   
    end if
    
    if(switch4==1)then
    do i=1,order
      if(i/=KPN*COE-1)then                
          PM(KPN*COE-1)=0                       !消除荷载向量上的对应Y向荷载（刚架）
          end if
     end do       
    end if
    
    if(switch5==1)then
     do i=1,order
      if(i/=KPN*COE)then                
          PM(KPN*COE)=0                       !消除荷载向量上的对应弯矩（刚架）
         
          end if
      end do  
    
    end if
end function
      
function switchtest(switchecho)              !消行控制开关
integer :: switchecho,i
integer :: switchtest(5)

do i=1,5
    switchtest(i)=0
end do

   if(CalType=='1')then                       !桁架消行  
      select case(switchecho)
      case(0)  
      
      case(1)                                !消除荷载向量上的对应Y向荷载
      switchtest(1)=1 
      case(2)                                !消除荷载向量上的对应Y向荷载
      switchtest(2)=1
      case(3)                                !消除荷载向量上的对应X向/Y向荷载
      switchtest(1)=1
      switchtest(2)=1
      end select
      
else if(CalType=='2')then                          !刚架消行
     select case(switchecho)
     case(0)  
     case(1)                                  !消除荷载向量上的对应X向荷载
     switchtest(3)=1   
     case(2)                                  !消除荷载向量上的对应Y向荷载
     switchtest(4)=1      
     case(3)                                  !消除荷载向量上的对应X向及Y向荷载
     switchtest(3)=1 
     switchtest(4)=1
     case(4)                                  !消除荷载向量上的对应弯矩
     switchtest(5)=1
     case(5)                                  !消除荷载向量上的对应X向及弯矩荷载
     switchtest(3)=1 
     switchtest(5)=1
     case(6)                                  !消除荷载向量上的对应Y向及弯矩荷载
     switchtest(4)=1 
     switchtest(5)=1
     case(7)                                  !消除全部荷载
     switchtest(3)=1
     switchtest(4)=1 
     switchtest(5)=1
    
      end select
end if

end function
  
    
    
function Jacobi(M,V,TOL,Order)                  !雅克比迭代法，用于解非齐次线性方程组  M：系数矩阵    V：向量
use Global
implicit none
integer :: Order
real(kind=8) :: M(Order,Order)                     
real(kind=8) :: V(Order)
real(kind=8),allocatable :: Jacobi(:)
real(kind=8),allocatable :: B(:,:)
real(kind=8),allocatable :: F(:)
real(kind=8),allocatable :: X0(:)
real(kind=8),allocatable :: X1(:)
real(kind=8):: error=100
real(kind=8):: D,L,TOL
integer:: i,j


allocate(Jacobi(Order))
allocate(B(Order,Order))
allocate(F(Order))
allocate(X0(Order))
allocate(X1(Order))

do i=1,Order                 !迭代矩阵求解
    do j=1,Order
    if(i==j)then
    B(i,j)=0
    else
    B(i,j)=-M(i,j)/M(i,i)
    end if    
        
    end do
    F(i)=V(i)/M(i,i)          !迭代向量求解
end do

do i=1,Order
    X0(i)=1
end do




do while (error>TOL)                     !当误差大于精度时，执行此循环
    
    
    X1=matmul(B,X0)                      !X1=B*X0+F
    X1=X1+F
 
    D=0
    L=0
    do i=1,Order                         !计算误差
        D=D+((X1(i)-X0(i))**2)
        L=L+(X1(i)**2)
    end do   
    error=sqrt(D)/sqrt(L)
    
    X0=X1
end do


Jacobi=X1

end function   

function InFor(ESM,U,ELN,order1,order2)    !桁架节点变形内力计算
use Global
implicit none
integer :: i,j,order1,order2,ELN
real(kind=8) :: ESM(order1,order1)
real(kind=8) :: U(order2)
real(kind=8) :: ER(order1)
real(kind=8) :: UE(order1)
real(kind=8),allocatable :: InFor(:)

allocate(InFor(order1))

!===============================位移向量降阶===============================

 UE(1)=U(NA(ELN)*2-1)
 UE(2)=U(NA(ELN)*2)   
 UE(3)=U(NB(ELN)*2-1)
 UE(4)=U(NB(ELN)*2) 


 
 ER=matmul(ESM,UE)                          !桁架节点变形内力=总体坐标系下的单元刚度矩阵*单元对应节点位移向量
 InFor=ER
 
end function

function EndReac(Inf,order1,order2)        !支座反力求解
 use Global
 implicit none

 integer :: i,j,k,order1,order2

 real(kind=8) :: Inf(order1,num2)          !桁架节点变形内力
 real(kind=8),allocatable :: EndR(:,:)
 real(kind=8),allocatable :: Vec(:)
 real(kind=8),allocatable :: MA1(:,:)
 real(kind=8),allocatable :: EndReac(:,:)
 real(kind=8) :: ConsKP(num1)              !被约束的节点号

 allocate(EndR(order2,num2))
 allocate(Vec(order2))

 do i=1,order2
     do j=1,num2
         EndR(i,j)=0
     end do
 end do
 

 do i=1,num2                            !从节点内力中读入2个U以及两个V
     EndR(NA(i)*2-1,i)=Inf(1,i)           
     EndR(NA(i)*2,i)=Inf(2,i)
     EndR(NB(i)*2-1,i)=Inf(3,i)
     EndR(NB(i)*2,i)=Inf(4,i)
 end do

 Vec=sum(EndR,2)                        !累加成一个向量
 k=1

 do i=1,num1
     if(cons(i)==0)then                 !若该处无约束，则消去对应节点的力
     Vec(i*2-1)=0
     Vec(i*2)=0
     else
     ConsKP(k)=i                        !若该有无约束，记录约束的节点
     k=k+1
     end if
 end do

 
 k=1

 allocate(MA1(3,num1))
 allocate(EndReac(3,num1))
 do i=1,3
     do j=1,num1
     MA1(i,j)=0
     end do
 end do
 do i=1,num1
     do j=1,num1
         if(i==ConsKP(j))then       !从Vec中读取约束节点所对应的力
         MA1(1,k)=Vec(i*2-1)        !PX
         MA1(2,k)=Vec(i*2)          !PY
         MA1(3,k)=i                 !对应节点号
         k=k+1
         end if
     end do
 end do
 
 EndReac=MA1                        !将MA1赋值给支座反力向量

end function

function SumForceCal(PE,LT)                    !刚架总内力求解
use Global
implicit none
real(kind=8) ::PE(6)                           !节点变形内力
real(kind=8) ::LT(6)                           !单元等效荷载
real(kind=8) ::SumForceCal(6)
integer :: i,j

do i=1,6
   SumForceCal(i)=PE(i)+(-1)*LT(i)
end do


end function


function ELMLoad(U,IESM,ELN)        !刚架节点变形内力
use Global
implicit none
real(kind=8) ::U(3*num1)         !位移向量
real(kind=8) ::IESM(6,6)         !局部坐标系下的单元刚度矩阵
real(kind=8) ::T(6,6)
real(kind=8) ::PE(6)
integer :: ELN,i
real(kind=8) ::UE(6)
real(kind=8) ::ELMLoad(6)        !刚架节点变形内力

UE(1)=U(NA(ELN)*3-2)             !降阶
UE(2)=U(NA(ELN)*3-1)
UE(3)=U(NA(ELN)*3)
UE(4)=U(NB(ELN)*3-2)
UE(5)=U(NB(ELN)*3-1)
UE(6)=U(NB(ELN)*3)

T=TransMatrix(ELN)
UE=matmul(transpose(T),UE)      !把单元位移向量由整体坐标系转回局部坐标系

PE=matmul(IESM,UE)              !节点变形内力
ELMLoad=PE
end function

    end module