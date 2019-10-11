module FunctionList
    use Global
    implicit none
    
    
    contains
function TransMatrix(ELN)                                 !��Ԫ����ת���������
    use Global                                          
    implicit none
    integer ELN,i,j 
    real(kind=8),allocatable ::T(:,:)                     !����ת������
    real(kind=8),allocatable ::TransMatrix(:,:)           !����ת������
    real(kind=8) ::AL                                     !��Ԫ����
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
    
    AL=sqrt((X(NA(ELN))-X(NB(ELN)))**2+(Y(NA(ELN))-Y(NB(ELN)))**2)  !��Ԫ���ȼ���
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
    

function LoadEQ(ELN)                !�ֲ�����ϵ�µ�Ԫ��Ч���ؼ���
    use Global                                          
    implicit none
    integer ELN,i,j    !��Ԫ��
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
    
    
    AL=sqrt((X(NA(ELN))-X(NB(ELN)))**2+(Y(NA(ELN))-Y(NB(ELN)))**2)  !���ȼ���
    C1=EMP(ELN)                                                     !��������λ����A�ڵ����
    C2=AL-EMP(ELN)                                                  !��������λ����B�ڵ����
    select case(EMLT(ELN))
    case(1)                                                         !������
        MA=-EML(ELN)*C1*(C2**2)/(AL**2)
        MB=EML(ELN)*C2*(C1**2)/(AL**2)
        FA=-EML(ELN)*(C2**2)*(1+2*C1/AL)/(AL**2)
        FB=-EML(ELN)*(C1**2)*(1+2*C2/AL)/(AL**2)
    case(2)                                                         !������ż
        MA=(2*C1-C2)*EML(ELN)*C2/(AL**2)
        MB=(2*C2-C1)*EML(ELN)*C1/(AL**2)
        FA=6*EML(ELN)*C1*C2/(AL**2)
        FB=-6*EML(ELN)*C1*C2/(AL**2)
    case(3)                                                         !�ֲ���
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
 
function LoadTrans(LoadEQ,ELN)                          !����Ԫ��Ч����ת������������ϵ������
    use Global                                          
    implicit none
    integer ELN,i,j    !��Ԫ��
    real(kind=8) :: T(6,6)
    real(kind=8):: Vecg(3*num1)
    real(kind=8):: Vecc(3*num1)
    real(kind=8):: LoadTrans(3*num1)
    real(kind=8):: LoadEQ(3*num1)
    T=TransMatrix(ELN)
    
    do i=1,3*num1
        Vecg(i)=0
    end do
    
    Vecc=matmul(T,LoadEQ)                     !����ת������������ϵ      
    
    Vecg(1+(NA(ELN)-1)*3)=Vecc(1)             !����
    Vecg(2+(NA(ELN)-1)*3)=Vecc(2)
    Vecg(3+(NA(ELN)-1)*3)=Vecc(3)
    Vecg(1+(NB(ELN)-1)*3)=Vecc(4)
    Vecg(2+(NB(ELN)-1)*3)=Vecc(5)
    Vecg(3+(NB(ELN)-1)*3)=Vecc(6)
    
    LoadTrans=Vecg
end function



function IESMC(ELN)                          !�ֲ�����ϵ�µĵ�Ԫ�նȾ�����㣬���ض�Ӧ�ĵ�Ԫ�նȾ���
    use Global                                          
    implicit none
    integer ELN,i,j    !��Ԫ��
    real(kind=8),allocatable ::IESM(:,:)                  !��ʼ�նȾ���
    real(kind=8),allocatable ::IESMC(:,:)                  !��Ԫ�նȾ���
    real(kind=8),allocatable ::T(:,:)                     !����ת������

    real(kind=8) ::AL                                     !��Ԫ����
    real(kind=8) ::C(4)                                   !��Ԫ�նȾ������Ԫ��
    
 
    
    
    AL=sqrt((X(NA(ELN))-X(NB(ELN)))**2+(Y(NA(ELN))-Y(NB(ELN)))**2)  !��Ԫ����
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
    
   
    if(CalType=='1')then                              !��ܵ�Ԫ�նȾ������
    C(1)=E*A(ELN)/AL
    IESM(1,1)=C(1)                          
    IESM(1,3)=-C(1)
    IESM(3,1)=-C(1)
    IESM(3,3)=C(1)
    
   
    
    
    
    else if(CalType=='2')then                   !�ռܵ�Ԫ�նȾ������
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
function ESMC(IESM,order,ELN)              !��������ϵ�µĵ�Ԫ�նȾ������
    use Global                                          
    implicit none
    integer order,ELN   !��Ԫ��
    real(kind=8) :: IESM(order,order)
    real(kind=8) :: ESMC(order,order)
    real(kind=8) :: T(order,order)
    T=TransMatrix(ELN)
    ESMC=matmul((matmul(T,IESM)),(transpose(T)))  

end function


function TSMC(ESM,order)           !����նȾ������
   use Global
   implicit none
   integer:: i,j,k,bu,order,COE
   real(kind=8) :: ESM(order,order,num2)
   real(kind=8),allocatable :: TSMC(:,:)
   real(kind=8),allocatable :: TM(:,:,:)
   
   
   
!========================��ʼ��=============================!
if(CalType=='1')then   
COE=2                               !�����鼯ϵ�����������������������ֲ�����ϵ����������ϵ�У���Ԫ�նȾ����������
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
     
!========================��Ԫ�նȾ�������=============================!
do k=1,num2                                   !kΪ��Ԫ��
    do i=1,COE                                !iΪ��Ԫ�նȾ�����к�
        do j=1,COE                            !jΪ��Ԫ�նȾ�����к�
            TM(i+(NA(k)-1)*COE,j+(NA(k)-1)*COE,k)=ESM(i,j,k)                             !���Ͻ�
        end do
        do j=COE+1,2*COE
            TM(i+(NA(k)-1)*COE,j+(NA(k)-1)*COE+(NB(k)-NA(k)-1)*COE,k)=ESM(i,j,k)         !���Ͻ�
        end do
    end do
    do i=COE+1,2*COE                                                                     
       do j=1,COE                                                                        !���½�
            TM(i+(NA(k)-1)*COE+(NB(k)-NA(k)-1)*COE,j+(NA(k)-1)*COE,k)=ESM(i,j,k)
        end do  
       do j=COE+1,2*COE                                                                  !���½�
            TM(i+(NA(k)-1)*COE+(NB(k)-NA(k)-1)*COE,j+(NA(k)-1)*COE+(NB(k)-NA(k)-1)*COE,k)=ESM(i,j,k)       
        end do 
    end do
end do

TSMC=sum(TM,3)                                 !�鼯����նȾ���
  
end function
function PC(KPN)                               !�ڵ���ؼ��㲢����
    use Global
    implicit none
    integer:: KPN,i,COE
    real(kind=8),allocatable :: PC(:)
    
    if(CalType=='1')then     !����������
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
   
   
   
function TSMM(TSM,KPN,consdir,order)           !����������  TSM:����նȾ���  KPN���ڵ��  consdir:Լ������
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
    SwitchRes=switchtest(consdir)           !���ƿ���
    switch1=SwitchRes(1)
    switch2=SwitchRes(2)
    switch3=SwitchRes(3)
    switch4=SwitchRes(4)
    switch5=SwitchRes(5)
   
 
 !==============================��������==================================   
if(switch1==1)then                        !��ȥ���Խ�Ԫ�ϵ�Ԫ��FX��أ���ܣ�
     do i=1,order
      if(i/=KPN*2-1)then                 
          TSMT(KPN*2-1,i)=0
          TSMT(i,KPN*2-1)=0
          end if
      end do
end if
if(switch2==1)then                        !��ȥ���Խ�Ԫ�ϵ�Ԫ��FX��أ���ܣ�
     do i=1,order
      if(i/=KPN*2)then                 
          TSMT(KPN*2,i)=0
          TSMT(i,KPN*2)=0
          end if
      end do
end if
if(switch3==1)then                             !��ȥ���Խ�Ԫ�ϵ�Ԫ��FX��أ��ռܣ�
do i=1,order  
      if(i/=KPN*COE-2)then                   
          TSMT(KPN*COE-2,i)=0
          TSMT(i,KPN*COE-2)=0
          
          end if
      end do    
end if
if(switch3==1)then
     do i=1,order 
         if(i/=KPN*COE-1)then                   !��ȥ���Խ�Ԫ�ϵ�Ԫ��FY��أ��ռܣ�
     
          TSMT(KPN*COE-1,i)=0 
          TSMT(i,KPN*COE-1)=0
          end if
          end do   
end if
if(switch5==1)then
     do i=1,order 
         if(i/=KPN*COE)then                        !��ȥ���Խ�Ԫ�ϵ�Ԫ��M��� ���ռܣ�
          TSMT(KPN*COE,i)=0 
          TSMT(i,KPN*COE)=0
          end if
          end do   
end if    
TSMM=TSMT
end function
      
function PM(P,KPN,consdir,order)               !����������  P:������������    KPN���ڵ��  consdir:Լ������
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
          PM(KPN*2-1)=0                      !�������������ϵĶ�ӦX����أ���ܣ�
          end if
      end do

end if

if(switch2==1)then
    do i=1,order
 if(i/=KPN*2)then                
          PM(KPN*2)=0                       !�������������ϵĶ�ӦY����أ���ܣ�
          end if
     end do        

end if
if(switch3==1)then
do i=1,order
      if(i/=KPN*COE-2)then                
          PM(KPN*COE-2)=0                       !�������������ϵĶ�ӦX����أ��ռܣ�
          end if
    end do   
    end if
    
    if(switch4==1)then
    do i=1,order
      if(i/=KPN*COE-1)then                
          PM(KPN*COE-1)=0                       !�������������ϵĶ�ӦY����أ��ռܣ�
          end if
     end do       
    end if
    
    if(switch5==1)then
     do i=1,order
      if(i/=KPN*COE)then                
          PM(KPN*COE)=0                       !�������������ϵĶ�Ӧ��أ��ռܣ�
         
          end if
      end do  
    
    end if
end function
      
function switchtest(switchecho)              !���п��ƿ���
integer :: switchecho,i
integer :: switchtest(5)

do i=1,5
    switchtest(i)=0
end do

   if(CalType=='1')then                       !�������  
      select case(switchecho)
      case(0)  
      
      case(1)                                !�������������ϵĶ�ӦY�����
      switchtest(1)=1 
      case(2)                                !�������������ϵĶ�ӦY�����
      switchtest(2)=1
      case(3)                                !�������������ϵĶ�ӦX��/Y�����
      switchtest(1)=1
      switchtest(2)=1
      end select
      
else if(CalType=='2')then                          !�ռ�����
     select case(switchecho)
     case(0)  
     case(1)                                  !�������������ϵĶ�ӦX�����
     switchtest(3)=1   
     case(2)                                  !�������������ϵĶ�ӦY�����
     switchtest(4)=1      
     case(3)                                  !�������������ϵĶ�ӦX��Y�����
     switchtest(3)=1 
     switchtest(4)=1
     case(4)                                  !�������������ϵĶ�Ӧ���
     switchtest(5)=1
     case(5)                                  !�������������ϵĶ�ӦX����غ���
     switchtest(3)=1 
     switchtest(5)=1
     case(6)                                  !�������������ϵĶ�ӦY����غ���
     switchtest(4)=1 
     switchtest(5)=1
     case(7)                                  !����ȫ������
     switchtest(3)=1
     switchtest(4)=1 
     switchtest(5)=1
    
      end select
end if

end function
  
    
    
function Jacobi(M,V,TOL,Order)                  !�ſ˱ȵ����������ڽ��������Է�����  M��ϵ������    V������
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

do i=1,Order                 !�����������
    do j=1,Order
    if(i==j)then
    B(i,j)=0
    else
    B(i,j)=-M(i,j)/M(i,i)
    end if    
        
    end do
    F(i)=V(i)/M(i,i)          !�����������
end do

do i=1,Order
    X0(i)=1
end do




do while (error>TOL)                     !�������ھ���ʱ��ִ�д�ѭ��
    
    
    X1=matmul(B,X0)                      !X1=B*X0+F
    X1=X1+F
 
    D=0
    L=0
    do i=1,Order                         !�������
        D=D+((X1(i)-X0(i))**2)
        L=L+(X1(i)**2)
    end do   
    error=sqrt(D)/sqrt(L)
    
    X0=X1
end do


Jacobi=X1

end function   

function InFor(ESM,U,ELN,order1,order2)    !��ܽڵ������������
use Global
implicit none
integer :: i,j,order1,order2,ELN
real(kind=8) :: ESM(order1,order1)
real(kind=8) :: U(order2)
real(kind=8) :: ER(order1)
real(kind=8) :: UE(order1)
real(kind=8),allocatable :: InFor(:)

allocate(InFor(order1))

!===============================λ����������===============================

 UE(1)=U(NA(ELN)*2-1)
 UE(2)=U(NA(ELN)*2)   
 UE(3)=U(NB(ELN)*2-1)
 UE(4)=U(NB(ELN)*2) 


 
 ER=matmul(ESM,UE)                          !��ܽڵ��������=��������ϵ�µĵ�Ԫ�նȾ���*��Ԫ��Ӧ�ڵ�λ������
 InFor=ER
 
end function

function EndReac(Inf,order1,order2)        !֧���������
 use Global
 implicit none

 integer :: i,j,k,order1,order2

 real(kind=8) :: Inf(order1,num2)          !��ܽڵ��������
 real(kind=8),allocatable :: EndR(:,:)
 real(kind=8),allocatable :: Vec(:)
 real(kind=8),allocatable :: MA1(:,:)
 real(kind=8),allocatable :: EndReac(:,:)
 real(kind=8) :: ConsKP(num1)              !��Լ���Ľڵ��

 allocate(EndR(order2,num2))
 allocate(Vec(order2))

 do i=1,order2
     do j=1,num2
         EndR(i,j)=0
     end do
 end do
 

 do i=1,num2                            !�ӽڵ������ж���2��U�Լ�����V
     EndR(NA(i)*2-1,i)=Inf(1,i)           
     EndR(NA(i)*2,i)=Inf(2,i)
     EndR(NB(i)*2-1,i)=Inf(3,i)
     EndR(NB(i)*2,i)=Inf(4,i)
 end do

 Vec=sum(EndR,2)                        !�ۼӳ�һ������
 k=1

 do i=1,num1
     if(cons(i)==0)then                 !���ô���Լ��������ȥ��Ӧ�ڵ����
     Vec(i*2-1)=0
     Vec(i*2)=0
     else
     ConsKP(k)=i                        !��������Լ������¼Լ���Ľڵ�
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
         if(i==ConsKP(j))then       !��Vec�ж�ȡԼ���ڵ�����Ӧ����
         MA1(1,k)=Vec(i*2-1)        !PX
         MA1(2,k)=Vec(i*2)          !PY
         MA1(3,k)=i                 !��Ӧ�ڵ��
         k=k+1
         end if
     end do
 end do
 
 EndReac=MA1                        !��MA1��ֵ��֧����������

end function

function SumForceCal(PE,LT)                    !�ռ����������
use Global
implicit none
real(kind=8) ::PE(6)                           !�ڵ��������
real(kind=8) ::LT(6)                           !��Ԫ��Ч����
real(kind=8) ::SumForceCal(6)
integer :: i,j

do i=1,6
   SumForceCal(i)=PE(i)+(-1)*LT(i)
end do


end function


function ELMLoad(U,IESM,ELN)        !�ռܽڵ��������
use Global
implicit none
real(kind=8) ::U(3*num1)         !λ������
real(kind=8) ::IESM(6,6)         !�ֲ�����ϵ�µĵ�Ԫ�նȾ���
real(kind=8) ::T(6,6)
real(kind=8) ::PE(6)
integer :: ELN,i
real(kind=8) ::UE(6)
real(kind=8) ::ELMLoad(6)        !�ռܽڵ��������

UE(1)=U(NA(ELN)*3-2)             !����
UE(2)=U(NA(ELN)*3-1)
UE(3)=U(NA(ELN)*3)
UE(4)=U(NB(ELN)*3-2)
UE(5)=U(NB(ELN)*3-1)
UE(6)=U(NB(ELN)*3)

T=TransMatrix(ELN)
UE=matmul(transpose(T),UE)      !�ѵ�Ԫλ����������������ϵת�ؾֲ�����ϵ

PE=matmul(IESM,UE)              !�ڵ��������
ELMLoad=PE
end function

    end module