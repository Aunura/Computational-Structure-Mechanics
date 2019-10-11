
include 'Global.f90'                                    !ȫ�ֱ���ģ��
include 'IO.f90'                                        !�������ģ��
include 'FunctionList.f90'                              !�Զ��庯��ģ��


 program MF
    use Global
    use IO
    use FunctionList
    
    !=================================������================================
    implicit none
    integer i,j
    real(kind=8),allocatable ::IESM(:,:,:)                !�ֲ�����ϵ�µĵ�Ԫ�նȾ���
    real(kind=8),allocatable ::ESM(:,:,:)                 !��������ϵ�µĵ�Ԫ�նȾ���
    real(kind=8),allocatable ::TSM(:,:)                   !����նȾ���
    real(kind=8),allocatable ::PK(:,:)                    !���ɺ�������
    real(kind=8),allocatable ::LT(:,:)                    !��ֲ��������
    real(kind=8),allocatable ::P(:)                       !��������
    real(kind=8),allocatable ::U(:)                       !λ������
    real(kind=8),allocatable :: InternalDeformationForce(:,:)!�ڵ�����
    real(kind=8),allocatable :: EndReaction(:,:)           !֧������
    real(kind=8),allocatable :: SumForce(:,:)              !�ڵ�������
    real(kind=8) :: TOL                                    !��������
    
    !=================================������================================
    write(*,*)
    write(*,*)'***********************************************************'
    write(*,*)'*                 �ṹ��ѧ�������                        *'
    write(*,*)'*                                                         *'
    write(*,*)'*                  Power by Matte                         *'
    write(*,*)'*                     2018.1.1                            *'
    write(*,*)'*              2017�ۿڡ���������������                   *'
    write(*,*)'*                                                         *'
    write(*,*)'***********************************************************'
    write(*,*)'*��������Խ������¼��㣺                                 *'
    write(*,*)'*1.ƽ�����  2.ƽ��ռ�                                   *'
    write(*,*)'*                                                         *'
    write(*,*)'***********************************************************'
    write(*,*)'��������Ҫ���еļ������ࣺ                                 '     
    read(*,*)CalType
    write(*,*)'�����뵯��ģ��(KN/m^2)��                                 ' 
    !=================================������================================
    read(*,*)E

    do while(CalType/='1'.AND.CalType/='2')
       write(*,*)'�����ָ������������������'  
       read(*,*)CalType
    end do
    !=================================������================================
    call InputKP()
    call InputEM()
    !================================��ʼ����================================
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
    !=================================������================================  
     do i=1,num2
        IESM(:,:,i)=IESMC(i)                                      !�ֲ�����ϵ�µĵ�Ԫ�նȾ������
        ESM(:,:,i)=ESMC(IESM(:,:,i),size(IESM,1),i)               !��������ϵ�µĵ�Ԫ�նȾ������
     end do
     
    TSM=TSMC(ESM,size(ESM,1))                                     !����նȾ������
    
    do i=1,size(P)                                                !��������������ʼ��
        P(i)=0
    end do
    if(CalType=='2')then 
    allocate(LT(6,num2))   
    do i=1,6
        do j=1,num2
            LT(i,j)=0                                             !��Ԫ��Ч���س�ʼ��            
        end do
    end do
    
    do i=1,num2
        LT(:,i)=LoadEQ(i)                                         !�ֲ�����ϵ�µ�Ԫ��Ч���ؼ���
        P=P+LoadTrans(LT(:,i),i)                                  !����Ԫ��Ч����ת������������ϵ������
    end do
    
    end if
    
    do i=1,num1
    PK(:,i)=PC(i)                                                 !�ڵ���ؼ��㲢����
    end do
    
    P=P+sum(PK,2)                                                 !�鼯������������
    
    

    do i=1,num1                                                   !����Լ�����ͽ�����������
      select case(cons(i))
      case(1)
          TSM=TSMM(TSM,i,CD(i),size(TSM,1))               !����նȾ�������������TSM����նȾ���i�ڵ��,CDΪԼ������
          P=PM(P,i,CD(i),size(P))                         !����������������  ��P����������i�ڵ��,CDΪԼ������
          
      case default
      continue
     end select 
    end do
   !=================================������================================   
   
    TOL=0.00001                                                                   !��������
    U=Jacobi(TSM,P,TOL,size(P))                                                   !�ſ˱ȵ����������Է�����
    
    !===============================������===============================
 
    allocate(InternalDeformationForce(size(ESM,1),num2))                          !�ڵ��������
    if(CalType=='1')then
    do i=1,num2
        InternalDeformationForce(:,i)=InFor(ESM(:,:,i),U,i,size(ESM,1),size(U))   !��ܽڵ��������
    end do
    
    allocate(EndReaction((size(U)/num1)+1,num1))
    EndReaction=EndReac(InternalDeformationForce,size(InternalDeformationForce,1),size(U))!֧���������
    
    else if(CalType=='2')then
    allocate(SumForce(6,num2))    
     
    do i=1,num2                                                                   !�ռܽڵ�����������
      InternalDeformationForce(:,i)=ELMLoad(U,IESM(:,:,i),i)              
    end do
    
    do i=1,num2
    SumForce(:,i)=SumForceCal(InternalDeformationForce(:,i),LT(:,i))              !�ռ����������
    end do
        
    end if    
    
    
    
    !===============================�����===============================
    Call Outputdis(U,size(U))                                                     !���λ��
    Call OutputIDF(InternalDeformationForce,size(InternalDeformationForce,1))     !����ڵ��������
    if(CalType=='1')then
    Call OutputENR(EndReaction,size(EndReaction,1))                               !��Ϊ��ܣ����֧������
    else if(CalType=='2')then
    Call OutputSUF(SumForce,size(SumForce,1))                                     !��Ϊ�ּܣ����������
        
    end if
    
    write(*,*)'�����д���Ŀ¼�е�Result�ļ���'
    write(*,*)'���������������...'
    read(*,*)
    
    stop
    end program MF
    
