Module IO
    use Global
    implicit none
    
    contains
    subroutine inputEM()                     !��Ԫ��Ϣ����
    use Global
    implicit none
    character(len=20) :: filename="EM.txt"   !�ڵ���Ϣ
   
    character(160):: cLine             !�ݴ����ݵ��ַ���
    integer:: i=0,j=0
    
    logical alive                            
    !�ж��ļ��Ƿ���ڵ�Flag
     inquire(file=filename,exist=alive)       
     !��ѯ�ļ��Ƿ����
        if(alive)then
        open(unit=11,file= filename,action='read')
        read(11,*)                              
        !���ô˷�������һ��
        
        if(CalType=='1')then
        do j=1,50   
        read(11,'(1Xa)',iostat=i)cline
        if(i/=0) exit
        read(cline,*)NA(j),NB(j),A(j)           
        !����ڵ��A���ڵ��B��������
        num2=j
        
        end do
        else
        do j=1,50   
        read(11,'(1Xa)',iostat=i)cline
        if(i/=0) exit
        read(cline,*)NA(j),NB(j),A(j),EMLT(j),EML(j),EMP(j),MOI(j)  
        !��ȡ�ڵ��A���ڵ��B������������Ԫ�ϵĺ������࣬��Ԫ�ϵĺ��ش�С����Ԫ�ϵĺ���λ�ã�������Ծ�
        num2=j
      
        end do    
            
        end if    
        close(11)
        else
        write(*,*)filename,"doesn't exist."
        stop
        end if 
       
    
    
    
    
    return
    end subroutine
    
    subroutine inputKP()                     !�ڵ���Ϣ����
    use Global
    implicit none
    character(len=20) :: filename="KP.txt"   !�ڵ���Ϣ
    
    integer:: i=0,j=0
    logical alive                             !�ж��ļ��Ƿ���ڵ�Flag
    character(80):: cLine                     !�ݴ����ݵ��ַ���
    
    
    inquire(file=filename,exist=alive)       !��ѯ�ڵ��ļ��Ƿ����
    if(alive)then
        open(unit=10,file= filename,action='read')
        read(10,*)                            !���ô˷�������һ��
        
        if(CalType=='1')then
        do j=1,50       !���ݶ�ȡ    
        read(10,'(1Xa)',iostat=i)cline
        if(i/=0) exit
        read(cline,*)X(j),Y(j),cons(j),CD(j),LX(j),LY(j) !��ȡ�ڵ�X���꣬Y���꣬Լ�����࣬Լ�����򣬽ڵ��ϵ���X��������ش�С���ڵ��ϵ���Y��������ش�С
        num1=j
        
        end do
       else 
         do j=1,50       !���ݶ�ȡ    
        read(10,'(1Xa)',iostat=i)cline
        if(i/=0) exit
        read(cline,*)X(j),Y(j),cons(j),CD(j),LX(j),LY(j),MZ(j)!��ȡ�ڵ�X���꣬Y���꣬Լ�����࣬Լ�����򣬽ڵ��ϵ���X��������ش�С���ڵ��ϵ���Y��������ش�С,�ڵ�����ش�С
        num1=j
        end do
        end if
        
        close(10)
    else
        write(*,*)filename,"doesn't exist."
        stop
        end if
       
       
    return
    end subroutine
    subroutine Outputdis(U,order)     
    use Global
    implicit none
    integer :: order,i
    character(len=80) :: C1="*******************************************************************************"   
    character(len=80) :: C2="*                                                                             *"   
    character(len=80) :: C3="*                                                                             *"  
    character(len=80) :: C4="*                              �ڵ�λ�Ƽ�����                               *"
    character(len=80) :: C5="*                                                                             *"
    character(len=80) :: C6="*                                                                             *"
    character(len=80) :: C7="*******************************************************************************"
    character(len=8) :: C8="�ڵ�ţ�"
    character(len=20) :: C9="�ڵ�λ��:"
    character(len=20) :: C10
    character(len=20) :: C11
    real(kind=8) :: U(order)
    call system('md Result\')         !����IVFϵͳ������ļ���
    open(12,file="Result\Displacement.txt")
    write(12,*)C1
    write(12,*)C2
    write(12,*)C3
    write(12,*)C4
    write(12,*)C5
    write(12,*)C6
    write(12,*)C7
    do i=1,order
        write(C10,'(I1)')i           !�ڵ�Ŷ����ַ���
        C11=trim(C8)//trim(C10)      !�ڵ��ƴ��
        write(12,*)C11,C9,U(i)       !���λ����Ϣ
    end do
    
    
    
    end subroutine
    
    subroutine OutputIDF(IDF,order)
    use Global
    implicit none
    integer :: order,i
    character(len=80) :: C1="*******************************************************************************"   
    character(len=80) :: C2="*                                                                             *"   
    character(len=80) :: C3="*                                                                             *"  
    character(len=80) :: C4="*                              �ڵ��������������                           *"
    character(len=80) :: C5="*                                                                             *"
    character(len=80) :: C6="*                                                                             *"
    character(len=80) :: C7="*******************************************************************************"
    character(len=20) :: C8="��Ԫ�ţ�"
    character(len=20) :: C9="U"
    character(len=20) :: C10="V"
    character(len=20) :: C11,C12,C13,C14,C15,C16,C17
    character(len=20) :: C18="M"
    character(len=20) :: C19,C20,C21,C22,C23,C24,C25
    
    real(kind=8) :: IDF(order,num2)
    
    open(12,file="Result\InternalDeformationForce.txt")
    write(12,*)C1
    write(12,*)C2
    write(12,*)C3
    write(12,*)C4
    write(12,*)C5
    write(12,*)C6
    write(12,*)C7
    
if(CalType=='1')then                   !��ܽڵ��������
    do i=1,num2
        write(C11,'(I1)')NA(i)           
        write(C12,'(I1)')NB(i)
        write(C13,'(I1)')i     
        C13=trim(C8)//trim(C13)        !�ڵ��ƴ��
        C14=trim(C9)//trim(C11)//':'   !Uƴ��
        C15=trim(C10)//trim(C11)//':'  !Vƴ��
        C16=trim(C9)//trim(C12)//':'   !Uƴ��
        C17=trim(C10)//trim(C12)//':'  !Vƴ��
        write(12,*)C13
        write(12,*)C14,IDF(1,i)        !����ڵ��������������
        write(12,*)C15,IDF(2,i)
        write(12,*)C16,IDF(3,i)
        write(12,*)C17,IDF(4,i)
    end do
    
else if(CalType=='2')then              !�ռܽڵ��������
    do i=1,num2
        write(C11,'(I1)')NA(i)
        write(C12,'(I1)')NB(i)
        write(C13,'(I1)')i
        C13=trim(C8)//trim(C13)        !�ڵ��ƴ��
        C14=trim(C9)//trim(C11)//':'   !Uƴ��
        C15=trim(C10)//trim(C11)//':'  !Vƴ��
        C16=trim(C9)//trim(C12)//':'   !Mƴ��
        C17=trim(C10)//trim(C12)//':'  !Uƴ��
        C19=trim(C18)//trim(C11)//':'  !Vƴ��
        C20=trim(C18)//trim(C12)//':'  !Mƴ��
        write(12,*)C13
        write(12,*)C14,IDF(1,i)
        write(12,*)C15,IDF(2,i)
        write(12,*)C19,IDF(3,i)
        write(12,*)C16,IDF(4,i)
        write(12,*)C17,IDF(5,i)
        write(12,*)C20,IDF(6,i)
    end do
    end if 
        
    end subroutine
    
    subroutine OutputENR(ENR,order)
    use Global
    implicit none
    integer :: order,i
    character(len=80) :: C1="*******************************************************************************"   
    character(len=80) :: C2="*                                                                             *"   
    character(len=80) :: C3="*                                                                             *"  
    character(len=80) :: C4="*                              ֧������������                               *"
    character(len=80) :: C5="*                                                                             *"
    character(len=80) :: C6="*                                                                             *"
    character(len=80) :: C7="*******************************************************************************"
    character(len=20) :: C8="�ڵ�ţ�"
    character(len=20) :: C9="PX:"
    character(len=20) :: C10="PY:"
    character(len=20) :: C12
    character(len=20) :: C13,C14,C15,C16,C17
    
    real(kind=8) :: ENR(order,num1)
    open(13,file="Result\EndReaction.txt")
    write(13,*)C1
    write(13,*)C2
    write(13,*)C3
    write(13,*)C4
    write(13,*)C5
    write(13,*)C6
    write(13,*)C7
    
    do i=1,num1
    if(ENR(3,i)==0)then
        exit
    else 
       write(C12,'(I1)')int(ENR(3,i))
       C13=trim(C8)//trim(C12)
       write(13,*)C13
       write(13,*)C9,ENR(1,i)
       write(13,*)C10,ENR(2,i)
    end if
    
    
    end do
    end subroutine
    
    subroutine OutputSUF(SUF,order)
    use Global
    implicit none
    integer :: order,i
    character(len=80) :: C1="*******************************************************************************"   
    character(len=80) :: C2="*                                                                             *"   
    character(len=80) :: C3="*                                                                             *"  
    character(len=80) :: C4="*                              �ڵ�������������                             *"
    character(len=80) :: C5="*                                                                             *"
    character(len=80) :: C6="*                                                                             *"
    character(len=80) :: C7="*******************************************************************************"
    character(len=20) :: C8="��Ԫ�ţ�"
    character(len=20) :: C9="N"
    character(len=20) :: C10="Q"
    character(len=20) :: C11,C12,C13,C14,C15,C16,C17
    character(len=20) :: C18="M"
    character(len=20) :: C19,C20,C21,C22,C23,C24,C25
    
    real(kind=8) :: SUF(order,num2)
    open(14,file="Result\SumForce.txt")
    write(14,*)C1
    write(14,*)C2
    write(14,*)C3
    write(14,*)C4
    write(14,*)C5
    write(14,*)C6
    write(14,*)C7
    
     do i=1,num2
        write(C11,'(I1)')NA(i)
        write(C12,'(I1)')NB(i)
        write(C13,'(I1)')i
        C13=trim(C8)//trim(C13)
        C14=trim(C9)//trim(C11)//':'
        C15=trim(C10)//trim(C11)//':'
        C16=trim(C9)//trim(C12)//':'
        C17=trim(C10)//trim(C12)//':'
        C19=trim(C18)//trim(C11)  
        C20=trim(C18)//trim(C12)  
        write(14,*)C13
        write(14,*)C14,SUF(1,i)
        write(14,*)C15,SUF(2,i)
        write(14,*)C19,SUF(3,i)
        write(14,*)C16,SUF(4,i)
        write(14,*)C17,SUF(5,i)
        write(14,*)C20,SUF(6,i)
    end do
    end subroutine
    end Module