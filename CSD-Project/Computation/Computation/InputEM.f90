subroutine inputEM()
    use global
    implicit none
    character(len=20) :: filename="EM.txt"   !�ڵ���Ϣ
    character(80):: cLine             !�ݴ����ݵ��ַ���
    integer(4) i,j
    
    logical alive                            !�ж��ļ��Ƿ���ڵ�Flag
     inquire(file=filename,exist=alive)       !��ѯ�ļ��Ƿ����
        if(alive)then
        open(unit=11,file= filename,action='read')
        read(11,*)                              !���ô˷�������һ��
        
        
        do j=1,50   
        read(11,'(1Xa)',iostat=i)cline
        read(cline,*)NA(j),NB(j),A(j)
        num2=j-1
        if(i/=0) exit
        end do
        close(11)
        else
        write(*,*)filename,"doesn't exist."
        stop
        end if 
       
    
    
    
    
    return
    end subroutine