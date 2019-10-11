subroutine inputKP()
    use global
    implicit none
    character(len=20) :: filename="KP.txt"   !�ڵ���Ϣ
    integer(4) i,j
    logical alive                             !�ж��ļ��Ƿ���ڵ�Flag
    character(80):: cLine                     !�ݴ����ݵ��ַ���
    
    
    inquire(file=filename,exist=alive)       !��ѯ�ڵ��ļ��Ƿ����
    if(alive)then
        open(unit=10,file= filename,action='read')
        read(10,*)                            !���ô˷�������һ��
        
                                              !���ݶ�ȡ
        do j=1,50
        read(10,'(1Xa)',iostat=i)cline
        read(cline,*)X(j),Y(j),cons(j),CD(j),LX(j),LY(j)
        num1=j-1
        if(i/=0) exit
        end do
        close(10)
        else
        write(*,*)filename,"doesn't exist."
        stop
        end if
        
       
    return
    end subroutine