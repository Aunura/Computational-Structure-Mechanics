subroutine inputKP()
    use global
    implicit none
    character(len=20) :: filename="KP.txt"   !节点信息
    integer(4) i,j
    logical alive                             !判断文件是否存在的Flag
    character(80):: cLine                     !暂存数据的字符串
    
    
    inquire(file=filename,exist=alive)       !查询节点文件是否存在
    if(alive)then
        open(unit=10,file= filename,action='read')
        read(10,*)                            !利用此法跳过第一行
        
                                              !数据读取
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