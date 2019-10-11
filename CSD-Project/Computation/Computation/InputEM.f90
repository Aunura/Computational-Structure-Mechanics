subroutine inputEM()
    use global
    implicit none
    character(len=20) :: filename="EM.txt"   !节点信息
    character(80):: cLine             !暂存数据的字符串
    integer(4) i,j
    
    logical alive                            !判断文件是否存在的Flag
     inquire(file=filename,exist=alive)       !查询文件是否存在
        if(alive)then
        open(unit=11,file= filename,action='read')
        read(11,*)                              !利用此法跳过第一行
        
        
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