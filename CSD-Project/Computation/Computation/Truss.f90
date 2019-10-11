module truss
    use global
    implicit none
    
    
    contains

real(kind=8) function ESMC(EM)
                                                   !EM为单元号
implicit none
integer i,j
real(kind=8) ::IESM(4,4)=0                         !初始刚度矩阵
real(kind=8) ::T(4,4)=0                            !坐标转换矩阵

real(kind=8) ::AL(EM) =0                          !单元长度
real(kind=8) ::C(EM)=0                            !EA/L
real(kind=8) ::st=0                                 !sin(theta)
real(kind=8) ::ct=0                                 !cos(theta)

    AL(i)=sqrt((X(NA(i))-X(NB(i)))**2+(Y(NA(i))-Y(NB(i)))**2)
    C(i)=E*A(i)/AL(i)
    IESM(1,1)=C
    IESM(1,3)=-C
    IESM(3,1)=-C
    IESM(3,3)=C
    
    st=(Y(NA(i))-Y(NB(i)))/AL(i)
    ct=(X(NA(i))-X(NB(i)))/AL(i)
    
    T(1,1)=ct
    T(1,2)=-st
    T(2,1)=st
    T(2,2)=ct
    T(3,3)=ct
    T(3,4)=-st
    T(4,3)=st
    T(4,4)=ct
    
    ESMC=matmul((matmul(T,IESM)),(transpose(T)))  !最终刚度矩阵等于坐标转换矩阵*初始刚度矩阵*坐标转换矩阵的转置
        





end function
    
    
    end module