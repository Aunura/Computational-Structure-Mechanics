module Global
    implicit none
    integer,parameter :: icol=50
    
    
    integer:: num1=0               !�ڵ���
    integer:: num2=0               !��Ԫ��
    real(kind=8):: X(icol)=0       !�ڵ��X����
    real(kind=8):: Y(icol)=0       !�ڵ��Y����
    real(kind=8):: LX(icol)=0      !�ڵ��ϵ���X��������ش�С����ֵ����X�ᷴ��
    real(kind=8):: LY(icol)=0      !�ڵ��ϵ���Y��������ش�С����ֵ����Y�ᷴ��
    real(kind=8):: MZ(icol)=0       !�ڵ��ϵ���غ��ش�С����ֵΪ��ʱ�룬��ֵΪ˳ʱ�룩
    integer:: cons(icol)=0         !�ڵ��ϵ�Լ������  0����Լ�� 1������֧�� 
    integer:: CD(icol)=0           !�ڵ���Լ������    0: ��Լ�� 1��X����Լ�� 2:Y����Լ�� 3:XY��Լ�� 4:���Լ��  5:X�����Լ�� 6:Y�����Լ�� 7:���Թ̶�(ȫԼ����
    integer:: NA(icol)=0           !��ԪA�˽ڵ���
    integer:: NB(icol)=0           !��ԪB�˽ڵ���
    integer:: EMLT(icol)=0         !��Ԫ�ϵĺ�������  0:�޺��� 1:������  2:������ż  3���ֲ���
    real:: EML(icol)=0             !��Ԫ�ϵĺ��ش�С���˴�ȡ�ֲ�����ϵ���ֲ�����ϵx����������ԶΪNAָ��NB��y��������Ϊx����������ʱ����ת90��,�������ʱ��Ϊ����
    real:: EMP(icol)=0             !��Ԫ�ϵĺ���λ�ã��ֲ�����ϵ�µ�X���꣩
    real(kind=8):: A(icol)=0       !��Ԫ������
    real(kind=8):: MOI(icol)=0     !������Ծ�
    real(kind=8):: E=1      !����ģ��
    character(len=1) :: CalType     !��������
    
    common  num1,num2,X,Y,LX,LY,MZ,cons,CD,NA,NB,EMLT,EML,EMP,A,MOI,E,CalType
    end module