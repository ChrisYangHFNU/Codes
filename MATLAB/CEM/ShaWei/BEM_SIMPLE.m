%   BEM_SIMPLE.m
%   �������ñ߽�Ԫ������������������ڵ�λ�ֲ�
%   �����    ɳ��(Wei Sha) ���մ�ѧ(Anhui University) ws108@ahu.edu.cn

clear;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  1.��������

a=6;          %  �����γ�
N=3;          %  ÿ�ߵ���
minstep=a/N;  %  ��С��ɢ����
TOTAL=N*4;    %  ���е���
C=1/2;        %  ��������
NN=100;       %  ������ɢ���� 
V_L=300;      %  ��֪��ѹ����
xx=a/2;       %  �����ڲ�����һ��X����
yy=a/2;       %  �����ڲ�����һ��Y����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2.���궨λ
%  �Է������½�Ϊ����ԭ�㽨������ϵ
%  ƥ��������ʱ�뷽���(0,a/6)��������α��
%  ������������XΪ����,������������YΪ����

value_b=-minstep/2;   %  �²��ֵ
value_r=-minstep/2;   %  �Ҳ��ֵ
value_t=a+minstep/2;  %  �ϲ��ֵ
value_l=a+minstep/2;  %  ����ֵ

for i=1:TOTAL;
    
    if (i>0 & i<N+1)          %  �²�
        value_b=value_b+minstep;
        point(1,i)=value_b;
        point(2,i)=0;
    
    elseif (i>N & i<2*N+1)    %  �Ҳ�
        value_r=value_r+minstep;
        point(1,i)=a;
        point(2,i)=value_r;
    
    elseif (i>2*N & i<3*N+1)  %  �ϲ�
        value_t=value_t-minstep;
        point(1,i)=value_t;
        point(2,i)=a;

    else                      %  ���
        value_l=value_l-minstep;
        point(1,i)=0;
        point(2,i)=value_l;
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3.H����h_stȷ��

for s=1:TOTAL      %  ����ѭ��
    for t=1:TOTAL  %  Դ��ѭ��  
        
        if (s==t)  %  ����㴦��
            h_st(s,t)=C;
        else
            fieldpoint_x=point(1,s);   %  ����X����
            currentpoint_x=point(1,t); %  Դ��X����
            fieldpoint_y=point(2,s);   %  ����Y����
            currentpoint_y=point(2,t); %  Դ��Y���� 
            
            current_x=linspace(currentpoint_x-minstep/2,currentpoint_x+minstep/2,NN); %  X���ֱ�����ɢ
            current_y=linspace(currentpoint_y-minstep/2,currentpoint_y+minstep/2,NN); %  Y���ֱ�����ɢ
            
            
            if (t>0 & t<N+1)|(t>2*N & t<3*N+1)  %  ���²�
                
                quad=abs(fieldpoint_y-currentpoint_y)./...
                    ((fieldpoint_x-current_x).^2+(currentpoint_y-fieldpoint_y).^2);
                h_st(s,t)=-(1/(2*pi))*trapz(current_x,quad);
                
                
            else    %  ���Ҳ�
                
                quad=abs(fieldpoint_x-currentpoint_x)./...
                    ((fieldpoint_x-currentpoint_x).^2+(current_y-fieldpoint_y).^2);
                h_st(s,t)=-(1/(2*pi))*trapz(current_y,quad);
                
            end;   
        end;
    end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  4.K����k_stȷ��

for s=1:TOTAL       %  ����ѭ��
    for t=1:TOTAL   %  Դ��ѭ��  
        
         if (s==t)  %  ����㴦��
             k_st(s,t)=-(log(minstep/2)-1)*minstep/(2*pi);
             
         else
            fieldpoint_x=point(1,s);   %  ����X����
            currentpoint_x=point(1,t); %  Դ��X����
            fieldpoint_y=point(2,s);   %  ����Y����
            currentpoint_y=point(2,t); %  Դ��Y���� 
            
            current_x=linspace(currentpoint_x-minstep/2,currentpoint_x+minstep/2,NN); %  X���ֱ�����ɢ
            current_y=linspace(currentpoint_y-minstep/2,currentpoint_y+minstep/2,NN); %  Y���ֱ�����ɢ
            
            
            if ((t>0 & t<N+1)|(t>2*N & t<3*N+1))   %  ���²�
                 
                quad=log( ( (fieldpoint_x-current_x).^2 + ...
                            (currentpoint_y-fieldpoint_y).^2 ).^(1/2) );
                k_st(s,t)=-(1/(2*pi))*trapz(current_x,quad);
                
            else  %  ���Ҳ�
                
                quad=log( ( (fieldpoint_x-currentpoint_x).^2 + ...
                            (current_y-fieldpoint_y).^2 ).^(1/2) );
                k_st(s,t)=-(1/(2*pi))*trapz(current_y,quad);
                
            end;
        end;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  5.��������
%   ���²��ɷֲ���֪
%   ���Ҳ��ѹ�ֲ���֪

H_K1=[h_st(:,[1:N]),-k_st(:,[N+1:2*N]),h_st(:,[2*N+1:3*N]),-k_st(:,[3*N+1:4*N])];
H_K2=[k_st(:,[1:N]),-h_st(:,[N+1:2*N]),k_st(:,[2*N+1:3*N]),-h_st(:,[3*N+1:4*N])];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  6.��֪���ֵ�ѹ�͵�ɾ���g_uȷ��

for u=1:TOTAL;
    
    if ( (u>3*N) & (u<4*N+1) )   %  �ϲ�
        g_u(u)=V_L;    
    else
        g_u(u)=0;
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7.���ʣ�µ�ɵ�ѹ�ֲ�����ʾ
charge_voltage=(H_K1)^(-1)*H_K2*g_u.';

disp('�²��λ������:');
disp(charge_voltage(1:N));
disp('�ϲ��λ������:')
disp(charge_voltage(3*N:-1:2*N+1));
disp('�Ҳ��ɴ��ϵ���:');
disp(charge_voltage(2*N:-1:N+1));
disp('����ɴ��ϵ���:');
disp(charge_voltage(3*N+1:4*N));

voltage=[charge_voltage(1:N);g_u(N+1:2*N)';charge_voltage(2*N+1:3*N);g_u(3*N+1:4*N)'];
charge= [g_u(1:N)';charge_voltage(N+1:2*N);g_u(:,[2*N+1:3*N])';charge_voltage(3*N+1:4*N)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  8.�����ڲ����Ե�H1����h_st1ȷ��


for t=1:TOTAL  %  Դ��ѭ��  
    
    fieldpoint_x=xx;           %  ����X����
    currentpoint_x=point(1,t); %  Դ��X����
    fieldpoint_y=yy;           %  ����Y����
    currentpoint_y=point(2,t); %  Դ��Y���� 
    
    current_x=linspace(currentpoint_x-minstep/2,currentpoint_x+minstep/2,NN); %  X���ֱ�����ɢ
    current_y=linspace(currentpoint_y-minstep/2,currentpoint_y+minstep/2,NN); %  Y���ֱ�����ɢ
    
    
    if (t>0 & t<N+1)|(t>2*N & t<3*N+1) %  ���²�
        
        quad=abs(fieldpoint_y-currentpoint_y)./...
            ((fieldpoint_x-current_x).^2+(currentpoint_y-fieldpoint_y).^2);
        h_st1(t)=-(1/(2*pi))*trapz(current_x,quad);
        
    else %  ���Ҳ�
        
        quad=abs(fieldpoint_x-currentpoint_x)./...
            ((fieldpoint_x-currentpoint_x).^2+(current_y-fieldpoint_y).^2);
        h_st1(t)=-(1/(2*pi))*trapz(current_y,quad);
        
    end;   
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  9.�����ڲ����Ե�K1����k_st1ȷ��

for t=1:TOTAL  %  Դ��ѭ��  
    
    fieldpoint_x=xx;           %  ����X����
    currentpoint_x=point(1,t); %  Դ��X����
    fieldpoint_y=yy;           %  ����Y����
    currentpoint_y=point(2,t); %  Դ��Y���� 
    
    current_x=linspace(currentpoint_x-minstep/2,currentpoint_x+minstep/2,NN); %  X���ֱ�����ɢ
    current_y=linspace(currentpoint_y-minstep/2,currentpoint_y+minstep/2,NN); %  Y���ֱ�����ɢ
    
    
    if ((t>0 & t<N+1)|(t>2*N & t<3*N+1))   %  ���²�
        
        quad=log( ( (fieldpoint_x-current_x).^2 + ...
            (currentpoint_y-fieldpoint_y).^2 ).^(1/2) );
        k_st1(t)=-(1/(2*pi))*trapz(current_x,quad);
        
    else  %  ���Ҳ�
        
        quad=log( ( (fieldpoint_x-currentpoint_x).^2 + ...
            (current_y-fieldpoint_y).^2 ).^(1/2) );
        k_st1(t)=-(1/(2*pi))*trapz(current_y,quad);
        
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  10.����ڲ����Ե��λ�������
resolve=k_st1*charge-h_st1*voltage;  %  ������ɢ�����ɹ�ʽ
show=[xx,yy];

disp('�ڷ����ڲ���λֵ');
disp('    x=    y=');
disp(show);
disp('BEM����Ϊ:');
disp(resolve);

analysis=V_L*(a-xx)/a;
disp('������Ϊ:')
disp(analysis)




