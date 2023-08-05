clc;clear

%  �ڵ��Ԫ����Ϣ
node=load('NODE_rectangular.txt');
element=load('ELEMENT_rectangular.txt');

%  cm-->m 
node(:,2)=node(:,2)*100;
node(:,3)=node(:,3)*100;
node(:,4)=node(:,4)*100;

%  ͳ�ƽڵ������������������߽������θ���
N_node=length(node);
N_tetra=sum(element(:,2)==1);
N_tri=sum(element(:,2)==101);

%  �������Ӧ�������
Index_edge=zeros(N_tetra,6);
edge1=[]; % ��߶�������1
edge2=[]; % ��߶�������2

%  ��¼�����Ϣ
for m=1:N_tetra  %  ������
    %  ÿ��Ԫ��Ӧ6�����
    e_edge1=[element(m,3),element(m,3),element(m,3),element(m,4),element(m,4),element(m,5)];
    e_edge2=[element(m,4),element(m,5),element(m,6),element(m,5),element(m,6),element(m,6)];
    %  6����ѭ��
    for n=1:6
        flag=0; %  edge1��edge2�Ҳ���
        for p=1:length(edge1)
            ref_sum=edge1(p)+edge2(p);
            ref_sub=abs(edge1(p)-edge2(p));
            curr_sum=e_edge1(n)+e_edge2(n);
            curr_sub=abs(e_edge1(n)-e_edge2(n));
            if (ref_sum==curr_sum & ref_sub==curr_sub)  %  edge���ҵ���
                flag=1;
                break;
            end
        end
        if (flag==0) %  edge1��edge2�Ҳ���(������߿�)
            edge1=[edge1,e_edge1(n)];
            edge2=[edge2,e_edge2(n)];
        end
    end
end

clear e_edge1
clear e_edge2

%  ȥ���߽����
for m=N_tetra+1:N_tetra+N_tri  %  �߽�������
    %  ÿ��Ԫ��Ӧ3�����
    e_edge1=[element(m,3),element(m,3),element(m,4)];
    e_edge2=[element(m,4),element(m,5),element(m,5)];
    %  3����ѭ��
    for n=1:3
        flag=0; %  edge1��edge2�Ҳ���
        for p=1:length(edge1)
            ref_sum=edge1(p)+edge2(p);
            ref_sub=abs(edge1(p)-edge2(p));
            curr_sum=e_edge1(n)+e_edge2(n);
            curr_sub=abs(e_edge1(n)-e_edge2(n));
            if (ref_sum==curr_sum & ref_sub==curr_sub)  %  edge���ҵ���
                flag=1;
                index_e=p; %  edge�е�λ��
                break;
            end
        end
        if (flag==1) %  edge1��edge2�ҵ�(ȥ����߿�)
            edge1(index_e)=[];
            edge2(index_e)=[];
        end
    end
end

clear e_edge1
clear e_edge2

%  ÿ���������Ӧ���������
for m=1:N_tetra  %  ������
    %  ÿ��Ԫ��Ӧ6�����
    e_edge1=[element(m,3),element(m,3),element(m,3),element(m,4),element(m,4),element(m,5)];
    e_edge2=[element(m,4),element(m,5),element(m,6),element(m,5),element(m,6),element(m,6)];
    %  6����ѭ��
    for n=1:6
        flag=0; %  edge1��edge2�Ҳ���
        for p=1:length(edge1)
            ref_sum=edge1(p)+edge2(p);
            ref_sub=abs(edge1(p)-edge2(p));
            curr_sum=e_edge1(n)+e_edge2(n);
            curr_sub=abs(e_edge1(n)-e_edge2(n));
            if (ref_sum==curr_sum & ref_sub==curr_sub)  %  edge���ҵ���
                flag=1;
                index_e=p; %  edge�е�λ��
                break;
            end
        end
        if (flag==0) %  edge1��edge2�Ҳ���
            Index_edge(m,n)=-1;
        else %  �ҵ�
            Index_edge(m,n)=index_e; %  edge������
        end
    end
end

%  �ն�������������
% I_index=[]; %  �У�Դ��
% J_index=[]; %  �У�����
% Value_A=[]; %  A����
% Value_B=[]; %  B����
I_index=ones(1,36*N_tetra); %  �У�Դ��
J_index=ones(1,36*N_tetra); %  �У�����
Value_A=zeros(1,36*N_tetra); %  A����
Value_B=zeros(1,36*N_tetra); %  B����

b=zeros(1,4);c=zeros(1,4);d=zeros(1,4);
edge_info=[1 2;1 3;1 4;2 3;2 4;3 4];  % 6����߱��
index_index=0;        
for m=1:N_tetra  %  ������ѭ��
    N_tetra-m
    v_m=[1 1 1 1;
         node(element(m,3),2),node(element(m,4),2),node(element(m,5),2),node(element(m,6),2);
         node(element(m,3),3),node(element(m,4),3),node(element(m,5),3),node(element(m,6),3);
         node(element(m,3),4),node(element(m,4),4),node(element(m,5),4),node(element(m,6),4)];
    v_e=1/6*det(v_m); %  ���������
    
    % �ڵ������ϵ��
    b(1)=-det([1 1 1;node(element(m,4),3),node(element(m,5),3),node(element(m,6),3);node(element(m,4),4),node(element(m,5),4),node(element(m,6),4)]);
    b(2)=det([1 1 1;node(element(m,3),3),node(element(m,5),3),node(element(m,6),3);node(element(m,3),4),node(element(m,5),4),node(element(m,6),4)]);
    b(3)=-det([1 1 1;node(element(m,3),3),node(element(m,4),3),node(element(m,6),3);node(element(m,3),4),node(element(m,4),4),node(element(m,6),4)]);
    b(4)=det([1 1 1;node(element(m,3),3),node(element(m,4),3),node(element(m,5),3);node(element(m,3),4),node(element(m,4),4),node(element(m,5),4)]);
    
    c(1)=det([1 1 1;node(element(m,4),2),node(element(m,5),2),node(element(m,6),2);node(element(m,4),4),node(element(m,5),4),node(element(m,6),4)]);
    c(2)=-det([1 1 1;node(element(m,3),2),node(element(m,5),2),node(element(m,6),2);node(element(m,3),4),node(element(m,5),4),node(element(m,6),4)]);
    c(3)=det([1 1 1;node(element(m,3),2),node(element(m,4),2),node(element(m,6),2);node(element(m,3),4),node(element(m,4),4),node(element(m,6),4)]);
    c(4)=-det([1 1 1;node(element(m,3),2),node(element(m,4),2),node(element(m,5),2);node(element(m,3),4),node(element(m,4),4),node(element(m,5),4)]);
    
    d(1)=-det([1 1 1;node(element(m,4),2),node(element(m,5),2),node(element(m,6),2);node(element(m,4),3),node(element(m,5),3),node(element(m,6),3)]);
    d(2)=det([1 1 1;node(element(m,3),2),node(element(m,5),2),node(element(m,6),2);node(element(m,3),3),node(element(m,5),3),node(element(m,6),3)]);
    d(3)=-det([1 1 1;node(element(m,3),2),node(element(m,4),2),node(element(m,6),2);node(element(m,3),3),node(element(m,4),3),node(element(m,6),3)]);
    d(4)=det([1 1 1;node(element(m,3),2),node(element(m,4),2),node(element(m,5),2);node(element(m,3),3),node(element(m,4),3),node(element(m,5),3)]);
    
    % �����Ϣ
    for p=1:6 % ���ѭ��������
        index_ef=Index_edge(m,p); %  edge�ľ����������ڼ�������ߣ�
        if (index_ef~=-1)  %  ��Ч�����
            
            sign_p=2*(edge1(index_ef)==element(m,edge_info(p,1)+2))-1;  %  ȫ������Ƿ�;ֲ����һ�£�1һ�£�-1��һ�£�
            
            xf1=node(edge1(index_ef),2);yf1=node(edge1(index_ef),3);zf1=node(edge1(index_ef),4);
            xf2=node(edge2(index_ef),2);yf2=node(edge2(index_ef),3);zf2=node(edge2(index_ef),4);
            edge_f=sqrt((xf1-xf2)^2+(yf1-yf2)^2+(zf1-zf2)^2);  %  ��߳���
            
            %  ���node������ϵ��
            bi1=b(edge_info(p,1));
            bi2=b(edge_info(p,2));
            ci1=c(edge_info(p,1));
            ci2=c(edge_info(p,2));
            di1=d(edge_info(p,1));
            di2=d(edge_info(p,2));
        end
        for q=1:6; % ���ѭ����Դ��
            index_index=index_index+1;
            index_es=Index_edge(m,q); %  edge�ľ����������ڼ�������ߣ�
            if (index_ef~=-1 & index_es~=-1)  %  ��Ч�����
                
                sign_q=2*(edge1(index_es)==element(m,edge_info(q,1)+2))-1;  %  ȫ������Ƿ�;ֲ����һ�£�1һ�£�-1��һ�£�
                
                xs1=node(edge1(index_es),2);ys1=node(edge1(index_es),3);zs1=node(edge1(index_es),4);
                xs2=node(edge2(index_es),2);ys2=node(edge2(index_es),3);zs2=node(edge2(index_es),4);
                edge_s=sqrt((xs1-xs2)^2+(ys1-ys2)^2+(zs1-zs2)^2);   %  ��߳���
                
                %  ���node������ϵ��
                bj1=b(edge_info(q,1));
                bj2=b(edge_info(q,2));
                cj1=c(edge_info(q,1));
                cj2=c(edge_info(q,2));
                dj1=d(edge_info(q,1));
                dj2=d(edge_info(q,2));

                %  ���к�
                %I_index=[I_index,index_es];
                %J_index=[J_index,index_ef];
                I_index(index_index)=index_es;
                J_index(index_index)=index_ef;
                %  A����
                value_a=4*edge_s*edge_f*v_e/((6*v_e)^4)*[(ci1*di2-di1*ci2)*(cj1*dj2-dj1*cj2)+...
                        (di1*bi2-bi1*di2)*(dj1*bj2-bj1*dj2)+(bi1*ci2-ci1*bi2)*(bj1*cj2-cj1*bj2)];
                %Value_A=[Value_A,value_a];
                Value_A(index_index)=sign_p*sign_q*value_a;
                
                %  B����
                fi2j2=bi2*bj2+ci2*cj2+di2*dj2;
                fi2j1=bi2*bj1+ci2*cj1+di2*dj1;
                fi1j2=bi1*bj2+ci1*cj2+di1*dj2;
                fi1j1=bi1*bj1+ci1*cj1+di1*dj1;
                
                sigma_i1j1=0;sigma_i1j2=0;sigma_i2j1=0;sigma_i2j2=0;
                if (edge_info(p,1)==edge_info(q,1))
                    sigma_i1j1=1;
                end
                if (edge_info(p,1)==edge_info(q,2))
                    sigma_i1j2=1;
                end
                if (edge_info(p,2)==edge_info(q,1))
                    sigma_i2j1=1;
                end
                if (edge_info(p,2)==edge_info(q,2))
                    sigma_i2j2=1;
                end
                value_b=edge_s*edge_f/((6*v_e)*120)*[(1+sigma_i1j1)*fi2j2-(1+sigma_i1j2)*fi2j1-...
                                                     (1+sigma_i2j1)*fi1j2+(1+sigma_i2j2)*fi1j1];
                %Value_B=[Value_B,value_b];
                Value_B(index_index)=sign_p*sign_q*value_b;
            end
        end
    end
end

%  ��������ֵ
A_matrix=sparse(I_index,J_index,Value_A);
B_matrix=sparse(I_index,J_index,Value_B);
result_d=eigs(A_matrix,B_matrix,8,7*7);
result_f=sqrt(result_d)

                
    
  




