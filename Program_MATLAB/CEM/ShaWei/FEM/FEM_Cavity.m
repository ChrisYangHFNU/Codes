clc;clear

%  节点和元的信息
node=load('NODE_rectangular.txt');
element=load('ELEMENT_rectangular.txt');

%  cm-->m 
node(:,2)=node(:,2)*100;
node(:,3)=node(:,3)*100;
node(:,4)=node(:,4)*100;

%  统计节点个数、四面体个数、边界三角形个数
N_node=length(node);
N_tetra=sum(element(:,2)==1);
N_tri=sum(element(:,2)==101);

%  四面体对应棱边索引
Index_edge=zeros(N_tetra,6);
edge1=[]; % 棱边顶点索引1
edge2=[]; % 棱边顶点索引2

%  记录棱边信息
for m=1:N_tetra  %  四面体
    %  每个元对应6条棱边
    e_edge1=[element(m,3),element(m,3),element(m,3),element(m,4),element(m,4),element(m,5)];
    e_edge2=[element(m,4),element(m,5),element(m,6),element(m,5),element(m,6),element(m,6)];
    %  6条边循环
    for n=1:6
        flag=0; %  edge1和edge2找不到
        for p=1:length(edge1)
            ref_sum=edge1(p)+edge2(p);
            ref_sub=abs(edge1(p)-edge2(p));
            curr_sum=e_edge1(n)+e_edge2(n);
            curr_sub=abs(e_edge1(n)-e_edge2(n));
            if (ref_sum==curr_sum & ref_sub==curr_sub)  %  edge库找到了
                flag=1;
                break;
            end
        end
        if (flag==0) %  edge1和edge2找不到(加入棱边库)
            edge1=[edge1,e_edge1(n)];
            edge2=[edge2,e_edge2(n)];
        end
    end
end

clear e_edge1
clear e_edge2

%  去掉边界棱边
for m=N_tetra+1:N_tetra+N_tri  %  边界三角形
    %  每个元对应3条棱边
    e_edge1=[element(m,3),element(m,3),element(m,4)];
    e_edge2=[element(m,4),element(m,5),element(m,5)];
    %  3条边循环
    for n=1:3
        flag=0; %  edge1和edge2找不到
        for p=1:length(edge1)
            ref_sum=edge1(p)+edge2(p);
            ref_sub=abs(edge1(p)-edge2(p));
            curr_sum=e_edge1(n)+e_edge2(n);
            curr_sub=abs(e_edge1(n)-e_edge2(n));
            if (ref_sum==curr_sum & ref_sub==curr_sub)  %  edge库找到了
                flag=1;
                index_e=p; %  edge中的位置
                break;
            end
        end
        if (flag==1) %  edge1和edge2找到(去掉棱边库)
            edge1(index_e)=[];
            edge2(index_e)=[];
        end
    end
end

clear e_edge1
clear e_edge2

%  每个四面体对应的棱边索引
for m=1:N_tetra  %  四面体
    %  每个元对应6条棱边
    e_edge1=[element(m,3),element(m,3),element(m,3),element(m,4),element(m,4),element(m,5)];
    e_edge2=[element(m,4),element(m,5),element(m,6),element(m,5),element(m,6),element(m,6)];
    %  6条边循环
    for n=1:6
        flag=0; %  edge1和edge2找不到
        for p=1:length(edge1)
            ref_sum=edge1(p)+edge2(p);
            ref_sub=abs(edge1(p)-edge2(p));
            curr_sum=e_edge1(n)+e_edge2(n);
            curr_sub=abs(e_edge1(n)-e_edge2(n));
            if (ref_sum==curr_sum & ref_sub==curr_sub)  %  edge库找到了
                flag=1;
                index_e=p; %  edge中的位置
                break;
            end
        end
        if (flag==0) %  edge1和edge2找不到
            Index_edge(m,n)=-1;
        else %  找到
            Index_edge(m,n)=index_e; %  edge的索引
        end
    end
end

%  刚度质量矩阵生成
% I_index=[]; %  行（源）
% J_index=[]; %  列（场）
% Value_A=[]; %  A矩阵
% Value_B=[]; %  B矩阵
I_index=ones(1,36*N_tetra); %  行（源）
J_index=ones(1,36*N_tetra); %  列（场）
Value_A=zeros(1,36*N_tetra); %  A矩阵
Value_B=zeros(1,36*N_tetra); %  B矩阵

b=zeros(1,4);c=zeros(1,4);d=zeros(1,4);
edge_info=[1 2;1 3;1 4;2 3;2 4;3 4];  % 6条棱边编号
index_index=0;        
for m=1:N_tetra  %  四面体循环
    N_tetra-m
    v_m=[1 1 1 1;
         node(element(m,3),2),node(element(m,4),2),node(element(m,5),2),node(element(m,6),2);
         node(element(m,3),3),node(element(m,4),3),node(element(m,5),3),node(element(m,6),3);
         node(element(m,3),4),node(element(m,4),4),node(element(m,5),4),node(element(m,6),4)];
    v_e=1/6*det(v_m); %  四面体体积
    
    % 节点基函数系数
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
    
    % 棱边信息
    for p=1:6 % 棱边循环（场）
        index_ef=Index_edge(m,p); %  edge的矩阵索引（第几条内棱边）
        if (index_ef~=-1)  %  有效内棱边
            
            sign_p=2*(edge1(index_ef)==element(m,edge_info(p,1)+2))-1;  %  全局棱边是否和局部棱边一致（1一致，-1不一致）
            
            xf1=node(edge1(index_ef),2);yf1=node(edge1(index_ef),3);zf1=node(edge1(index_ef),4);
            xf2=node(edge2(index_ef),2);yf2=node(edge2(index_ef),3);zf2=node(edge2(index_ef),4);
            edge_f=sqrt((xf1-xf2)^2+(yf1-yf2)^2+(zf1-zf2)^2);  %  棱边长度
            
            %  棱边node基函数系数
            bi1=b(edge_info(p,1));
            bi2=b(edge_info(p,2));
            ci1=c(edge_info(p,1));
            ci2=c(edge_info(p,2));
            di1=d(edge_info(p,1));
            di2=d(edge_info(p,2));
        end
        for q=1:6; % 棱边循环（源）
            index_index=index_index+1;
            index_es=Index_edge(m,q); %  edge的矩阵索引（第几条内棱边）
            if (index_ef~=-1 & index_es~=-1)  %  有效内棱边
                
                sign_q=2*(edge1(index_es)==element(m,edge_info(q,1)+2))-1;  %  全局棱边是否和局部棱边一致（1一致，-1不一致）
                
                xs1=node(edge1(index_es),2);ys1=node(edge1(index_es),3);zs1=node(edge1(index_es),4);
                xs2=node(edge2(index_es),2);ys2=node(edge2(index_es),3);zs2=node(edge2(index_es),4);
                edge_s=sqrt((xs1-xs2)^2+(ys1-ys2)^2+(zs1-zs2)^2);   %  棱边长度
                
                %  棱边node基函数系数
                bj1=b(edge_info(q,1));
                bj2=b(edge_info(q,2));
                cj1=c(edge_info(q,1));
                cj2=c(edge_info(q,2));
                dj1=d(edge_info(q,1));
                dj2=d(edge_info(q,2));

                %  行列号
                %I_index=[I_index,index_es];
                %J_index=[J_index,index_ef];
                I_index(index_index)=index_es;
                J_index(index_index)=index_ef;
                %  A矩阵
                value_a=4*edge_s*edge_f*v_e/((6*v_e)^4)*[(ci1*di2-di1*ci2)*(cj1*dj2-dj1*cj2)+...
                        (di1*bi2-bi1*di2)*(dj1*bj2-bj1*dj2)+(bi1*ci2-ci1*bi2)*(bj1*cj2-cj1*bj2)];
                %Value_A=[Value_A,value_a];
                Value_A(index_index)=sign_p*sign_q*value_a;
                
                %  B矩阵
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

%  广义特征值
A_matrix=sparse(I_index,J_index,Value_A);
B_matrix=sparse(I_index,J_index,Value_B);
result_d=eigs(A_matrix,B_matrix,8,7*7);
result_f=sqrt(result_d)

                
    
  




