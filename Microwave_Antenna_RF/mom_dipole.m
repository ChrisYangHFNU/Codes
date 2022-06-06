% Ӧ�÷�Ӧ���ַ���(1)�����÷ֶ����һ�����٤�ɽ𷨻�ȡ�Գ����ӵ������迹�����ߵ����ֲ��뷽��ͼ
clear;clc;
global beta a zzm zzn Deltz;
C = 3.0e8;
f = 3.0e8;
beta = 2*pi*f/C; % ��λ����
a = 0.0005; % ���߰뾶
L = [0.1:.01:2.0]; % �Գ����ӵ糤�� L=2l ==>0.1*lambda~2*lambda�����0.01*lambda
len = length(L);
ZA = zeros(1,len); % �����迹������ʼֵ����
for num=1:len
    if L(1,num)<=1.0  % ���L<=1.0*lambda
       N = 11; % ���Ӿ��ȷָ�Ϊ12���ӶΣ�11���ֶ����һ�������Խ��Щ�Ӷ�
    else
       N = 23;  % ���Ӿ��ȷָ�Ϊ24���ӶΣ�23���ֶ����һ�������Խ��Щ�Ӷ�
    end
    Deltz = L(1,num)/(N+1); %�Ӷγ���
    zm = [(N-1)/2:-1:-(N-1)/2]*Deltz; % ��������λ������ԭ�㣬��z����ã��ֶνڵ����꣨�����˵㣩
    % =======================================================
    % ��ȡ���廥�迹����[Zmn]
    Zmn = zeros(N,N); % [Zmn]��ʼ��
    zzm = zm(1,1); % �����˵㣬���϶˷ֶνڵ�����
    for n=1:N
        zzn = zm(1,n); % �����˵㣬���ϵ��£����ֶνڵ�����
        Zmn(1,n) = -j*30/(sin(beta*Deltz)^2)*(quad('Fun_Zmn01',zzn,zzn+Deltz)+quad('Fun_Zmn02',zzn-Deltz,zzn));
    end
    Zmn = toeplitz(Zmn(1,:),Zmn(1,:)); % [Zmn]ΪToeplitz����
    % ============================
    % ��ȡ�����ѹ����[Vm]
    VA = 1.0; % ��Ӽ�����ѹ
    Vm = zeros(N,1); % [Vm]��ʼ��
    % ��϶�ļ���ģ��(slice generator)
    Vm((N+1)/2,1) = -1/sin(beta*Deltz)*(VA/Deltz)*(quad('Fun_Vm01',0,Deltz)+quad('Fun_Vm02',-Deltz,0));
    % ===============================
    % ��ȡ����������� [In]
    In = zeros(N,1); % [In]��ʼ��
    In = Zmn\Vm;
    % ===============================
    ZA(1,num) = VA/In((N+1)/2,1);
end

figure(1)
plot(L,real(ZA),'r.-'); % �������������Գ����ӵ糤��L�ı仯����
axis([0 2.0 0 2500]);
xlabel('2l/\lambda');
ylabel('�������\Omega');
hold on;
figure(2)
plot(L,imag(ZA),'b.-'); % ��������翹��Գ����ӵ糤��L�ı仯����
axis([0 2.0 -1500 1200]);
xlabel('2l/\lambda');
ylabel('����翹\Omega');
hold on;
% =================================
% �糤��L=2.0ʱ�����ƶԳ��������ߵ����ֲ�ͼ
% ���ƹ����������ߵ����ֲ��ķֶ����һ�����ͼ
pw = 50; % ���Ӿ��ȷָ�ΪN+1�ӶΣ����Ӷξ���ϸ��Ϊ50�ֶ�
Sn = zeros((N+1)*pw,N); % Ϊ���ڻ���������ߵ����ֲ����ع���ֶ����ҵ������������󣬳�ֵ����
for n=1:N
    cn = zm(1,n); % ��n�ڵ�����Zn
    zp = linspace(cn-Deltz,cn+Deltz,2*pw)'; % ���˵㣬���µ��ϣ����Ӷ�ϸ�ֽڵ���������
    % �����������������n��
    Sn((n-1)*pw+1:(n+1)*pw,n) = In(n,1)*sin(beta*(Deltz-abs(zp-cn)))/sin(beta*Deltz);
end
% ====================================
figure(3)  % ���Ʒֶ����ҵ���������ͼ
for n=1:N
    plot([1:(N+1)*pw]*(Deltz/pw),abs(Sn(:,n)),'k-'); % ���Ƶ�n�ֶ����ҵ���������ͼ
    hold on;
end
I = zeros((N+1)*pw,1); % �Գ��������ߵ����ֲ���������ֵ����
I = sum(Sn,2); % �ɷֶ����ҵ������������󣬻�ȡ�Գ��������ߵ����ֲ�����
plot([pw:N*pw]*(Deltz/pw),abs(I(pw:N*pw)),'r.-'); % ���ƶԳ��������ߵ����ֲ�ͼ
set(gca,'XDir','reverse'); % �����귴ת
% ====================================
% �糤��L=2.0ʱ�����ƶԳ����ӷ���ͼ��ͨ���Ƚϣ�˵��������ģ�ͣ�����פ�������ֲ����ǶԳ�����Զ�����Է���������ͼ�������ý���
Theta = 0:pi/180:2*pi; % Ϊ����E�淽��ͼ��ȡ����� 0<theta<2*pi
Len = length(Theta);
F = zeros(1,Len); % ��������������������ֵ����
Fp = zeros(1,Len); % ������ģ�ͷ���ͼ������������ֵ����
% ��ȡ������������������Ϊ�ˣ���ԭ�Գ�������Ϊ�ɷֶ����ҵ����������̶Գ����ӹ��ɵĲ�����ֱ�����У��䷽��������ȡ�μ��̲�
for n=1:N % �����n�ֶ����ҵ����������̶Գ����ӵķ�����ֵ
    cn = zm(1,n); % ��n�ֶ����ҵ����������̶Գ����ӵ��е㣨��Ϊ����㣩����
    for m=1:Len % �����n�ֶ����ҵ����������̶Գ�������Theta���ķ�����ֵ
        theta = (m-1)*pi/180;
        if sin(theta)~=0 % ��ʽ������ ��sin(theta)=0,��F(theta)=0
            F(1,m) = F(1,m)+j*(60/sin(beta*Deltz))*((cos(beta*Deltz*cos(theta))-cos(beta*Deltz))/sin(theta))*In(n,1)*exp(j*beta*cn*cos(theta));
        end
    end
end
% ====== ��ȡ������ģ�ͷ������������μ��̲�(l=L/2=1.0lambda) =======
for m=1:Len
    theta = (m-1)*pi/180;
    if sin(theta)~=0
       Fp(1,m) = Fp(1,m)+(cos(beta*cos(theta))-cos(beta))/sin(theta);
    end
end
figure(4)
polar(Theta,abs(Fp)/max(abs(Fp)),'ko'); % ���ƴ�����ģ�͹�һ������ͼ
hold on;
polar(Theta(1,1:3:Len),abs(F(1,1:3:Len))/max(abs(F)),'r-'); % ���ƾ���������ͼ
hold on;

