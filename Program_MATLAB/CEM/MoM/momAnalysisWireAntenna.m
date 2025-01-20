% Matlab code for Pocklington Integral Equation
% Method of moment Analysis of wire antenna of length L
% March 05, 2019
% Ajeet Kumar IIST, India

clc
close all
clear all
fr = 5.7*1e9;       % Operating frequency in Hz
v0 = 3*1e8;         % velocity of light in air
lambda = v0/fr;
% L = 1;  % Length of wire antenna
L = 0.5*lambda;
a = 0.001;       % Radius of wire
N = 500;                % Number of sections
M = N;                  % Number of testing function

Ezi = 1;                % incident Electric field
ep_0 = 8.854*1e-12;     % Permittivity of free space
ep_r = 1;               % relative Permittivity
mu_0 = 4*pi*1e-7;       % Permeability of free space
mu_r = 1;               % relative Permeability
lim1 = -L/2;            % Lower limit of wire antenna placed in coordinate
lim2 = L/2;             % Upper limit of wire antenna placed in coordinate
k = 2*pi*fr*sqrt(ep_0*ep_r*mu_0*mu_r);      % Propagation constant

hz = (lim2-lim1)/(N);   % Step width of each section
fs = 1000;              % Sampling rate
z = lim1:1/fs:lim2;     % x vector along the line
z1 = lim1+hz:hz:lim2;
dz = hz;

Zmn = zeros(M,N);       % Solution matris initialization


for n = 1:N
    for m = 1:N
        zn = lim1+(n-0.5)*dz;
        za = zn-(dz/2);
        zb = zn+(dz/2);
        zm = lim1+(m-0.5)*dz;
        R = sqrt(a^2+(zm-zn)^2);
        R1n = sqrt(a^2+(zm-zn+(dz/2))^2);
        R2n = sqrt(a^2+(zm-zn-(dz/2))^2);
        t1 = (zm-zn+dz/2)*(2+1j*k*R1n)/(4*pi*R1n^2)*exp(-1j*k*R1n);
        t2 = (zm-zn-dz/2)*(2+1j*k*R2n)/(4*pi*R2n^2)*exp(-1j*k*R2n);
        phi_n = t1-t2;
        if m == n
%             if dz >= 1000*a             % If dz >> a
            Zmn(m,n) = ((k^2/(4*pi))*(2*log(dz/a)-1j*k*dz)) + phi_n;
%             else
%                 Zmn(m,n) = (k^2/(4*pi))*(log(((dz/2)+sqrt(a^2+(dz/2)^2))/((dz/2)+sqrt(a^2+(dz/2)^2)))-1j*k*dz) + phi_n;
%             end
        else
            Zmn(m,n) = (((k^2*dz)/(4*pi*R))*exp(-1j*k*R))+phi_n;
        end
           
        
        
        
    end
end

bm = -1j*2*pi*fr*ep_0*ep_r*Ezi*ones(M,1);           % Constant Matrix  (N x 1)
an = Zmn\bm;                        % Solution of unknown constant (N x 1)
Ln = length(an);
an1= abs(an);
f1 = 0;                             % Current Distribution claculation vector
fvar = 0;                           % Variable for Current distribution claulation
for n = 1:N
    zn = lim1+(n-0.5)*dz;           % nth shifted pulse basis function
    fvar = an(n)*rectpuls(z-zn,dz);
    f1 = f1+fvar;
end
Curr = abs(f1);           % ABsolute value of current Distribution

figure()
plot(z1,abs(an),'k--','LineWidth', 2);
grid on
title('Plot of solution coefficient a_n');
xlabel('dipole length (m)');
ylabel('Value a_n(n)');
legend(['N = ',num2str(N)]);

figure();
plot(z,Curr,'k--','LineWidth', 2);
grid on
hold on
title('Plot of current distribution');
xlabel('dipole length (m)');
ylabel('Current (A)');
plot(z,imag(f1),'LineWidth', 2);
plot(z,real(f1),'LineWidth', 2);
% legend(['N = ',num2str(N)]);
legend(['|I|,N =',num2str(N)],'Im[I]','Re[I]');

