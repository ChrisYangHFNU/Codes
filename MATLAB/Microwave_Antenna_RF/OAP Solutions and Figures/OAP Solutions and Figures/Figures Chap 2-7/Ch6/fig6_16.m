%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.16
% Xiaomin Lu 6/26/00
% updated by K. Bell 11/19/00, 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
Vm = ones(N,1);
u = -1:1/1000:1;

k1 = 1;
for INR = 10.^([10 20]/10)
   for ui = [0.0433 0.02]
      Vi = exp(j*n*pi*ui);
      Sn = INR*Vi*Vi' + eye(N);
      W = inv(Sn)*Vm/(Vm'*inv(Sn)*Vm);
      beam(k1,:) = W'*exp(j*n*pi*u);
      k1 = k1 + 1;
   end
end

beam = 20*log10(abs(beam));
subplot(2,2,1)
plot(u,beam(1,:))
xlabel('{\itu}')
ylabel('Beam pattern (dB)')
%title('\sigma_I^2 = 10dB, ui = 0.0433')
axis([-0.2 0.2 -30 20])
line([0.0433 0.0433],[-30 15])
grid

subplot(2,2,2)
plot(u,beam(2,:))
xlabel('{\itu}')
ylabel('Beam pattern (dB)')
%title('\sigma_I^2 = 10dB, ui = 0.02')
axis([-0.2 0.2 -30 20])
line([0.02 0.02],[-30 15])
grid

subplot(2,2,3)
plot(u,beam(3,:))
xlabel('{\itu}')
ylabel('Beam pattern (dB)')
%title('\sigma_I^2 = 20dB, ui = 0.0433')
axis([-0.2 0.2 -30 20])
grid
line([0.0433 0.0433],[-30 15])

subplot(2,2,4)
plot(u,beam(4,:))
xlabel('{\itu}')
ylabel('Beam pattern (dB)')
%title('\sigma_I^2 = 20dB, ui = 0.02')
axis([-0.2 0.2 -30 20])
grid
line([0.02 0.02],[-30 15])

         




