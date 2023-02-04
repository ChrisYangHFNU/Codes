% Chapter1_GIFanimation.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees = pi/180;

% DASHBOARD GRID
lam0  = 0.1;
theta = 30*degrees;
Sx    = 1;
Sy    = 1;
Nx    = 100;
Ny    = 100;

NFRAMES = 40;
gif_name = 'FDFD_animation.gif';

% CALCULATE MESHGRID
dx = Sx/Nx;
xa = [1:Nx]*dx;
xa = xa - mean(xa);
dy = Sy/Ny;
ya = [1:Ny]*dy;
ya = ya - mean(ya);
[Y,X] = meshgrid(ya,xa);

% CALCULATE A WAVE (REPLACES FDFD SIMULATION)
k0 = 2*pi/lam0;
kx = k0*sin(theta);
ky = k0*cos(theta);
f  = exp(-1i*(kx*X + ky*Y));

% VISUALIZE USING IMAGESC
for nframe = 1 : NFRAMES
    
    % Add Phase
    phase = 2*pi*nframe/NFRAMES;
    fi = f*exp(1i*phase);
    
    % Draw Field
    pcolor(xa,ya,real(fi).');
    shading interp;
    axis equal tight off;
    set(gca,'YDir','reverse');
    drawnow;
    
    % Capture Frame
    F = getframe(gca);
    F = frame2im(F);
    [ind,cmap] = rgb2ind(F,256,'nodither');
    
    % Add Frame to GIF
    if nframe == 1
        imwrite(ind,cmap,gif_name,'gif',...
                'DelayTime',0,'Loopcount',inf);
    else
        imwrite(ind,cmap,gif_name,'gif',...
                'DelayTime',0,'WriteMode','append');
    end
end