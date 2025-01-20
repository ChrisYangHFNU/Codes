%% Microstrip line-fed rectangular antenna analysis using 3D FDTD code
%
% FDTD routine used to simulate the propagation of ultra wide-band pulse
% through line-fed rectangular antenna to compute return loss (S_{11} 
% scattering coefficient). 
%
function FDTD_3D_Microstrip_Antenna
close all; clear; clc;
%% Physical constants
   epsilon0 = 8.85418782e-12; mu0 = 1.25663706e-6;
   c = 1.0/sqrt(mu0*epsilon0);

%% Gaussian half-width
   t_half = 15.0e-12;

%% Microstrip fed parameters
   lineW = 1.6604e-3; 
   lineH = 0.4000e-3;
   % Roger's 5880 Duroid parameters
   lineEr = 2.2;      % eps_r 
   lineTan = 0.00041; % loss tangent 
   
%% End time
   t_end = 2.0e-9;

%% Total mesh dimensions and grid cells sizes (without PML)
   nx = 70; ny = 100; nz = 6;
   dx = 0.2372e-3; dy = 0.2265e-3; dz = 0.1000e-3;

%% Number of PML layers (ten or more!)
   PML = 10;

% Add PML layers
   nx = nx + 2*PML; ny = ny + 2*PML; nz = nz + 2*PML;
% Calculate dt, number of FDTD iterations and Z0 microstrip line matched load    
   dt = (1.0/c/sqrt( 1.0/(dx^2) + 1.0/(dy^2) + 1.0/(dz^2)))*0.999;
   number_of_iterations = ceil(t_end/dt);
   eps_eff = 0.5*(lineEr+1) + 0.5*(lineEr-1)/sqrt(1+12*lineH/lineW);
   Z0 = 377/sqrt(eps_eff)*1/(lineW/lineH + 1.393 + 0.667*log(lineW/lineH) + 1.444);

%% Matrix of material's constants
   number_of_materials = 4;
   % For material of number x = 1,2,3... :
   % Material(x,1) - relative permittivity, Material(x,2) - relative permeability,
   % Material(x,3) - specific conductivity
   % Vacuum
   Material(1,1) = 1.0;   Material(1,2) = 1.0;   Material(1,3) = 0.0;
   % Metal (Copper)
   Material(2,1) = 1.0;   Material(2,2) = 1.0;   Material(2,3) = 5.88e+7;
   % Substrate material (RT/Duroid 5880)
   Material(3,1) = lineEr;   Material(3,2) = 1.0;
   % Calculate conductivity of Duroid at 20 GHz from loss tangent and eps_r
   Material(3,3) = 2*pi*10e9*lineTan*lineEr*epsilon0;
   % Matched load material is calculated from transmission line parameters
   Material(4,1) = 1.0;   Material(4,2) = 1.0;   Material(4,3) = lineH/(Z0*lineW*dy);

%% 3D array for geometry
   Index = ones(nx, ny, nz);

%% Define microstrip antenna geometry
   % Ground plane 
   Index((1+PML):(nx-PML), (1+PML):(ny-PML), PML+1) = 2;
   % Rectangular patch with feed (thickness is equal to one cell)
   Index((nx/2-25):(nx/2+25), (ny-PML-50):(ny-PML-10), PML+6) = 2;
   Index((nx/2-4):(nx/2+4), (ny-PML-50):(ny-PML-45), PML+6) = 1; 
   Index((nx/2-3):(nx/2+3), (PML+1):(ny-PML-40), PML+6) = 2;
   % Dielectric substrate between ground plane and filter patch
   Index((1+PML):(nx-PML), (1+PML):(ny-PML), (PML+2):(PML+5)) = 3;
   % Matched load before port 1
   Index((nx/2-3):(nx/2+3), PML+1, (PML+2):(PML+5)) = 4;
   
%% 3D FDTD physical (fields) and additional arrays are defined as 'single' 
%% to increase performance
   Ex = zeros(nx, ny+1, nz+1, 'single'); 
   Gx = zeros(nx, ny+1, nz+1, 'single'); 
   Fx = zeros(nx, ny+1, nz+1, 'single');  
   Ey = zeros(nx+1, ny, nz+1, 'single'); 
   Gy = zeros(nx+1, ny, nz+1, 'single'); 
   Fy = zeros(nx+1, ny, nz+1, 'single');
   Ez = zeros(nx+1, ny+1, nz, 'single'); 
   Gz = zeros(nx+1, ny+1, nz, 'single'); 
   Fz = zeros(nx+1, ny+1, nz, 'single');
   Hx = zeros(nx+1, ny, nz, 'single'); 
   Bx = zeros(nx+1, ny, nz, 'single'); 
   Hy = zeros(nx, ny+1, nz, 'single');
   By = zeros(nx, ny+1, nz, 'single'); 
   Hz = zeros(nx, ny, nz+1, 'single'); 
   Bz = zeros(nx, ny, nz+1, 'single');

%% FDTD PML coefficients arrays (also 'singles'). 
%% They are already filled with values corresponding to free space 
   m = 4; ka_max = 1.0; R_err = 1.0e-16;
   eta = sqrt(mu0/epsilon0*Material(1,1)/Material(1,2));
   k_Ex_c = ones(nx, ny, nz, 'single')*2.0*epsilon0;  
   k_Ex_d = ones(nx, ny, nz, 'single')*(-2.0*epsilon0);
   k_Ey_a = ones(nx+1, ny, nz, 'single');
   k_Ey_b = ones(nx+1, ny, nz, 'single')/(2.0*epsilon0);
   k_Gz_a = ones(nx+1, ny, nz, 'single');
   k_Gz_b = ones(nx+1, ny, nz, 'single');
   k_Hy_a = ones(nx, ny, nz, 'single'); 
   k_Hy_b = ones(nx, ny, nz, 'single')/(2.0*epsilon0);
   k_Hx_c = ones(nx+1, ny, nz, 'single')*2.0*epsilon0/mu0;
   k_Hx_d = ones(nx+1, ny, nz, 'single')*(-2.0*epsilon0/mu0);
   k_Bz_a = ones(nx, ny, nz, 'single');
   k_Bz_b = ones(nx, ny, nz, 'single')*dt;
   k_Gx_a = ones(nx, ny+1, nz, 'single');
   k_Gx_b = ones(nx, ny+1, nz, 'single');
   k_Ey_c = ones(nx, ny, nz, 'single')*2.0*epsilon0; 
   k_Ey_d = ones(nx, ny, nz, 'single')*(-2.0*epsilon0);
   k_Ez_a = ones(nx, ny+1, nz, 'single'); 
   k_Ez_b = ones(nx, ny+1, nz, 'single')/(2.0*epsilon0);
   k_Bx_a = ones(nx, ny, nz, 'single'); 
   k_Bx_b = ones(nx, ny, nz, 'single')*dt;
   k_Hy_c = ones(nx, ny+1, nz, 'single')*2.0*epsilon0/mu0;
   k_Hy_d = ones(nx, ny+1, nz, 'single')*(-2.0*epsilon0/mu0);
   k_Hz_a = ones(nx, ny, nz, 'single');
   k_Hz_b = ones(nx, ny, nz, 'single')/(2.0*epsilon0);
   k_Ex_a = ones(nx, ny, nz+1, 'single');
   k_Ex_b = ones(nx, ny, nz+1, 'single')/(2.0*epsilon0);
   k_Gy_a = ones(nx, ny, nz+1, 'single'); 
   k_Gy_b = ones(nx, ny, nz+1, 'single');
   k_Ez_c = ones(nx, ny, nz, 'single')*2.0*epsilon0;
   k_Ez_d = ones(nx, ny, nz, 'single')*(-2.0*epsilon0);
   k_Hx_a = ones(nx, ny, nz, 'single');
   k_Hx_b = ones(nx, ny, nz, 'single')/(2.0*epsilon0);
   k_By_a = ones(nx, ny, nz, 'single');  
   k_By_b = ones(nx, ny, nz, 'single')*dt;
   k_Hz_c = ones(nx, ny, nz+1, 'single')*2.0*epsilon0/mu0;   
   k_Hz_d = ones(nx, ny, nz+1, 'single')*(-2.0*epsilon0/mu0);

%% General FDTD coefficients 
   I = 1:number_of_materials;
   K_a(I) = (2.0*epsilon0*Material(I,1) - Material(I,3)*dt)./...
            (2.0*epsilon0*Material(I,1) + Material(I,3)*dt);
   K_b(I) = 2.0*dt./(2.0*epsilon0*Material(I,1) + Material(I,3)*dt);
   K_c(I) = Material(I,2);
   Ka = single(K_a(Index)); 
   Kb = single(K_b(Index)); 
   Kc = single(K_c(Index));  
   
%% PML coefficients along x-axis
   sigma_max = -(m + 1.0)*log(R_err)/(2.0*eta*PML*dx);
   for I=0:(PML-1)
        sigma_x = sigma_max*((PML - I)/PML)^m;
        ka_x = 1.0 + (ka_max - 1.0)*((PML - I)/PML)^m;
        k_Ey_a(I+1,:,:) = (2.0*epsilon0*ka_x - sigma_x*dt)/...
                          (2.0*epsilon0*ka_x + sigma_x*dt);
        k_Ey_a(nx-I,:,:) = k_Ey_a(I+1,:,:);
        k_Ey_b(I+1,:,:) = 1.0/(2.0*epsilon0*ka_x + sigma_x*dt);
        k_Ey_b(nx-I,:,:) = k_Ey_b(I+1,:,:);
        k_Gz_a(I+1,:,:) = (2.0*epsilon0*ka_x - sigma_x*dt)/...
                          (2.0*epsilon0*ka_x + sigma_x*dt);
        k_Gz_a(nx-I,:,:) = k_Gz_a(I+1,:,:);
        k_Gz_b(I+1,:,:) = 2.0*epsilon0/(2.0*epsilon0*ka_x + sigma_x*dt);
        k_Gz_b(nx-I,:,:) = k_Gz_b(I+1,:,:);
        k_Hx_c(I+1,:,:) = (2.0*epsilon0*ka_x + sigma_x*dt)/mu0;
        k_Hx_c(nx-I,:,:) = k_Hx_c(I+1,:,:);
        k_Hx_d(I+1,:,:) = -(2.0*epsilon0*ka_x - sigma_x*dt)/mu0;
        k_Hx_d(nx-I,:,:) = k_Hx_d(I+1,:,:);

        sigma_x = sigma_max*((PML - I - 0.5)/PML)^m;
        ka_x = 1.0 + (ka_max - 1.0)*((PML - I - 0.5)/PML)^m;
        k_Ex_c(I+1,:,:) = 2.0*epsilon0*ka_x + sigma_x*dt;
        k_Ex_c(nx-I-1,:,:) = k_Ex_c(I+1,:,:);
        k_Ex_d(I+1,:,:) = -(2.0*epsilon0*ka_x - sigma_x*dt);
        k_Ex_d(nx-I-1,:,:) = k_Ex_d(I+1,:,:);
        k_Hy_a(I+1,:,:) = (2.0*epsilon0*ka_x - sigma_x*dt)/...
                          (2.0*epsilon0*ka_x + sigma_x*dt);
        k_Hy_a(nx-I-1,:,:) = k_Hy_a(I+1,:,:);
        k_Hy_b(I+1,:,:) = 1.0/(2.0*epsilon0*ka_x + sigma_x*dt);
        k_Hy_b(nx-I-1,:,:) = k_Hy_b(I+1,:,:);
        k_Bz_a(I+1,:,:) = (2.0*epsilon0*ka_x - sigma_x*dt)/...
                          (2.0*epsilon0*ka_x + sigma_x*dt);
        k_Bz_a(nx-I-1,:,:) = k_Bz_a(I+1,:,:);
        k_Bz_b(I+1,:,:) = 2.0*epsilon0*dt/(2.0*epsilon0*ka_x + sigma_x*dt);
        k_Bz_b(nx-I-1,:,:) = k_Bz_b(I+1,:,:);
   end

%% PML coefficients along y-axis
   sigma_max = -(m + 1.0)*log(R_err)/(2.0*eta*PML*dy);
   for J=0:(PML-1)
        sigma_y = sigma_max*((PML - J)/PML)^m;
        ka_y = 1.0 + (ka_max - 1.0)*((PML - J)/PML)^m;
        k_Gx_a(:,J+1,:) = (2.0*epsilon0*ka_y - sigma_y*dt)/...
                          (2.0*epsilon0*ka_y + sigma_y*dt);
        k_Gx_a(:,ny-J,:) = k_Gx_a(:,J+1,:);
        k_Gx_b(:,J+1,:) = 2.0*epsilon0/(2.0*epsilon0*ka_y + sigma_y*dt);
        k_Gx_b(:,ny-J,:) = k_Gx_b(:,J+1,:);
        k_Ez_a(:,J+1,:) = (2.0*epsilon0*ka_y - sigma_y*dt)/...
                          (2.0*epsilon0*ka_y + sigma_y*dt);
        k_Ez_a(:,ny-J,:) = k_Ez_a(:,J+1,:);
        k_Ez_b(:,J+1,:) = 1.0/(2.0*epsilon0*ka_y + sigma_y*dt);
        k_Ez_b(:,ny-J,:) = k_Ez_b(:,J+1,:);
        k_Hy_c(:,J+1,:) = (2.0*epsilon0*ka_y + sigma_y*dt)/mu0;
        k_Hy_c(:,ny-J,:) = k_Hy_c(:,J+1,:);
        k_Hy_d(:,J+1,:) = -(2.0*epsilon0*ka_y - sigma_y*dt)/mu0;
        k_Hy_d(:,ny-J,:) = k_Hy_d(:,J+1,:);

        sigma_y = sigma_max*((PML - J - 0.5)/PML)^m;
        ka_y = 1.0 + (ka_max - 1.0)*((PML - J - 0.5)/PML)^m;
        k_Ey_c(:,J+1,:) = 2.0*epsilon0*ka_y+sigma_y*dt;
        k_Ey_c(:,ny-J-1,:) = k_Ey_c(:,J+1,:);
        k_Ey_d(:,J+1,:) = -(2.0*epsilon0*ka_y-sigma_y*dt);
        k_Ey_d(:,ny-J-1,:) = k_Ey_d(:,J+1,:);
        k_Bx_a(:,J+1,:) = (2.0*epsilon0*ka_y-sigma_y*dt)/...
                          (2.0*epsilon0*ka_y+sigma_y*dt);
        k_Bx_a(:,ny-J-1,:) = k_Bx_a(:,J+1,:);
        k_Bx_b(:,J+1,:) = 2.0*epsilon0*dt/(2.0*epsilon0*ka_y+sigma_y*dt);
        k_Bx_b(:,ny-J-1,:) = k_Bx_b(:,J+1,:);
        k_Hz_a(:,J+1,:) = (2.0*epsilon0*ka_y-sigma_y*dt)/...
                          (2.0*epsilon0*ka_y+sigma_y*dt);
        k_Hz_a(:,ny-J-1,:) = k_Hz_a(:,J+1,:);
        k_Hz_b(:,J+1,:) = 1.0/(2.0*epsilon0*ka_y+sigma_y*dt);
        k_Hz_b(:,ny-J-1,:) = k_Hz_b(:,J+1,:);
   end

%% PML coefficients along z-axis 
   sigma_max = -(m + 1.0)*log(R_err)/(2.0*eta*PML*dz);
   for K=0:(PML-1)
        sigma_z = sigma_max*((PML - K)/PML)^m;
        ka_z = 1.0 + (ka_max - 1.0)*((PML-K)/PML)^m;
        k_Ex_a(:,:,K+1) = (2.0*epsilon0*ka_z - sigma_z*dt)/...
                          (2.0*epsilon0*ka_z + sigma_z*dt);
        k_Ex_a(:,:,nz-K) = k_Ex_a(:,:,K+1);
        k_Ex_b(:,:,K+1) = 1.0/(2.0*epsilon0*ka_z + sigma_z*dt);
        k_Ex_b(:,:,nz-K) = k_Ex_b(:,:,K+1);
        k_Gy_a(:,:,K+1) = (2.0*epsilon0*ka_z - sigma_z*dt)/...
                          (2.0*epsilon0*ka_z + sigma_z*dt);
        k_Gy_a(:,:,nz-K) = k_Gy_a(:,:,K+1);
        k_Gy_b(:,:,K+1) = 2.0*epsilon0/(2.0*epsilon0*ka_z + sigma_z*dt);
        k_Gy_b(:,:,nz-K) = k_Gy_b(:,:,K+1);
        k_Hz_c(:,:,K+1) = (2.0*epsilon0*ka_z + sigma_z*dt)/mu0;
        k_Hz_c(:,:,nz-K) = k_Hz_c(:,:,K+1);
        k_Hz_d(:,:,K+1) = -(2.0*epsilon0*ka_z - sigma_z*dt)/mu0;
        k_Hz_d(:,:,nz-K) = k_Hz_d(:,:,K+1);

        sigma_z = sigma_max*((PML - K - 0.5)/PML)^m;
        ka_z = 1.0 + (ka_max - 1.0)*((PML - K - 0.5)/PML)^m;
        k_Ez_c(:,:,K+1) = 2.0*epsilon0*ka_z + sigma_z*dt;
        k_Ez_c(:,:,nz-K-1) = k_Ez_c(:,:,K+1);
        k_Ez_d(:,:,K+1) = -(2.0*epsilon0*ka_z - sigma_z*dt);
        k_Ez_d(:,:,nz-K-1) = k_Ez_d(:,:,K+1);
        k_Hx_a(:,:,K+1) = (2.0*epsilon0*ka_z - sigma_z*dt)/...
                          (2.0*epsilon0*ka_z + sigma_z*dt);
        k_Hx_a(:,:,nz-K-1) = k_Hx_a(:,:,K+1);
        k_Hx_b(:,:,K+1) = 1.0/(2.0*epsilon0*ka_z + sigma_z*dt);
        k_Hx_b(:,:,nz-K-1) = k_Hx_b(:,:,K+1);
        k_By_a(:,:,K+1) = (2.0*epsilon0*ka_z - sigma_z*dt)/...
                          (2.0*epsilon0*ka_z + sigma_z*dt);
        k_By_a(:,:,nz-K-1) = k_By_a(:,:,K+1);
        k_By_b(:,:,K+1) = 2.0*epsilon0*dt/(2.0*epsilon0*ka_z + sigma_z*dt);
        k_By_b(:,:,nz-K-1) = k_By_b(:,:,K+1);
   end
 
%% Variables for stored values
   Time = []; Incident = []; Port_1 = [];
   
%% Main 3D FDTD+UPML routine (operates with 'singles' to increase speed) 
   figure('units','normalized','outerposition',[0 0 1 1]);
   set(gcf, 'doublebuffer', 'on');    
   tic;
   for T=0:(number_of_iterations-1)
        %% Calculate Fx -> Gx -> Ex
        I = 1:nx; J = 2:ny; K = 2:nz;
        Fx_r = Fx(I,J,K);
        Fx(I,J,K) = Ka(I,J,K).*Fx(I,J,K) + Kb(I,J,K).*...
                    ((Hz(I,J,K) - Hz(I,J-1,K))/dy - (Hy(I,J,K) - Hy(I,J,K-1))/dz);
        Gx_r = Gx(I,J,K);
        Gx(I,J,K) = k_Gx_a(I,J,K).*Gx(I,J,K) + k_Gx_b(I,J,K).*(Fx(I,J,K) - Fx_r);
        Ex(I,J,K) = k_Ex_a(I,J,K).*Ex(I,J,K) + k_Ex_b(I,J,K).*...
                    (k_Ex_c(I,J,K).*Gx(I,J,K) + k_Ex_d(I,J,K).*Gx_r);

        %% Calculate Fy -> Gy -> Ey
        I = 2:nx; J = 1:ny; K = 2:nz;
        Fy_r = Fy(I,J,K);
        Fy(I,J,K) = Ka(I,J,K).*Fy(I,J,K) + Kb(I,J,K).*...
                    ((Hx(I,J,K) - Hx(I,J,K-1))/dz - (Hz(I,J,K) - Hz(I-1,J,K))/dx);
        Gy_r = Gy(I,J,K);
        Gy(I,J,K) = k_Gy_a(I,J,K).*Gy(I,J,K) + k_Gy_b(I,J,K).*(Fy(I,J,K) - Fy_r);
        Ey(I,J,K) = k_Ey_a(I,J,K).*Ey(I,J,K) + k_Ey_b(I,J,K).*...
                    (k_Ey_c(I,J,K).*Gy(I,J,K) + k_Ey_d(I,J,K).*Gy_r);

        %% Calculate Fz -> Gz -> Ez
        I = 2:nx; J = 2:ny; K = 1:nz;
        Fz_r = Fz(I,J,K);
        Fz(I,J,K) = Ka(I,J,K).*Fz(I,J,K) + Kb(I,J,K).*...
                    ((Hy(I,J,K) - Hy(I-1,J,K))/dx - (Hx(I,J,K) - Hx(I,J-1,K))/dy);
        Gz_r = Gz(I,J,K);
        Gz(I,J,K) = k_Gz_a(I,J,K).*Gz(I,J,K) + k_Gz_b(I,J,K).*(Fz(I,J,K) - Fz_r);
        Ez(I,J,K) = k_Ez_a(I,J,K).*Ez(I,J,K) + k_Ez_b(I,J,K).*...
                    (k_Ez_c(I,J,K).*Gz(I,J,K) + k_Ez_d(I,J,K).*Gz_r);
       
        %% Source of vertical electric field Ez applied to the whole face 
        %% plane at the port 1 during short time interval
        if (T*dt<=8.0*t_half)
            Ez((nx/2-3):(nx/2+3), PML+2, (PML+2):(PML+5)) = Source(t_half, 1.0, T*dt);
        end
        
        %% Save reflected Ez at port 1
        if ( mod(T,4) == 0 )
            Time = [Time, T*dt];
            Incident = [Incident, Source(t_half, 1.0, T*dt)]; 
            % Remove incident signal from Port 1
            Port_1 = [Port_1, Ez(nx/2, PML+2, PML+3)];
        end
        
        %% Calculate Bx -> Hx
        I = 2:nx; J = 1:ny; K = 1:nz;
        Bx_r = Bx(I,J,K);
        Bx(I,J,K) = k_Bx_a(I,J,K).*Bx(I,J,K) + k_Bx_b(I,J,K).*...
                    ((Ey(I,J,K+1) - Ey(I,J,K))/dz - (Ez(I,J+1,K) - Ez(I,J,K))/dy);
        Hx(I,J,K) = k_Hx_a(I,J,K).*Hx(I,J,K) + k_Hx_b(I,J,K).*...
                    (k_Hx_c(I,J,K).*Bx(I,J,K) + k_Hx_d(I,J,K).*Bx_r)./Kc(I,J,K);
 
        %% Calculate By -> Hy
        I = 1:nx; J = 2:ny; K = 1:nz;
        By_r = By(I,J,K);
        By(I,J,K) = k_By_a(I,J,K).*By(I,J,K) + k_By_b(I,J,K).*...
                    ((Ez(I+1,J,K) - Ez(I,J,K))/dx - (Ex(I,J,K+1) - Ex(I,J,K))/dz);
        Hy(I,J,K) = k_Hy_a(I,J,K).*Hy(I,J,K) + k_Hy_b(I,J,K).*...
                    (k_Hy_c(I,J,K).*By(I,J,K) + k_Hy_d(I,J,K).*By_r)./Kc(I,J,K);

        %% Calculate Bz -> Hz
        I = 1:nx; J = 1:ny; K = 2:nz;
        Bz_r = Bz(I,J,K);
        Bz(I,J,K) = k_Bz_a(I,J,K).*Bz(I,J,K) + k_Bz_b(I,J,K).*...
                    ((Ex(I,J+1,K) - Ex(I,J,K))/dy - (Ey(I+1,J,K) - Ey(I,J,K))/dx);
        Hz(I,J,K) = k_Hz_a(I,J,K).*Hz(I,J,K) + k_Hz_b(I,J,K).*...
                    (k_Hz_c(I,J,K).*Bz(I,J,K) + k_Hz_d(I,J,K).*Bz_r)./Kc(I,J,K);
                
        %% Life graphics
        if ( mod(T, fix(number_of_iterations/50)) == 0 || T==number_of_iterations )
            subplot(2,2,[2,4]);
            contour(1e3*dx*(0:(nx-2*PML-1)), 1e3*dy*(0:(ny-2*PML-1)), Index((1+PML):(nx-PML), (1+PML):(ny-PML), PML+6)',1);
            hold on;
            pcolor(1e3*dx*(0:(nx-2*PML-1)), 1e3*dy*(0:(ny-2*PML-1)), double(Ez((1+PML):(nx-PML), (1+PML):(ny-PML), PML+3)'));
            shading interp;
            xlabel('[mm]',  'FontSize', 16);
            ylabel('[mm]',  'FontSize', 16);
            caxis([-1 1]);
            colorbar;
            axis image;
            title(['E_z at t = ',num2str(ceil(T*dt*1e12)),' ps'], 'FontSize', 18);
            drawnow;
        end
    end
    toc
    
%% Plot signals and return loss
   % Input and reflected (at port 1) signal
     subplot(2,2,1);
     plot(Time, Incident, 'r-', Time, Port_1, 'k-', 'LineWidth', 2);
     title('Input and Port 1 signals', 'FontSize', 18);
     xlabel('Time, [ns]', 'FontSize', 16);
     ylabel('Strength, [V/cm]', 'FontSize', 16);
     legend('Input' ,'Port 1');
     grid on;
   % |S_{11}|
     subplot(2,2,3);
     NFFT = 2^nextpow2(length(Incident));
     S_11 = abs(fft(Port_1-Incident, NFFT))./abs(fft(Incident, NFFT));
     f = (0.5/(Time(2)-Time(1)))*linspace(0, 1, NFFT/2);
     plot(f*1e-9, 20*log10(S_11(1:NFFT/2)), 'k-', 'LineWidth', 2);
     title('Return loss from antenna - |S_{11}|', 'FontSize', 18);
     xlabel('Frequency, [GHz]', 'FontSize', 16);
     ylabel('Magnitude, [dB]', 'FontSize', 16);
     xlim([0 40]);
     grid on;
     drawnow;

%% Gauss function for voltage source
function [res] = Source(t_half, amplitude, t)
% Pulse delay
t0 = 3.0*t_half;
res = amplitude*exp(-0.5*((t-t0)/(t_half/2.35482))^2.0);
