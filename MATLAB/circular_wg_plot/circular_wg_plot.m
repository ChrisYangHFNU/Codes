close all
lambda = 1;
v = physconst('LightSpeed');
k = 2*pi/lambda;
eps0 = 8.854187817e-12;
mu0 = 4*pi*1e-7;
w = 2*pi*v/lambda; 
k0 = sqrt((w^2)*eps0*mu0);% this is the same as 2*pi/lambda
Zw = sqrt(mu0/eps0);

a = 1; %waveguide radius (not important for plotting)
m_max = 1;%mode list max azumuthal order
n_max = 2;%mode list max radial order
mode_nr = 1; %mode number to plot from the list sorted according to cuttof wavelenghts  


p_TM = zeros(m_max*n_max, 5);
p_TE = zeros(m_max*n_max, 5);

n = 1;
for m = 0:m_max
    m_indexes = ones(n_max,1)*m;
    n_indexes = (1:n_max)';
    p_TM(n:n+n_max-1,1) = m_indexes;
    p_TM(n:n+n_max-1,2) = n_indexes;
    p_TM(n:n+n_max-1,4) = zerobess('J',m,n_max);

    p_TE(n:n+n_max-1,1) = m_indexes;
    p_TE(n:n+n_max-1,2) = n_indexes;
    p_TE(n:n+n_max-1,3) = ones(n_max,1);
    if m == 0
        dj_roots = zerobess('DJ',m,n_max+1);
        p_TE(n:n+n_max-1,4) = dj_roots(2:end);
    else
        p_TE(n:n+n_max-1,4) = zerobess('DJ',m,n_max);
    end
    
    n = n+n_max;
end

p_array = cat(1, p_TM, p_TE);
[~,idx] = sort(p_array(:,4)); %sortging accroding to cuttof wavelenght
p_array = p_array(idx,:);
lambda_cutoff = 2*pi*a./p_array(:, 4);
p_array(:,5) = (1e-9)*v./lambda_cutoff;%cutoff freq, GHz

a_cuttoff_lambdas = p_array(:, 4)/(2*pi);

fprintf('# mode Cuttoff_lambda root\n');
for i=1:length(p_array)
    fprintf('%d. ', i);
    
    if p_array(i, 3) == 1
        fprintf('TE')
    else
        fprintf('TM')
    end


    fprintf('%d%d %0.4f %0.4f\n', p_array(i, 1), p_array(i, 2), a_cuttoff_lambdas(i), p_array(i, 4));
end

n_points = 10; %plot resolution
[xy_array, r_phi_array] = get_circular_r_phi(a, n_points);

m = p_array(mode_nr,1);
pq = p_array(mode_nr,4);
TE_TM = p_array(mode_nr,3);

beta = sqrt(k0^2 - (pq./a).^2);
prop_c = 1j*beta;

if TE_TM == 1%TE
    Z = w*mu0./beta;
elseif TE_TM == 0%TM
    Z = beta./(w*eps0);
end
    
[ E, H ] = cwg_field_3( r_phi_array(:, 1), r_phi_array(:, 2), a, m, pq, TE_TM, -1);

E = sqrt(Z).*E;
H = sqrt(1/Z).*H;

A = get_c2p_matrix(r_phi_array(:,2));
Ecart = zeros(size(E));
Hcart = zeros(size(H));
for i = 1:length(E)
    Ecart(i,:) = E(i,:)*A(:,:,i);
    Hcart(i,:) = H(i,:)*A(:,:,i);
end

x = xy_array(:,1);
y = xy_array(:,2);
ampl = sqrt(abs(E(:, 1)).^2 + abs(E(:, 2)).^2);

[xi_p,yi_p] = meshgrid(linspace(-a,a,n_points), linspace(-a,a,n_points));
pattern = griddata(x,y,ampl,xi_p,yi_p);

figure(1)
% pcolor(xi_p, yi_p, pattern)
% colormap(jet)
% shading flat
% colorbar

hold on
quiver(x,y,Ecart(:,1),Ecart(:,2))
hold on
quiver(x,y,Hcart(:,1),Hcart(:,2)) %uncomment to plot H field vectors

if TE_TM == 1%TE
    plot_str = ['TE', num2str(p_array(mode_nr,1)), num2str(p_array(mode_nr,2))];
elseif TE_TM == 0%TM
    plot_str = ['TM', num2str(p_array(mode_nr,1)), num2str(p_array(mode_nr,2))];
end
title(plot_str);

saveas(gcf,'circ_wg_plot.png')
