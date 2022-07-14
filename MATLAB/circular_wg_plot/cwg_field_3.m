%from waveguide handbook
function [ E, H] = cwg_field_3(r, phi, a, m, pq, TE_TM, degen)
%CWG_FIELD Summary of this function goes here
%   Detailed explanation goes here
%returns only transverse components,z set to 0

% z_unit = cat(3, zeros(size(r)), zeros(size(r)), ones(size(r)));
zeros_array = zeros(size(r));
k = pq./a; %wavenumber

if TE_TM == 1 %TE (H modes), q
    
    if degen == -1
        grad_psi_r = k.*besselj_deriv(m, k.*r).*cos(m.*phi);
        grad_psi_phi = (1./r).*besselj(m, k.*r).*(-m.*sin(m.*phi));
    elseif degen == 1
        grad_psi_r = k.*besselj_deriv(m, k.*r).*sin(m.*phi);
        grad_psi_phi = (1./r).*besselj(m, k.*r).*m.*cos(m.*phi);
    else
        error('degen must be -1 or +1')
    end

    N_te = get_N_te(m, pq);
    grad_psi_r = N_te.*grad_psi_r;
    grad_psi_phi = N_te.*grad_psi_phi;    

%     E = cross(z_unit, cat(3, grad_psi_r, grad_psi_phi, zeros(size(r))), 3);
%     H = cross(z_unit, E, 3);    
    E = cat(3, -grad_psi_phi, grad_psi_r, zeros_array);
    H = cat(3, -grad_psi_r, -grad_psi_phi, zeros_array);
    
elseif TE_TM == 0 %TM (E modes), p
    
    if degen == -1
        grad_psi_r = k.*besselj_deriv(m, k.*r).*sin(m.*phi);
        grad_psi_phi = (1./r).*besselj(m, k.*r).*m.*cos(m.*phi);
    elseif degen == 1
        grad_psi_r = k.*besselj_deriv(m, k.*r).*cos(m.*phi);
        grad_psi_phi = (1./r).*besselj(m, k.*r).*(-m.*sin(m.*phi));
    else
        error('degen must be -1 or +1')
    end
    
    N_tm = get_N_tm(m, pq);
    grad_psi_r = N_tm.*grad_psi_r;
    grad_psi_phi = N_tm.*grad_psi_phi;    

%     E = -cat(3, grad_psi_r, grad_psi_phi, zeros(size(r)));
%     H = cross(z_unit, E, 3);
    E = cat(3, -grad_psi_r, -grad_psi_phi, zeros_array);
    H = cat(3, grad_psi_phi, -grad_psi_r, zeros_array);
    
end
end