function [N] = get_N_te(m, q)
    N = sqrt(em(m)/pi)./(sqrt(q.^2 - m^2).*besselj(m, q));
end