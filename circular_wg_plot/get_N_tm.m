function [N] = get_N_tm(m, p)
    N = sqrt(em(m)/pi)./(p.*besselj(m+1, p));
end