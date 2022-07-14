function [ deriv_z ] = besselj_deriv( nu, z )
%BESSELJ_DERIV Summary of this function goes here
%   Detailed explanation goes here
if nu == 0
    deriv_z = -besselj(1,z);
else
    deriv_z = 0.5*(besselj(nu-1,z)-besselj(nu+1,z));
end

end

