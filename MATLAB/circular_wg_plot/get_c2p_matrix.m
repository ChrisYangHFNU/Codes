function [ A ] = get_c2p_matrix( phi )
%GET_C2P_MATRIX Summary of this function goes here
%   Detailed explanation goes here
A = zeros(3,3,length(phi));

A(1,1,:) = cos(phi);
A(1,2,:) = sin(phi);
A(2,1,:) = -sin(phi);
A(2,2,:) = cos(phi);
A(3,3,:) = ones(length(phi),1);

end

