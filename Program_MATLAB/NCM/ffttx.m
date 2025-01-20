function y = ffttx(x)
%FFTTX  Textbook Fast Finite Fourier Transform.
%    FFTTX(X) computes the same finite Fourier transform as FFT(X).
%    The code uses a recursive divide and conquer algorithm for
%    even order and matrix-vector multiplication for odd order.
%    If length(X) is m*p where m is odd and p is a power of 2, the
%    computational complexity of this approach is O(m^2)*O(p*log2(p)).

%   Copyright 2012 Cleve Moler and The MathWorks, Inc.
%   Revised by 张志涌（Zhiyong Zhang）,2022

[mm,nn]=size(x);
if min(mm,nn)==1
    x = x(:);
    n = length(x);
    omega = exp(-2*pi*1i/n);

    if rem(n,2) == 0
       % Recursive divide and conquer
       k = (0:n/2-1)';
       w = omega .^ k;
       u = ffttx(x(1:2:n-1));
       v = w.*ffttx(x(2:2:n));
       y = [u+v; u-v];
    else
       % The Fourier matrix.
       j = 0:n-1;
       k = j';
       F = omega .^ (k*j);
       y = F*x;
    end
else
    y=zeros(mm,nn);
    for kk=1:nn
        xc=x(:,kk);
        y(:,kk)=ffttx(xc);
    end
end
    