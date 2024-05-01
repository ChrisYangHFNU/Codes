function [A,b] = qrsteps(A,b)
%QRSTEPS  Orthogonal-triangular decomposition.
%   Demonstrates M-file version of built-in QR function.
%   R = QRSTEPS(A) is the upper trapezoidal matrix R that
%   results from the orthogonal transformation, R = Q'*A.
%   With no output argument, QRSTEPS(A) shows the steps in
%   the computation of R.  Press <enter> key after each step.
%   [R,bout]  = QRSTEPS(A,b) also applies the transformation to b.

%   Copyright 2012 Cleve Moler and The MathWorks, Inc.
%   Revised by 张志涌（Zhiyong Zhang）,2022

[m,n] = size(A);
if nargin < 2, b = zeros(m,0); end

% Householder reduction to triangular form.

if nargout == 0
   clc
   A
   if ~isempty(b), b, end
   disp('以上是输入的设计矩阵和观察向量。(在命令窗中，按空格键可继续!)')
   pause
end

for k = 1:min(m-1,n)
   % Introduce zeros below the diagonal in the k-th column.
   % Use Householder transformation, I - rho*u*u'.
   i = k:m;
   u = A(i,k);
   sigma = norm(u);

   % Skip transformation if column is already zero.
   if sigma ~= 0

      if u(1) ~= 0, sigma = sign(u(1))*sigma; end
      u(1) = u(1) + sigma;
      rho = 1/(conj(sigma)*u(1));

      % Update the k-th column.
      A(i,k) = 0;
      A(k,k) = -sigma;

      % Apply the transformation to remaining columns of A.
      j = k+1:n;
      v = rho*(u'*A(i,j));
      A(i,j) = A(i,j) - u*v;

      % Apply the transformation to b.
      if ~isempty(b)
         tau = rho*(u'*b(i));
         b(i) = b(i) - tau*u;
      end
   end
   if nargout == 0
      clc
      disp(['以下是第 ',int2str(k),' 次反射变换结果。(在命令窗中，按空格键可继续！)'])
      A
      if ~isempty(b), b, end
      pause
   end
end
