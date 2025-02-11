function u = fmintx(F,a,b,tol,varargin)
%FMINTX  Textbook version of FMINBND
%   x = FMINTX(F,a,b) finds a local minimizer x of the function F
%   in the interval a <= x <= b. F accepts scalar input x and returns 
%   a scalar function value, F(x).
%
%   x = FMINTX(F,a,b,tol) uses stopping tolerance tol instead of 1.e-6.
%
%   x = FMINTX(F,a,b,tol,p1,p2,...) provides for additional
%   arguments, which are passed to the objective function, F(x,p1,p2,...).
%   (Use [] as a place holder for tol to get the default tolerance.)
%
%   Example
%      x = fmintx(@cos,3,4)
%      computes pi to six decimal places.
%
%   See also FMINBND, FMINSEARCH, FZERO, @, INLINE.

%   Copyright 2012 Cleve Moler and The MathWorks, Inc.
%   Reference: "Computer Methods for Mathematical Computations",
%   Forsythe, Malcolm, and Moler, Prentice-Hall, 1976.
%   Revised by 张志涌（Zhiyong Zhang）,2022 

% Initialization

if nargin < 4 || isempty(tol)
   tol = 1.e-6; 
end
clf
fplot(F,[a,b],'b'),xlabel('x'),ylabel('f(x)'),shg
phi = (1 + sqrt(5))/2;
rho = 2 - phi;
u = a + rho*(b-a);
v = u; w = u; x = u;
fu = F(u,varargin{:}); 
fv = fu; fw = fu; fx = fu;
xm = 0.5*(a+b);
d = 0.0;
e = 0.0;
sstr='init';
fprintf('%s\n',[' step',blanks(14),'x',blanks(14),'f(x)'])
fprintf('%s	%.10f	%.10f\n',['0  ',sstr,':'],u,fu)
line(u,fu,'Marker','.','MarkerSize',10,'Color','r')
text(u,fu,'0 '),shg

% Main loop
k=0;
while abs(x-xm) > tol
   k=k+1;
   % Is a parabolic fit possible?
   para = abs(e) > tol;
   if para
      % Try parabolic fit.
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      s = 2*(q-r);
      if s > 0.0, p = -p; end
      s = abs(s);
      % Is the parabola acceptable?
      para = ( (abs(p)<abs(0.5*s*e)) & (p>s*(a-x)) & (p<s*(b-x)) );
      if para
         % Parabolic interpolation step
         e = d;
         d = p/s;
      end
      sstr='para';
   end
   if ~para
      % Golden-section step
      if x >= xm
         e = a-x;
      else
         e = b-x;
      end
      d = rho*e;
      sstr='gold';
   end
   
   u = x + d;
   fu = F(u,varargin{:}); 
   kstr=int2str(k);
   fprintf('%s	%.10f	%.10f\n',[kstr,'  ',sstr,':'],u,fu)
   line(u,fu,'Marker','.','MarkerSize',10,'Color','r')
   text(u,fu,kstr),shg
 
   % Update a, b, x, v, w, xm
   if fu <= fx
      if u >= x, a = x; else, b = x; end
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
   else
      if u < x, a = u; else, b = u; end
      if ( (fu <= fw) || (w == x) )
         v = w; fv = fw;
         w = u; fw = fu;
      elseif ( (fu <= fv) || (v == x) || (v == w) )
         v = u; fv = fu;
      end
   end
   xm = 0.5*(a+b);
end

