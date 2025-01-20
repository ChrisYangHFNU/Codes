function x = pagerank(U,G,p)
% PAGERANK  Google's PageRank
% pagerank(U,G,p) uses the URLs and adjacency matrix produced by SURFER,
% together with a damping factory p, (default is .85), to compute and plot
% a bar graph of page rank, and print the dominant URLs in page rank order.
% x = pagerank(U,G,p) returns the page ranks instead of printing.
% See also SURFER, SPY.

%   Revised by 张志涌（Zhiyong Zhang）,2022

if nargin < 3, p = .85; end

% Eliminate any self-referential links

G = G - diag(diag(G));
  
% c = out-degree, r = in-degree

[~,n] = size(G);
c = full(sum(G,1));     %revised
r = full(sum(G,2));     %revised

% Scale column sums to be 1 (or 0 where there are no out links).

k = find(c~=0);
D = sparse(k,k,1./c(k),n,n);

% Solve (I - p*G*D)*x = e

e = ones(n,1);
I = speye(n,n);
x = (I - p*G*D)\e;

% Normalize so that sum(x) == 1.

x = x/sum(x);

% Bar graph of page rank.

shg
bar(x,'b')
axis([-10,510,0,0.02])
title('Page Rank')

% Print URLs in page rank order.

if nargout < 1
   [~,q] = sort(-x);
   disp('     page-rank  in  out  url')
   k = 1;
   while (k <= n) && (x(q(k)) >= .005)
      j = q(k);
      fprintf(' %3.0f %8.4f %4.0f %4.0f  %s\n', ...
         j,x(j),r(j),c(j),U{j})         %revised
      k = k+1;
   end
end
