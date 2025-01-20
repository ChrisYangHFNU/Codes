clear all

M = rand(2000,2000);
tic
[A1,B1] = eig(M);
t1 = toc

tic
M = gpuArray(M);
[A2,B2] = eig(M);
A2 = gather(A2);
t2 = toc