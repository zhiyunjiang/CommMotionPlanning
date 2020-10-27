function out=polyfit2D(d,x)
M=length(x);
A=[ones(M,1) -10*log10(d)];
out=lsqr(A,x(:),1e-6);







