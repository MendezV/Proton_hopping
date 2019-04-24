function out=Li32(a)

N=10000;


U=sqrt(2*log(2^1023))-1;
u=[0:U/N:U];

out=sum(u.^2./(a*exp(0.5*u.^2)+1),2)*U/N;

