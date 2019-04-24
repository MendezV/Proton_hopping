function out=Li32(a)

N=1000;


Na=length(a); a=reshape(a,Na,1)*ones(1,N+1);

U=sqrt(2*log(2^1023))-1;
u=ones(Na,1)*[0:U/N:U];

out=sum(u.^2./(a.*exp(0.5*u.^2)+1),2)*U/N;

