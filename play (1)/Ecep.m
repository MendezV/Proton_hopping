function out=Ecep(sig,Sig)

a=2.35; b=2.4; c=3.2;

pre=1/(2*pi*sig*Sig)^3;
al=0.5*(1/sig^2+1/Sig^2);

N=10000;
U=sqrt(2*log(2^1023))-1;
u=[0:U/N:U];

%out=4*pi*pre*sum( (u.^2.*exp(-al*u.^2))./(a-b*sqrt(pre*exp(-al*u.^2))+c*pre*exp(-al*u.^2)) )*U/N;

out=4*pi*pre*al^(-1.5)*sum( (u.^2.*exp(-u.^2))./(a-b*sqrt(pre*exp(-u.^2))+c*pre*exp(-u.^2)) )*U/N;

return
