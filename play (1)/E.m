function out=E(sig,Sig)

Minv=5.55e-4;
sqrtwoopi=0.797884560802865;

a=2.35; b=2.4; c=3.2;
S=sqrt(1/(1/sig^2+1/Sig^2)); alpha=a/c*(2*pi)^3*sig^3*Sig^3;

%out=0.375*(1/sig^2+Minv/Sig^2)-sqrtwoopi/sqrt(sig^2+Sig^2)-4*pi/c*S^3*Li32(alpha);
%out=0.5*(0.2*3^2)*Sig^2+0.375*(1/sig^2+Minv/Sig^2)-sqrtwoopi/sqrt(sig^2+Sig^2)-4*pi/c*S^3*Li32(alpha);
%out=0.5*(0.2*3^2)*Sig^2+0.375*Minv/Sig^2;

%neo-hf
%out=0.375*(1/sig^2+Minv/Sig^2)-sqrtwoopi/sqrt(sig^2+Sig^2);

%neo-dft
%out=0.375*(1/sig^2+Minv/Sig^2)-sqrtwoopi/sqrt(sig^2+Sig^2)-Ecep(sig,Sig);

%neo-sho1
%out=0.5*(0.2)*Sig^2+0.375*(1/sig^2+Minv/Sig^2)-sqrtwoopi/sqrt(sig^2+Sig^2)-Ecep(sig,Sig);

%neo-sho10
out=0.5*(0.2*10)*Sig^2+0.375*(1/sig^2+Minv/Sig^2)-sqrtwoopi/sqrt(sig^2+Sig^2)-Ecep(sig,Sig);

%out=0.5*(0.2*10)*Sig^2+0.375*(Minv/Sig^2);





return