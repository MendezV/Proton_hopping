sig=1
Sig=0.01

ec=Ecep(sig,Sig)

a=2.35; b=2.4; c=3.2;
S=sqrt(1/(1/sig^2+1/Sig^2)); alpha=a/c*(2*pi)^3*sig^3*Sig^3;
4*pi/c*S^3*Li32(alpha)

4*pi/c*S^3*Li32(alpha)/ec


