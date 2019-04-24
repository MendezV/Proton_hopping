function [a,b]=fdgradE(sig,Sig)

del=1e-5;

a=(E(sig*(1+del),Sig)-E(sig*(1-del),Sig))/(2*del*sig);
b=(E(sig,Sig*(1+del))-E(sig,Sig*(1-del)))/(2*del*Sig);

return
