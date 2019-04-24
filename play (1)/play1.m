more off
format long

sig=1;
Sig=1;

step=0.1;

while(1)
E(sig,Sig)
[ds,dS]=fdgradE(sig,Sig);
sig-=step*ds;
Sig-=step*dS;
end

