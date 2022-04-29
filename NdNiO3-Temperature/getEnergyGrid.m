function [ Ngrid,Ecor,Etot,A,B,C ] = getEnergyGrid( N,Q,g,k )
%Electronic energy fit based on Delta N, Q, g for Nickelate system 

%symmetrize Q and Delta N to make sure no residual odd terms appear in the
%polynomial
Q=[-flip(Q) ; Q]
N=[-flip(N) ; N]
pE=polyfit(N,g/2*Q,6)
A=pE(2);
B=pE(4);
C=pE(6);

%Generate Electronic Energy E(Delta N)
Ngrid=[-1.8:0.001:1.8];
Ecor=(A/6)*Ngrid.^6+(B/4)*Ngrid.^4+(C/2)*Ngrid.^2;
%Generate total energy, optimized over Q. First term is obtained via direct
%substitution.
Etot=-1/8*g^2/k*Ngrid.^2+Ecor;

end

