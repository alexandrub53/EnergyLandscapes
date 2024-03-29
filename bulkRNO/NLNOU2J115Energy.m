%Data and Fits for NNO and Nd_0.75La_0.25NiO3, U=2.1eV, J=1.15eV, as used in
%the solid solution work on NLNO

clear all
bulk=[
0	0
0.01	0.043767
0.02	0.090553
0.03	0.132531
0.04	0.178382
0.05	0.22091
0.06	0.271343
0.07	0.317688
0.08	0.368391
0.09	0.422357
0.1	0.479379
0.11	0.541599
0.12	0.617438
0.13	0.709913
0.14	0.82761
0.15	0.914203
0.16	0.969272
0.17	1.006483
0.18	1.033624
0.19	1.059321
0.2	1.079804
0.21	1.101465
0.22	1.12087
0.23	1.134559
0.24	1.154913
0.25	1.170698
0.26	1.183213
0.27	1.196465
0.28	1.207281
0.29	1.224756
0.3	1.23756
0.31	1.251913
0.32	1.261235
0.33	1.271748
0.34	1.283685
0.35	1.295006
0.36	1.30539
0.37	1.314554
0.38	1.324145
0.39	1.336595
0.4	1.342392
0.41	1.352256
0.42	1.357852
0.43	1.368671
0.44	1.376227]





bulk25=[0	0
0.01	0.039846
0.02	0.082272
0.03	0.122777
0.04	0.164574
0.05	0.207349
0.06	0.250562
0.07	0.29607
0.08	0.339906
0.09	0.38815
0.1	0.440795
0.11	0.49192
0.12	0.556833
0.13	0.624229
0.14	0.709888
0.15	0.815234
0.16	0.905492
0.17	0.961482
0.18	0.99918
0.19	1.031263
0.2	1.054508
0.21	1.073619
0.22	1.096508
0.23	1.114797
0.24	1.132288
0.25	1.15244
0.26	1.163253
0.27	1.179902
0.28	1.193986
0.29	1.206602
0.3	1.218626
0.31	1.231211
0.32	1.244952
0.33	1.254542
0.34	1.267047
0.35	1.276829
0.36	1.289142
0.37	1.301712
0.38	1.305257
0.39	1.320473
0.4	1.331336
0.41	1.338217
0.42	1.345795
0.43	1.355344
0.44	1.363176]

bulk(:,2)=bulk(:,2);
%Values From DFT+U
%chi0b=1.164
%chi0b25=1.168
%Values from DFT
chi0b=1.09681
chi0b25=1.08061

U=2.1
J=1.15
DC=U-5/3*J
gb0=2.8865
gb025=2.896
gb=gb0*(1+chi0b*DC)
gb25=gb025*(1+chi0b25*DC)
kb=2*13.4303+1/2*gb0*chi0b*gb0
kb25=2*13.75+1/2*gb025*chi0b25*gb025
Q=bulk(:,1)/gb
N=bulk(:,2)
Q25=bulk25(:,1)/gb25
N25=bulk25(:,2)
x=[0:0.001:0.07]
by=2*kb/gb*x
ry=2*kb25/gb25*x
plot(Q,N,'-ob',x,by,'-b')
hold on
plot(Q25,N25,'-or',x,ry,'-r')
set( gca , 'FontSize'   , 37 );
legend('NNO','NNOLine','25%','25%')
figure



[ Ngrid,Ecorb,Etot,A,B,C ] = getEnergyGrid( N,Q,gb,kb )
[ Ngrid,Ecorb25,Etot25,A25,B25,C25 ] = getEnergyGrid( N25,Q25,gb25,kb25 )

plot(Ngrid,Etot,'-ob',Ngrid,Etot25,'-or')

set( gca , 'FontSize'   , 37 );
hold on
[  Ngrid,Ecor2,Etot2,A2,B2,C2,Constant ] = ConstrainedEnergyGrid(  N,Q,gb,kb, 4 )
[  Ngrid,Ecor225,Etot225,A225,B225,C225,Constant ] = ConstrainedEnergyGrid(  N25,Q25,gb25,kb25, 4 )
%plot(Ngrid,Etot2,'--b',Ngrid,Etot225,'--r')

legend('NdNiO3','25% La')
ylabel('E(N) (eV)')
xlabel('N')
grid on
figure
plot(Q,2/gb*(A*N.^5+B*N.^3+C*N),'-bo',Q,Q,'LineWidth',5)
grid on
hold on
plot(Q25,2/gb25*(A25*N25.^5+B25*N25.^3+C25*N25),'-ro','LineWidth',5)
legend('NdNiO3','Line','25 % La')

set( gca , 'FontSize'   , 37 );
[kb gb A B C]
