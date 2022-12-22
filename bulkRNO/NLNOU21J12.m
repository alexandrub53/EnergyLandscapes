%Data and Fits for NNO and Nd_0.75La_0.25NiO3, U=2.1eV, J=1.2eV, as used in
%the solid solution work on NLNO

clear all
bulk=[
0.22	1.209106
0.21	1.194028
0.2	      1.176952
0.19	1.162163
0.18	1.142579
0.17	1.116733
0.16	1.094949
0.15	1.068016
0.14	1.037892
0.13	1.004513
0.12	0.942834
0.11	0.840647
0.1	0.650089
0.09	0.543565
0.08	0.462801
0.07	0.391415
0.06	0.328137
0.05	0.26928
0.04	0.213294
0.03	0.157391
0.02	0.100942
0.01	0.053838
0	0.000481]
bulk=flip(bulk)





bulk25=[0	0.00000
0.01	0.047937
0.02	0.09636
0.03	0.147635
0.04	0.197679
0.05	0.249989
0.06	0.303743
0.07	0.36014
0.08	0.419442
0.09	0.486355
0.1	0.560207
0.11	0.657668
0.12	0.803783
0.13	0.941051
0.14	1.000352
0.15	1.03668
0.16	1.067247
0.17	1.093627
0.18	1.116931
0.19	1.135956
0.2	1.15656
0.21	1.175226
0.22	1.188325
0.23	1.203312
0.24	1.218717
0.25	1.230643]

bulk(:,2)=bulk(:,2);
%Values From DFT+U
%chi0b=1.164
%chi0b25=1.168
%Values from DFT
chi0b=1.159592851
chi0b25=1.08061

U=2.1
J=1.2
DC=U-5/3*J
gb0=2.8865
gb025=2.896
gb=gb0*(1+chi0b*DC)
gb25=gb025*(1+chi0b25*DC)
kb=2*13.4303+1/2*gb0*chi0b*gb0
kb25=2*13.75+1/2*gb025*chi0b25*gb025
%kb25=kb
Q=bulk(:,1)/gb
%Q2=bulkup(:,1)/gb
N=bulk(:,2)
%N2=bulkup(:,2)
%plot(Q,N,'-o',Q2,N2,'-o')
Q25=bulk25(:,1)/gb25
N25=bulk25(:,2)


x=[0:0.001:0.07]
by=2*kb/gb*x
ry=2*kb25/gb25*x

plot(Q,N,'-ob',x,by,'-b','MarkerSize',10,'LineWidth',5)
hold on

plot(Q25,N25,'-or',x,ry,'-r','MarkerSize',10,'LineWidth',5)

set( gca , 'FontSize'   , 37 );
legend('NNO','NNOLine','25%','25% Line')


%plot(Q,N,'o',Q25,N25,'o')

[ Ngrid,Ecorb,Etot,A,B,C ] = getEnergyGrid( N,Q,gb,kb )
[ Ngrid,Ecorb25,Etot25,A25,B25,C25 ] = getEnergyGrid( N25,Q25,gb25,kb25 )
%plot(2/gb*(A*Ngrid.^5+B*Ngrid.^3+C*Ngrid.^1),Ngrid,'-b')
%plot(2/gb25*(A25*Ngrid.^5+B25*Ngrid.^3+C25*Ngrid.^1),Ngrid,'-r')
%plot(Ngrid,Etot,'-ob',Ngrid,Etot25,'-or')

[  Ngrid,Ecor2,Etot2,A2,B2,C2,Constant ] = ConstrainedEnergyGrid(  N,Q,gb,kb, 4 )
[  Ngrid,Ecor225,Etot225,A225,B225,C225,Constant ] = ConstrainedEnergyGrid(  N25,Q25,gb25,kb25, 4 )
%plot(2/gb*(A*Ngrid.^5+B*Ngrid.^3+C*Ngrid.^1),Ngrid,'--b','LineWidth',5)
%plot(2/gb25*(A25*Ngrid.^5+B25*Ngrid.^3+C25*Ngrid.^1),Ngrid,'--r','LineWidth',5)

set( gca , 'FontSize'   , 37 );
hold on
title('U=2.1 J=1.2')
%[  Ngrid,Ecor2,Etot2,A2,B2,C2,Constant ] = ConstrainedEnergyGrid(  N,Q,gb,kb, 4 )





%legend('NdNiO3','25% La','NNO line','25% Line')
ylabel('N')
xlabel('Q(Angstrom)')
grid on
figure
plot(Q,2/gb*(A*N.^5+B*N.^3+C*N),'-bo',Q,Q,'LineWidth',5)
grid on
hold on
plot(Q25,2/gb25*(A25*N25.^5+B25*N25.^3+C25*N25),'-ro','LineWidth',5)
legend('NdNiO3','Line','25 % La')

set( gca , 'FontSize'   , 37 );
figure
plot(Ngrid,Etot,'-bo','LineWidth',5)
grid on
hold on
plot(Ngrid,Etot25,'-ro','LineWidth',5)

[kb,gb,A,B,C]
[kb25,gb25,A25,B25,C25]
