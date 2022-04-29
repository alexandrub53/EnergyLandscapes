%Main function to Generate Temperature-Dependent Plots for NdNiO3, shown in
% https://arxiv.org/abs/2105.02271 ;


figure('Position',[0 0 1080 1080])
clear all

%Load the Delta N vs applied on-site field from another file
LoadBeta

% Interaction parameters used in https://www.pnas.org/doi/abs/10.1073/pnas.1818728116   
U=2.1
J=1.2
%double counting term
DC=U-5/3*J
x=[0:0.01:0.1]

%Electronic linear response from DFT.
chi0b=1.159592851
% g0 obtained directly from DFT
gb0=2.8865
% g corrected for Double Counting. See https://www.pnas.org/doi/abs/10.1073/pnas.1818728116
gb=gb0*(1+chi0b*DC)
% lattice stiffess correctd for double counting. See https://www.pnas.org/doi/abs/10.1073/pnas.1818728116
kbulk=2*13.4303+1/2*gb0*chi0b*gb0
by=2*kbulk/gb*x

%Transform Delta_S into its corresponding lattice distortion Q.

Nb=bulk(: , 2);
Qb=bulk(: , 1)/gb;
Nb30=bulk30(:,2);
Qb30=bulk30(:,1)/gb;
Nb35=bulk35(:,2);
Qb35=bulk35(:,1)/gb;
Nb20=bulk202(:,2);
Qb20=bulk202(:,1)/gb;
Nb12=bulk12(:,2);
Qb12=bulk12(:,1)/gb;
Nb10=bulk10(:,2);
Qb10=bulk10(:,1)/gb;
Nb1025=bulk1025(:,2);
Qb1025=bulk1025(:,1)/gb;
Nb105=bulk105(:,2);
Qb105=bulk105(:,1)/gb;
Nb1075=bulk1075(:,2);
Qb1075=bulk1075(:,1)/gb;
Nb115=bulk115(:,2);
Qb115=bulk115(:,1)/gb;
Nb15=bulk15(:,2);
Qb15=bulk15(:,1)/gb;
Nb25=bulk25(:,2);
Qb25=bulk25(:,1)/gb;
Nb11=bulk11(:,2);
Qb11=bulk11(:,1)/gb;

%Fits
[ Ngrid,Ecorb,Etotb,Ab,Bb,Cb ] = getEnergyGrid( Nb,Qb,gb,kbulk )
[ Ngrid,Ecorb30,Etotb30,Ab30,Bb30,Cb30 ] = getEnergyGrid( Nb30,Qb30,gb,kbulk )
[ Ngrid,Ecorb20,Etotb20,Ab20,Bb20,Cb20 ] = getEnergyGrid( Nb20,Qb20,gb,kbulk )
[ Ngrid,Ecorb12,Etotb12,Ab12,Bb12,Cb12 ] = getEnergyGrid( Nb12,Qb12,gb,kbulk )
[ Ngrid,Ecorb10,Etotb10,Ab10,Bb10,Cb10 ] = getEnergyGrid( Nb10,Qb10,gb,kbulk )
[ Ngrid,Ecorb11,Etotb11,Ab11,Bb11,Cb11 ] = getEnergyGrid( Nb11,Qb11,gb,kbulk )
[ Ngrid,Ecorb15,Etotb15,Ab15,Bb15,Cb15 ] = getEnergyGrid( Nb15,Qb15,gb,kbulk )
[ Ngrid,Ecorb25,Etotb25,Ab25,Bb25,Cb25 ] = getEnergyGrid( Nb25,Qb25,gb,kbulk )
[ Ngrid,Ecorb35,Etotb35,Ab35,Bb35,Cb35 ] = getEnergyGrid( Nb35,Qb35,gb,kbulk )
[ Ngrid,Ecorb1025,Etotb1025,Ab1025,Bb1025,Cb1025 ] = getEnergyGrid( Nb1025,Qb1025,gb,kbulk )
[ Ngrid,Ecorb105,Etotb105,Ab105,Bb105,Cb105 ] = getEnergyGrid( Nb105,Qb105,gb,kbulk )
[ Ngrid,Ecorb1075,Etotb1075,Ab1075,Bb1075,Cb1075 ] = getEnergyGrid( Nb1075,Qb1075,gb,kbulk )
[ Ngrid,Ecorb115,Etotb115,Ab115,Bb115,Cb115 ] = getEnergyGrid( Nb115,Qb115,gb,kbulk )
set(gca,'FontSize',30)
grid on


Qgrid=[0:0.001:0.14]
Ngridb=2*kbulk/gb*Qgrid


%Time to plot!

plot(Qb,Nb,'+','Color',[0 1.0 1.0],'MarkerSize',15,'LineWidth',3)
hold on
plot(Qb35,Nb35,'d','Color',[0.05 0.8 0.8],'MarkerSize',15,'LineWidth',3)
plot(Qb30,Nb30,'*','Color',[0.1  0.6 0.6],'MarkerSize',15,'LineWidth',3)
plot(Qb25,Nb25,'o','Color',[0.15 0.5 0.5],'MarkerSize',15,'LineWidth',3)
plot(Qb20,Nb20,'+','Color',[0.2 0.4 0.4],'MarkerSize',15,'LineWidth',3)
plot(Qb15,Nb15,'d','Color',[0.3 0.3 0.3],'MarkerSize',15,'LineWidth',3)
plot(Qb12,Nb12,'*','Color',[0.4 0.2 0.2],'MarkerSize',15,'LineWidth',3)
plot(Qb115,Nb115,'o','Color',[0.5 0.1 0.1],'MarkerSize',15,'LineWidth',3)
plot(Qb11,Nb11,'+','Color',[0.6 0.1 0.1],'MarkerSize',15,'LineWidth',3)
plot(Qb1075,Nb1075,'d','Color',[0.7 0.0 0.0],'MarkerSize',15,'LineWidth',3)
plot(Qb105,Nb105,'*','Color',[0.8 0.0 0.0],'MarkerSize',15,'LineWidth',3)
plot(Qb1025,Nb1025,'o','Color',[0.9 0.0 0.0],'MarkerSize',15,'LineWidth',3)
plot(Qb10,Nb10,'+','Color',[1.0 0.0 0.0],'MarkerSize',15,'LineWidth',3)


plot((Ab*Ngrid.^5+(Bb)*Ngrid.^3+(Cb)*Ngrid)*2/gb,Ngrid,'Color',[0 1.0 1.0],'LineWidth',2)
plot((Ab35*Ngrid.^5+(Bb35)*Ngrid.^3+(Cb35)*Ngrid)*2/gb,Ngrid,'Color',[0.05 0.8 0.8],'LineWidth',2)
plot((Ab30*Ngrid.^5+(Bb30)*Ngrid.^3+(Cb30)*Ngrid)*2/gb,Ngrid,'Color',[0.1  0.6 0.6],'LineWidth',2)
plot((Ab25*Ngrid.^5+(Bb25)*Ngrid.^3+(Cb25)*Ngrid)*2/gb,Ngrid,'Color',[0.15 0.5 0.5],'LineWidth',2)
plot((Ab20*Ngrid.^5+(Bb20)*Ngrid.^3+(Cb20)*Ngrid)*2/gb,Ngrid,'Color',[0.2 0.4 0.4],'LineWidth',2)
plot((Ab15*Ngrid.^5+(Bb15)*Ngrid.^3+(Cb15)*Ngrid)*2/gb,Ngrid,'Color',[0.3 0.3 0.3],'LineWidth',2)
plot((Ab12*Ngrid.^5+(Bb12)*Ngrid.^3+(Cb12)*Ngrid)*2/gb,Ngrid,'Color',[0.4 0.2 0.2],'LineWidth',2)
plot((Ab115*Ngrid.^5+(Bb115)*Ngrid.^3+(Cb115)*Ngrid)*2/gb,Ngrid,'Color',[0.5 0.1 0.1],'LineWidth',2)
plot((Ab11*Ngrid.^5+(Bb11)*Ngrid.^3+(Cb11)*Ngrid)*2/gb,Ngrid,'Color',[0.6 0.1 0.1],'LineWidth',2)
plot((Ab1075*Ngrid.^5+(Bb1075)*Ngrid.^3+(Cb1075)*Ngrid)*2/gb,Ngrid,'Color',[0.7 0.0 0.0],'LineWidth',2)
plot((Ab105*Ngrid.^5+(Bb105)*Ngrid.^3+(Cb105)*Ngrid)*2/gb,Ngrid,'Color',[0.8 0.0 0.0],'LineWidth',2)
plot((Ab1025*Ngrid.^5+(Bb1025)*Ngrid.^3+(Cb1025)*Ngrid)*2/gb,Ngrid,'Color',[0.9 0.0 0.0],'LineWidth',2)
plot((Ab10*Ngrid.^5+(Bb10)*Ngrid.^3+(Cb10)*Ngrid)*2/gb,Ngrid,'Color',[1.0 0.0 0.0],'LineWidth',2)
xlabel('Q(Angstrom)')
ylabel('\Delta N')
legend('T=290K','T=330K','T=386K','T=464K','T=580K','T=773K','T=966K','T=1000K','T=1050K','T=1080K','T=1105K','T=1131K','T=1160K')
axis([0 0.08 0 1.5])




figure 
plot(Ngrid,Etotb,'Color',[0 1.0 1.0],'LineWidth',5)
hold on
plot(Ngrid,Etotb35,'Color',[0.05 0.8 0.8],'LineWidth',5)
plot(Ngrid,Etotb30,'Color',[0.1  0.6 0.6],'LineWidth',5)
plot(Ngrid,Etotb25,'Color',[0.15 0.5 0.5],'LineWidth',5)
plot(Ngrid,Etotb20,'Color',[0.2 0.4 0.4],'LineWidth',5)
plot(Ngrid,Etotb15,'Color',[0.3 0.3 0.3],'LineWidth',5)
plot(Ngrid,Etotb12,'Color',[0.4 0.2 0.2],'LineWidth',5)
plot(Ngrid,Etotb115,'Color',[0.5 0.1 0.1],'LineWidth',5)
plot(Ngrid,Etotb11,'Color',[0.6 0.1 0.1],'LineWidth',5)
plot(Ngrid,Etotb1075,'Color',[0.7 0.0 0.0],'LineWidth',5)
plot(Ngrid,Etotb105,'Color',[0.8 0.0 0.0],'LineWidth',5)
plot(Ngrid,Etotb1025,'Color',[0.9 0.0 0.0],'LineWidth',5)
plot(Ngrid,Etotb10,'Color',[1.0 0.0 0.0],'LineWidth',5)


xlabel('Q(Angstrom)')
ylabel('\Delta N (electrons)')


figure
plot(Ngrid,Etotb,'Color',[0 1.0 1.0],'LineWidth',5)
hold on
plot(Ngrid,Etotb35,'Color',[0.05 0.8 0.8],'LineWidth',5)
plot(Ngrid,Etotb30,'Color',[0.1  0.6 0.6],'LineWidth',5)
plot(Ngrid,Etotb25,'Color',[0.15 0.5 0.5],'LineWidth',5)
plot(Ngrid,Etotb20,'Color',[0.2 0.4 0.4],'LineWidth',5)
plot(Ngrid,Etotb15,'Color',[0.3 0.3 0.3],'LineWidth',5)
plot(Ngrid,Etotb12,'Color',[0.4 0.2 0.2],'LineWidth',5)
plot(Ngrid,Etotb115,'Color',[0.5 0.1 0.1],'LineWidth',5)
plot(Ngrid,Etotb11,'Color',[0.6 0.1 0.1],'LineWidth',5)
plot(Ngrid,Etotb1075,'Color',[0.7 0.0 0.0],'LineWidth',5)
plot(Ngrid,Etotb105,'Color',[0.8 0.0 0.0],'LineWidth',5)
plot(Ngrid,Etotb1025,'Color',[0.9 0.0 0.0],'LineWidth',5)
plot(Ngrid,Etotb10,'Color',[1.0 0.0 0.0],'LineWidth',5)

legend('290K','330K','386K','464K','580K','773K','966K','1000K','1050K','1080K','1105K','1131K','1160K','290K Fit','330K Fit','386K Fit','464K Fit','580K Fit','773K Fit','966K Fit','1000K Fit','1050K Fit','1080K Fit','1105K Fit','1131K Fit','1160K Fit')


grid on
set(gca,'FontSize',50)
xlabel('\Delta N')
ylabel('E_{tot}(eV)')
[V,I]=min(Etotb)
Emin40=Etotb(I)
Nmin40=abs(Ngrid(I))

[V,I]=min(Etotb35)
Emin35=Etotb35(I)
Nmin35=abs(Ngrid(I))

[V,I]=min(Etotb30)
Emin30=Etotb30(I)
Nmin30=abs(Ngrid(I))

[V,I]=min(Etotb25)
Emin25=Etotb25(I)
Nmin25=abs(Ngrid(I))

[V,I]=min(Etotb20)
Emin20=Etotb20(I)
Nmin20=abs(Ngrid(I))

[V,I]=min(Etotb15)
Emin15=Etotb15(I)
Nmin15=abs(Ngrid(I))


[V,I]=min(Etotb12)
Emin12=Etotb12(I)
Nmin12=abs(Ngrid(I))

[V,I]=min(Etotb115)
Emin115=Etotb115(I)
Nmin115=abs(Ngrid(I))

[V,I]=min(Etotb11)
Emin11=Etotb11(I)
Nmin11=abs(Ngrid(I))

[V,I]=min(Etotb1075)
Emin1075=Etotb1075(I)
Nmin1075=abs(Ngrid(I))

[V,I]=min(Etotb105)
Emin105=Etotb105(I)
Nmin105=abs(Ngrid(I))

[V,I]=min(Etotb1025)
Emin1025=Etotb1025(I)
Nmin1025=abs(Ngrid(I))

[V,I]=min(Etotb10)
Emin10=Etotb10(I)
Nmin10=abs(Ngrid(I))

figure
plot([300 342 400 480 600 800 1000 1043 1090 1116 1143 1170 1200]*29/30,[Emin40 Emin35 Emin30 Emin25 Emin20 Emin15 Emin12 Emin115 Emin11 Emin1075 Emin105 Emin1025 Emin10],'-o')
set(gca,'FontSize',30)
grid on
xlabel('T(K)')
ylabel('E Ground State')

figure
plot([300 342 400 480 600 800 1000 1043 1090 1116 1143 1170 1200]*29/30,[Nmin40 Nmin35 Nmin30 Nmin25 Nmin20 Nmin15 Nmin12 Nmin115 Nmin11 Nmin1075 Nmin105 Nmin1025 Nmin10]*gb/(2*kbulk),'-o')
set(gca,'FontSize',30)
grid on
xlabel('T(K)')
ylabel('N Equilibrium')




