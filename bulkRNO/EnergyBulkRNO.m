clear all

BulkData

[ Ngrid,EcorL,EtotL,AL,BL,CL ] = getEnergyGrid( NL,QL,gL,kL );

plot(Ngrid,EtotL,'r','LineWidth',5)

hold on

[ Ngrid,EcorS,EtotS,AS,BS,CS ] = getEnergyGrid( NS,QS,gS,kS );

plot(Ngrid,EtotS,'b','LineWidth',7)



[ Ngrid,EcorP,EtotP,AP,BP,CP ] = getEnergyGrid( NP,QP,gP,kP );

plot(Ngrid,EtotP,'g-','LineWidth',9)


set(gca,'FontSize',30)
grid on
axis([-2 2 -0.09 0.05])
xlabel('\Delta N(electrons)')
ylabel('E_{tot}(eV)')


legend('LuNiO3','SmNiO3','PrNiO3')
set(gcf,'position',[0,0,1000,1000])

% Print Useful Constants
[CL, BL, AL,-1/4*gL^2/kL]
[CS, BS, AS,-1/4*gS^2/kS]
[CP, BP, AP,-1/4*gP^2/kP]
