% Script containing data used to generate total energies for NdNiO3, and
% NdNiO3/NdAlO3 heterostructures

% Variables containing bulk, or ending in b correspond to bulk NdNiO3; het
% or h to a heterostructure containing of 1 layer NdNiO3 laternating with 3 layers of NdAlO3 
% while het22 or h22 corresponds to a 2/2 NdNiO3/NdAlO3 layered structure.
% The data here used in our paper, was originally generated by us for:
% https://www.pnas.org/doi/abs/10.1073/pnas.1818728116 
% 
clear all


het=[
0.22	1.386688    
0.21	1.376951
0.2	1.360838
0.19	1.347948
0.18	1.334966
0.17	1.313047
0.16	1.295357
0.15	1.271522
0.14	1.251725
0.13	1.22334
0.12	1.187036
0.11	1.141044
0.1	1.069734
0.09	0.887644
0.08	0.679954
0.07	0.55088
0.06	0.453436
0.05	0.36733
0.04	0.291844
0.03	0.214917
0.02	0.147542
0.01	0.081025
0.0  	0.009997
]

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

het22=[
0.22	1.286193
0.21	1.271699 
0.2	1.254759
0.19	1.243578
0.18	1.227261
0.17	1.207059
0.16	1.187521
0.15	1.167506
0.14	1.150411
0.13	1.122942
0.12	1.09475
0.11	1.061225
0.1	1.020128
0.09	0.952842
0.08	0.80315
0.07	0.56011
0.06	0.442974
0.05	0.352573
0.04	0.2746
0.03	0.202595
0.02	0.136836
0.01	0.064419
0	0.005059
    ]



U=2.1
J=1.2
DC=U-5/3*J
x=[0:0.01:0.1]


chi0b=1.159592851
chih220= 1.246190802
chih0=1.392567967

gb0=2.8865
gh0=2.965066369
gh220=2.975670521



gb=gb0*(1+chi0b*DC)
gh=gh0*(1+chih0*DC)
gh22=gh220*(1+chih220*DC)





kbulk=2*13.4303+1/2*gb0*chi0b*gb0
by=2*kbulk/gb*x

kh=2*17.1769+1/2*gh0*chih0*gh0
hy=2*kh/gh*x

kh22=2*15.0548+1/2*gh220*chih220*gh220
h22=2*kh22/gh22*x






Nb=bulk(: , 2);
Qb=bulk(: , 1)/gb;
Nh=het(: , 2);
Qh=het(: , 1)/gh;
Nh22=het22(: , 2);
Qh22=het22(: , 1)/gh22;


[ Ngrid,Ecorb,Etotb,Ab,Bb,Cb ] = getEnergyGrid( Nb,Qb,gb,kbulk )
[ Ngrid,Ecorh,Etoth,Ah,Bh,Ch ] = getEnergyGrid( Nh,Qh,gh,kh )
[ Ngrid,Ecorh22,Etoth22,Ah22,Bh22,Ch22 ] = getEnergyGrid( Nh22,Qh22,gh22,kh22 )

figure('Position',[0 0 1080 1080])
set(gca,'FontSize',30)
grid on
plot(Ngrid,Ecorb*1000,'LineWidth',9)
hold on
plot(Ngrid,Ecorh22*1000,'LineWidth',7)
plot(Ngrid,Ecorh*1000,'LineWidth',5)
set(gca,'FontSize',30)
xlabel('\Delta N(electrons)');
ylabel('E_{el}(meV)')
grid on
axis([-1.7 1.7 -0 50])
legend('Bulk NdNiO_3','(NdNiO_3)_2/(NdAlO_3)_2','(NdNiO_3)_1/(NdAlO_3)_3')

% Print out useful variables
[Cb,Bb,Ab,-gb^2/(4*kbulk)]
[Ch22,Bh22,Ah22,-gh22^2/(4*kh22)]
[Ch,Bh,Ah,-gh^2/(4*kh)]
box on
figure

% Plot total Energies, in meV
plot(Ngrid,Etotb*1000,'LineWidth',9)
hold on
box on
plot(Ngrid,Etoth22*1000,'LineWidth',7)
plot(Ngrid,Etoth*1000,'LineWidth',5)
xlabel('\Delta N(electrons)');
ylabel('E_{tot}(meV)')
grid on
axis([-1.7 1.7 -17 5])
set(gca,'FontSize',30)
box on

