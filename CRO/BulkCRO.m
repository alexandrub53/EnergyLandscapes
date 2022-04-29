%Script Used To Plot Ca2RuO4, bulk, unstrained data and perform fit.
%Notation follows notation in paper.
clear all


Delta=[0.07608346709470304
0.10698234349919741    
0.1394863563402889
0.1687800963081862
0.1964686998394863
0.2237560192616371
0.2490369181380417]

QlQ=[0.017718737601442727
0.0018524436429215169
-0.013203606853020772
-0.028753074842200202
-0.061188133453561784
-0.07774456266907126]

eta=[ -0.053997569550496616
    0.007031407386226984 
    0.09479662856363047
    0.22139705356736716
      0.9166680746810847
       0.963067565196483  ]




QQiang=QlQ(:,1);
NQiang=eta(:,1);
lambda0=0.45;
F3=-2.8;

objectiveFunc=@(p)sum((F3*QQiang'-p(1)-p(2)*NQiang'-p(3)*NQiang'.^2-p(4)*NQiang'.^3).^2)
p0 = [ 1 1 1 1 ];
p = fminsearch( objectiveFunc, p0 );
NgridPlot=[-0.3:0.001:1.2]
s0=p(1)
a=p(2);
b=p(3);
c=p(4);
d=0;
figure

hold on
plot((d*NgridPlot.^5+c*NgridPlot.^3+b*NgridPlot.^2+a*NgridPlot.^1+s0)/(F3),NgridPlot,'k','LineWidth',7)
hold on
plot(QQiang,NQiang,'dk','MarkerSize',20,'LineWidth',10)

set(gca,'FontSize',30)
axis([-0.1 0.02 -0.1 1.0 ])
legend('Polynomial Fit Ca_2RuO_4','DFT+DMFT Data Ca_2RuO_4')
grid on
xlabel('Q_3-\lambda_0 Q_0 (Angstrom)')
ylabel('\Delta N(electrons)')




K=[17.69 7.63
   7.63   46.16]
k33=K(1,1)
k03=K(1,2)
k00=K(2,2)



CorrD=-F3^2*[ 1 -lambda0]*inv(K)*[1 ; -lambda0];


ap=a+CorrD;
figure
NgridPlot=[-0.3:0.0001:1.0]

plot(NgridPlot,(d*NgridPlot.^6)/6+c*(NgridPlot.^4)/4+b*(NgridPlot.^3)/3+a*(NgridPlot.^2)/2+s0*NgridPlot,'-o','LineWidth',4,'MarkerSize',4)
hold on
plot(NgridPlot,(d*NgridPlot.^6)/6+c*(NgridPlot.^4)/4+b*(NgridPlot.^3)/3+ap*(NgridPlot.^2)/2+s0*NgridPlot,'-d','LineWidth',7,'MarkerSize',7)
axis([-0.3 1.0 -0.22 0.15])
legend('Ca_2RuO_4 F_{el}', 'Ca_2RuO_4 F_{tot} ' )
xlabel('\Delta N')
ylabel('Energy (eV)')
set(gca,'FontSize',30)
grid on
