%Script to Generate CRO Total Energy W Applied Epitaxial Strain
clear all
K=[17.69 7.63
   7.63   46.16]
k33=K(1,1)
k03=K(1,2)
k00=K(2,2)

QlQ=[0.017718737601442727
0.0018524436429215169
-0.013203606853020772
-0.028753074842200202
-0.061188133453561784
-0.07774456266907126]

Delta=[0.07608346709470304
0.10698234349919741    
0.1394863563402889
0.1687800963081862
0.1964686998394863
0.2237560192616371
0.2490369181380417]

eta=[ -0.053997569550496616
    0.007031407386226984 
    0.09479662856363047
    0.22139705356736716
      0.9166680746810847
       0.963067565196483]


%plot(QlQ,eta,'-o')





%plot(QiangCROStrain(:,1),QiangCROStrain2(:,1),'-o')

QQiang=QlQ(:,1);
NQiang=eta(:,1);
lambda0=0.45
%plot(QQiang,NQiang,'-o')

%F3=2.8;
F3=-2.8

objectiveFunc=@(p)sum((F3*QQiang'-p(1)-p(2)*NQiang'-p(3)*NQiang'.^2-p(4)*NQiang'.^3).^2)
p0 = [ 1 1 1 1  ];
p = fminsearch( objectiveFunc, p0 );

s0=p(1)
a=p(2);
b=p(3);
c=p(4);
d=0;

maxZ=0.13;
maxN=1.2;
gridtot=1000;
Zy=zeros(gridtot,gridtot);
Nx=zeros(gridtot,gridtot);
for k=1:gridtot
    for j=1:gridtot
        Zy(j,k)=(k-gridtot/2)*maxZ/(gridtot/2);
        Nx(j,k)=(j-gridtot/2)*maxN/(gridtot/2);
    end
end

NgridPlot=Nx(:,1)';
ZgridPlot=Zy(1,:);

x0p0grid=-0.08

for p=[1:length(x0p0grid)]
x0p0=x0p0grid(p);
lambda=lambda0;
for k=1:gridtot
    for j=1:gridtot
        Zy(j,k)=(k-gridtot/2)*maxZ/(gridtot/2);
        Nx(j,k)=(j-gridtot/2)*maxN/(gridtot/2);
        ECOR(j,k)=(d*Nx(j,k).^6)/6+c*(Nx(j,k).^4)/4+b*(Nx(j,k).^3)/3+a*(Nx(j,k).^2)/2+s0*Nx(j,k);
        ELIN=-Nx(j,k)*F3/sqrt(6)*(((2-lambda*sqrt(2))*Zy(j,k)-(1+sqrt(2)*lambda)*x0p0));
        ESQUARE=(Zy(j,k)^2)/3*(k33+sqrt(2)*k03+1/2*k00)-Zy(j,k)*x0p0/3*(k33-k00-k03/sqrt(2));
        ECRO(j,k)=ECOR(j,k)+ELIN+ESQUARE;
    end
end
EvsN=(min(ECRO'));
end


plot(Nx(:,1),EvsN,'Linewidth',3)
hold on

x0p0grid=0.016

for p=[1:length(x0p0grid)]
x0p0=x0p0grid(p);
lambda=lambda0;
for k=1:gridtot
    for j=1:gridtot
        Zy(j,k)=(k-gridtot/2)*maxZ/(gridtot/2);
        Nx(j,k)=(j-gridtot/2)*maxN/(gridtot/2);
        ECOR(j,k)=(d*Nx(j,k).^6)/6+c*(Nx(j,k).^4)/4+b*(Nx(j,k).^3)/3+a*(Nx(j,k).^2)/2+s0*Nx(j,k);
        ELIN=-Nx(j,k)*F3/sqrt(6)*(((2-lambda*sqrt(2))*Zy(j,k)-(1+sqrt(2)*lambda)*x0p0));
        ESQUARE=(Zy(j,k)^2)/3*(k33+sqrt(2)*k03+1/2*k00)-Zy(j,k)*x0p0/3*(k33-k00-k03/sqrt(2));
        ECRO(j,k)=ECOR(j,k)+ELIN+ESQUARE;
    end
end
EvsN=(min(ECRO'));
end


plot(Nx(:,1),EvsN,'Linewidth',5)
hold on


x0p0grid=0.06

for p=[1:length(x0p0grid)]
x0p0=x0p0grid(p);
lambda=lambda0;
for k=1:gridtot
    for j=1:gridtot
        Zy(j,k)=(k-gridtot/2)*maxZ/(gridtot/2);
        Nx(j,k)=(j-gridtot/2)*maxN/(gridtot/2);
        ECOR(j,k)=(d*Nx(j,k).^6)/6+c*(Nx(j,k).^4)/4+b*(Nx(j,k).^3)/3+a*(Nx(j,k).^2)/2+s0*Nx(j,k);
        ELIN=-Nx(j,k)*F3/sqrt(6)*(((2-lambda*sqrt(2))*Zy(j,k)-(1+sqrt(2)*lambda)*x0p0));
        ESQUARE=(Zy(j,k)^2)/3*(k33+sqrt(2)*k03+1/2*k00)-Zy(j,k)*x0p0/3*(k33-k00-k03/sqrt(2));
        ECRO(j,k)=ECOR(j,k)+ELIN+ESQUARE;
    end
end
EvsN=(min(ECRO'));
end


plot(Nx(:,1),EvsN,'Linewidth',7)
hold on

x0p0grid=0.092

for p=[1:length(x0p0grid)]
x0p0=x0p0grid(p);
lambda=lambda0;
for k=1:gridtot
    for j=1:gridtot
        Zy(j,k)=(k-gridtot/2)*maxZ/(gridtot/2);
        Nx(j,k)=(j-gridtot/2)*maxN/(gridtot/2);
        ECOR(j,k)=(d*Nx(j,k).^6)/6+c*(Nx(j,k).^4)/4+b*(Nx(j,k).^3)/3+a*(Nx(j,k).^2)/2+s0*Nx(j,k);
        ELIN=-Nx(j,k)*F3/sqrt(6)*(((2-lambda*sqrt(2))*Zy(j,k)-(1+sqrt(2)*lambda)*x0p0));
        ESQUARE=(Zy(j,k)^2)/3*(k33+sqrt(2)*k03+1/2*k00)-Zy(j,k)*x0p0/3*(k33-k00-k03/sqrt(2));
        ECRO(j,k)=ECOR(j,k)+ELIN+ESQUARE;
    end
end
EvsN=(min(ECRO'));
end


plot(Nx(:,1),EvsN,'Linewidth',7)
hold on



legend('NdAlO_3','LaAlO_3','NSAT','NdGaO3')
xlabel('Delta N')
ylabel('Total Energy(eV)')
set(gca,'FontSize',30)
grid on
%plot(ZgridPlot,min(ECRO)')