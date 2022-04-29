%Constants and data extracted from https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.245127

LuNiO3=[0.0 0.0
0.001327433628318591 0.0482625482625485
0.002389380530973456 0.0791505791505791
0.003982300884955756 0.12548262548262556
0.005044247787610624 0.1563706563706566
0.006460176991150448 0.2065637065637067
0.014955752212389383 1.333976833976834
0.01743362831858408 1.3494208494208495
0.022566371681415936 1.3725868725868726
0.027699115044247793 1.3996138996138998
0.03000000000000001 1.4034749034749037
0.03495575221238939 1.422779922779923
0.04008849557522125 1.4382239382239383
0.044867256637168146 1.44980694980695
0.05 1.4652509652509653
0.054955752212389394 1.480694980694981
0.05991150442477877 1.4961389961389964
0.06486725663716815 1.507722007722008
0.07 1.5193050193050195
0.07495575221238938 1.530888030888031
0.07991150442477878 1.5424710424710426
0.08486725663716815 1.5540540540540542
0.09 1.5617760617760619
0.09495575221238939 1.5810810810810811
0.10982300884955752 1.6003861003861006];

SmNiO3=[0.0  0.0
0.005044247787610624  0.09073359073359089
0.010000000000000002  0.20270270270270285
%0.014955752212389383  0.3262548262548264
%0.020088495575221247  0.5077220077220079
0.024867256637168152  1.1911196911196913
0.03000000000000001  1.2374517374517375
0.03477876106194691  1.2683397683397684
0.03991150442477877  1.2876447876447878
0.044867256637168146  1.310810810810811
0.05  1.333976833976834
0.054955752212389394  1.3571428571428572
0.05991150442477877  1.3764478764478767
0.06504424778761063  1.391891891891892
0.07  1.4073359073359075
0.07991150442477878  1.4420849420849422
0.09  1.4691119691119692
0.09991150442477877  1.4922779922779923
0.10964601769911504  1.5154440154440156];

PrNiO3=[0.0 0.0
0.0006194690265486774 0.013513513513513598
0.004867256637168145 0.05984555984555984
0.010000000000000002 0.13320463320463327
0.015132743362831862 0.20270270270270285
0.020088495575221247 0.2799227799227799
0.025044247787610625 0.35328185328185335
0.03017699115044249 0.43436293436293427
0.03513274336283187 0.5193050193050193
0.05 1.1756756756756759
0.054955752212389394 1.2142857142857144
0.06008849557522125 1.2451737451737452
0.06486725663716815 1.2722007722007724
0.07 1.2953667953667956
0.0747787610619469 1.3185328185328187
0.08008849557522124 1.337837837837838
0.08486725663716815 1.3571428571428572
0.09 1.3764478764478767
0.09495575221238939 1.391891891891892
0.10946902654867256 1.4343629343629345];



NL=LuNiO3(:,2);
QL=LuNiO3(:,1);

gL=3.75;
kL=39.29;

NS=SmNiO3(:,2);
QS=SmNiO3(:,1);

gS=4.02;
kS=41.45;

NP=PrNiO3(:,2);
QP=PrNiO3(:,1);

gP=4.24;
kP=44.47;

