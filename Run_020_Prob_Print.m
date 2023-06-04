## A Convective Fractional & Probabilistic Inelastic Damage
## Constitutive Model for Double-Leaf Stone Masonry Under Cyclic Loads
##  Master’s Thesis University of Basilicata, Italy June 2018
##  Structural Analysis of Monuments and Historical Constructions - Heritage & Intervention Design
##  © Copyright 2018 by Fernando Maximo Valdivia-Vilca All Rights Reserved
##  cip73651@gmail.com

  clear all; close all; clc
  graphics_toolkit qt
load Magenes_StressStrain.txt;  % load experimental data
eps =  [0.000000  0.000660  0.000335  0.001040  0.000610  0.001585  0.001040  0.002340  0.001650  0.003660  0.002750  0.004175  0.003000];
d_eps = 0.000005;
epsL01 =   eps(1) : d_eps :  eps(2) ;  epsU01 =   eps(2) :-d_eps :  eps(3) ;
epsL02 =   eps(3) : d_eps :  eps(4) ;  epsU02 =   eps(4) :-d_eps :  eps(5) ;
epsL03 =   eps(5) : d_eps :  eps(6) ;  epsU03 =   eps(6) :-d_eps :  eps(7) ;
epsL04 =   eps(7) : d_eps :  eps(8) ;  epsU04 =   eps(8) :-d_eps :  eps(9) ;
epsL05 =   eps(9) : d_eps : eps(10) ;  epsU05 =  eps(10) :-d_eps : eps(11) ;
epsL06 =  eps(11) : d_eps : eps(12) ;  epsU06 =  eps(12) :-d_eps : eps(13) ;

#Pushover(       f2, f4, f6, eps, eps0, sig0, epsY, sigY, epsU, sigU, epsT, sigT )
%epsP =  0.00:0.0001:0.0045;  ###                      0              Y                 U                T
%sigP  =  Pushover( 280,	170,	375,  epsP,  0.00,0.00,   0.00094,1.79,   0.00415,3.0925,   0.010,0.25 );

E0 = 3445;  %  degradation vector
epsL = [ epsL01(1)  epsL02(1)  epsL03(1)  epsL04(1)  epsL05(1)  epsL06(1) ] ;
epsU = [ epsL01(end)  epsL02(end)  epsL03(end)  epsL04(end)  epsL05(end)  epsL06(end) ] ;

cum_epsL = cumsum([ 0 (eps(2) - eps(1))  (eps(4) - eps(3))  (eps(6) - eps(5))   (eps(8) - eps(7))  (eps(10) - eps(9))    ]) ;  % loading ranges only
cum_epsU = cumsum([ (eps(2) - eps(1))  (eps(4) - eps(3))  (eps(6) - eps(5))   (eps(8) - eps(7))  (eps(10) - eps(9))   (eps(12) - eps(11))  ]) ;  % Unloading ranges only


%epsL = [ epsU01(end)  epsU02(end)  epsU03(end)  epsU04(end)  epsU05(end)  epsU06(end) ] ;

sigU = Pushover( 280, 170, 375,  epsU,  0.00,0.00,   0.00094,1.79,   0.00415,3.0925,   0.010,0.25 );

f1 = [  0.9234	1.1494	1.3999	1.6198	1.774	1.8148 ];

f3 = [   2.794	  3.06	  3.502	  4.187	  5.172	  5.652  ] ;

D = [   3.6547	3.6073	  3.495	  3.2852	  3.0659	  3.0276  ] ;




%%%%%%%%%%%  completed ********************************************************
Ampx =   [  19113	  10762	  7298	  4750	  3844	  3040  ];
Ampy =   [  0.6535	1.0935	1.4915	1.909	2.1924	2.476  ] ;
nN   =   [  7        6        5       3        3        3       ] ;

a    = [0.47844   0.48378   0.48947   0.49735   0.50786   0.52410 ] ;
b    = [  0.07	0.104	0.12	0.129 	0.138	  0.151 ] ;
c    = [0.684	      0.68	    0.67	    0.64	0.62	    0.615 ];
d    = [ 1.0358   1.0766   1.1201   1.1803   1.2606   1.3847 ] ;
%%%%%%%%%%%  completed ********************************************************


%        Sig_Load( Ampx, Ampy, xx, nN, a, b, c, d, eps_00, sig_00 )
%        Sig_Unload( f1,     f3,   E0,   d,    eps,    eps0,   sig0)
eps_0 = 0; sig_0 = 0;
sigL01 =   Sig_Load( Ampx(1), Ampy(1), epsL01, nN(1), a(1), b(1), c(1), d(1), eps_0, sig_0 );   eps_0=epsL01(end); sig_0=sigL01(end);
sigU01 = Sig_Unload( f1(1), f3(1), E0,  D(1), epsU01, eps_0, sig_0 );                   eps_0=epsU01(end); sig_0=sigU01(end);
sigL02 =   Sig_Load( Ampx(2), Ampy(2), epsL02, nN(2), a(2), b(2), c(2), d(2), eps_0, sig_0 );   eps_0=epsL02(end); sig_0=sigL02(end);
sigU02 = Sig_Unload( f1(2), f3(2), E0,  D(2), epsU02, eps_0, sig_0 );                   eps_0=epsU02(end); sig_0=sigU02(end);
sigL03 =   Sig_Load( Ampx(3), Ampy(3), epsL03, nN(3), a(3), b(3), c(3), d(3), eps_0, sig_0 );   eps_0=epsL03(end); sig_0=sigL03(end);
sigU03 = Sig_Unload( f1(3), f3(3), E0,  D(3), epsU03, eps_0, sig_0 );                   eps_0=epsU03(end); sig_0=sigU03(end);
sigL04 =   Sig_Load( Ampx(4), Ampy(4), epsL04, nN(4), a(4), b(4), c(4), d(4), eps_0, sig_0 );   eps_0=epsL04(end); sig_0=sigL04(end);
sigU04 = Sig_Unload( f1(4), f3(4), E0,  D(4), epsU04, eps_0, sig_0 );                   eps_0=epsU04(end); sig_0=sigU04(end);
sigL05 =   Sig_Load( Ampx(5), Ampy(5), epsL05, nN(5), a(5), b(5), c(5), d(5), eps_0, sig_0 );   eps_0=epsL05(end); sig_0=sigL05(end);
sigU05 = Sig_Unload( f1(5), f3(5), E0,  D(5), epsU05, eps_0, sig_0 );                   eps_0=epsU05(end); sig_0=sigU05(end);
sigL06 =   Sig_Load( Ampx(6), Ampy(6), epsL06, nN(6), a(6), b(6), c(6), d(6), eps_0, sig_0 );   eps_0=epsL06(end); sig_0=sigL06(end);
sigU06 = Sig_Unload( f1(6), f3(6), E0,  D(6), epsU06, eps_0, sig_0 );                   eps_0=epsU06(end); sig_0=sigU06(end);


figure
plot( Magenes_StressStrain(:,1), Magenes_StressStrain(:,2), '.m', 'MarkerSize', 6 ); hold on
plot( epsL01, sigL01, '--b','LineWidth',2);
plot( epsU01, sigU01, '--b','LineWidth',2);
plot( epsL02, sigL02, '--b','LineWidth',2);
plot( epsU02, sigU02, '--b','LineWidth',2);
plot( epsL03, sigL03, '--b','LineWidth',2);
plot( epsU03, sigU03, '--b','LineWidth',2);
plot( epsL04, sigL04, '--b','LineWidth',2);
plot( epsU04, sigU04, '--b','LineWidth',2);
plot( epsL05, sigL05, '--b','LineWidth',2);
plot( epsU05, sigU05, '--b','LineWidth',2);
plot( epsL06, sigL06, '--b','LineWidth',2);
plot( epsU06, sigU06, '--b','LineWidth',2);
xlim([0.0  0.005]);
ylim([0  4]); grid on
title('Vertical Compression Stress-Strain')
legend('Experimental','Proposed model')
xlabel('Strain') % x-axis label
ylabel('Stress (MPa)') % y-axis label

 % plot parameters for loading  Ampx Ampy nN a b c d
#figure
#subplot(3,1,1);

figure
plot(cum_epsL, Ampx, 'r--x', 'MarkerSize', 10)
title('Local strain factor AmpX calibration')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('AmpX') % y-axis label

#subplot(3,1,2);

figure
plot(cum_epsL, Ampy, 'b--x', 'MarkerSize', 10)
title('Global stress factor AmpY calibration')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('AmpY') % y-axis label

#subplot(3,1,3);


figure
plot(cum_epsL,  nN, '--x', 'MarkerSize', 10)
title('Number of fractional contributing N calibration')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('N') % y-axis label

figure
plot(cum_epsL, a, 'r--x', 'MarkerSize', 10)
title(' fractional coefficient a ')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('a') % y-axis label

figure
plot(cum_epsL, b, 'r--x', 'MarkerSize', 10)
title(' Post yielding pseudo hardening b calibration ')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('b') % y-axis label

figure
plot(cum_epsL, c, 'r--x', 'MarkerSize', 10)   # k replaces c for plots
title('Diffusion function coefficient  k calibration')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('k') % y-axis label

figure
plot(cum_epsL, d, 'r--x', 'MarkerSize', 10)
title(' Diffusion equation number d calibration')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('d') % y-axis label



% plots for Unloading f1 f3 D
#figure
#subplot(2,2,1);

figure
plot( cum_epsU , f1, 'r--x', 'MarkerSize', 10)
title('Boundary function f1 - calibration')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('f1') % y-axis label

#subplot(2,2,3);

figure
plot( cum_epsU, f3, 'b--x', 'MarkerSize', 10)
title(' Boundary function f3 - calibration')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('f3') % y-axis label

#subplot(2,2,[2,4]);

figure
plot( cum_epsU,  D  , '--x', 'MarkerSize', 10)
title(' Local degradation function D - calibration ')
xlabel('Cumulative Loading Strain') % x-axis label
ylabel('D') % y-axis label









