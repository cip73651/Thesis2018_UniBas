## A Convective Fractional & Probabilistic Inelastic Damage
## Constitutive Model for Double-Leaf Stone Masonry Under Cyclic Loads
##  Master’s Thesis University of Basilicata, Italy June 2018
##  Structural Analysis of Monuments and Historical Constructions - Heritage & Intervention Design
##  © Copyright 2018 by Fernando Maximo Valdivia-Vilca All Rights Reserved
##  cip73651@gmail.com


function [ Push_over ] = Sig_Unload( f1, f3, E0, d, eps, eps0, sig0) %, epsE, sigE, epsY, sigY, epsU, sigU, epsT, sigT )
%UNTITLED2 Summary of this function goes here

%E0   = 3445;    d = 3.2; %df = 0.085;# b    = .25;
#E0   = 3475; d    = 1; #

f2 = (E0 * exp(sig0 / f1) * d * eps0 * LambertW(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) * log(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) / 0.2e1) * f1 * f3 * (-0.1e1 + exp(sig0 / f1)) * (0.1e1 + exp(sig0 / f1)) / E0 / d / eps0 * exp((f1 ^ 2 * log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) * (-0.1e1 + exp(sig0 / f1)) * (0.1e1 + exp(sig0 / f1)) * log(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) / 0.2e1) - 0.2e1 * sig0 * E0 * exp(sig0 / f1) * d * eps0) / f1 / E0 / exp(sig0 / f1) / d / eps0 / 0.2e1) / 0.2e1) - log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) * f1 * log(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) / 0.2e1) * (-0.1e1 + exp(sig0 / f1)) * (0.1e1 + exp(sig0 / f1)) / 0.2e1) * f3 / E0 / exp(sig0 / f1) / d / eps0 ^ 2 / LambertW(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) * log(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) / 0.2e1) * f1 * f3 * (-0.1e1 + exp(sig0 / f1)) * (0.1e1 + exp(sig0 / f1)) / E0 / d / eps0 * exp((f1 ^ 2 * log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) * (-0.1e1 + exp(sig0 / f1)) * (0.1e1 + exp(sig0 / f1)) * log(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) / 0.2e1) - 0.2e1 * sig0 * E0 * exp(sig0 / f1) * d * eps0) / f1 / E0 / exp(sig0 / f1) / d / eps0 / 0.2e1) / 0.2e1);
f4 = 0.2e1 * eps0 * log(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) / 0.2e1) * exp(sig0 / f1) * E0 * d / (0.2e1 * E0 * exp(sig0 / f1) * d * eps0 * LambertW(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) * log(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) / 0.2e1) * f1 * f3 * (-0.1e1 + exp(sig0 / f1)) * (0.1e1 + exp(sig0 / f1)) / E0 / d / eps0 * exp((f1 ^ 2 * log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) * (-0.1e1 + exp(sig0 / f1)) * (0.1e1 + exp(sig0 / f1)) * log(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) / 0.2e1) - 0.2e1 * sig0 * E0 * exp(sig0 / f1) * d * eps0) / f1 / E0 / exp(sig0 / f1) / d / eps0 / 0.2e1) / 0.2e1) - log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) * f1 * log(log((0.1e1 + exp(sig0 / f1)) / (-0.1e1 + exp(sig0 / f1))) / 0.2e1) * (-0.1e1 + exp(sig0 / f1)) * (0.1e1 + exp(sig0 / f1)));

Push_over = -f1 .* log(tanh((-f2 .* eps + f3) .^ f4));

end







