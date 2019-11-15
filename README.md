# SEEC
SEEC model library

This library contains the script need to use the SEEC evolutionary model as described in:

"Epistatic contributions promote the unification of incompatible models of neutral molecular evolution" Jose Alberto de la Paz, Charisse M. Nartey, Monisha Yuvaraj, Faruck Morcos. 

Scripts are intended to be read using MATLAB. The library contains:

Generalhamiltonian.m
General Potts Hamiltonian

Neoheterotachy.m
Calculation of rates and heterotachy degree from 1000 evolutionary simulations

Probevolution.m
Code to generate evolutionary trajectories

Rates.m
Inference of the individual site rates. Input generated with RatesEvolution.m

RatesEvolution.m
Code to generate evolutionary trajectories (specific for evolutionary rates measurement)

Stokes.m
Stokes Shift calculation

effAlphabet.m
Effective alphabet measure from the single-site probability distribution along trajectory

siteprobdistribution.m
Site probability distribution for evolution scripts


All of these codes are based on MATLAB. Below is one example for the basic analysis of lambda SEQRS data. 


Any publication resulting from applications of SEEC should cite:

Jose Alberto de la Paz, Charisse M. Nartey, Monisha Yuvaraj, Faruck Morcos, Epistatic contributions promote the unification of incompatible models of neutral molecular evolution. Under revision.

