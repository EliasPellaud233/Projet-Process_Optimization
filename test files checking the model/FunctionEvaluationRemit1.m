clear all;clc
%inputs
x(1) = 750 + 273.15; %stack temperature [K]
x(2) = 1.1;          %fuel flowrate [L/min]
x(3) = 40;           %air flowrate [L/min]
x(4) = 17.5;         %electric current[A]

%function evaluation
f = SOFC_model_remit1(x);

%outputs
Pel   = f(1) %electric power [W]
Ucell = f(2) %cell voltage [V]
FU    = f(3) %fuel utilization [-]
lair  = f(4) %air-excess ratio [-]
eff   = f(5) %efficiency [%]


