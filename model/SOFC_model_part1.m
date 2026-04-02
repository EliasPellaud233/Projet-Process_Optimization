function y = SOFC_model_part1(u)

% Inputs
qH2_in  = u(1);   % [L/min]
qAir_in = u(2);   % [L/min]
I       = u(3);   % [A]
T       = u(4);   % [K]

%% Parameters
Ncells = 5;
ne = 2;
F = 96485;
R = 8.31451;
AAct = 50e-4;
he = 10e-6;
RMIC1 = 4.9361e-9;
RMIC2 = 4.9361e-9;
p0 = 1; % [bar]

Eact_c = 153260.5;
k0_c = 4.1031e12;
Edis_c = 1.4895e5;
R0_c = 9.225228e-14;
Ee = 7.6223e4;
sigma0_e = 1.63e4;

Tref = 298.15;
E0_ref = 1.184;
dE0dT = -230e-6;

LHV_H2 = 242000;

%% Inlet molar flow rates [mol/s]
nH2_an_in = (1/22.4) * (1/60) * qH2_in;
nH2O_an_in = (3/97) * nH2_an_in;

nO2_cath_in = (1/22.4) * (1/60) * (21/100) * qAir_in;
nN2_cath_in = (79/21) * nO2_cath_in;

%% Fuel utilization and air excess ratio
FU = (Ncells/(2*F)) * I / nH2_an_in;
lambda_air = 2 * nO2_cath_in / nH2_an_in;

%% Reaction rates [mol/s]
nH2_r  = -Ncells/(2*F) * I;
nO2_r  = -(1/2)*Ncells/(2*F) * I;
nH2O_r =  Ncells/(2*F) * I;
nN2_r  = 0;

%% Outlet molar flow rates [mol/s]
nH2_an_out  = nH2_an_in  + nH2_r;
nH2O_an_out = nH2O_an_in + nH2O_r;

nO2_cath_out = nO2_cath_in + nO2_r;
nN2_cath_out = nN2_cath_in + nN2_r;

%% Total outlet flow rates
ntot_an_out = nH2_an_out + nH2O_an_out;
ntot_cath_out = nO2_cath_out + nN2_cath_out;

%% Outlet mole fractions
xH2_an_out = nH2_an_out / ntot_an_out;
xH2O_an_out = nH2O_an_out / ntot_an_out;

xO2_cath_out = nO2_cath_out / ntot_cath_out;
xN2_cath_out = nN2_cath_out / ntot_cath_out;

%% Inlet mole fractions
ntot_an_in = nH2_an_in + nH2O_an_in;
ntot_cath_in = nO2_cath_in + nN2_cath_in;

xH2_an_in = nH2_an_in / ntot_an_in;
xH2O_an_in = nH2O_an_in / ntot_an_in;

xO2_cath_in = nO2_cath_in / ntot_cath_in;
xN2_cath_in = nN2_cath_in / ntot_cath_in;

%% Partial pressures [bar]
pH2_an_out = xH2_an_out * p0;
pH2O_an_out = xH2O_an_out * p0;
pO2_cath_out = xO2_cath_out * p0;
pO2_cath_in = xO2_cath_in * p0;

%% Standard potential
E0 = E0_ref + dE0dT * (T - Tref);

%% Redox potential
E = E0 - (R*T/(ne*F)) * log( pH2O_an_out / (pH2_an_out * sqrt(pO2_cath_out)) );

%% Current density
i = I / AAct;

%% Cathode activation loss
i0_c = (R*T/F) * k0_c * exp(-Eact_c/(R*T));
Uact_c = (R*T/F) * asinh(i/(2*i0_c));

%% Electrolyte ohmic loss
sigma_i_e = sigma0_e * exp(-Ee/(R*T));
Ui_e = i * he / sigma_i_e;

%% Cathode dissociation loss
Udis_c = R0_c * sqrt(p0/pO2_cath_in) * exp(-Edis_c/(R*T)) * i;

%% Diffusion losses
Udif_a = -(R*T/(2*F)) * log(1 - FU);
Udif_c = -(R*T/(2*F)) * log(1 - FU/lambda_air);

%% MIC loss
UMIC = (RMIC1 + RMIC2) * I;

%% Cell voltage
Ucell = E - Uact_c - Ui_e - Udis_c - Udif_a - Udif_c - UMIC;

%% Electrical power and efficiency
Pel = Ncells * Ucell * I;
eta_el = Pel / (nH2_an_in * LHV_H2) * 100;

%% Output vector
y = [Pel; Ucell; FU; lambda_air; eta_el];

end