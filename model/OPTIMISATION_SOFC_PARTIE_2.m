%% ============================================================
%  OPTIMISATION SOFC - PARTIE 2
%
%  Le modèle SOFC_model_part1.m doit être dans le même dossier.
%
%  Variables optimisées :
%       x = [I, qH2, qAir]
%
%  Objectif :
%       maximiser le rendement eta_el
%
%  Contraintes :
%       Pel = Ptarget
%       FU = 0.8
%       750°C <= T <= 800°C
%       Ucell >= 0.7 V
%       4 <= lambda_air <= 15
%% ============================================================

clear; clc; close all;

%% Puissances à atteindre
Ptargets = [50 75 100 130];

%% Bornes des variables : x = [I, qH2, qAir]
lb = [0,    0.01,  1];
ub = [50,   2,     250];

%% Options du solveur
options = optimoptions('fmincon', ...
    'Display', 'off', ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 1000, ...
    'ConstraintTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6);

%% Tableau de résultats
results = zeros(length(Ptargets), 12);

%% Optimisation pour chaque puissance cible
for k = 1:length(Ptargets)

    Ptarget = Ptargets(k);

    fprintf('\n=====================================\n');
    fprintf('Optimisation pour Ptarget = %.0f W\n', Ptarget);
    fprintf('=====================================\n');

    %% Point initial cohérent
    x0 = initial_guess(Ptarget, lb, ub);

    fprintf('Point de depart : I0 = %.2f A, qH2_0 = %.3f L/min, qAir_0 = %.2f L/min\n', ...
        x0(1), x0(2), x0(3));

    %% Optimisation
    [x_opt, fval, exitflag, output] = fmincon( ...
        @(x) objectif(x), ...
        x0, ...
        [], [], [], [], ...
        lb, ub, ...
        @(x) contraintes(x, Ptarget), ...
        options);

    %% Simulation finale
    [T, Pel, Ucell, FU, lambda_air, eta_el] = simulate_SOFC(x_opt);

    %% Vérification des contraintes
    [c, ceq] = contraintes(x_opt, Ptarget);

    %% Stockage
    results(k,:) = [Ptarget, x_opt(1), x_opt(2), x_opt(3), ...
                    T, T-273.15, Pel, Ucell, FU, lambda_air, eta_el, exitflag];

    %% Affichage
    fprintf('\n--- Résultat pour %.0f W ---\n', Ptarget);
    fprintf('I          = %.3f A\n', x_opt(1));
    fprintf('qH2        = %.3f L/min\n', x_opt(2));
    fprintf('qAir       = %.3f L/min\n', x_opt(3));
    fprintf('T          = %.2f K / %.2f °C\n', T, T-273.15);
    fprintf('Pel        = %.2f W\n', Pel);
    fprintf('Ucell      = %.3f V\n', Ucell);
    fprintf('FU         = %.3f\n', FU);
    fprintf('lambda_air = %.3f\n', lambda_air);
    fprintf('eta        = %.2f %%\n', eta_el);
    fprintf('exitflag   = %d\n', exitflag);
    fprintf('Erreur puissance Pel - Ptarget = %.4e W\n', ceq(1));
    fprintf('Erreur FU - 0.8 = %.4e\n', ceq(2));
    fprintf('Contrainte max c = %.4e\n', max(c));

end

%% Tableau final
ResultsTable = array2table(results, ...
    'VariableNames', {'Ptarget_W','I_A','qH2_Lmin','qAir_Lmin', ...
                      'T_K','T_C','Pel_W','Ucell_V','FU','lambda_air', ...
                      'eta_percent','exitflag'});

disp(' ');
disp('=========== RESULTATS POUR TOUTES LES PUISSANCES ===========');
disp(ResultsTable);


%% ============================================================
% POINT DE DEPART
%% ============================================================
function x0 = initial_guess(Ptarget, lb, ub)

    Ncells = 5;
    Ucell_est = 0.80;
    F = 96485;

    % Estimation initiale du courant
    I0 = Ptarget / (Ncells * Ucell_est);

    % Ici on impose déjà FU initial proche de la contrainte finale
    FU0 = 0.8;

    % Débit H2 estimé pour obtenir FU ≈ 0.8
    nH2_in0 = Ncells * I0 / (2 * F * FU0);
    qH2_0 = nH2_in0 * 22.4 * 60;

    % Lambda initial au milieu de la plage admissible
    lambda0 = 8;

    % Débit air estimé
    qAir_0 = lambda0 * qH2_0 / 0.42;

    % Saturation dans les bornes
    x0 = [I0, qH2_0, qAir_0];
    x0 = max(x0, lb);
    x0 = min(x0, ub);

end


%% ============================================================
% OBJECTIF
%% ============================================================
function J = objectif(x)

    [~, ~, ~, ~, ~, eta_el] = simulate_SOFC(x);
    J = -eta_el;

end


%% ============================================================
% CONTRAINTES
%% ============================================================
function [c, ceq] = contraintes(x, Ptarget)

    [T, Pel, Ucell, FU, lambda_air, ~] = simulate_SOFC(x);

    Tmin = 750 + 273.15;
    Tmax = 800 + 273.15;

    % Contraintes d'inégalité : c <= 0
    c = [
        Tmin - T
        T - Tmax
        0.7 - Ucell
        4 - lambda_air
        lambda_air - 15
    ];

    % Contraintes d'égalité : ceq = 0
    ceq = [
        Pel - Ptarget
        FU - 0.8
    ];

end


%% ============================================================
% SIMULATION COMPLETE
%% ============================================================
function [T_end, Pel, Ucell, FU, lambda_air, eta_el] = simulate_SOFC(x)

    I       = x(1);
    qH2_in  = x(2);
    qAir_in = x(3);

    T0 = 1048.15;
    tspan = [0 600];

    ode_options = odeset('RelTol',1e-6,'AbsTol',1e-8);

    [~, T_profile] = ode45(@(t,T) energy_balance(T, qH2_in, qAir_in, I), ...
                           tspan, T0, ode_options);

    T_end = T_profile(end);

    y = SOFC_model_part1([qH2_in, qAir_in, I, T_end]);

    Pel        = y(1);
    Ucell      = y(2);
    FU         = y(3);
    lambda_air = y(4);
    eta_el     = y(5);

end


%% ============================================================
% BILAN ENERGIE
%% ============================================================
function dTdt = energy_balance(T, qH2_in, qAir_in, I)

    Ncells = 5;
    F = 96485;

    mcp = 2500;
    A_stack = 0.0235;
    sigma_SB = 5.6703e-8;
    alpha = 2/3;

    T_furn = 780 + 273.15;
    Tfuel_in = 600 + 273.15;
    Tair_in  = 600 + 273.15;

    %% Débits molaires entrants
    nH2_in  = (1/22.4) * (1/60) * qH2_in;
    nH2O_in = (3/97) * nH2_in;

    nO2_in = (1/22.4) * (1/60) * 0.21 * qAir_in;
    nN2_in = (79/21) * nO2_in;

    %% Réaction H2 + 1/2 O2 -> H2O
    nH2_r  = -Ncells/(2*F) * I;
    nO2_r  = -0.5 * Ncells/(2*F) * I;
    nH2O_r =  Ncells/(2*F) * I;

    %% Débits sortants
    nH2_out  = nH2_in  + nH2_r;
    nH2O_out = nH2O_in + nH2O_r;
    nO2_out  = nO2_in  + nO2_r;
    nN2_out  = nN2_in;

    % Protection numérique
    nH2_out  = max(nH2_out,  1e-12);
    nH2O_out = max(nH2O_out, 1e-12);
    nO2_out  = max(nO2_out,  1e-12);
    nN2_out  = max(nN2_out,  1e-12);

    %% Enthalpies entrantes/sortantes
    h_in  = stack_enthalpy(Tfuel_in, Tair_in);
    h_out = stack_enthalpy(T, T);

    %% Puissance électrique à la température T
    y = SOFC_model_part1([qH2_in, qAir_in, I, T]);
    Pel = y(1);

    %% Pertes radiatives
    Qloss = A_stack * alpha * sigma_SB * (T^4 - T_furn^4);

    %% Bilan énergétique
    H_in = ...
        nH2_in  * h_in(1) + ...
        nH2O_in * h_in(5) + ...
        nO2_in  * h_in(3) + ...
        nN2_in  * h_in(4);

    H_out = ...
        nH2_out  * h_out(1) + ...
        nH2O_out * h_out(5) + ...
        nO2_out  * h_out(3) + ...
        nN2_out  * h_out(4);

    dTdt = (H_in - H_out - Pel - Qloss) / mcp;

end


%% ============================================================
% FONCTION D'ENTHALPIE
%% ============================================================
function h = stack_enthalpy(tempF,tempA)

    hH2O = -241826 ...
        + 143.05.*(tempF-298.15) ...
        - 46.432.*((tempF.^1.25)-(298.15.^1.25)) ...
        + 8.2751./1.5.*((tempF.^1.5)-(298.15.^1.5)) ...
        - 0.0184945.*((tempF.^2)-(298.15.^2));

    hH2 = 56.505.*(tempF-298.15) ...
        - 88890.4.*((tempF.^0.25)-(298.15.^0.25)) ...
        + 116500.*log(tempF./298.15) ...
        + 1121400.*((tempF.^(-0.5))-(298.15.^(-0.5)));

    hO2 = 37.432.*(tempA-298.15) ...
        + 8.0408e-6.*((tempA.^2.5)-(298.15.^2.5)) ...
        + 357140.*((tempA.^(-0.5))-(298.15.^(-0.5))) ...
        - 2368800.*((tempA.^(-1.0))-(298.15.^(-1.0)));

    hN2 = 39.060.*(tempA-298.15) ...
        + 1025580.*((tempA.^(-0.5))-(298.15.^(-0.5))) ...
        - 10727e3.*((tempA.^(-1.0))-(298.15.^(-1.0))) ...
        + 4.102e8.*((tempA.^(-2.0))-(298.15.^(-2.0)));

    h = [hH2 0 hO2 hN2 hH2O];

end
%% ============================================================
% GRAPHIQUES POUR LE RAPPORT
%% ============================================================

Ptarget_vec = ResultsTable.Ptarget_W;

%% 1) Rendement optimal en fonction de la puissance
figure;
plot(Ptarget_vec, ResultsTable.eta_percent, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Puissance cible P_{target} [W]');
ylabel('Rendement électrique \eta_{el} [%]');
title('Rendement optimal en fonction de la puissance cible');

%% 2) Courant optimal en fonction de la puissance
figure;
plot(Ptarget_vec, ResultsTable.I_A, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Puissance cible P_{target} [W]');
ylabel('Courant optimal I [A]');
title('Courant optimal en fonction de la puissance cible');

%% 3) Débits optimaux H2 et air
figure;
plot(Ptarget_vec, ResultsTable.qH2_Lmin, '-o', 'LineWidth', 1.5);
hold on;
plot(Ptarget_vec, ResultsTable.qAir_Lmin, '-s', 'LineWidth', 1.5);
grid on;
xlabel('Puissance cible P_{target} [W]');
ylabel('Débit [L/min]');
title('Débits optimaux en fonction de la puissance cible');
legend('q_{H2}', 'q_{air}', 'Location', 'best');

%% 4) Température finale du stack
figure;
plot(Ptarget_vec, ResultsTable.T_C, '-o', 'LineWidth', 1.5);
hold on;
yline(750, '--', 'Limite basse 750°C');
yline(800, '--', 'Limite haute 800°C');
grid on;
xlabel('Puissance cible P_{target} [W]');
ylabel('Température finale [°C]');
title('Température finale du stack en fonction de la puissance cible');
legend('T finale', 'Limite basse', 'Limite haute', 'Location', 'best');

%% 5) Tension cellule
figure;
plot(Ptarget_vec, ResultsTable.Ucell_V, '-o', 'LineWidth', 1.5);
hold on;
yline(0.7, '--', 'Limite U_{cell} = 0.7 V');
grid on;
xlabel('Puissance cible P_{target} [W]');
ylabel('Tension cellule U_{cell} [V]');
title('Tension cellule en fonction de la puissance cible');
legend('U_{cell}', 'Limite basse', 'Location', 'best');


%% 6) Excès d’air lambda
figure;
plot(Ptarget_vec, ResultsTable.lambda_air, '-o', 'LineWidth', 1.5);
hold on;
yline(4, '--', '\lambda_{air} min = 4');
yline(15, '--', '\lambda_{air} max = 15');
grid on;
xlabel('Puissance cible P_{target} [W]');
ylabel('\lambda_{air} [-]');
title('Taux d''excès d''air en fonction de la puissance cible');
legend('\lambda_{air}', 'Limite basse', 'Limite haute', 'Location', 'best');