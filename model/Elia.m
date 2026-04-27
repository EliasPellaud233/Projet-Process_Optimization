clear all;clc
u = [1.1; 40; 0; 1023.15];
% Inputs
% qH2_in  = u(1);   % [L/min]
% qAir_in = u(2);   % [L/min]
% I       = u(3);   % [A]
% T       = u(4);   % [K]

%50.0000   → Pel
%0.9650    → Ucell
%0.3281    → FU
%15.2727   → lambda_air
%25.2442   → eta_el


%% 50 W
Ptarget = 50;
fun = @(I) get_power(I, u) - Ptarget;
I50 = fzero(fun, [0.01 50]);

u50 = u;
u50(3) = I50;
y50 = SOFC_model_part1(u50);

%% 75 W
Ptarget = 75;
fun = @(I) get_power(I, u) - Ptarget;
I75 = fzero(fun, [0.01 50]);

u75 = u;
u75(3) = I75;
y75 = SOFC_model_part1(u75);

%% 100 W
Ptarget = 100;
fun = @(I) get_power(I, u) - Ptarget;
I100 = fzero(fun, [0.01 50]);

u100 = u;
u100(3) = I100;
y100 = SOFC_model_part1(u100);

%% affichage
disp('--- 50 W ---')
disp(I50) 
disp(y50)

disp('--- 75 W ---')
disp(I75)
disp(y75)

disp('--- 100 W ---')
disp(I100)
disp(y100)

%% fonction locale
function P = get_power(I, u)
    u(3) = I;
    y = SOFC_model_part1(u);
    P = y(1);
end



