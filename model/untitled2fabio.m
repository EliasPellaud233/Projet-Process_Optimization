u = [1.1; 40; 17.5; 1023.15];
y = SOFC_model_part1(u)
%Pel = y(1)
%Ucell = y(2)
%FU = y(3)
lair = y(4)
%eff = y(5)


% Ex 2
%------------------------------------
%Leistung const (75W)
u = [0.00806; 600; 20; 1023.15];
y = SOFC_model_part1(u)
Pel = y(1)

%Point nominal
%u = [1.1; 40; 16.7; 1023.15]

%Piont 1
%u = [3.3; 40; 15; 1023.15];

%Piont 2
%u = [350; 600; 13.65; 1023.15];

%Point 3
%u = [28; 600; 13.9; 1023.15];

%Piont 4
%u = [0.806; 600; 19; 1023.15];
Q_h2= [350; 28; 3.3; 1.1; 0.806;];
I_2 = [13.65; 13.9; 15; 16.7; 19];

plot(I_2, Q_h2, 'O-')

scatter(I_2, Q_h2, 'filled')
xlabel('Cell current [A]')
ylabel('H2 flow rate [l/min]')
title('Constant electrical power 75 [W]')
grid on
