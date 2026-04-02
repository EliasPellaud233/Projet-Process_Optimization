u = [1.1; 40; 17.5; 1023.15];
y = SOFC_model_part1(u)
Pel = y(1)
Ucell = y(2)
FU = y(3)
lair = y(4)
eff = y(5)