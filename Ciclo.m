
%%FLUSSI SEPARATI:

R=287;   %[J/kg*K] costante dei gas caratteristica dell'aria

cp=1004; %[J/kg*K] calore specifico a pressione costante per aria secca

gam=1.4;  %rapporto tra cp e cv di un gas reale biatomico (circa aria)

de=(gam-1)/2;  %delta, grandezza di comodo

Hi=43000000; %[J/kg] potere calorifico inferiore

v0=M_inf*sqrt(gam*R*Ta);  %[m/s] velocità di volo

if M_inf <= 1   %se volo subsonico

    M_PD=M_inf;

else             %se volo supersonico (onda normale)

    %M_PD=...

end

%PRESA DINAMICA
T01=Ta*(1+de*M_inf^2);  %temperatura totale in ingresso

p01=pa*(1+de*M_inf^2)^(gam/(gam-1)); %pressione totale imbocco presa dinamica

T02=T01;

p02=eps_PD*p01;   %pressione totale all'ingresso del Fan