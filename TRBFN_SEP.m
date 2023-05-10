function [f, BPRmax, T, I, TSFC, ETA_p, ETA_th, ETA_o] = TRBFN_SEP(M_inf, Ta, pa, ma, BPR, eps_PD, beta_F, eta_F, beta_C, eta_C, T_max, eps_CB, eta_HPT, eta_LPT, eta_N1, eta_N2)
%questa funzione permette, dato il set di input richiesti, di ricavare
%rendimenti e prestazioni utili di un propulsore di tipo turbofan
%tutte le unità di misura sono espresse secondo SI

% M_inf   è il mach di volo
% Ta      è la temperatura ambientale
% pa      è la pressione ambientale
% ma      è la portata trattata dal gruppo turbogas
% BPR     è il rapporto di bypass
% eps_PD  è il rendimento pneumatico della presa dinamica (supersonica)
% beta_F  è il rapporto di compressore del Fan
% eta_F   è il rendimento adiabatico del Fan
% beta_C  è il rapporto di compressore del Compressore
% eta_C   è il rendimento adiabatico del Compressore
% T_max   è la temperatura massima raggiungibile in camera di combustione
% eps_CB  è il rendimento pneumatico del combustore
% eta_HPT è il rendimento adiabatico della turbina HPT
% eta_LPT è il rendimento adiabatico della turbina LPT
% eta_N1  è il rendimento adiabatico dell'ugello primario
% eta_N2  è il rendimento adiabatico dell'ugello secondario

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

%FAN
T03id=T02*beta_F^((gam-1)/gam);   %temperatura ideale in uscita dal Fan e ingresso compressore
    
T03=T02+(T03id-T02)/eta_F;        %temperatura reale in uscita dal Fan e ingresso compressore
    
p03=beta_F*p02;   %pressione in uscita dal Fan e ingresso compressore
    
%L_F=ma*(1+BPR)*cp*(T03-T02); %potenza assorbita dal Fan
    
%COMPRESSORE
T04id=T03*beta_C^((gam-1)/gam);   %temperatura ideale in uscita dal compressore e ingresso combustore
    
T04=T03+(T04id-T03)/eta_C;        %temperatura reale in uscita dal compressore e ingresso combustore
    
p04=beta_C*p03;   %pressione in uscita dal compressore e ingresso combustore
    
L_C=ma*cp*(T04-T03); %potenza assorbita dal compressore
    
%COMBUSTORE:
T05=T_max;          %temperatura in uscita dal combustore
    
f=(cp*(T05-T04)) / (Hi-cp*T05);  %rapporto tra mf e ma;
    
if f >= 0.067
    
    warning('VIENE IMMESSO TROPPO CARBURANTE');
    
end
    
p05=eps_CB*p04;
    
%TURBINA HPT
T06=T05-(T04-T03)/(1+f); %temperatura in uscita da HPT: il suo lavoro muove il compressore
    
beta_HPT=(T05/(T05-(T05-T06)/eta_HPT))^(gam/(gam-1)); %rapporto di compressione della HPT
    
p06=p05/beta_HPT;  %pressione in uscita da HPT

%TURBINA LPT

T07min=T06*(1-eta_LPT+eta_LPT*(p01/p06)^((gam-1)/gam));  %temperatura totale minima a cui posso
% espandere facendo in modo di avere una velocità di efflusso uguale
% alla velocità di volo (considerando ugello di scarico primario perfetto)

BPRmax=-1+( (1+f)*(T05-T07min)+T03-T04 )/( T03-T02 );    %massimo bpr utilizzabile

if BPR >= BPRmax

    warning('la velocità di efflusso dello scarico principale è inferiore alla velocità di volo');

end

T07=T06- ( (1+BPR)*(T03-T02) )/( 1+f ); %temperatura in uscita da LPT: il suo lavoro muove il Fan
    
p07=p06*( 1-(T06-T07)/(T06*eta_LPT) )^( gam/(gam-1) );  %pressione in uscita da LPT

beta_LPT=p06/p07;    %rapporto di compressione della LPT

%UGELLO PRIMARIO
T8id=T07*(pa/p07)^((gam-1)/gam);   %temperatura ideale in uscita dall'ugello primario

ve1=sqrt(2*cp*eta_N1*(T07-T8id));  %velocità in efflusso da ugello primario

%UGELLO SECONDARIO
T9id=T03*(pa/p03)^((gam-1)/gam);   %temperatura ideale in uscita dall'ugello secondario

ve2=sqrt(2*cp*eta_N2*(T03-T9id));  %velocità in efflusso da ugello secondario

%SPINTA E SPINTA SPECIFICA 
T=ma*( ((1+f)*ve1-v0)+BPR*(ve2-v0) );  %[N] spinta

I=T/(ma*(1+BPR)); %[m/s] spinta specifica

%TSFC
TSFC=f/ ( ((1+f)*ve1-v0)+BPR*(ve2-v0) ); %[kg/N*s] tsfc

%RENDIMENTI

ETA_p=(2*v0*( ((1+f)*ve1-v0)+BPR*(ve2-v0) )) / ( (1+f)*ve1^2-v0^2 + BPR*(ve2^2-v0^2) ); %rendimento propulsivo

ETA_th=( (1+f)*ve1^2-v0^2 + BPR*(ve2^2-v0^2) ) / (2*f*Hi); %rendimento termodinamico

ETA_o=ETA_th*ETA_p;   %rendimento globale

t=[Ta, T02, T03, T04, T05, T06, T07, T8id];

p=[pa, p02, p03, p04, p05, p06, p07, pa];

%Piano_Ts(t, p);

%T01

%T02

%p01

%p02

end



