function [f, beta_F, T, I, TSFC, ETA_p, ETA_th, ETA_o] = TRBFN_AS(M_inf, Ta, pa, ma, BPR, eps_PD, eta_F, beta_C, eta_C, T_max, eps_CB, eta_HPT, eta_LPT, eta_N, T_AB, eps_AB)
%questa funzione permette, dato il set di input richiesti, di ricavare
%rendimenti e prestazioni utili di un propulsore di tipo turbofan
%tutte le unità di misura sono espresse secondo SI

% M_inf   è il mach di volo
% Ta      è la temperatura ambientale
% pa      è la pressione ambientale
% ma      è la portata trattata dal gruppo turbogas
% BPR     è il rapporto di bypass
% eps_PD  è il rendimento pneumatico della presa dinamica (supersonica)
% eta_F   è il rendimento adiabatico del Fan
% beta_C  è il rapporto di compressore del Compressore
% eta_C   è il rendimento adiabatico del Compressore
% T_max   è la temperatura massima raggiungibile in camera di combustione
% eps_CB  è il rendimento pneumatico del combustore
% eta_HPT è il rendimento adiabatico della turbina HPT
% eta_LPT è il rendimento adiabatico della turbina LPT
% eta_N   è il rendimento adiabatico dell'ugello primario
% T_AB    è la temperatura massima raggiunta nel postbruciatore (se c'è)
% eps_AB  è il rendimento pneumatico del postbruciatore (se presente)


R=287;   %[J/kg*K] costante dei gas caratteristica dell'aria

cp=1004; %[J/kg*K] calore specifico a pressione costante per aria secca

gam=1.4;  %rapporto tra cp e cv di un gas reale biatomico (circa aria)

de=(gam-1)/2;  %delta, grandezza di comodo

Hi=43000000; %[J/kg] potere calorifico inferiore

g=(gam-1)/gam;     %costante di comodo

v0=M_inf*sqrt(gam*R*Ta);  %[m/s] velocità di volo

if M_inf <= 1   %se volo subsonico

    M_PD=M_inf;

else             %se volo supersonico (onda normale)

    %M_PD=...

    error('la function tratta solo il volo subsonico!');

end

%PRESA DINAMICA
T01=Ta*(1+de*M_inf^2);  %temperatura totale in ingresso

p01=pa*(1+de*M_inf^2)^(1/g); %pressione totale imbocco presa dinamica

T02=T01;

p02=eps_PD*p01;   %pressione totale all'ingresso del Fan

%BISOGNA ORA RICAVARE IL RAPPORTO DI COMPRESSIONE DEL FAN, IL RAPPORTO DI MISCELA f E LE ALTRE GRANDEZZE
%CHE PERMETTONO DI REALIZZARE IL CICLO: L'EQUAZIONE NON SI RISOLVE IN FORMA CHIUSA, USO METODO ITERATIVO:

T05=T_max;   %questo dato è fissato

F=[0.005:1e-6:0.067];

M=length(F);

C=zeros(M, 1);

C1=eps_CB*beta_C;           %rapporto di pressioni costante

c=1+((beta_C^g) -1)/eta_C;  %costante di comodo

k=1;

T3=zeros(M, 1);

T4=zeros(M, 1);

T6=zeros(M, 1);

T7=zeros(M, 1);

p3=zeros(M, 1);

p4=zeros(M, 1);

p5=zeros(M, 1);

b_F=zeros(M, 1);

b_HPT=zeros(M, 1);

b_LPT=zeros(M, 1);

while k <= M

    T4(k)=(1+F(k))*T05-F(k)*Hi/cp;  %dal bilancio in camera di combustione

    T3(k)=T4(k)/c;   %dalla relazione adiabatica considerando il rendimento del compressore

    b_F(k)=(1+eta_F*(-1+(T3(k)/T02)))^(1/g);   %dalla relazione adiabatica considerando il rendimento del Fan
    
    p3(k)=b_F(k)*p02;    %dal rapporto di compressione
    
    p4(k)=beta_C*p3(k);    %dal rapporto di compressione
    
    p5(k)=eps_CB*p4(k);    %dall'efficienza pneumatica del combustore
    
    T6(k)=T05-(T4(k)-T3(k))/(1+F(k)); %per l'equilibrio dei lavori all'asse del compressore
    
    T7(k)=T6(k)-(T3(k)-T02)*(1+BPR)/(1+F(k)); %per l'equilibrio dei lavori all'asse del Fan
    
    b_HPT(k)=(1+((T6(k)/T05)-1)/eta_HPT)^(-1/g);
    
    b_LPT(k)=(1+((T7(k)/T6(k))-1)/eta_LPT)^(-1/g);
    
    C2=b_LPT(k)*b_HPT(k);

    C(k)=abs(C2-C1);

    k=k+1;

end

J=find(C==min(C)); %posizione, nel vettore C, in cui è minore la discrepanza tra il prodotto b_LPT*b_HPT e il prodotto eps_CB*beta_C

%W=find(b_F==real(min(b_F)));

%b_F(W)

f1=F(J);

T03=T3(J);

T04=T4(J);

T06=T6(J);

T07=T7(J);

beta_F=b_F(J);

beta_HPT=b_HPT(J);

beta_LPT=b_LPT(J);

p03=p3(J);

p04=p4(J);

p05=p5(J);

p06=p05/beta_HPT;
    
p07=p06/beta_LPT;
    
T07min=T6*(1-eta_LPT+eta_LPT*(p01/p06)^(g));  %temperatura totale minima a cui posso
        % espandere facendo in modo di avere una velocità di efflusso uguale
        % alla velocità di volo (considerando ugello di scarico primario perfetto)
    
BPRmax=-1+( (1+f1)*(T05-T07min)+T03-T04 )/( T03-T02 );    %massimo bpr utilizzabile
    
if BPR >= BPRmax
    
    warning('la velocità di efflusso dello scarico principale è inferiore alla velocità di volo');
    
end

%f_1=(cp*(T05-T04)) / (Hi-cp*T05);

L_F=ma*(1+BPR)*(T03-T02);

L_LPT=ma*(1+f1)*(T06-T07);

if abs(L_F-L_LPT)/L_F > 0.001

    warning('ci sono errori di roundoff nel ciclo!');

end

L_C=ma*(T04-T03);

L_HPT=ma*(1+f1)*(T05-T06);

if abs(L_C-L_HPT)/L_C > 0.001

    warning('ci sono errori di roundoff nel ciclo!');

end


%MIXER
T08=(T03*BPR+T07*(1+f1))/(1+f1+BPR);

p08=p03; %consideriamo nulle le perdite di carico nel mixer, per semplicità

if nargin > 14  %c'è il postbruciatore

    T09=T_AB;

    f2=cp*(1+f1+BPR)*(T09-T08)/(Hi-cp*T09);

    f=f1+f2;  %carburante totale utulizzato

    if f > 0.067  %saturazione di carburante

        warning('la combustione completa avviene con eccesso di carburante');

    end

    p09=eps_AB*p03; %p08=p03=p07

    t=[Ta, T02, T03, T04, T05, T06, T07, T08, T09, T09*(pa/p09)^g]; %vettore temperature per avviare la funzione Piano_Ts

    p=[pa, p02, p03, p04, p05, p06, p07, p08, p09, pa]; %vettore pressioni per avviare la funzione Piano_Ts

else   %no postbruciatore
    
    T09=T08;

    f=f1;

    p09=p03;

    t=[Ta, T02, T03, T04, T05, T06, T07, T08, T09*(pa/p09)^g]; %vettore temperature per avviare la funzione Piano_Ts

    p=[pa, p02, p03, p04, p05, p06, p07, p08, pa]; %vettore pressioni per avviare la funzione Piano_Ts

end

%UGELLO
T10id=T09*(pa/p09)^g;

ve=sqrt(2*cp*eta_N*(T09-T10id));

%SPINTA E SPINTA SPECIFICA 
T=ma*((1+f+BPR)*ve-(1+BPR)*v0);  %[N] spinta

I=T/(ma*(1+BPR)); %[m/s] spinta specifica

%TSFC
TSFC=f/((1+f+BPR)*ve-(1+BPR)*v0); %[kg/N*s] tsfc

%RENDIMENTI

ETA_p=2*v0*(((1+f+BPR)*ve-(1+BPR)*v0))/((1+f+BPR)*ve^2-(1+BPR)*v0^2); %rendimento propulsivo

ETA_th=((1+f+BPR)*ve^2-(1+BPR)*v0^2)/(2*f*Hi); %rendimento termodinamico

ETA_o=ETA_th*ETA_p;   %rendimento globale

Piano_Ts(t, p)

%T01

%T02

%p01

%p02

end

