function [outputArg1,outputArg2] = Piano_Ts(T, p)
%tutte le unità di misura sono espresse secondo SI

%T, p vettori di temperature e pressioni
%N è la risoluzione desiderata (intesa come densità di puntini uniti nel
%grafico)



R=287;

cp=1004;

M=length(T);

s=zeros(M, 1);

ds=zeros(M, 1);

for k=1:M-1 %questo loop costruisce il vettore delle M entropie relative ai punti caratteristici del ciclo

    ds(k+1)=cp*log(T(k+1)/T(k))-R*log(p(k+1)/p(k));

    s(k+1)=s(k)+ds(k+1); %vettore delle M entropie nei punti caratteristici del ciclo

end

N=100; %qtà di punti per il plotting di ogni trasformazione

for k=1:M-1

    y=log(p(k+1)/p(k))/log(T(k+1)/T(k));  %esponente della politropica

    pk=linspace(p(k), p(k+1), N); %variazione lineare di pressione nella k_esima trasformazione

    tk(1)=T(k); %inizializzazione

    sk(1)=s(k); %inizializzazione

    if y==0  %se k_esima trasformazione isobara

        tk=linspace(T(k), T(k+1), N); %variazione lineare della temperatura

        for i=1:N-1

            sk(i+1)=sk(1)+cp*log(tk(i+1)/tk(1));  %i+1_esima entropia nella k_esima trasformazione

        end

    else  %se k_esima trasformazione non isobara

        for i=1:N-1

        tk(i+1)=tk(1)*(pk(i+1)/pk(1))^(1/y);  %i+1_esima temperatura nella k_esima trasformazione

        sk(i+1)=sk(1)+cp*log(tk(i+1)/tk(1))-R*log(pk(i+1)/pk(1)); %i+1_esima entropia nella k_esima trasformazione

        end

    end

    plot(sk, tk);

    hold on

end

xlabel('s [J/Kg*K]');
ylabel('T [K]');
axis equal;

for n=1:M   %disegno delle isobare

    S=linspace( s(n)-200, s(n)+200, 100);

    T_= @(S) T(n)*exp((S-s(n))/cp);

    plot(S, T_(S), "--k");

end

grid on

end

