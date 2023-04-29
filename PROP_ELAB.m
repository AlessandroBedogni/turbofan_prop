%SESSOPAZZO
clc;
clear;
%flussi associati
[f, beta_F, T, I, TSFC, ETA_p, ETA_th, ETA_o]=TRBFN_AS(0.7, 270, 0.7e5, 42, 1.2, 0.99, 0.95, 8, 0.9, 1550, 0.95, 0.95, 0.95, 0.99)



%% flussi separati
clc;
clear;


[f, BPRmax, T, I, TSFC, ETA_p, ETA_th, ETA_o]=TRBFN_SEP(0.7, 270, 1e5, 42, 1.2, 0.99, 3.3277, 0.95, 27, 0.9, 1550, 0.95, 0.95, 0.95, 0.99, 0.99)


%%
%if nargin==2

%    N=100;

%end

M=length(T);

%s=zeros(M, 1);

TT=zeros(M-1, N);

pp=zeros(M-1, N);

ss=zeros(M-1, N);

dss=zeros(M-1, N);

for k=1:M-1

    TT(k, :)=linspace(T(k), T(k+1), N);

    pp(k, :)=linspace(p(k), p(k+1), N);

    for i=1:M-1

        dss(k, i+1)=cp*log(TT(k, i+1)/TT(k, i))-R*log(pp(k, i+1)/pp(k, i));

    end

    %ds(k)=sum(dss(k, :));   %incremento di entropia dopo ogni trasformazione (riga)

end

for m=1:M-2

    for n=1:N

        ss(m+1, n)=sum(dss(m, :))+dss(m+1, n);

    end

    plot(ss(m, :), TT(m, :))

    hold on

    grid on

end

%%

for i=1:M-1

    fun= @(x) x^(s(i+1))-x^(s(i))+T(i+1)-T(i);

    dfun= @(y) s(i+1)*y^(s(i+1)-1)-s(i)*y^(s(i)-1);

    m(i)=biseznewton(-1e2, 1e2, 1e6, 1e3, 1e-4, fun, dfun);

    q(i)=T(i)-m(i)^(s(i));

    F= @(t) q(i)+m(i)^t;

    X=linspace(s(i), s(i+1), 100);

    plot(X, F(X));

    grid on

    hold on

end

%%

N=100;
t=linspace(T(5), T(6), N);
P=linspace(p(5), p(6), N);
S=zeros(N, 1);
dS=zeros(N, 1);

for k=1:N-1

    dS(k+1)=cp*log(t(k+1)/t(k))-R*log(P(k+1)/P(k));

    S(k+1)=S(k)+dS(k+1);

end


plot(S, t)

%%

%F= @(b) (1+ ( (b(k)^g) -1 )/eta_F)*T02*(BPR+c)+T02*(1+BPR)+(1+ ((T_max-T02*c*(1+ ( (b(k)^g) -1 )/eta_F)) / (E-T02*c*(1+ ( (b(k)^g) -1 )/eta_F))))*((T_max- (T02*(c-1)*(1+ ( (b(k)^g) -1 )/eta_F)) / (1+((T_max-T02*c*(1+ ( (b(k)^g) -1 )/eta_F)) / (E-T02*c*(1+ ( (b(k)^g) -1 )/eta_F)))))*(1-eta_LPT+eta_LPT*((b(k)*p02/p06(k))^g))-T_max);



%t=[Ta:1:T02:1:T03:1:T04:1:T05:1:T06:1:T07:1:T08:1:T09*(pa/p09)^g];

 %   p=[pa:10:p02:10:p03:10:p04:10:p05:10:p06:10:p07:10:p08:10:pa];

 %%
 %G= @(y) y.^2;

%for i=1:M-1

%    m(i)=(T(i+1)-T(i))/(G(s(i+1))-G(s(i)));

%    q(i)=T(i)-m(i)*G(s(i));

%    F= @(x) q(i)+m(i)*G(x);

%    X=linspace(s(i), s(i+1), 100);

%    plot(X, F(X));

%    grid on

%    hold on

%end
