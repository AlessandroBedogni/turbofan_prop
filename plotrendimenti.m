%% Plot rendimenti
clc;clear all;close all;

Ta=270;
pa=0.7e5;
ma=42;
eps_PD=0.99;
eta_F=0.95;
beta_C=8;
eta_C=0.9;
T_max=1550;
eps_CB=0.95;
eta_HPT=0.95;
eta_LPT=0.95;
eta_N=0.99;
M_inf=[0.2,0.3,0.4,0.5,0.6,0.8,1];
BPR=[0:0.1:30];

eta_p=zeros(length(M_inf),length(BPR));
eta_th=zeros(length(M_inf),length(BPR));
eta_o=zeros(length(M_inf),length(BPR));
tsfc=zeros(length(M_inf),length(BPR));

for i=1:length(M_inf)
    for j=1:length(BPR)
        [f, beta_F, T, I, TSFC, ETA_p, ETA_th, ETA_o] = TRBFN_AS(M_inf(i), Ta, pa, ma, BPR(j), eps_PD, eta_F, beta_C, eta_C, T_max, eps_CB, eta_HPT, eta_LPT, eta_N);
        eta_p(i,j)=ETA_p;
        eta_th(i,j)=ETA_th;
        eta_o(i,j)=ETA_o;
        tsfc(i,j)=TSFC;
    end
end

%%
figure(1)
for i=1:length(M_inf)
    plot(BPR(:),eta_p(i,:),LineWidth=1.5)
    hold on;grid on,grid minor
end
xlabel("BPR");
ylabel("\eta_p")
legend("M=0.2","M=0.3","M=0.4","M=0.5","M=0.6","M=0.8","M=1")

figure(2)
for i=1:length(M_inf)
    plot(BPR(:),eta_th(i,:),LineWidth=1.5)
    hold on;grid on,grid minor
end
xlabel("BPR");
ylabel("\eta_t_h")
legend("M=0.2","M=0.3","M=0.4","M=0.5","M=0.6","M=0.8","M=1")

figure(3)
for i=1:length(M_inf)
    plot(BPR(:),eta_o(i,:),LineWidth=1.5)
    hold on;grid on,grid minor
end
xlabel("BPR");
ylabel("\eta_o")
legend("M=0.2","M=0.3","M=0.4","M=0.5","M=0.6","M=0.8","M=1")

figure(4)
for i=1:length(M_inf)
    plot(BPR(:),tsfc(i,:),LineWidth=1.5)
    hold on;grid on,grid minor
end
xlabel("BPR");
ylabel("TSFC [Kg/N]")
legend("M=0.2","M=0.3","M=0.4","M=0.5","M=0.6","M=0.8","M=1")
