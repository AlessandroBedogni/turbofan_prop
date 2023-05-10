%% PLOT POST BRUCIATORE CONSUMI

clc;clear all;close all;

Ta=270;
pa=0.7e5;
ma=42;
eps_PD=0.99;
eta_F=0.95;
beta_C=8;
eta_C=0.9;
T_max=1300;
eps_CB=0.95;
eta_HPT=0.95;
eta_LPT=0.95;
eta_N=0.99;
M_inf=[0.3,0.8];
BPR=[0:0.1:1];

eta_p=zeros(length(M_inf),length(BPR));
eta_th=zeros(length(M_inf),length(BPR));
eta_o=zeros(length(M_inf),length(BPR));
tsfc=zeros(length(M_inf),length(BPR));

eta_p_ab=zeros(length(M_inf),length(BPR));
eta_th_ab=zeros(length(M_inf),length(BPR));
eta_o_ab=zeros(length(M_inf),length(BPR));
tsfc_ab=zeros(length(M_inf),length(BPR));

for i=1:length(M_inf)
    for j=1:length(BPR)
        [f, beta_F, T, I, TSFC, ETA_p, ETA_th, ETA_o] = TRBFN_AS(M_inf(i), Ta, pa, ma, BPR(j), eps_PD, eta_F, beta_C, eta_C, T_max, eps_CB, eta_HPT, eta_LPT, eta_N);
        eta_p(i,j)=ETA_p;
        eta_th(i,j)=ETA_th;
        eta_o(i,j)=ETA_o;
        tsfc(i,j)=TSFC;
    end
end

for i=1:length(M_inf)
    for j=1:length(BPR)
        [f, beta_F, T, I, TSFC, ETA_p, ETA_th, ETA_o] = TRBFN_AS(M_inf(i), Ta, pa, ma, BPR(j), eps_PD, eta_F, beta_C, eta_C, T_max, eps_CB, eta_HPT, eta_LPT, eta_N,1550,0.99);
        eta_p_ab(i,j)=ETA_p;
        eta_th_ab(i,j)=ETA_th;
        eta_o_ab(i,j)=ETA_o;
        tsfc_ab(i,j)=TSFC;
    end
end

%%

figure(1)
for i=1:length(M_inf)
    plot(BPR(:),eta_p(i,:),LineWidth=1.5)
    hold on;grid on,grid minor
    plot(BPR(:),eta_p_ab(i,:),LineWidth=1.5)
end
xlabel("BPR");
ylabel("\eta_p")
legend("M=0.3","M_a_b=0.3","M=0.8","M_a_b=0.8")

figure(2)
for i=1:length(M_inf)
    plot(BPR(:),eta_th(i,:),LineWidth=1.5)
    hold on;grid on,grid minor
    plot(BPR(:),eta_th_ab(i,:),LineWidth=1.5)
end
xlabel("BPR");
ylabel("\eta_t_h")
legend("M=0.3","M_a_b=0.3","M=0.8","M_a_b=0.8")

figure(3)
for i=1:length(M_inf)
    plot(BPR(:),eta_o(i,:),LineWidth=1.5)
    hold on;grid on,grid minor
    plot(BPR(:),eta_o_ab(i,:),LineWidth=1.5)
end
xlabel("BPR");
ylabel("\eta_o")
legend("M=0.3","M_a_b=0.3","M=0.8","M_a_b=0.8")

figure(4)
for i=1:length(M_inf)
    plot(BPR(:),tsfc(i,:),LineWidth=1.5)
    hold on;grid on,grid minor
    plot(BPR(:),tsfc_ab(i,:),LineWidth=1.5)
end
xlabel("BPR");
ylabel("TSFC [Kg/N]")
legend("M=0.3","M_a_b=0.3","M=0.8","M_a_b=0.8")


