clear; close all;

%% System definition
s = tf('s');
Ts = 0.0001; % 10 kHz
omega = logspace(0, log10(2*pi*10^3), 20001);

% gaurav's plant
Gyy= -8.38e6*(s^2+21.12*s+2.28e4)*(s^2+11.27*s+1.58e6)*...
(s^2+3.25*s+4.72e6)*(s^2+20.72*s+6.71e6)*(s-7500)/...
    ((s^2+5.77*s+1.164e4)*(s^2+24.43*s+3.14e4)*...
    (s^2+17.50*s+1.54e6)*(s^2+6.99*s+4.49e6)*...
    (s^2+14.76*s+7.04e6)*(s+7500));
% gaurav's compensator
Cyy=260*(s+100)^2/(s*(s+2000)*(s+3000));

CG_ss = ss(Cyy*Gyy); % critical to convert to ss before c2d
CG_ssd = c2d(CG_ss,Ts,'zoh');
T_ssd = feedback(CG_ssd,1);

% 8th order plant:
sysg = frd(Gyy, omega);
sys = fitfrd(sysg,8);
sys_d = c2d(sys,Ts);

% Gaurav's in SS representation
CG_ss = ss(Cyy*Gyy); % feedforward TF
CG_ssd = c2d(CG_ss,Ts,'zoh'); % discretize
T_ssd = feedback(CG_ssd,1);

% Testing for NMPZ
% bode(sys);
% title('Continuous plant')

% sys_d = c2d(sys,1e-5);
% figure
% bode(sys_d);
% title('Discrete plant, 100 Hz')

%% LQR
Q = eye(8).*[1 1 1 1 1 1 1 1]';
ind = 8;
%   Q(ind,ind) = 100;
Q(1,1) = 1000;
Q(2,2) = 750;
Q(3,3) = 1;
Q(4,4) = 1000;
Q(5,5) = 1000;
Q(6,6) = 1;
Q(7,7) = 1;
Q(8,8) = 1;
R = 200; 

%Continuous LQR
[A,B,C,D] = ssdata(sys);
[K,P] = lqr(A,B,Q,R);
G = 1/((C*(-A+B*K)^(-1)*B)+D); %precomp
TsfG = ss(A-B*K,B*G,C,D);

%Discretizing
[Ad,Bd,Cd,Dd] = ssdata(c2d(sys,Ts));
[Kd,Pd] = dlqr(Ad,Bd,Q,R);
Ad_sf = Ad-Bd*Kd;
Bd_sf = Bd;
Cd_sf = Cd; 
Dd_sf = Dd;
Gd = 1; %55.3166;%.01287;%.031; % precompensator
Tsf_d = ss(Ad_sf,Bd,Cd,Dd,Ts);% for sens
% n = rscale(c2d(sys,Ts),K)
Gd = dcgain(Tsf_d/Gd);

%sensitivity trial
Lsf_d=ss(Ad,Bd,Kd,Dd,Ts);
S_sf_d=1/(1+Lsf_d);
w = logspace(-2,5,200);
[mS_sf_d,pS_sf_d] = bode(1-Tsf_d/Gd,w);
mS_sf_d = squeeze(mS_sf_d);pS_sf = squeeze(pS_sf_d);

figure
bode(Tsf_d/Gd,w)
title('LQR bode plot')

[mag_lqr,ph_lqr] = bode(Tsf_d/Gd,w);
mag_lqr = squeeze(mag_lqr);ph_lqr = squeeze(ph_lqr);

figure
k=bodeplot(1-Tsf_d/Gd);
title('LQR Sensitivity');
setoptions(k,'FreqUnits','Hz','Xlim',[1 1000],'MagLowerLim',-40,'PhaseVisible','off');

sens_max = max(db(mS_sf_d))
[Gm,Pm] = margin(Tsf_d/Gd);
Gm
% finding phase margin
ind_ph = find(abs(db(mag_lqr)+0.1)<0.05,1);
Pm = 180+(ph_lqr(ind_ph(1))-720)
bw = bandwidth(Tsf_d/Gd)/(2*pi)

% control effort
control_eff = ss(Ad_sf,Bd,-K,0,Ts);

figure
step(control_eff)
title('LQR Control Effort')

