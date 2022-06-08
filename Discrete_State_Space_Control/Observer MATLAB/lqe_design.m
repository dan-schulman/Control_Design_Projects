clear; close all; clc;

%% Transfer function from (3.1)

s = tf('s');

Ps = 1.28e10*(s^2+5.63*s+3.34e5)/((s+333.1)*(s^2+150.5*s+3.31e4)*...
     (s^2 + 12.43*s + 3.87e5));
Cs = 1.57e4*(s+141.5)*(s^2+159.5*s+5.01e4)/(s*(s+4000)*(s^2+6700*s+1.92e7));
Fs = tf([3.807, 3350], [1, 3350]);

Ss = 1/(1+Ps*Cs);
Ts = (Ps*Cs)/(1+Ps*Cs);

Gyy = (-8.38e6/(s^2+5.77*s+1.164e4))*((s^2 + 21.12*s+2.28e4)/(s^2+24.43*s+3.14e4)) ...
      *((s^2+11.27*s+1.58e6)/(s^2+17.50*s+1.54e6))*((s^2+3.25*s+4.72e6)/(s^2+6.99*s+4.49e6)) ...
      *((s^2+20.72*s+6.71e6)/(s^2+14.76*s+7.04e6))*(s-7500)/(s+7500);
Gxx = Gyy;
  
omega = logspace(0, log10(2*pi*10^3), 20001);

[mag, phase] = bode(Gyy, omega);

sysg = frd(Gyy, omega);

mag = squeeze(mag);
phase = squeeze(phase);


% Fit different dimensions of state space models
max_dim = 15;
ss_cells = cell(max_dim, 1);
for dim = 1:max_dim
    sys = fitfrd(sysg, dim);
    ss_cells{dim} = sys;
    [new_mag, new_phase] = bode(sys, omega);

    new_mag = squeeze(new_mag);
    new_phase = squeeze(new_phase);
    
    norm_diff_mag(dim) = norm(new_mag - mag, 2);
    actual_dim(dim) = length(sys.A);
    
    %figure;
    subplot(2, 1, 1)
    hold on
    plot(omega./(2*pi), 20*log10(mag), 'LineWidth', 1.5)
    plot(omega./(2*pi), 20*log10(new_mag), ':', 'LineWidth', 1.5)
    set(gca, 'XScale', 'log')
    box off
    grid on
    subplot(2, 1, 2)
    hold on
    plot(omega./(2*pi), unwrap(wrapTo180(squeeze(phase))), 'LineWidth', 1.5);
    plot(omega./(2*pi), unwrap(wrapTo180(squeeze(new_phase))), ':', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log')
    xlabel('Frequency (Hz)')
    legend('Original TF', sprintf('%0.0f DOF SS Fit', length(sys.A)))
    box off
    grid on
end
    
sys1 = ss_cells{8};
sys2 = ss_cells{11};
sys3 = ss_cells{4};

save('StateSpaceExamples.mat', 'sys1', 'sys2', 'Gxx', 'ss_cells')

[A,B,C,D] = ssdata(sys1);

%LQR Design

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


%LQE Design
W=eye(8);
W(1,1) = 200;
W(2,2) = 300;
W(3,3) = 1000;
W(4,4) = 1000;
W(5,5) = 100;
W(6,6) = 100;
W(7,7) = 150;
W(8,8) = 100;
W=W*10^(-3);
V=15e3;

G1=eye(8);

%Discretizing
sysd = c2d(sys1,1e-4);
[Ad,Bd,Cd,Dd] = ssdata(c2d(sys1,1e-4));
[Kd,Pd] = lqrd(A,B,Q,R,1e-4);
[Ld,Sd] = lqed(A,G1,C,W,V,1e-4);
%Gd = 1/((Cd*(-Ad+Bd*Kd)^(-1)*Bd)+Dd);
w = logspace(-2,5,200);
t = linspace(0,6);


%Discrete observer Compensator
Ad_obs = [Ad -Bd*Kd; Ld*Cd Ad-Bd*Kd-Ld*Cd]; 
Bd_obs = [Bd; Bd];
Cd_obs = [Cd zeros(1,8)]; 
Dd_obs = Dd;
Tref_obs_d = ss(Ad_obs,Bd_obs,Cd_obs,Dd_obs,1e-4);

Gd = dcgain(Tref_obs_d);

%Observer compensator
% state space of Cobs, the observer compensator
A_C = Ad-Bd*Kd-Ld*Cd;
B_C = Ld;
C_C = Kd;
D_C = Dd;
Cobs_d = ss(A_C,B_C,C_C,Dd,1e-4);


Lobs_d = Cobs_d*sysd;
Sobs_d = 1/(1+Lobs_d);


figure
k=bodeplot(Sobs_d);
setoptions(k,'FreqUnits','Hz','Xlim',[0.01 1000],'MagLowerLim',-40,'PhaseVisible','off');
title('')

figure
k=bodeplot(Tref_obs_d/Gd);
setoptions(k,'FreqUnits','Hz','Xlim',[0.01 1000],'MagLowerLim',-40,'PhaseVisible','off');
title('')

figure
k=bodeplot(Lobs_d);
setoptions(k,'FreqUnits','Hz','Xlim',[0.01 1000],'MagLowerLim',-40);
title('')

figure
step(Tref_obs_d/Gd);
xlim([0 1]);
title('')
