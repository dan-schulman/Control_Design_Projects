clear; close all; clc;

load('StateSpaceExamples.mat')

%% Discretize the model first
Ts = 1/10000;
Gxx_ss = ss(Gxx);
Gxxd = c2d(Gxx_ss, Ts);

%% Coarse search to see if something even exists
test_crossover = (50:50:5000)./(2*pi);
K_cells = cell(length(test_crossover), 1);

% Storage for the intermediate weight matrices
W1_cells = cell(length(test_crossover), 1);
W2_cells = W1_cells;
W3_cells = W1_cells;

for i = 1:length(test_crossover)
    fprintf('%0.0f\n', i)
    
    W1 = makeweight(db2mag(120),  test_crossover(i), db2mag(-6), Ts);
    W2 = makeweight(db2mag(-20), test_crossover(i), db2mag(20), Ts);
    W3 = W1^-1;
    
    W1_cells{i} = W1;
    W2_cells{i} = W2;
    W3_cells{i} = W3;
    
    [K, CL, gamma, info] = mixsyn(Gxxd, W1, W2, W3);
    K_cells{i} = K;
    gammas_d(i) = gamma;
    
    % Calculate the wanted metrics
    L = Gxxd*K;
    I = eye(size(L));
    S = feedback(I, L);
    
    test_omega = logspace(-3, 4, 20001);
    [mag, phase, w] = bode(S, test_omega);
    
    bandwidth(i) = w(find(mag2db(mag) < -3, 1, 'last'))/(2*pi);
    max_peak(i) = hinfnorm(S);
    
    [Gm, Pm, Wcg, Wcp] = margin(L);
    
    Gms(i) = Gm;
    Pms(i) = Pm;
    Wcgs(i) = Wcg;
    Wcps(i) = Wcp;
    
    step_response = step(1-S, 0:Ts:2);
    ss_error(i) = 1 - step_response(end);
end

% Find optimal controller to meet bounds
bandwidth_goal = 45;
peak_goal = 6.75;            % dB
gain_goal = 10;              % dB
min_goal  = 6;               % dB
phase_goal = 30;             % deg

bandwidth_good = bandwidth > bandwidth_goal;
peak_good = mag2db(max_peak) < peak_goal;
gain_good = mag2db(Gms) > gain_goal;
phase_good = Pms > phase_goal;

all_good = bandwidth_good & peak_good & gain_good & phase_good;

% Get best bandwidth index
good_bandwidths = bandwidth(all_good);
ind = find(good_bandwidths == max(good_bandwidths), 1);

actual_ind = find(all_good == 1, 1, 'first') + ind - 1;

K = K_cells{actual_ind};

L = Gxxd*K;
I = eye(size(L));
S = feedback(I, L);
T = feedback(L, I);

test_omega = logspace(-2, 4, 20001);
[S_mag, S_phase] = bode(S, test_omega);
[T_mag, T_phase] = bode(T, test_omega);

time = 0:Ts:0.25;
step_response = step(1-S, time);

% on longer time frame, calculate error
time = 0:Ts:1;
step_response = step(1-S, time);
error = step_response(end) - 1;

% Print the sensitivity
S_mag = mag2db(squeeze(S_mag));
S_phase = squeeze(S_phase);

T_mag = mag2db(squeeze(T_mag));
T_phase = squeeze(T_phase);

W1 = makeweight(db2mag(120),  test_crossover(actual_ind), db2mag(-6), Ts);
W3 = W1^-1;

% Actual wanted shape is the inverse of the shaping function
W1 = W1^-1;
W3 = W3^-1;

[W1_mag, W1_phase] = bode(W1, test_omega);
[W3_mag, W3_phase] = bode(W3, test_omega);

W1_mag = mag2db(squeeze(W1_mag));
W1_phase = squeeze(W1_phase);
W3_mag = mag2db(squeeze(W3_mag));
W3_phase = squeeze(W3_phase);

%% Metrics
crossover = test_crossover(actual_ind)
gain_margin = mag2db(Gms(actual_ind))
phase_margin = Pms(actual_ind)
final_bandwidth = bandwidth(actual_ind)
S_hinf = mag2db(hinfnorm(S))

controller = K
system = Gxxd*K

%% Noise
Delta = ultidyn('Delta', [1, 1]);       %SISO system
Sunc = feedback(1, Gxxd*(1+Delta)*K);



%% Plots and Figures
figure;
subplot(2, 1, 1)
hold on
plot(test_omega./(2*pi), S_mag, 'LineWidth', 1.5)
%plot(test_omega./(2*pi), W1_mag, 'LineWidth', 1.5)
ylabel('|S| (dB)')
set(gca, 'XScale', 'log')
grid on
axis([10^-2 10^3 -inf 6])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
subplot(2, 1, 2)
hold on
plot(test_omega./(2*pi), S_phase, 'LineWidth', 1.5)
%plot(test_omega./(2*pi), W1_phase, 'LineWidth', 1.5)
grid on
set(gca, 'XScale', 'log')
%legend('System', 'W1^{-1}', 'Location', 'SouthWest')
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
axis([10^-2 10^3 -inf inf])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print -dpdf -fillpage SensitivityMixedBode

figure;
subplot(2, 1, 1)
hold on
plot(test_omega./(2*pi), T_mag, 'LineWidth', 1.5)
%plot(test_omega./(2*pi), W3_mag, 'LineWidth', 1.5)
ylabel('|T| (dB)')
set(gca, 'XScale', 'log')
axis([10^-2 10^3 -inf 1])
grid on
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

subplot(2, 1, 2)
hold on
plot(test_omega./(2*pi), T_phase, 'LineWidth', 1.5)
%plot(test_omega./(2*pi), W3_phase, 'LineWidth', 1.5)
grid on
set(gca, 'XScale', 'log')
%legend('System', 'W3^{-1}', 'Location', 'SouthWest')
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
axis([10^-2 10^3 -inf inf])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print -dpdf -fillpage CompSensitivityMixedBode

figure;
plot(time, step_response, '-', 'LineWidth', 1.5)
xlabel('Time (sec)')
ylabel('Position')
grid on
box off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print -dpdf -fillpage StepResponseMixedController

figure;
hold on
plot(test_crossover, bandwidth, '.', 'MarkerSize', 15)
plot(test_crossover, bandwidth_goal.*ones(size(Gms)), ':', 'LineWidth', 1.5)
axis([0 inf 0 inf])
ylabel('Bandwidth (Hz)')
xlabel('Crossover Frequency (rad)')
grid on
box off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print -dpdf -fillpage BandwidthMixed

figure;
hold on
plot(test_crossover, mag2db(max_peak), '.', 'MarkerSize', 15)
plot(test_crossover, peak_goal.*ones(size(Gms)), ':', 'LineWidth', 1.5)
axis([0 inf 0 inf])
ylabel('|S|_\infty (dB)')
xlabel('Crossover Frequency (rad)')
axis([-inf inf -inf 7.5])
grid on
box off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print -dpdf -fillpage SensitivityMixed

figure;
hold on
plot(test_crossover, mag2db(Gms), '.', 'MarkerSize', 15);
plot(test_crossover, gain_goal.*ones(size(Gms)), ':', 'LineWidth', 1.5)
axis([0 inf 0 inf])
ylabel('Gain Margin (dB)')
xlabel('Crossover Frequency (rad)')
grid on
box off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print -dpdf -fillpage GainMixed

figure;
hold on
plot(test_crossover, Pms, '.', 'MarkerSize', 15);
plot(test_crossover, phase_goal.*ones(size(Pms)), ':', 'LineWidth', 1.5)
axis([0 inf 0 inf])
ylabel('Phase Margin (deg)')
xlabel('Crossover Frequency (rad)')
grid on
box off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print -dpdf -fillpage PhaseMixed

figure;
hold on
for i = 1:length(W1_cells)
    [mag_W1, phase_W1] = bode(W1_cells{i}, test_omega);
    plot(test_omega, squeeze(mag2db(mag_W1)), 'Color', [72 157 255]./255, 'LineWidth', 1.5);
end
set(gca, 'XScale', 'log')
ylabel('W_1 (dB)')
xlabel('Frequency (rad/s)')
grid on
box off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print -dpdf -fillpage W1

figure;
hold on
for i = 1:length(W1_cells)
    [mag_W2, phase_W2] = bode(W2_cells{i}, test_omega);
    plot(test_omega, squeeze(mag2db(mag_W2)), 'Color', [255 170 72]./255, 'LineWidth', 1.5);
end
set(gca, 'XScale', 'log')
ylabel('W_2 (dB)')
xlabel('Frequency (rad/s)')
grid on
box off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print -dpdf -fillpage W2

figure;
hold on
for i = 1:length(W1_cells)
    [mag_W3, phase_W3] = bode(W3_cells{i}, test_omega);
    plot(test_omega, squeeze(mag2db(mag_W3)), 'Color', [170 72 255]./255, 'LineWidth', 1.5);
end
set(gca, 'XScale', 'log')
ylabel('W_3 (dB)')
xlabel('Frequency (rad/s)')
grid on
box off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print -dpdf -fillpage W3