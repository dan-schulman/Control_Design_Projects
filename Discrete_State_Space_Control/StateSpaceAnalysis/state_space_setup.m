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
    
    figure;
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

figure;
plot(1:max_dim, norm_diff_mag, 'o')
set(gca, 'YScale', 'log')

figure;
plot(actual_dim, norm_diff_mag, 'o')
set(gca, 'YScale', 'log')

sys1 = ss_cells{8};
sys2 = ss_cells{11};

save('StateSpaceExamples.mat', 'sys1', 'sys2', 'Gxx', 'ss_cells')