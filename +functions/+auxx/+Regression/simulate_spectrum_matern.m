% Parameters
A = 1;
B = 0.1;
C = 2.5;
f0 = 10;
fs = 1000;        % sampling frequency (Hz)
T = 20;           % duration in seconds
dt = 1/fs;
t = 0:dt:(T-dt);
N = length(t);
f = linspace(-fs/2, fs/2, N);

% Build theoretical spectrum S(f) (piecewise mirror peaks)
S_theoretical = zeros(1, N);
for k = 1:N
    if f(k) >= 0
        S_theoretical(k) = A^2 / (1 + B*(f(k) - f0)^2)^C;
    else
        S_theoretical(k) = A^2 / (1 + B*(f(k) + f0)^2)^C;
    end
end

% Impulse response h_GM(t) (causal, truncated)
t_h = 0:dt:2;   % truncate to finite causal kernel
Nh = length(t_h);
h_GM = zeros(1, Nh);
parfor n = 1:Nh
    tau = t_h(n);
    integrand_pos = @(u) exp(2*pi*1i*u*tau) ./ (1 + B*u.^2).^(C/2);
    integrand_neg = @(u) exp(2*pi*1i*u*tau) ./ (1 + B*u.^2).^(C/2);
    % Integrate piecewise
    I_pos = integral(integrand_pos, -f0, Inf, 'ArrayValued', true);
    I_neg = integral(integrand_neg, -Inf, f0, 'ArrayValued', true);
    h_GM(n) = 2*A * real(exp(2*pi*1i*f0*tau) * I_pos + exp(-2*pi*1i*f0*tau) * I_neg);
end

% Normalize kernel
h_GM = h_GM / norm(h_GM);

% Generate white noise and filter it
n = randn(1, N);
u = conv(n, h_GM, 'same');

% Estimate empirical power spectrum of u(t)
[Sp_empirical, f_emp] = periodogram(u, [], N, fs, 'centered');

% Plot results
figure;
plot(f, 10*log10(S_theoretical), 'LineWidth', 2); hold on;
xlim([-20 20])
plot(f_emp, 10*log10(Sp_empirical), 'r--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend('Theoretical Spectrum', 'Empirical Spectrum');
title('Generalized Matérn Process: Simulated vs. Theoretical Spectrum');
grid on;
