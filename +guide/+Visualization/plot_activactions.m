
function plot_activactions(x,parameters);
% Figure 2 Show the initial solution and the final solution

%% Plot initial solution
x0= generateRandomSample(parameters,0.1);
[e0,a0,s0] = x2v(x0);
% Create a single figure for all subplots
figure;
% Alpha Activation Plot
subplot(2, 3, 1); % 3 rows, 1 column, 1st subplot
J = a0(:,1); % Pass the current axes handle to esi_plot
esi_plot(gca,J);
title('Alpha Activation');

% Xi Activation Plot
subplot(2, 3, 2); % 3 rows, 1 column, 2nd subplot
J = e0(:,1);
esi_plot(gca,J); % Pass the current axes handle to esi_plot
title('Xi Activation');

% Alpha Peak Frequency
subplot(2, 3, 3); % 3 rows, 1 column, 3rd subplot
J = a0(:,4).*a0(:,1)/max(a0(:,1));
esi_plot(gca,J); % Pass the current axes handle to esi_plot
title('Alpha Peak Freq');

%% Plot estimate initial solution
[e,a,s2] = x2v(x);

% Alpha Activation Plot
subplot(2, 3, 4); % 3 rows, 1 column, 1st subplot
J = a(:,1); % Pass the current axes handle to esi_plot
esi_plot(gca,J);


% Xi Activation Plot
subplot(2, 3, 5); % 3 rows, 1 column, 2nd subplot
J = e(:,1);
esi_plot(gca,J); % Pass the current axes handle to esi_plot


% Alpha Peak Frequency
subplot(2, 3, 6); % 3 rows, 1 column, 3rd subplot
J = a(:,4);
esi_plot(gca,J); % Pass the current axes handle to esi_plot
end