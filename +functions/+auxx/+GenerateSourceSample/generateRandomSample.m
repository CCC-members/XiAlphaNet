function x = generateRandomSample(x0, var)

import functions.*
import functions.auxx.*
import functions.auxx.GenerateSourceSample.*
import functions.auxx.ModelVectorization.*

% %Nv = parameters.Dimensions.Nv;
% t = activation_threshold(S,T);
% %[me,se,ma,sa] = activation_threshold(parameters)
% [mf,stdf,mr,~] = init_freq(Ne,S,freq);
% % Mean values
% meanE1 = (4/5)*t;     % Target mean for e(:,1)
% meanE2 = 0.001; % Target mean for e(:,2)
% meanE3 = 2.5;   % Target mean for e(:,3)
% 
% meanA1 = (1/5)*t;     % Target mean for a(:,1), assuming you want the range 0.5 to 2
% meanA2 = 0.008; % Target mean for a(:,2)
% meanA3 = 15;    % Target mean for a(:,3), assuming you want the range 10 to 20
% meanA4 = mf;     % Midpoint of the range [5, 13] for a(:,4)
% 
% % Adjusting gamma distribution parameters for e and a (except e(:,2) and a(:,2))
% gammaShapeE1 = meanE1^2 / var;
% gammaRateE1 = meanE1 / var;
% 
% gammaShapeE3 = meanE3^2 / var;
% gammaRateE3 = meanE3 / var;
% 
% gammaShapeA1 = meanA1^2 / var;
% gammaRateA1 = meanA1 / var;
% 
% gammaShapeA3 = meanA3^2 / var;
% gammaRateA3 = meanA3 / var;
% 
% % Generate sigma^2 using a gamma distribution
% sigma2 = gamrnd(10, 1);  % Shape=10, Rate=1, Mean=10
% 
% % Generate matrix e (Nr x 3)
% e1 = abs(normrnd(meanE1, sqrt(var)/10, Nv, 1));  % e(:,1)
% e2 = abs(normrnd(meanE2, sqrt(var)/1000, Nv, 1));       % e(:,2) using normal distribution
% e3 = gamrnd(gammaShapeE3, 1/gammaRateE3, Nv, 1);  % e(:,3)
% e = [e1, e2, e3];
% 
% % Generate matrix a (Nr x 4)
% a1 = abs(normrnd(meanA1, sqrt(var)/1000, Nv, 1));  % a(:,1)
% a2 = abs(normrnd(meanA2, sqrt(var)/1000, Nv, 1));       % a(:,2) using normal distribution
% a3 = gamrnd(gammaShapeA3, 1/gammaRateA3, Nv, 1);  % a(:,3)
% a4 = normrnd(meanA4, sqrt(stdf), Nv, 1);          % a(:,4)
% a4 = max(min(a4, 13), 5);  %// Clipping to [7, 13]
% 
% a = [a1, a2, a3, a4];
%%
% --- Step 1: Decompose ---
     % --- Step 1: Decompose x0 into components ---
    [e, a, s2] = x2v(x0);  % Decompose x0 into components

% --- Step 2: Add Gaussian noise to e ---
e_noisy = e;
for col = 1:size(e, 2)
    std_col = std(e(:, col));
    noise = var * std_col * randn(size(e, 1), 1);  % Gaussian noise
    e_col_noisy = e(:, col) + noise;
    e_noisy(:, col) = max(e_col_noisy, 0);  % Ensure positivity
end

% --- Step 3: Add Gaussian noise to a (except a(:,4)) ---
a_noisy = a;
for col = 1:size(a, 2)
    std_col = std(a(:, col));
    noise = var * std_col * randn(size(a, 1), 1);  % Gaussian noise

    if col == 4
        % Preserve [7, 13] constraint via Gaussian + clip
        a_col_noisy = min(max(a(:, col) + noise, 7), 13);
    else
        a_col_noisy = max(a(:, col) + noise, 0);  % Ensure positivity
    end

    a_noisy(:, col) = a_col_noisy;
end

% --- Step 4: Sample s2 from Gaussian distribution and ensure positivity ---
s2_sample = max(normrnd(10, var), 0);  % Mean = 10, ensure positive

% --- Step 5: Reconstruct noisy parameter vector ---
x = v2x(e_noisy, a_noisy, s2_sample);


% Reshape matrices e and a into vectors and concatenate with sigma^2
%x = v2x(e, a, sigma2);
end
