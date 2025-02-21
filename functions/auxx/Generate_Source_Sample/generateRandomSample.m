function x = generateRandomSample(Ne,Nv,S,T,freq, var)
    %Nv = parameters.Dimensions.Nv;
    t = activation_threshold(S,T);
    %[me,se,ma,sa] = activation_threshold(parameters)
    [mf,stdf,mr,~] = init_freq(Ne,S,freq);
    % Mean values
    meanE1 = (4/5)*t;     % Target mean for e(:,1)
    meanE2 = 0.001; % Target mean for e(:,2)
    meanE3 = 2.5;   % Target mean for e(:,3)
    
    meanA1 = (1/5)*t;     % Target mean for a(:,1), assuming you want the range 0.5 to 2
    meanA2 = 0.008; % Target mean for a(:,2)
    meanA3 = 15;    % Target mean for a(:,3), assuming you want the range 10 to 20
    meanA4 = mf;     % Midpoint of the range [5, 13] for a(:,4)

    % Adjusting gamma distribution parameters for e and a (except e(:,2) and a(:,2))
    gammaShapeE1 = meanE1^2 / var;
    gammaRateE1 = meanE1 / var;
    
    gammaShapeE3 = meanE3^2 / var;
    gammaRateE3 = meanE3 / var;

    gammaShapeA1 = meanA1^2 / var;
    gammaRateA1 = meanA1 / var;
    
    gammaShapeA3 = meanA3^2 / var;
    gammaRateA3 = meanA3 / var;

    % Generate sigma^2 using a gamma distribution
    sigma2 = gamrnd(10, 1);  % Shape=10, Rate=1, Mean=10

    % Generate matrix e (Nr x 3)
    e1 = abs(normrnd(meanE1, sqrt(var)/10, Nv, 1));  % e(:,1)
    e2 = abs(normrnd(meanE2, sqrt(var)/1000, Nv, 1));       % e(:,2) using normal distribution
    e3 = gamrnd(gammaShapeE3, 1/gammaRateE3, Nv, 1);  % e(:,3)
    e = [e1, e2, e3];

    % Generate matrix a (Nr x 4)
    a1 = abs(normrnd(meanA1, sqrt(var)/1000, Nv, 1));  % a(:,1)
    a2 = abs(normrnd(meanA2, sqrt(var)/1000, Nv, 1));       % a(:,2) using normal distribution
    a3 = gamrnd(gammaShapeA3, 1/gammaRateA3, Nv, 1);  % a(:,3)
    a4 = normrnd(meanA4, sqrt(stdf), Nv, 1);          % a(:,4)
    a4 = max(min(a4, 13), 5);  %// Clipping to [7, 13]

    a = [a1, a2, a3, a4];

    % Reshape matrices e and a into vectors and concatenate with sigma^2
    x = v2x(e, a, sigma2);
end
