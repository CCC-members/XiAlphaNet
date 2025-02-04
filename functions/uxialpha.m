%x = generateRandomSample(parameters,0.01);
[e,a,s2] = x2v(x0);
os = parameters.Data.freq;
xi_omega = @(o) e(:,1) ./ (1 + e(:,2) .* o.^2).^e(:,3);
alpha_omega = @(o) a(:,1) ./ (1 + a(:,2) .* (o- a(:,4)).^2).^a(:,3);
xi = xi_omega(os);
al = alpha_omega(os);
plot(os,mean(log(1+xi'+al'),2));