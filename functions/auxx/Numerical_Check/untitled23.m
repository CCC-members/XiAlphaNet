% Define the small example
N=10;
x= rand(N,1);
%x = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12; 13; 14;15];

% Use the x2v function to decompose x into e, a, and sigma2
[e, a, sigma2] = x2v(x);

% Display the results
% disp('e:');
% disp(e);
% disp('a:');
% disp(a);
% disp('sigma2:');
% disp(sigma2);

% Use the v2x function to convert e, a, and sigma2 back to x
x_reconstructed = v2x(e, a, sigma2);

% % Display the reconstructed x
% disp('Reconstructed x:');
% disp(x_reconstructed);

mean(x - x_reconstructed)