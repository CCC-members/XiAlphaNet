function [Log_Spec,freq]= log_spectrum(A,parameters);
freq = parameters.Data.freq;
Nf  = length(freq);
for f = 1:Nf
    ls(:,f) = real(log(real(diag(A(:,:,f))))); 
end
Log_Spec = ls;
end
