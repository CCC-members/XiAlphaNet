function [mf,stdf,mr,stdr] = init_freq(Ne,S,freq)

import functions.*
import functions.auxx.*
import functions.auxx.CrossSpectrum.*

alpha_band = [7,13];
spec = exp(log_spectrum(S, freq));
freq_list = [];
ratio = [];
for j=1:Ne
    alpha_indices = find(freq>=alpha_band(1) & freq<=alpha_band(2));
    alpha_spect = spec(j,alpha_indices);
    alpha_freq = freq(alpha_indices);
    [~,positions] = max(alpha_spect);
    freq_list = [freq_list,alpha_freq(positions)];
    ratio = [ratio, spec(j,1)/max(alpha_spect)];
end
mf = mean(freq_list);
stdf = std(freq_list);
mr = mean(ratio);
stdr = std(ratio);
end