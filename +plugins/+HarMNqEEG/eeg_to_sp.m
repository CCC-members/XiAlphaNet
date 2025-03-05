function [Spec, fmin, freqres, fmax, CrossM, ffteeg] = eeg_to_sp(EEG, sp)
% EEG: EEG scalp data matrix, organized as nd x nt x ne, where
% nd : number of channels
% nt : epoch size (# of instants of times in an epoch)
% ne : number of epochs
% sp : sampling period in seconds. Eg: 0.005 means 5 millisec

[nd, nt, ne] = size(EEG);
NFFT = nt;
freqres = 1 ./ (sp*NFFT);
fmin = freqres;
fmax = 1 ./ (2*sp);

ncoefs = floor(nt./2);
ffteeg = zeros(nd, ncoefs, ne);

% wind = reshape(hanning(nt), 1, nt); %multiply the EEG signal by Hanning windows

CrossM = zeros(nd, nd, ncoefs); % allocated matrix for the cross spectrum
for k=1:ne
  % A=EEG(:,:,k) .* repmat(wind, nd, 1);  %Applying Hanning window
  A=EEG(:,:,k);  %In this version we don't use windows to be compatible with the cross-spectra of the Cuban Normative data 1990
  f = fft(A,NFFT,2);
  f = f(:, 2:ncoefs+1); %Nyquist
  ffteeg(:,:,k) = f;
  
  for h=1:ncoefs
      CrossM(:,:,h) = CrossM(:,:,h) + f(:,h) * f(:,h)';
  end

end
CrossM = CrossM ./ ne;
Spec = zeros(nd, ncoefs); % allocated matrix for the spectrum
for k=1:size(CrossM,3)
    Spec(:,k) = diag(CrossM(:,:,k));
end
Spec = real(Spec);



