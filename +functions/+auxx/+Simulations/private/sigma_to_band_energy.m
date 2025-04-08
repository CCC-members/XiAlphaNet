function energy=sigma_to_band_energy(sigma,band,deltaf)
% Authors:
% - Ying Wang

% Date: Nov 24, 2019

f=[round(band(1)/deltaf)+1:1:round(band(2)/deltaf)+1];
[m,~,s]=size(sigma);
energy=zeros(m,length(f));
% energy=zeros(m,1);
for i=1:m
energy(i,:)=sqrt(squeeze(sigma(i,i,f)));
end
end


