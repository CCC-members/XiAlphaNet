x=x0(:,1);
S1 = reconstructSDM(x,parameters);
x=x0(:,2);
S2 = reconstructSDM(x,parameters);
[Nd,~,Nf] =size(S);
for f = 1:Nf
    ls(:,f) = real(log(real(diag(S1(:,:,f))))); 
    lsr(:,f) = real(log(real(diag(S2(:,:,f))))); 
end
fs =49;
figure(1)
subplot(2,1,1)
plot(parameters.Data.freq(1:fs),ls(:,1:fs)'-lsr(:,1:fs)')
title('Estimated Scalp LogSpectra')
subplot(2,1,2)
plot(parameters.Data.freq(1:fs),lsr(:,1:fs)')
title('Scalp LogSpectra')

  figure(2)
    subplot(1,2,1)
    surf(abs(mean(S1,3)),'EdgeColor','none')
     subplot(1,2,2)
    surf(abs(mean(S2,3)),'EdgeColor','none')
    colormap('hot')

     d = w_distance(S2(:,:,1),S1(:,:,1))
    d2  = norm(S1(:,:,1)-S2(:,:,1),'fro')
    norm(ls-lsr,'fro')