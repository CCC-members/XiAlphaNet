
parameters.Stochastic.stoch = 0; 
parameters = sample_frequencies(parameters);
%x = generateRandomSample(parameters,0.0001);
S = reconstructSDM(x,parameters);
Sr = parameters.Data.Cross;
[Nd,~,Nf] =size(S);
    for f = 1:Nf
        ls(:,f) = real(log(real(diag(S(:,:,f))))); 
        lsr(:,f) = real(log(real(diag(Sr(:,:,f))))); 
    end
    fs =49;
    figure(4)
    subplot(2,1,1)
    plot(parameters.Data.freq(1:fs),ls(:,1:fs)')
    title('Estimated Scalp LogSpectra')
    subplot(2,1,2)
    plot(parameters.Data.freq(1:fs),lsr(:,1:fs)')
    title('Scalp LogSpectra')
    d = w_distance(S(:,:,1),Sr(:,:,1))
    d2  = norm(S(:,:,1)-Sr(:,:,1),'fro')
    figure(2)
    subplot(1,2,1)
    surf(abs(mean(S,3)),'EdgeColor','none')
     subplot(1,2,2)
    surf(abs(mean(Sr,3)),'EdgeColor','none')
    colormap('hot')