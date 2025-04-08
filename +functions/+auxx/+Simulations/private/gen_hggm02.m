function [Sigma,Theta] = gen_hggm02(q,options)
% q: number of active nodes(include extensions)(sum(options.extensions))


% Pedro Valdes-Sosa, Oct 2017
% Deirel Paz Linares, Oct 2017
% Eduardo Gonzalez-Moreira, Oct 2017
% Ying Wang, Oct 2019
%%
extensions     = options.extensions;
connections    = options.connections;
%%
nblocks        = size(extensions,1);
nconnections   = size(connections,1);
Index_cell     = cell(nblocks,1);
index          = 0;
Theta              = eye(q);
%Cycle by nblocks to fill the matrix X
for i = 1:nblocks
    %Infimum and maximum size of blocks
    ext        = extensions(i);
    blockRe    = randn(ext);
    blockIm    = randn(ext);
    block      = blockRe + 1i*blockIm;
    index      = index(end)+1:index(end) + ext;
    Index_cell{i} = index;
    Theta(index,index)    = block;
end
%% Cycle by nconnections to fill the matrix X
for i = 1:nconnections
    pair       = sort(connections(i,:));
    index_row  = Index_cell{pair(1)};
    index_col  = Index_cell{pair(2)};
    blockRe    = randn(length(index_row),length(index_col)); %2*(rand(length(index_row),length(index_col)) - 0.5);
    blockIm    = randn(length(index_row),length(index_col)); %2*(rand(length(index_row),length(index_col)) - 0.5);
    block      = blockRe + 1i*blockIm;
    Theta(index_row,index_col)    = block;
    Theta(index_col,index_row)    = block';
end
%% Positive Definiteness
th                       = 1.5;
% th                       = 0.1;
Theta                    = (Theta + Theta')/2;
Xndiag                   = Theta - diag(diag(Theta));
Theta(abs(Xndiag) < th)  = 0;
Theta(abs(Xndiag) >= th) = Theta(abs(Xndiag) >= th)./abs(Theta(abs(Xndiag) >= th));
Theta                    = (Theta + Theta')/2;
dmin                     = min(eig(Theta));
if dmin < 0
    Theta  = Theta + abs(dmin)*eye(q) + eye(q);
else
    Theta  = Theta - abs(dmin)*eye(q) + eye(q);
end
%% Applying isomorphism and generating data
Sigma           = eye(q)/Theta;


end









