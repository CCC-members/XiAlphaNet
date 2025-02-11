function [S,Data,X] = gen_hggm2(m,q,options)
%% Hermitian Gaussian Graphical Model generator 
%  the gaussian engine samples come from a random binary precision matrix
%  the precision matrix is defined as a block structure with nonoverlaping
%  diagonal blocks with size defined in options and given off-diagonal nonzero 
%  blocks defined in options 
%%
% Pedro Valdes-Sosa, Oct 2017
% Deirel Paz Linares, Oct 2017
% Eduardo Gonzalez-Moreira, Oct 2017
%%
extensions     = options.extensions;
connections    = options.connections;
%%
nblocks        = size(extensions,1);
nconnections   = size(connections,1);
Index_cell     = cell(nblocks,1);
index          = 0;
X              = eye(q);
%Cycle by nblocks to fill the matrix X
for cont1 = 1:nblocks
    %Infimum and maximum size of blocks
    ext        = extensions(cont1);
    blockRe    = randn(ext);
    blockIm    = randn(ext);
    block      = blockRe + 1i*blockIm;
    index      = index(end)+1:index(end) + ext; 
    Index_cell{cont1} = index;
    X(index,index)    = block;
end
%% Cycle by nconnections to fill the matrix X
for cont2 = 1:nconnections
    pair       = sort(connections(cont2,:));
    index_row  = Index_cell{pair(1)};
    index_col  = Index_cell{pair(2)};
    blockRe    = randn(length(index_row),length(index_col)); %2*(rand(length(index_row),length(index_col)) - 0.5);
    blockIm    = randn(length(index_row),length(index_col)); %2*(rand(length(index_row),length(index_col)) - 0.5);
    block      = blockRe + 1i*blockIm;
    X(index_row,index_col)    = block;
    X(index_col,index_row)    = block';
end
%% Positive Definiteness
th                   = 1.5; 
X                    = (X + X')/2; 
Xndiag               = X - diag(diag(X));
X(abs(Xndiag) < th)  = 0; 
X(abs(Xndiag) >= th) = X(abs(Xndiag) >= th)./abs(X(abs(Xndiag) >= th));
X                    = (X + X')/2; 
dmin                 = min(eig(X));
if dmin < 0
    X  = X + abs(dmin)*eye(q) + eye(q);
else
    X  = X - abs(dmin)*eye(q) + eye(q);
end
%% Applying isomorphism and generating data
W           = eye(q)/X;
Wisomph     = [real(W) -imag(W); imag(W) real(W)];
Data_isomph = mvnrnd(zeros(1,2*q),Wisomph,m);
DataRe      = Data_isomph(:,1:q);
DataIm      = Data_isomph(:,q+1:2*q);
Data        = DataRe + 1i*DataIm;
Data        = transpose(Data);
S           = (1/m)*(Data*Data');