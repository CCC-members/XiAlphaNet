function [S,Data,Theta] = gen_hggm1(param)
%% Hermitian Gaussian Graphical Model generator 
%  the gaussian engine samples come from a random binary precision matrix
%  the precision matrix is defined as a block structure with overlaping
%  diagonal blocks with random sizes 
%%
% Pedro Valdes-Sosa, Oct 2017
% Deirel Paz Linares, Oct 2017
% Eduardo Gonzalez-Moreira, Oct 2017
%%
m        = param.m;
q        = param.q;
nblocks  = param.nblocks;
if param.config == 1
    %% Non overlapping blocks
    X                  = eye(q); 
    %Infimum and maximum size of blocks
    infsize            = ceil(q/(2*nblocks));
    maxsize            = ceil(q/nblocks) - 1;
    %Cycle by nblocks to fill the matrix X 
    for i = 1:nblocks
        size           = randi([infsize maxsize]);
        blockRe        = 2*(rand(size) - 0.5);
        blockRe        = (blockRe + blockRe')/2;
        if param.var == 1
            block = blockRe;
        elseif param.var == 2
            blockIm    = 2*(rand(size) - 0.5);
            blockIm    = (blockIm - blockIm')/2;
            block      = blockRe + 1i*blockIm;
        end
        block          = block - diag(diag(block)) + diag(abs(diag(block)));
        dmin           = min(eig(block));
        if dmin < 0
            block      = block + abs(dmin)*eye(size) + eye(size);
        else
            block      = block - abs(dmin)*eye(size) + eye(size);
        end
        index          = randi([(i-1)*maxsize + 1 i*maxsize],size,1);
        X(index,index) = block;
    end
end
if param.config == 2
    %% Overlapping blocks 
    X        = eye(q); 
    %Infimum and maximum size of blocks
    infsize  = ceil(q/(2*nblocks));
    maxsize  = ceil(q/nblocks) - 1;
    %Cycle by nblocks to fill the matrix X
    for i = 1:nblocks
        size           = randi([infsize maxsize]);
        blockRe        = 2*(rand(size) - 0.5);
        blockRe        = (blockRe + blockRe')/2;
        if param.var == 1
            block = blockRe;
        elseif param.var == 2
            blockIm    = 2*(rand(size) - 0.5);
            blockIm    = (blockIm - blockIm')/2;
            block      = blockRe + 1i*blockIm;
        end      
        if i == nblocks
            index      = randi([(i-1)*maxsize + 1 i*maxsize],size,1);
        else
            index      = randi([(i-1)*maxsize + 1 ceil((i + 1/2)*maxsize)],size,1);
        end
        X(index,index) = block; 
    end
    X        = X - diag(diag(X)) + diag(abs(diag(X)));
    dmin     = min(eig(X));
    if dmin < 0
        X    = X + abs(dmin)*eye(q) + eye(q);
    else
        X    = X - abs(dmin)*eye(q) + eye(q);
    end    
end
%% Applying isomorphism and generating data
W               = eye(q)/X;
W               = (W + W')/2;
if param.var == 1
    Data        = mvnrnd(zeros(1,q),W,m);
elseif param.var == 2
    Wisomph     = [real(W) -imag(W); imag(W) real(W)];
    Data_isomph = mvnrnd(zeros(1,2*q),Wisomph,m);
    DataRe      = Data_isomph(:,1:q);
    DataIm      = Data_isomph(:,q+1:2*q);
    Data        = DataRe + 1i*DataIm;
end
Theta       = eigendecomposition(X,param);
S           = eigendecomposition((1/m)*(Data'*Data),param);