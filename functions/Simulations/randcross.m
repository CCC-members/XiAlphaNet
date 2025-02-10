function svv = randcross(Ne,Nf);
svv=zeros(Ne,Ne,Nf);
for j=1:Nf
    svv(:,:,j) =rand(Ne,Ne);
    svv(:,:,j) = regularization(svv(:,:,j));
end
end

