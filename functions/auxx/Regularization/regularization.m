function [Xreg] = regularization(X)

X=(X+X')/2;

emin=min(eig(X));
p=0;
delta=0.0001;
[d,~]=size(X);
L=trace(X)/d;
Xreg=X;

%------------------------------------
while emin<0&& p<1
    Xreg=(1-p)*X+p*L*eye(size(X));
    emin=min(eig(Xreg));
    p=p+delta;
end
Xreg =Xreg;
%------------------------------------
% if lplot==1
%     figure(1)
%       surf(abs(X));
%     figure(2)
%       surf(abs(Xreg));
% end

end