function [a,l] = extract_non_nan(A,L)
    % A: Input matrix that may contain NaN values
    [N,M] =size(A);
    a=[];
    l=[];
    for i=1:N
        for j=1:i
            if ~isnan(A(i,j))
                a=[a;A(i,j)];
                l=[l;L(i,j)];
            end
        end
    end         
   
end