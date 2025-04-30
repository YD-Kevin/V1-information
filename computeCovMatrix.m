function covMatrix = computeCovMatrix(data, T)

    n=size(data,1);
    
    TMatrix=zeros(n, numel(T_func(data(1,:)))); 
    for i=1:n
        TMatrix(i,:)=T(data(i,:)); 
    end
    covMatrix=cov(TMatrix);
end