function misRecogSample=stocSim(n,grid,sample,sigmaRFs,noiseLevel,pos)

    sigmaNoise=(noiseLevel/2)*(1/(sigmaRFs*sqrt(2*pi)));
    
    gap=1/(n+1);
    peaks=linspace(gap,1-gap,n);
    responses=normpdf(pos,peaks,sigmaRFs);
    
    refPoints=linspace(0,1,grid);
    refResponses=zeros(grid,n);
    for i=1:grid
        refResponses(i,:)=normpdf(refPoints(i),peaks,sigmaRFs);
    end
    
    points=zeros(sample,n);
    for i=1:n
        points(:,i)=responses(i)+sigmaNoise*randn(1,sample);
    end
    
    exPoints=repmat(reshape(points, [1,sample, n]), [grid,1,1]);
    exRefPoints=repmat(reshape(refResponses, [grid, 1, n]), [1,sample,1]);
    
    dis=exPoints-exRefPoints;
    dis=squeeze(sum(dis.^2,[3]));
    dis=normpdf(dis,0,sigmaNoise);
    misRecogSample=zeros(sample,1);
    for i=1:sample
        misRecogSample(i)=weightedRandomSelection(squeeze(dis(:,i)));
    end
    
    edges=linspace(0,1,grid);
    misRecogSample=misRecogSample/grid;
    % figure,histogram(misRecogSample,edges);
    
    function idx = weightedRandomSelection(weights)

        cumWeights = cumsum(weights);
        
        totalWeight = cumWeights(end);
        r = rand() * totalWeight;
        
        idx = find(cumWeights >= r, 1, 'first');

        if(norm(weights)==0)
            idx=randi(length(weights)); 
        end
    end
end
