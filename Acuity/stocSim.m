n=3;
sample=1000;
refPointsNumber=25;
sigmaRFs=0.1;
noiseLevel=0.2;
pos=0.5;

sigmaNoise=(noiseLevel/2)*(1/(sigmaRFs*sqrt(2*pi)));

gap=1/(n+1);
peaks=linspace(gap,1-gap,n);
responses=normpdf(pos,peaks,sigmaRFs);

refPoints=linspace(0,1,refPointsNumber);
refResponses=zeros(refPointsNumber,n);
for i=1:refPointsNumber
    refResponses(i,:)=normpdf(refPoints(i),peaks,sigmaRFs);
end

points=zeros(sample,n);
for i=1:n
    points(:,i)=responses(i)+sigmaNoise*randn(1,sample);
end

exPoints=repmat(reshape(points, [1,sample, n]), [refPointsNumber,1,1]);
exRefPoints=repmat(reshape(refResponses, [refPointsNumber, 1, n]), [1,sample,1]);

dis=exPoints-exRefPoints;
 dis=squeeze(sum(dis.^2,[3]));
 dis=normpdf(dis,0,sigmaNoise);
misRecogSample=zeros(sample,1);
for i=1:sample
    misRecogSample(i)=weightedRandomSelection(squeeze(dis(:,i)));
end

edges=linspace(0,1,refPointsNumber);
misRecogSample=misRecogSample/refPointsNumber;
figure,histogram(misRecogSample,edges);

function idx=weightedRandomSelection(weights)

    cumWeights=cumsum(weights);
    
    
    totalWeight=cumWeights(end);
    r=rand()*totalWeight;
    
    idx=find(cumWeights >= r, 1, 'first');
    
    if isempty(idx)
        idx=randi(length(weights)); 
    end
end
