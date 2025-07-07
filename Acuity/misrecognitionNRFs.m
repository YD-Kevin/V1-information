function misRec=misrecognitionNRFs(n,sample,sigmaRFs,noiseLevel,pos)


sigmaNoise=(noiseLevel/2)*(1/(sigmaRFs*sqrt(2*pi)));
GTresponses=posToResponses(pos,n,sigmaRFs);

sampleResponses=zeros(sample,n);
for i=1:sample
    sampleResponses(i,:)=posToResponses((i-1)*(1/(sample-1)),n,sigmaRFs);
end

extGTresponses=repmat(reshape(GTresponses, [1, n]), [sample,1]);

distance=sum((sampleResponses-extGTresponses).^2,[2]);

misRec=1./exp((distance.^2)/(sigmaNoise^2));
misRec=sample*misRec./(sum(misRec,"all"));

    function responses=posToResponses(pos,n,sigmaRFs)
        gap=1/(n+1);
        peaks=linspace(gap,1-gap,n);
        responses=normpdf(pos,peaks,sigmaRFs);
    end
end
