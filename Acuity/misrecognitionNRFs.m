function misRec=misrecognitionNRFs(n,sample,FWHM,noiseLevel,pos)

    sigmaRFs=FWHM/2*sqrt(2*log(2)); 
    sigmaNoise=(noiseLevel/2)*(1/(sigmaRFs*sqrt(2*pi)));
    GTresponses=posToResponses(pos,n,sigmaRFs);
    gVal=zeros(sample,1);


    sampleResponses=zeros(sample,n);
    for i=1:sample
        temSampleResponses=posToResponses((i-1)*(1/(sample-1)),n,sigmaRFs);
        l=sqrt(sum((temSampleResponses-GTresponses).^2,"all"));
        sampleResponses(i,:)=temSampleResponses;
        integrand = @(r1) (1 ./ (sqrt(2*pi))) .* sech(l.*r1./(2.*sigmaNoise)) .* exp(-r1.^2./2);
        gVal(i)=integral(integrand, -Inf, Inf);
    end
    
    extGTresponses=repmat(reshape(GTresponses, [1, n]), [sample,1]);
    
    distanceSQ=sum((sampleResponses-extGTresponses).^2,[2]);
    
    B=(0.5*exp(-distanceSQ/(8*sigmaNoise^2))).*gVal;
    misRec=B./(1-B);
    
    misRec=sample*misRec./(sum(misRec,"all"));
    
        function responses=posToResponses(pos,n,sigmaRFs)
            gap=1/(n+1);
            peaks=linspace(gap,1-gap,n);
            responses=normpdf(pos,peaks,sigmaRFs);
        end
end
