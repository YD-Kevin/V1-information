function misRec=misrecognitionNRFs2D(n,m,sample,sigmaRFsX,sigmaRFsY,rho,noiseLevel,posX,posY)


    sigmaNoise=(noiseLevel/2)*(1/((sigmaRFsX*sigmaRFsY)*2*pi*sqrt(1-rho^2)));
    GTresponses=posToResponses(posX,posY,n,m,sigmaRFsX,sigmaRFsY,rho);
    
    sampleResponses=zeros(sample,sample,n,m);
    for i=1:sample
        for j=1:sample
            sampleResponses(i,j,:,:)=posToResponses((i-1)*(1/(sample-1)),(j-1)*(1/(sample-1)),n,m,sigmaRFsX,sigmaRFsY,rho);
        end
    end
    
    extGTresponses=repmat(reshape(GTresponses, [1, 1, n, m]), [sample,sample,1,1]);
    
    distanceSQ=sum((sampleResponses-extGTresponses).^2,[3,4]);
    
    misRec=(0.5*exp(-distanceSQ/(8*sigmaNoise^2)))./(1-0.5*exp(-distanceSQ/(8*sigmaNoise^2)));
    
    misRec=sample*misRec./(sum(misRec,"all"));
    
        function responses=posToResponses(posX,posY,n,m,sigmaRFsX,sigmaRFsY,rho)
            gapX=1/(n+1);
            gapY=1/(m+1);
            peaksX=linspace(gapX,1-gapX,n);
            peaksY=linspace(gapY,1-gapY,m);
            [X,Y]=meshgrid(peaksX,peaksY);

            responsesX=(normpdf(posX,X,sigmaRFsX)')*(sqrt(2*pi)*sigmaRFsX);
            responsesY=normpdf(posY,Y,sigmaRFsY)'*(sqrt(2*pi)*sigmaRFsY);
            responsesCor=exp((rho*(posX-X').*(posY-Y'))/((1-rho^2)*sigmaRFsX*sigmaRFsY));
            responses=responsesX.^(1/(1-rho^2)).*responsesY.^(1/(1-rho^2)).*responsesCor;
            responses=responses*(1/((sigmaRFsX*sigmaRFsY)*2*pi*sqrt(1-rho^2)));
        end
end