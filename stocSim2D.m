function misRecogSample=stocSim2D(n,m,grid,sample,sigmaRFsX,sigmaRFsY,noiseLevel,posX,posY)

    sigmaNoise=(noiseLevel/2)*(1/((sigmaRFsX*sigmaRFsY)*2*pi));
    
    gapX=1/(n+1);
    gapY=1/(m+1);
    peaksX=linspace(gapX,1-gapX,n);
    peaksY=linspace(gapY,1-gapY,m);
    [X,Y]=meshgrid(peaksX,peaksY);
    responses=normpdf(posX,X,sigmaRFsX).*normpdf(posY,Y,sigmaRFsY);
    

    
    refPoints=linspace(0,1,grid);
    refResponsesX=normpdf(0,peaksX-refPoints',sigmaRFsX);
    refResponsesY=normpdf(0,peaksY-refPoints',sigmaRFsY);
    refResponsesX=repmat(reshape(refResponsesX, [grid, 1, n, 1]), [1,grid,1,m]);
    refResponsesY=repmat(reshape(refResponsesY, [1, grid, 1, m]), [grid,1,n,1]);
    refResponses=refResponsesX.*refResponsesY;

    
    
    points=repmat(reshape(responses, [1, n, m]), [sample,1,1]);
    points=points+sigmaNoise*randn(sample,n,m);

    
    exPoints=repmat(reshape(points, [1,1,sample, n,m]), [grid,grid,1,1,1]);
    exRefPoints=repmat(reshape(refResponses, [grid,grid, 1, n,m]), [1,1,sample,1,1]);
    
    dis=exPoints-exRefPoints;
    dis=sqrt(squeeze(sum(dis.^2,[4,5])));
    pdf=normpdf(dis,0,sigmaNoise);
    pdfS=sum(pdf,[1,2])+1e-6;
    pdfNorm=(pdf+1e-6)./pdfS;
    misRecogSample=sum(pdfNorm,[3]);
    misRecogSample=misRecogSample/sum(misRecogSample,'all');
end
