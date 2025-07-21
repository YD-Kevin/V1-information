function misRecogList=stocSim2DFast(n,m,grid,sample,sigmaRFsX,sigmaRFsY,noiseLevel,posListX,posListY,posGrid)

    sigmaNoise=(noiseLevel/2)*(1/((sigmaRFsX*sigmaRFsY)*2*pi));
    
    gapX=1/(n+1);
    gapY=1/(m+1);
    peaksX=linspace(gapX,1-gapX,n);
    peaksY=linspace(gapY,1-gapY,n);
    [X,Y]=meshgrid(peaksX,peaksY);
    refPoints=linspace(0,1,grid);
    refResponsesX=normpdf(0,peaksX-refPoints',sigmaRFsX);
    refResponsesY=normpdf(0,peaksY-refPoints',sigmaRFsY);
    refResponsesX=repmat(reshape(refResponsesX, [grid, 1, n, 1]), [1,grid,1,m]);
    refResponsesY=repmat(reshape(refResponsesY, [1, grid, 1, m]), [grid,1,n,1]);
    refResponses=refResponsesX.*refResponsesY;
    exRefPoints=repmat(reshape(refResponses, [grid,grid, 1, n,m]), [1,1,sample,1,1]);
    misRecogList=zeros(posGrid,posGrid,grid,grid);


    for i=1:posGrid
        for j=1:posGrid
            responses=normpdf(posListX(i),X,sigmaRFsX).*normpdf(posListY(j),Y,sigmaRFsY);
            points=repmat(reshape(responses, [1, n, m]), [sample,1,1]);
            points=points+sigmaNoise*randn(sample,n,m);
            exPoints=repmat(reshape(points, [1,1,sample, n,m]), [grid,grid,1,1,1]);
            dis=exPoints-exRefPoints;
            dis=squeeze(sum(dis.^2,[4,5]));
            dis=normpdf(dis,0,sigmaNoise);
            misRecogSample=sum(dis,[3]);
            misRecogSample=misRecogSample./sum(misRecogSample,"all");
            misRecogList(i,j,:,:)=misRecogSample;
        end
    end
  
end
