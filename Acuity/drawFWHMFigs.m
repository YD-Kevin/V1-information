nGrid=60;
sizeGrid=30;
noiseGrid=1;
posGrid=10;
count=0;
sample=1000;
posList=rand(1,10);
compSdAlbum=zeros(nGrid,sizeGrid,noiseGrid,posGrid);
for i=1:nGrid
    for j=1:sizeGrid
        for k=1:noiseGrid
            for l=1:posGrid
                compPDF=misrecognitionNRFs(i*10,sample,j*0.025,0.37,posList(l));
                compSdAlbum(i,j,k,l)=pdfStats(linspace(0,1,1000), compPDF,posList(l));
            end
        end
        count=((i-1)*30+j)/1800
    end
end

FWHM=mean(compSdAlbum,[4])*2*sqrt(2*log(2));


figure,
imagesc([0.25,0.75],[200,600],FWHM(20:60,10:30)),colorbar,colormap jet;
xlabel("RF size");
ylabel("RF number");
title("FWHM for Decoding Distributions (Macaque, NoiseLevel 0.37)");