sample=100;
sigmaRFs=0.12;
noiseLevel=0.1;
pos=0.5;

misrecognition3RFs(sample,sigmaRFs,noiseLevel,pos);
xcoord=linspace(1/sample,1,sample);

figure,plot(xcoord,misRec);


misRecOverNoise=zeros(10,100);
for i=1:10
    misRecOverNoise(i,:)=misrecognition3RFs(100,0.12,0.1*i,0.5);
end
figure,surf(misRecOverNoise);
