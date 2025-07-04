function misRec=misrecognition3RFs(sample,sigmaRFs,noiseLevel,pos)

    %Input:
    %sample:number of sample points taken in [0,1]
    %sigmaRFs:SD of the tuning curve of RFs
    %noiseLevel:a value in [0,1] arranging the noise intensity, (SD of noise)/(max of the tuning curve/2)
    %pos:ground truth we are interested in, position in [0,1]
    
    
    
    t=0.25;
    sigmaNoise=(noiseLevel/2)*(1/(sigmaRFs*sqrt(2*pi)));
    responseMax=(1/(sigmaRFs*sqrt(2*pi)))+(2*sigmaNoise);
    pos=round(pos*sample);
    
    [I,J,K]=meshgrid(1:sample,1:sample,1:sample);
    coreA=zeros(sample,sample,sample,3);
    coreA(:,:,:,1)=(1/sample)*I;
    coreA(:,:,:,2)=(1/sample)*J;
    coreA(:,:,:,3)=(1/sample)*K;
    coreA=permute(coreA,[2,1,3,4]);
    
    
    A=repmat(reshape(coreA, [1, sample, sample, sample,3]), [sample, 1, 1, 1, 1])*responseMax;
    %A: the tensor that records coordinates in neural space
    
    
    aVec=(1/sample:1/sample:1)';
    B_base=[aVec-t,aVec-2*t,aVec-3*t]; 
    
    
    B=repmat(B_base,[1,1,sample,sample,sample]); 
    B=permute(B,[1,3,4,5,2]);          
    
    
    FB=normpdf(B,0,sigmaRFs);
    %FB:the tensor that records ideal neural acts of RFs in response to stimuli
    
    
    dif=FB-A;
    pdf3RFs=normpdf(dif,0,sigmaNoise);
    pdf3RFs=squeeze(pdf3RFs(:,:,:,:,1).*pdf3RFs(:,:,:,:,2).*pdf3RFs(:,:,:,:,3));
    %pdf3RFs:pdf for joint distribution p(pos,act1,act2,act3)
    
    pdfPos=squeeze(pdf3RFs(pos,:,:,:));
    [n1, n2, n3] = size(pdfPos);
    pdfPosExt = repmat(reshape(pdfPos, [1, n1, n2, n3]), [n1, 1, 1, 1]);
    
    
    misRec=pdfPosExt.*pdf3RFs;
    misRec=sum(misRec,[2,3,4]);
    misRec=(misRec*sample)/sum(misRec,"all");
    %xcoord=linspace(1/sample,1,sample);
    %figure,plot(xcoord,misRec);
end
