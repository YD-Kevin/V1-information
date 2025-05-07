function blobs=buildColorBlobs_1(rho,lcRad,n,threshold,majorAxis,minorAxis,noiseLevel)

map=buildMap(rho,lcRad*2+n,0,0);
shiftedMap=buildMap(rho,lcRad*2+n,0.5,0.5);
Fx=cos(map*2-pi);
Fy=sin(map*2-pi);

[dfydx, ~]=gradient(Fy); % dFy/dx
[~, dfxdy]=gradient(Fx); % dFx/dy

curl =dfydx-dfxdy;
binaryCurl=double(abs(curl) >= threshold);

    

[nRows, nCols]=size(binaryCurl);
finalMask=false(nRows, nCols);


[yCoords, xCoords]=find(binaryCurl);


theta=linspace(0, 2*pi, 100); 
baseR=(majorAxis * minorAxis) ./ sqrt((minorAxis*cos(theta)).^2+(majorAxis*sin(theta)).^2);


for i=1:length(xCoords)
   
    noise=1 + noiseLevel * randn(size(theta));
    distortedR=baseR .* noise;
    
   
    x=xCoords(i) + distortedR .* cos(theta);
    y=yCoords(i) + distortedR .* sin(theta);
    
    
    mask=poly2mask(x, y, nRows, nCols);
    finalMask=finalMask | mask;
end


blobs=shiftedMap .* finalMask;

end

