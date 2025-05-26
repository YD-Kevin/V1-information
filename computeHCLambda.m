function Lambda = computeHCLambda(map,n,nDist,sigma)


%sigma = 1.7 cycles/mm


imgFFT = fft2(img - mean(img(:)));
powerSpectrum = abs(fftshift(imgFFT)).^2;

[M, N] = size(img);
[KX, KY] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
k = sqrt(KX.^2 + KY.^2) * (1/(pixelSize*n)); 
k = k(:);
power = powerSpectrum(:);


kEdges = 0:0.1:max(k); 
kCenters = (kEdges(1:end-1) + kEdges(2:end))/2;
P = zeros(size(kCenters));
for i = 1:length(kCenters)
    mask = k >= kEdges(i) & k < kEdges(i+1);
    P(i) = mean(power(mask));
end

P_smooth = imgaussfilt(P, sigma/(kEdges(2)-kEdges(1)));


[~, peakIdx] = max(P_smooth);
q = kCenters(peakIdx); % cycles/mm


Lambda = 2*pi / q; % mm
end