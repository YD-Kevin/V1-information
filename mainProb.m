%% Biological Factors

%macaque

rf=1; %receptive field size by degree
nd=3000; %2-D neural density in V1 by cell/mm^2 (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
mf=2; %cortical magnification factor degree/mm (Roger B. H. Tootell et.alca., 1988)
lcScale=0.55; %a reference scale factor of lateral connection by mm
refAct=50;
%Cat

rf=3; %receptive field size by degree (Gilbert CD and Wiesel 1989, Fuyuki Karube and Zoltan F. Kisvarday, 2010)
nd=1800; %2-D neural density in V1 by cell/mm^2 (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
mf=5; %cortical magnification factor degree/mm (Gilbert CD and Wiesel 1989)
lcScale=0.55; %a reference scale factor of lateral connection by mm (Fuyuki Karube and Zoltan F. Kisvarday, 2010)
refAct=35;
% %Tree Shrew

% rf=3; %receptive field size by degree (Francois Mooser et.al 2004)
% nd=2000; %2-D neural density in V1 by cell/mm^2 (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
% mf=20; %cortical magnification factor degree/mm (Francois Mooser et.al 2004)
% lcScale=0.2; %a reference scale factor of lateral connection by mm (Francois Mooser et.al 2004)
%refAct=20;
%Mouse

rf=10; %receptive field size by degree (Jiakun Fu et.al, 2024)
nd=2200; %2-D neural density in V1 by cell/mm^2 (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
mf=120; %cortical magnification factor degree/mm (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
lcScale=0.03; %a reference scale factor of lateral connection by mm (Peijia Yu et.al. 2025)
refAct=10;
%% Parameters

n=64; %128*128 cells, number of cells in this patch of V1 cortex


m=128; %image input is 256*256 unit pixels
nDist=sqrt(1/nd); %(sqrt(n^2/nd)/n), convert the unit distance between cells to mm
visWindow=n*nDist*mf; %the correspondent vision range of this patch by degree
rfRad=round((m/visWindow)*rf/2); %the radius of receptive field by unit pixels
simFreq=6; %
lcSigma=lcScale/nDist;
lcRad=round(lcSigma*1.8);
aParam=1;

candidateFuncs={
    @(x) x,    % T1(x)=x
    @(x) x.^2,    % T2(x)=x^2
    @(x) x.^3,    % T3(x)=x^3
    @(x) log(abs(x)+1e-3), % T4(x)=log(|x|)
    @(x) abs(x)    % T5(x)=|x|
};

maxk=5;

%% MainProb


nCells=64;
mPixels=128;
aParam=10;
tParam=12;
rParam=8;

mapAlbum=zeros(10,lcRad*2+nCells,lcRad*2+nCells);

mapIndex=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1,1.5];

for i=1:1:aParam
    mapAlbum(i,:,:)=buildMap(mapIndex(i),lcRad*2+nCells,0,0);
end
actSetTheta=zeros(aParam,tParam,rParam,simFreq,nCells,nCells);
spkSetTheta=zeros(aParam,tParam,rParam,simFreq,nCells,nCells);

for i=1:aParam
    map=squeeze(mapAlbum(i,:,:));
    imgList=zeros(simFreq,mPixels,mPixels);
    album=buildAlbum(lcRad,nCells,lcSigma,map,i);
    
    for thetaVal=0:tParam-1
        for j=1:simFreq
            imgList(j,:,:)=200*buildRaster(thetaVal*15,j*0.2,9,mPixels);
        end
        
        for k=1:rParam
            IList=input(rfRad,lcRad,nCells,mPixels,simFreq,imgList,map)*refAct;
            [act,spk]=evo(lcRad,nCells,simFreq,album,IList,(i-1)*tParam*rParam+thetaVal*rParam+k);
            actSetTheta(i,thetaVal+1,k,:,:,:)=act;
            spkSetTheta(i,thetaVal+1,k,:,:,:)=spk;
        end
    end
end





%% 建立互相关场图册

function album=buildAlbum(lcRad,nCells,lcSigma,map,aParam)

album=zeros(nCells,nCells,lcRad*2+nCells,lcRad*2+nCells);
subMap=zeros(lcRad*2+1);
countVal=0;
%根据偏好图绘制各个神经元对周边神经元的作用分布图
for x=1:1:nCells
    for y=1:1:nCells
        subMap=map(x:x+lcRad*2,y:y+lcRad*2);
        subMap=buildLCProfile(lcRad,lcSigma,subMap);
        album(x,y,x:x+2*lcRad,y:y+2*lcRad)=subMap;
    end
end

end

%% 计算前馈刺激

function IList=input(rfRad,lcRad,nCells,mPixels,simFreq,imgList,map)
% [X, Y]=meshgrid(1:256);
% SX=sin(X/2)*100;
% SY=sin(Y/2)*100;
%
% temImg=ones(256)*30;
% I=buildInput(SY,map,lambda,nCells);

% randseed=randn(nCells);
% randseed=randseed./abs(randseed);

noiseLevel=mean(abs(imgList),'all')/5;


for i=1:1:simFreq
    imgList(i,:,:)=squeeze(imgList(i,:,:))+noiseLevel*randn(mPixels);
end

%根据图片输入，绘制前馈刺激
rfAlbum=buildRFAlbum(map,rfRad,lcRad,nCells);

IList=zeros(simFreq,nCells,nCells);
for i=1:1:simFreq
    IList(i,:,:)=buildInput(squeeze(imgList(i,:,:)),rfAlbum,rfRad,nCells);
    IList(i,:,:)=squeeze(IList(i,:,:))+randn(nCells).*squeeze(IList(i,:,:))*0.1;
end

end

%% 计算神经活动演化
function [act,spk]=evo(lcRad,nCells,simFreq,album,IList,aParam)
    act=zeros(simFreq,nCells,nCells);
    act(1,:,:)=ones(nCells);
    spk=zeros(simFreq,nCells,nCells);
    for t=2:1:simFreq
        prevAct=squeeze(act(t-1,:,:));
        actualSpk=polssrnd(prevAct/simFreq);
        spk(t,:,:)=actualSpk;



        bingo=squeeze(sum(album .* (actualSpk*simFreq), [1 2]));
        bingo=bingo(lcRad+1:lcRad+nCells,lcRad+1:lcRad+nCells);
        
        act(t,:,:)=relu(-0.3*prevAct+bingo+squeeze(IList(t,:,:)));
        

    end

end

function pos=relu(rate)
pos=(abs(rate)+rate)/2;
end

function cs=buildRaster(theta,posX,wid,m)
thetaRad=deg2rad(theta);
[H,S]=meshgrid(1:m,1:m);
hTheta=H*cos(thetaRad)+S*sin(thetaRad);
cs=cos(2*pi*hTheta/(wid*2)+posX*2*pi);
end
