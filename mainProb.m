%% Biological Factors
% macaque
rf=1; %receptive field size by degree 
nd=3000; %2-D neural density in V1 by cell/mm^2 (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
mf=2; %cortical magnification factor degree/mm (Roger B. H. Tootell et.al., 1988)
lcScale=0.55; %a reference scale factor of lateral connection by mm 

%% Parameters

n=128; %128*128 cells, number of cells in this patch of V1 cortex
m=256; %image input is 256*256 unit pixels
nDist=sqrt(1/nd); %(sqrt(n^2/nd)/n), convert the unit distence between cells to mm
vWindow=n*nDist*mf; %the corespondent vision range of this patch by degree 
rfRad=round((m/vWindow)*rf/2); %the radius of receptive field by unit pixels
simFreq=40; %
lcSigma=lcScale/nDist;
lcRad=round(lcSigma*1.8);
a=1;


%% Test
 % map=buildMap(0.2,lcRad*2+n);
 % 
 % album=buildAlbum(lcRad,n,lcSigma,map,a);
 % 
 % img=buildRaster(0,0,15,m);
 % 
 % IList=input(rfRad,lcRad,n,m,simFreq,img,map);
 %
 % [act,Spk]=evo(lcRad,n,simFreq,album,IList,1);



mapAlbum=zeros(10,lcRad*2+n,lcRad*2+n);

for a=1:1:10
   mapAlbum(a,:,:)=buildMap(3/a-0.29,lcRad*2+n);
end
% album of reference maps with different density.

inputImg = double(imread('cameraman.tif'));
actAlbum=zeros(10,simFreq,n,n);
SpkAlbum=zeros(10,simFreq,n,n);


for a=1:1:10
    map=squeeze(mapAlbum(a,:,:)); 
    album=buildAlbum(lcRad,n,lcSigma,map,a);
    IList=input(rfRad,lcRad,n,m,simFreq,inputImg,map);
    [act,Spk]=evo(lcRad,n,simFreq,album,IList,a);
    actAlbum(a,:,:,:)=act;
    SpkAlbum(a,:,:,:)=Spk;
end



%% 建立互相关场图册

function album=buildAlbum(lcRad,n,lcSigma,map,a)

    album=zeros(n,n,lcRad*2+n,lcRad*2+n); 
    subMap=zeros(lcRad*2+1);
    count=0;
    %根据偏好图绘制各个神经元对周边神经元的作用分布图
    for x=1:1:n
        for y=1:1:n
    
            subMap=map(x:x+lcRad*2,y:y+lcRad*2);
            subMap=buildLCProfile(lcRad,lcSigma,subMap);
   
            album(x,y,x:x+2*lcRad,y:y+2*lcRad)=subMap;
            
        end
        count=((a-1)*128+x)/n
    end

end

%% 计算前馈刺激

 function IList=input(rfRad,lcRad,n,m,simFreq,inputImg,map)
% [X, Y]=meshgrid(1:256);
% SX=sin(X/2)*100;
% SY=sin(Y/2)*100;
% 
%  temImg=ones(256)*30;
%   I=buildInput(SY,map,lambda,n);
    

    randseed=randn(n);
    randseed=randseed./abs(randseed);

    imgList=zeros(simFreq,m,m);

    randIndex=mean(abs(inputImg),'all')/5;


    for i=1:1:simFreq
        imgList(i,:,:)=inputImg+randIndex*randn(m);
    end



    %根据图片输入，绘制前馈刺激
    IList=zeros(simFreq,n,n);
    for i=1:1:simFreq
        IList(i,:,:)=(buildInput(squeeze(imgList(i,:,:)),map,rfRad,lcRad,n).*randseed)/30;
        i=i
    end

     % I=buildInput(inputImg,map,lambda,n);
     % I=I.*randseed;

end


%% 计算神经活动演化
function [act,Spk]=evo(lcRad,n,simFreq,album,IList,a)
 count=0;
 act=zeros(simFreq,n,n);
 act(1,:,:)=rand(n);
 Spk=zeros(simFreq,n,n);
 for t=2:1:simFreq
     prevAct=squeeze(act(t-1,:,:));
     actualSpk=poissrnd(prevAct/simFreq);
     Spk(t,:,:)=actualSpk;
     % bingo=zeros(n+2*lcRad);
     % for x=1:n
     %     for y=1:n
     %         bingo(x:x+2*lcRad,y:y+2*lcRad)=bingo(x:x+2*lcRad,y:y+2*lcRad)+(actualSpk(x,y)*simFreq*squeeze(album(x,y,:,:)));
     %     end
     %     x=x*count
     % end
     bingo=squeeze(sum(album .* (actualSpk*simFreq), [1 2]));
     bingo=bingo(lcRad+1:lcRad+n,lcRad+1:lcRad+n);
     act(t,:,:)=relu(-0.6*prevAct+bingo+squeeze(IList(t,:,:)));
     count=((a-1)*40+t)/40
 end

end

function pos=relu(rate)
    pos=(abs(rate)+rate)/2;
end

function cs=buildRaster(theta,posX,wid,m)
theta_rad = deg2rad(theta);
[H,S] = meshgrid(1:m,1:m);
h_theta =  H * cos(theta_rad) + S * sin(theta_rad);
cs=cos(2*pi*h_theta/(wid*2)+posX);
end
