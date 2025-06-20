%% Biological Factors
% macaque

rf=1; %receptive field size by degree 
nd=3000; %2-D neural density in V1 by cell/mm^2 (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
mf=2; %cortical magnification factor degree/mm (Roger B. H. Tootell et.al., 1988)
lcScale=0.55; %a reference scale factor of lateral connection by mm 

%Cat
rf=3; %receptive field size by degree (Gilbert CD and Wiesel 1989, Fuyuki  Karube and  Zoltan F. Kisvarday, 2010)
nd=1800; %2-D neural density in V1 by cell/mm^2 (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
mf=5; %cortical magnification factor degree/mm (Gilbert CD and Wiesel 1989)
lcScale=0.55; %a reference scale factor of lateral connection by mm (Fuyuki  Karube and  Zoltan F. Kisvarday, 2010)

%Tree Shrew
rf=3; %receptive field size by degree (François Mooser et.al 2004)
nd=2000; %2-D neural density in V1 by cell/mm^2 (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
mf=20; %cortical magnification factor degree/mm  (François Mooser et.al 2004)
lcScale=0.2; %a reference scale factor of lateral connection by mm (François Mooser et.al 2004)

%Mouse
rf=10; %receptive field size by degree (Jiakun Fu et.al, 2024)
nd=2200; %2-D neural density in V1 by cell/mm^2 (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
mf=150; %cortical magnification factor degree/mm (Jaeson Jang, Min Song, Se-Bum Paik, 2020)
lcScale=0.03; %a reference scale factor of lateral connection by mm (Peijia Yu et.al. 2025)

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

candidateFuncs = {
    @(x) x,             % T1(x)=x
    @(x) x.^2,          % T2(x)=x^2
    @(x) x.^3,            % T3(x)=x^3
    @(x) log(abs(x)+1e-3), % T4(x)=log(|x|)
    @(x) abs(x)         % T5(x)=|x|
};
maxk=5;

%% Main Prob



%% LinearFIActual
n=64;
m=128;
t=12;

mapIndex=[ 0.05,  0.1 ,  0.2 ,  0.3 ,  0.4 ,  0.5 ,  0.6 ,  0.8 ,  1  ,1.5 ];
a=size(mapIndex);

mapAlbum=zeros(a,lcRad*2+n,lcRad*2+n);



for i=1:1:a
    mapAlbum(i,:,:)=buildMap(mapIndex(i),lcRad*2+n,0,0);
end
 
actSetTheta=zeros(a,t,simFreq,n,n);
SpkSetTheta=zeros(a,t,simFreq,n,n);

for i=1:a
    
    map=squeeze(mapAlbum(i,:,:));
    imgList=zeros(simFreq,m,m);
    
    album=buildAlbum(lcRad,n,lcSigma,map,i);
    
    for theta=0:t-1
        for j=1:simFreq
            imgList(j,:,:)=200*buildRaster(theta*15,j*1.5,9,m);
        end
        IList=input(rfRad,lcRad,n,m,simFreq,imgList,map)*10;
        [act,Spk]=evo(lcRad,n,simFreq,album,IList,(i-1)*t+theta+1);
        actSetTheta(i,theta+1,:,:,:)=act;
        SpkSetTheta(i,theta+1,:,:,:)=Spk;
    end
end
LFIAveEachDensity=LinearFIActual(actSetTheta,a,t,n);







%% LinearFISimple
% mapAlbum=zeros(10,lcRad*2+n,lcRad*2+n);
%  mapFixAlbum=zeros(10,n,n);
%  mapIndex=[ 0.05,  0.1 ,  0.2 ,  0.3 ,  0.4 ,  0.5 ,  0.6 ,  0.8 ,  1  ,1.5 ];
% 
%  actAlbum=zeros(10,simFreq,n,n);
%  SpkAlbum=zeros(10,simFreq,n,n);
% 
%  for a=1:1:10
%      mapAlbum(a,:,:)=buildMap(mapIndex(a),lcRad*2+n,0,0);
%  end
% 
% 
%  imgList=zeros(simFreq,m,m);
% for i=1:simFreq
%     imgList(i,:,:)=200*buildRaster(0,i*1.5,9,m);
% end
% 
% 
% 
% for a=1:1:10
%     map=squeeze(mapAlbum(a,:,:)); %生成偏好图
%     album=buildAlbum(lcRad,n,lcSigma,map,a);
% 
%     IList=input(rfRad,lcRad,n,m,simFreq,imgList,map);
% 
%     [act,Spk]=evo(lcRad,n,simFreq,album,IList,a);
%     actAlbum(a,:,:,:)=act;
%     SpkAlbum(a,:,:,:)=Spk;
% end
% 
% LFI=zeros(10,n,n);
% for a=1:10
%     act=squeeze(actAlbum(a,:,:,:));
%     map=squeeze(mapAlbum(a,:,:));
%     LFI(a,:,:)=LinearFISimple(act,map,lcRad,n);
% end
% plot(mean(LFI,[2,3]));


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

function IList=input(rfRad,lcRad,n,m,simFreq,imgList,map)


    randIndex=mean(abs(imgList),'all')/5;


    for i=1:1:simFreq
        imgList(i,:,:)=squeeze(imgList(i,:,:))+randIndex*randn(m);
    end


    %根据图片输入，绘制前馈刺激
    rfAlbum=buildRFAlbum(map,rfRad,lcRad,n);

    IList=zeros(simFreq,n,n);
    for i=1:1:simFreq
        IList(i,:,:)=(buildInput(squeeze(imgList(i,:,:)),rfAlbum,rfRad,n));
        i=i
    end
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
