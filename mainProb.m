 lambda=5;
 sigma=5;
 n=128;
 time=40;

%% 
% mapAlbum=zeros(10,lambda*8+n,lambda*8+n);
% 
% for a=1:1:10
%     mapAlbum(a,:,:)=buildMap(3/a-0.29,lambda*8+n);
% end
% 
% actAlbum=zeros(10,40,n,n);
% 
% 
% inputImg = double(imread('cameraman.tif'));%图片输入
% actAlbum=zeros(10,time,n,n);
read1=read;
readCount1=readCount;

map=squeeze(mapAlbum(5,:,:));

IList=input(lambda,n,time,inputImg,map);

[read,readCount]=evo(lambda,n,time,album,IList,a);

% for a=1:1:10
%     map=squeeze(mapAlbum(a,:,:)); %生成偏好图
%     album=buildAlbum(lambda,n,sigma,map,a);
%     IList=input(lambda,n,time,inputImg,map);
%     act=evo(lambda,n,time,album,IList,a);
%     actAlbum(a,:,:,:)=act;
% end



%% 建立互相关场图册

function album=buildAlbum(lambda,n,sigma,map,a)

    album=zeros(n,n,lambda*8+n,lambda*8+n); 
    subMap=zeros(lambda*8+1);
    
    %根据偏好图绘制各个神经元对周边神经元的作用分布图
    for x=1:1:n
        for y=1:1:n
    
            subMap=map(x:x+lambda*8,y:y+lambda*8);
            subMap=buildCoProfile(lambda,sigma,subMap);
    
            chart=zeros(lambda*8+n);
            chart(x:x+lambda*8,y:y+lambda*8)=subMap;
    
            album(x,y,:,:)=chart;
            count=((a-1)*168+x)/1680
        end
    end

end

%% 计算前馈刺激
 
 function IList=input(lambda,n,time,inputImg,map)
% [X, Y]=meshgrid(1:256);
% SX=sin(X/2)*100;
% SY=sin(Y/2)*100;
% 
%  temImg=ones(256)*30;
%   I=buildInput(SY,map,lambda,n);

    randseed=randn(n);
    randseed=randseed./abs(randseed);
    
    imgList=zeros(time,256,256);
    
    randIndex=sum(abs(inputImg),'all')/(256*256*5);
    
    
    for i=1:1:time
        imgList(i,:,:)=inputImg+randIndex*randn(256);
    end
    
    
    
    %根据图片输入，绘制前馈刺激
    IList=10*ones(time,n,n);
    for i=1:1:time
        IList(i,:,:)=(buildInput(squeeze(imgList(i,:,:)),map,lambda,n).*randseed)/30;
    end
    
     % I=buildInput(inputImg,map,lambda,n);
     % I=I.*randseed;

end


%% 计算神经活动演化
function [act,actCount]=evo(lambda,n,time,album,IList,a)

 act=zeros(time,n,n);
 actCount=zeros(n,n);
 act(1,:,:)=rand(n);
 for t=2:1:time
     tem=squeeze(act(t-1,:,:));
     musk=poissrnd(tem/time);
     actCount=actCount+musk;

     bingo=squeeze(sum(album .* (musk*time), [1 2]));
     bingo=bingo(4*lambda+1:4*lambda+n,4*lambda+1:4*lambda+n);
     act(t,:,:)=relu(-0.6*tem+bingo+squeeze(IList(t,:,:)));
     count=((a-1)*168+128+t)/1680
 end
end

function pos=relu(rate)
    pos=(abs(rate)+rate)/2;
end

function cs=buildRaster(theta,posX)
theta_rad = deg2rad(theta);
[H,S] = meshgrid(1:256,1:256);
h_theta =  H * cos(theta_rad) + S * sin(theta_rad);
cs=cos(2*pi*h_theta/8+posX);
end