function I=buildInput(img,map,lambda,n)

%根据视觉刺激和方向偏好分布，生成前馈刺激
%输入：
%img- 视觉刺激
%map- 方向偏好图
%lambda- 感受野尺度参数
%n- V1细胞数量参数（边长）

paddedImg = padarray(img, [lambda lambda], 'replicate');
rfAlbum=zeros(n,n,2*lambda+1,2*lambda+1);



theta=0;
for x=1:1:n
    for y=1:1:n
        theta=180*map(x+4*lambda,y+4*lambda);
        rfAlbum(x,y,:,:)=generateGabor(lambda, theta);
    end
end

I=zeros(n);
for x=1:1:n
    for y=1:1:n
        neighbor=paddedImg(x*2:x*2+2*lambda,y*2:y*2+2*lambda);
        I(x,y)=sum(neighbor.* squeeze(rfAlbum(x,y,:,:)),'all');
    end
end
