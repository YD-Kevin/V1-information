function I=buildInput(img,map,rfRad,lcRad,n)

%根据视觉刺激和方向偏好分布，生成前馈刺激
%输入：
%img- 视觉刺激
%map- 方向偏好图
%rfRad- recptive field radius
%n- V1细胞数量参数（边长）

paddedImg = padarray(img, [rfRad rfRad], 'replicate');
rfAlbum=zeros(n,n,2*rfRad+1,2*rfRad+1);



theta=0;
for x=1:1:n
    for y=1:1:n
        theta=180*map(x+lcRad,y+lcRad);
        rfAlbum(x,y,:,:)=generateGabor(rfRad, theta);
    end
end

I=zeros(n);
for x=1:1:n
    for y=1:1:n
        neighbor=paddedImg(x*2:x*2+2*rfRad,y*2:y*2+2*rfRad);
        I(x,y)=sum(neighbor.* squeeze(rfAlbum(x,y,:,:)),'all');
    end
end

