function rfAlbum=buildRFAlbum(map,rfRad,lcRad,n)
rfAlbum=zeros(n,n,2*rfRad+1,2*rfRad+1);


for x=1:1:n
    for y=1:1:n
        theta=180*map(x+lcRad,y+lcRad);
        rfAlbum(x,y,:,:)=generateGabor(rfRad, theta);
    end
end

end