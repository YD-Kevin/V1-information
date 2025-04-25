function profileTotal=buildCoProfile(lambda,sigma,map)

profileSE=buildProfileSE(lambda,sigma);
profileSI=buildProfileSI(lambda,sigma);
profileLE=buildProfileLE(lambda,sigma,map);

profileTotal=(profileSE+profileLE)/2-(profileSI);

%% 短程兴奋性作用
function profileSE=buildProfileSE(lambda,sigma)
    n=8*lambda+1;
    c=4*lambda+1; 
    [X, Y]=meshgrid(1:n); 
    A=sqrt((X-c).^2+(Y-c).^2); 
    A(c,c)=0; 

    profileSE=normpdf(A,0,sigma/4);
    profileSE=profileSE/sum(profileSE,"all");
end


%% 短程抑制性作用
function profileSI=buildProfileSI(lambda,sigma)
    n=8*lambda+1;
    c=4*lambda+1; 
    [X, Y]=meshgrid(1:n); 
    A=sqrt((X-c).^2+(Y-c).^2); 
    A(c,c)=0; 

    profileSI=normpdf(A,0,sigma/2);
    profileSI=profileSI/sum(profileSI,"all");
end

%% 长程兴奋性作用
function profileLE=buildProfileLE(lambda,sigma,map)
    n=8*lambda+1;
    c=4*lambda+1; 
    [X, Y]=meshgrid(1:n); 
    A=sqrt((X-c).^2+(Y-c).^2); 
    A(c,c)=0; 

    profileLE=normpdf(A,sigma,sigma/6);

    clock=abs(map-map(c,c));
    cotClock=abs(map-pi-map(c,c));
    thetaDif=min(clock,cotClock);

    profileLE=thetaDif.*profileLE;
    profileLE=profileLE/sum(profileLE,"all");

end


end
