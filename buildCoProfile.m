function profileTotal=buildLCProfile(lcRad,lcSigma,map)

profileSE=buildProfileSE(lcRad,lcSigma);
profileSI=buildProfileSI(lcRad,lcSigma);
profileLE=buildProfileLE(lcRad,lcSigma,map);

profileTotal=(profileSE+profileLE)/2-(profileSI);

%% 短程兴奋性作用
function profileSE=buildProfileSE(lcRad,lcSigma)
    n=2*lcRad+1;
    c=lcRad+1; 
    [X, Y]=meshgrid(1:n); 
    A=sqrt((X-c).^2+(Y-c).^2); 
    A(c,c)=0; 

    profileSE=normpdf(A,0,lcSigma/4);
    profileSE=profileSE/sum(profileSE,"all");
end


%% 短程抑制性作用
function profileSI=buildProfileSI(lcRad,lcSigma)
    n=2*lcRad+1;
    c=lcRad+1; 
    [X, Y]=meshgrid(1:n); 
    A=sqrt((X-c).^2+(Y-c).^2); 
    A(c,c)=0; 

    profileSI=normpdf(A,0,lcSigma/2);
    profileSI=profileSI/sum(profileSI,"all");
end

%% 长程兴奋性作用
function profileLE=buildProfileLE(lcRad,lcSigma,map)
    n=2*lcRad+1;
    c=lcRad+1; 
    [X, Y]=meshgrid(1:n); 
    A=sqrt((X-c).^2+(Y-c).^2); 
    A(c,c)=0; 

    profileLE=normpdf(A,lcSigma,lcSigma/6);

    clock=abs(map-map(c,c));
    cotClock=abs(map-pi-map(c,c));
    thetaDif=min(clock,cotClock);

    profileLE=thetaDif.*profileLE;
    profileLE=profileLE/sum(profileLE,"all");

end

end
