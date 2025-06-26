function [testPoints,testPointsAct,disTestPoints,corTestPointsAct] = visAcNorm(m,d,r,l)

xBase=linspace(0,m,d);
yBase=linspace(0,m,d);
[X, Y] = meshgrid(xBase, yBase);


testPoints=rand(l,2)*m;


sigmaSq=r;
normConst=1/(2*pi*sigmaSq);


diffX=reshape(testPoints(:,1),[l,1,1])-reshape(X, [1,d,d]);
diffY=reshape(testPoints(:,2),[l,1,1])-reshape(Y, [1,d,d]);
exponent=-(diffX.^2 + diffY.^2)/(2*sigmaSq);
testPointsAct=normConst * exp(exponent);
testPointsAct=reshape(testPointsAct,[l,d*d]);

permIndex = randperm(l);
testPointsPerm=testPoints(permIndex, :);
diffSqTestPoints= (testPoints-testPointsPerm).^2;
disTestPoints=sqrt(sum(diffSqTestPoints, 2));

testPointsActPerm=testPointsAct(permIndex, :);

corTestPointsAct=sum(testPointsAct.*testPointsActPerm,[2]);
end