----------------------------------
-- Unlabeled Multiview Variety --
----------------------------------

-- n = number of cameras
n=2

-- m = number of points
m=2

load "computerVision.m2"

World=worldProductRing2Points(n);
Image=imageProductRing2Points(n);
ChowRing=QQ[a00,a01,a02,a11,a12,a22]**QQ[b00,b01,b02,b11,b12,b22]
MultipleRing=QQ[l0,l1,k0,k1]

Q=(X_0*Y_3-X_3*Y_0)^2+(X_1*Y_3-X_3*Y_1)^2+(X_2*Y_3-X_3*Y_2)^2-X_3^2*Y_3^2
dehomTrafoQ=transformQ(Q,2);
dehomTrafoQ:=substitute(Q, {X_3=>1,Y_3=>1})

A_1=matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0}}
A_2=matrix{{1,0,0,0},{0,1,0,0},{0,0,0,1}}

M={A_1,A_2};

(R,I)=multiviewMapGraph(M,m);
F=map(R,World);
I=I+F dehomTrafoQ;

S=ChowRing**Image ** MultipleRing

time J=eliminateWorld(I,World,S);

ChowCoordA=l0*matrix{{a00,a01,a02},
    {0,a11,a12},
    {0,0,a22}}-matrix{
    {2*u_(1,0)*v_(1,0),u_(1,1)*v_(1,0)+u_(1,0)*v_(1,1),u_(1,2)*v_(1,0)+u_(1,0)*v_(1,2)},
    {0,2*u_(1,1)*v_(1,1),u_(1,2)*v_(1,1)+u_(1,1)*v_(1,2)},
    {0,0,2*u_(1,2)*v_(1,2)}}
ChowCoordB= k0*matrix{{b00,b01,b02},
    {0,b11,b12},
    {0,0,b22}}-matrix{
    {2*u_(2,0)*v_(2,0),u_(2,1)*v_(2,0)+u_(2,0)*v_(2,1),u_(2,2)*v_(2,0)+u_(2,0)*v_(2,2)},
    {0,2*u_(2,1)*v_(2,1),u_(2,2)*v_(2,1)+u_(2,1)*v_(2,2)},
    {0,0,2*u_(2,2)*v_(2,2)}}

K=ideal(J,ChowCoordA,ChowCoordB,l0*l1-1,k0*k1-1);

time K=eliminateWorld(K,Image ** MultipleRing,ChowRing);


