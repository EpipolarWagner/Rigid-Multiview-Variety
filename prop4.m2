-------------------
-- Proposition 4 --
-------------------

-- n = number of cameras
n=4

-- m = number of points
m=2

load "computerVision.m2"

World=worldProductRing2Points(n);
Image=imageProductRing2Points(n);

Q=(X_0*Y_3-X_3*Y_0)^2+(X_1*Y_3-X_3*Y_1)^2+(X_2*Y_3-X_3*Y_2)^2-X_3^2*Y_3^2
dehomTrafoQ=transformQ(Q,10);

A_1=matrix{{1,0,0,0},{0,1,0,0},{0,0,1,0}}
A_2=matrix{{1,0,0,0},{0,1,0,0},{0,0,0,1}}
A_3=matrix{{1,0,0,0},{0,0,1,0},{0,0,0,1}}
A_4=matrix{{0,1,0,0},{0,0,1,0},{0,0,0,1}}
M={A_1,A_2,A_3,A_4};

(R,I)=multiviewMapGraph(M,m);
F=map(R,World);
I=I+F dehomTrafoQ;


time J=eliminateWorld(I,World,Image);

betti J -- print first column of the betti tally of J
