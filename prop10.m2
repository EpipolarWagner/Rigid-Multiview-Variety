--------------------
-- Proposition 10 --
--------------------

-- n = number of cameras
n=2

-- m = number of points
m=4

load "computerVision.m2"

World=worldProductRing(m,n);
Image=imageProductRing(m,n);

(R,I)=multiviewMapGraph(genericCameras(2,3),4);

-- planarity constraint
I=I+ideal(det(
	(genericMatrix(R,X_(1,0),3,1)||1)|
	(genericMatrix(R,X_(2,0),3,1)||1)|
	(genericMatrix(R,X_(3,0),3,1)||1)|
	(genericMatrix(R,X_(4,0),3,1)||1)));

time J=eliminateWorld(I,World,Image);
betti J -- print first column of the betti tally of J

