----------------------------------
-- Proposition 1 (+betti tally) --
----------------------------------

-- n number of cameras
n=2

load "computerVision.m2"

World=worldProductRing2Points(n);
Image=imageProductRing2Points(n);

(R,I)=rigidMultiviewMapGraph(genericCameras(2,7));
J=eliminateWorld(I,World,Image);

betti J     -- print generators by degree
betti res J -- print entire Betti tally
