worldProductRing2Points = (n) -> (
    --creates the ring (PP^2xPP^2)^n2 of Image coordinates
    -- n ia thew number of cameras
    -- m is the number of points
    --  the variable names of our paper [JKSW]
    World:=QQ[X_0..X_3]**QQ[Y_0..Y_3]
    )

worldProductRing = (m,n) -> (
    --creates the ring (PP^2xPP^2)^n2 of Image coordinates
    -- n ia thew number of cameras
    -- m is the number of points
    --  the variable names of our paper [JKSW]
    if (m==2) then (
    	World:=worldProductRing2Points (2)
    	)
    else (
	World=QQ[X_(1,0)..X_(1,3)];
    	for i from 2 to m do(
	    World=World**QQ[X_(i,0)..X_(i,3)];
    	    );
    	);
    World
    )

imageProductRing2Points = (n) -> (
    --creates the ring (PP^2xPP^2)^n2 of Image coordinates
    -- n ia thew number of cameras
    -- m is the number of points
    --  the variable names of our paper [JKSW]
    Image:=QQ[u_(1,0)..u_(1,2)];
    for i from 2 to n do(
    	Image=Image**QQ[u_(i,0)..u_(i,2)]
    	);
    for i from 1 to n do(
    	Image=Image**QQ[v_(i,0)..v_(i,2)]
    	);
    Image
    )

imageProductRing = (m,n) -> (
    --creates the ring (PP^2xPP^2)^nm of Image coordinates
    -- n ia thew number of cameras
    -- m is the number of points
    -- if m==2 then  the variable names of our paper [JKSW]
    -- else u_(m,n,i) is used
    if (m==2) then (
    	Image := imageProductRing2Points(n)
    	)
    else (
    	for i from 1 to m do(
	    for j from 1 to n do(
	    	if (i==1 and j==1) then (Image = QQ[u_(1,1,0)..u_(1,1,2)];)
    	    	else Image = Image ** QQ[u_(i,j,0)..u_(i,j,2)];
	    	);
    	    );
    	);
    Image
    )

genericCameras = (n,h) -> (
    --creats n random cameras with integer entries bounded by h
    for i from 1 to n do(Camera_i=random(ZZ^3,ZZ^4,Height=>h););
    -- check if the cameras are generic 
    MList:={};
    for i from 1 to n do(MList=MList|entries Camera_i;);
    while (numgens minors(4,transpose matrix MList)=!=binomial(n*3,4)) do(
    	for i from 1 to n do(Camera_i=random(ZZ^3,ZZ^4,Height=>h););
    	MList={};
    	for i from 1 to n do(MList=MList|entries Camera_i;);
    	);
    A:={};
    for i from 1 to n do(A=A|{Camera_i};);
    A
    )

multiviewMapGraph = (A,m) -> (
    -- create the Graph of the Multiview Map and it Ring
    -- A is the List of input matrices 
    -- m number of points
    -- n number of cameras
    n:=length A;
    if (m==2) then(
	Image := imageProductRing2Points(n);
	World := worldProductRing2Points(n);
	-- creat the ring in which the multiview map graph lives
    	R = World ** Image;
    	--creates the ideal of dehomogenized equation AX=lu
    	I:=ideal();
    	for j from 1 to n do(
    	    I=I+minors(2,A_(j-1)*(genericMatrix(R,X_0,3,1)||matrix{{1}})|genericMatrix(R,u_(j,0),3,1));
	    I=I+minors(2,A_(j-1)*(genericMatrix(R,Y_0,3,1)||matrix{{1}})|genericMatrix(R,v_(j,0),3,1));
    	    );
 	)
    else (
    	Image = imageProductRing(m,n);
    	World = worldProductRing(m,n);
    	R = World ** Image;
    	I = ideal();
    	for i from 1 to m do(
	    for j from 1 to n do(
    	    	I=I+minors(2,(A_(j-1))*(genericMatrix(R,X_(i,0),3,1)||matrix{{1}})|genericMatrix(R,u_(i,j,0),3,1));
    	    	);
    	    );
    	);
    (R,I)
    )


transformQ = (Q,h)->(
    --transforms the the constraint Q(X,Y) by a random projective 
    --tranfomation with entry height h "M=projectiveTrafo"(ZZ^(4x4)) to Q(MX,MY) and
    -- dehomogenizes 
    -- Input: Q polynomial in World coordinates, h random entry height
    --Output: tranfrom dehomogenized Q
    projectiveTrafo=random(ZZ^4,ZZ^4,Height=>h);
    trafoX:=matrix{{X_0,X_1,X_2,X_3}}*projectiveTrafo;
    trafoY:=matrix{{Y_0,Y_1,Y_2,Y_3}}*projectiveTrafo;
    trafoQ:=substitute(Q, {
	    X_0=>trafoX_0_0,X_1=>trafoX_1_0,X_2=>trafoX_2_0,X_3=>trafoX_3_0,
	    Y_0=>trafoY_0_0,Y_1=>trafoY_1_0,Y_2=>trafoY_2_0,Y_3=>trafoY_3_0});
    dehomTrafoQ:=substitute(trafoQ, {X_3=>1,Y_3=>1})
)


rigidMultiviewMapGraph = (A) -> (
    --creates the rigid multiview variety of 2 points and n cameras from [JKSW]
    -- A is the List of input matrices 
    -- n number of cameras 
    (R,I):=multiviewMapGraph (A,2);
    I=I+ideal((X_0-Y_0)^2+(X_1-Y_1)^2+(X_2-Y_2)^2-1);
    (R,I) 
    )

eliminateWorld = (I,World,Image) -> (
    --eliminate the World coordinates and map to image ring
    R := ring I;
    WorldToR:= map(R,World);
    J:=eliminate(apply((entries vars World)_0,i -> WorldToR i),I);
    RToImage:=map(Image,R);
    RToImage J
    )

cramer = method()
cramer(ZZ, Matrix) := Matrix => (i, M) -> (
transpose matrix{{det submatrix'(M,{i-1},{0}), - det submatrix'(M,{i-1},{1}), det submatrix'(M,{i-1},{2}), - det submatrix'(M,{i-1},{3})}}    
)    

makeQuadrilinear = method()
makeQuadrilinear RingElement := Thing => Q -> (
local T;
T = symbol T;
for i from 1 to 4 do(
    for j from 1 to 4 do(
	for k from 1 to 4 do(
	    for l from 1 to 4 do( 
		T_(i,j,k,l) = (1/#set{i,j})*(1/#set{k,l})*coefficient(X_(i-1)*X_(j-1)*Y_(k-1)*Y_(l-1), Q)
		);	
	);
    );
    
);    
T
)


contractQuadrilinear = method()
contractQuadrilinear(IndexedVariableTable, List) := RingElement => (T,L) -> (
    local sum;
    sum = 0;
    for i from 1 to 4 do(
	for j from 1 to 4 do(
	    for k from 1 to 4 do(
		for l from 1 to 4 do(
		sum = sum + (T_(i,j,k,l))*((L#0)_(i-1,0))*((L#1)_(j-1,0))*((L#2)_(k-1,0))*((L#3)_(l-1,0));
		    );
		);
	    );
	);
    sum
    )

 
