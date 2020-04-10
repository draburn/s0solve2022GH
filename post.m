iterToViz = 3
studyPt_prm.vecXSecret = funcPrm.x0;
tic();
this_funchF = grootDatOut.funchF;
this_vecX = grootDatOut.iterDat(iterToViz).vecX;
this_vecV = eye(sizeX,sizeX);
this_vecW = prm.funchJ(this_vecX);
p = studyPt( this_funchF, this_vecX, this_vecV, this_vecW, studyPt_prm );
toc();
