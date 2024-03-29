function [ funchFG, vecX0, solverPrm, datOut ] = z__init0112( initPrm=[] )
	datOut = [];
	datOut.prngstates = setprngstates(0);
	%
	sizeX = 1E1
	datOut.sizeX = sizeX;
	%
	switch ( 0 )
	case 0
		sizeL = 1
	case 10
		sizeL = min([ sizeX, round(0.1*sqrt(sizeX)) ])
	otherwise
		error( "Invalid case." );
	endswitch
	datOut.sizeL = sizeL;
	%
	switch ( 10 )
	case 0
		cVals = [ 1.0, 0.0, 0.0, 0.0 ]
	case 10
		cVals = [ 1.0, 1.0E-2, 1.0E-4, 1.0E-4 ]
	otherwise
		error( "Invalid case." );
	endswitch
	datOut.cVals = cVals;
	%
	switch ( 10 )
	case 0
		noisePrm = [ 0.0, 0.0; 0.0, 0.0; 0.0, 0.0 ]
	case 10
		noisePrm = [ 1.0E-10, 1.0E-4; 1.0E-6, 1.0E-6; 1.0E-6, 1.0E-6 ]
	otherwise
		error( "Invalid case." );
	endswitch
	datOut.noisePrm = noisePrm;
	%
	%
	msgnnl( __FILE__, __LINE__, "Generating function..." );
	tic();
	vecXCrit = randn(sizeX,1);
	fCrit = 0.0;
	matAS = ...
	   cVals(1)*sparse(eye(sizeX,sizeX)) ...
	 + cVals(2)*sparse(diag(randn(sizeX,1))) ...
	 + cVals(3)*sprandn(sizeX,sizeX,sizeL*1.0/sizeX);
	matAW = cVals(4)*randn(sizeL,sizeX);
	toc();
	datOut.vecXCrit = vecXCrit;
	datOut.fCrit = fCrit;
	%
	funchFG = @(x) funcQuad1230( x, vecXCrit, fCrit, matAS, matAW, noisePrm );
	funchFG_noiseless = @(x) funcQuad1230( x, vecXCrit, fCrit, matAS, matAW, zeros(3,2) );
	datOut.funchFG_noiseless = funchFG_noiseless;
	vecX0 = zeros(sizeX,1);
	%
	%
	numStudyPts = 100
	msgnnl( __FILE__, __LINE__, "Studying function... " );
	tic();
	[ rvecFStudy, matGStudy ] = funchFG( vecX0 + zeros(sizeX,numStudyPts) );
	toc();
	fAvg = sum( rvecFStudy ) / numStudyPts;
	rvecFDevi = rvecFStudy - fAvg; % Scalar auto-broadcast.
	fVar = sqrt( sum( rvecFDevi.^2 ) / numStudyPts )
	vecGAvg = sum( matGStudy, 2 ) / numStudyPts;
	matGDevi = matGStudy - vecGAvg; % Auto-broadcast.
	vecGVar = sqrt( sum( matGDevi.^2, 2 ) / numStudyPts );
	gVar = norm( vecGVar )
	datOut.study_fAvg = fAvg;
	datOut.study_fVar = fVar;
	datOut.study_vecGAvg = vecGAvg;
	datOut.study_vecGVar = vecGVar;
	datOut.study_gVar = gVar;
	%
	%
	solverPrm = [];
	%
	solverPrm.learningRate = 0.1;
	solverPrm.momentumFactor = 0.9;
	%
	solverPrm.numFevalPerSuperPt = numStudyPts;
	solverPrm.xTol = (eps^0.6)*norm(vecX0) + (eps^0.8)*fAvg/(norm(vecGAvg)+gVar);
	solverPrm.fTol = (eps^0.4)*fVar + (eps^0.6)*fAvg;
	solverPrm.gTol = (eps^0.4)*gVar + (eps^0.6)*norm(vecGAvg);
	solverPrm.fBail = max(rvecFStudy)/eps;
	solverPrm.fevalLimit = -1;
	solverPrm.iterLimit = 20000;
	solverPrm.timeLimit = 600.0;
	solverPrm.stopSignalCheckInterval = 1.0;
	solverPrm.progressReportInterval = 1.0;
	%
	solverPrm.maxNumRecords = 50;
	solverPrm.useQNJ = false;
return;
endfunction
