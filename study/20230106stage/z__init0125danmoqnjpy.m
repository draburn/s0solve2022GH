function [ funchFG, vecX0, solverPrm, datOut ] = z__init0113( initPrm=[] )
	datOut = [];
	datOut.prngstates = setprngstates(0);
	%%%datOut.prngstates = setprngstates();
	%
	sizeX = 5
	datOut.sizeX = sizeX;
	%
	sizeL = sizeX
	datOut.sizeL = sizeL;
	%
	msgnnl( __FILE__, __LINE__, "Generating function..." );
	tic();
	vecXCrit = ones( sizeX, 1 );
	fCrit = 10.0
	matAS = [];
	matAW = diag([1:sizeX]);
	%noisePrm = [ 1.0E-6, 0.0; 0.0, 0.0; 0.0, 0.0 ]
	noisePrm = [ 0.0, 0.0; 0.0, 0.0; 0.0, 0.0 ]
	toc();
	datOut.vecXCrit = vecXCrit;
	datOut.fCrit = fCrit;
	%
	funchFG = @(x) funcQuad1230( x, vecXCrit, fCrit, matAS, matAW, noisePrm );
	funchFG_noiseless = @(x) funcQuad1230( x, vecXCrit, fCrit, matAS, matAW, zeros(3,2) );
	datOut.funchFG_noiseless = funchFG_noiseless;
	vecX0 = zeros(sizeX,1);
	
	if (0)
		matH = full( matAS'*matAS + matAW'*matAW );
		eig(matH)'
		error( "HALT!" );
	endif
	
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
	solverPrm.learningRate = 0.1; % Unless...
	solverPrm.learningRate = 0.01; % Overwrite?
	solverPrm.momentumFactor = 0.9;
	%
	solverPrm.numFevalPerSuperPt = numStudyPts;
	solverPrm.xTol = (eps^0.6)*norm(vecX0) + (eps^0.8)*fAvg/(norm(vecGAvg)+gVar);
	solverPrm.fTol = (eps^0.4)*fVar + (eps^0.6)*fAvg;
	solverPrm.gTol = (eps^0.4)*gVar + (eps^0.6)*norm(vecGAvg);
	solverPrm.fBail = max(rvecFStudy)/eps;
	solverPrm.fevalLimit = 1000000;
	solverPrm.iterLimit = 500;
	solverPrm.timeLimit = 600.0;
	solverPrm.stopSignalCheckInterval = 1.0;
	solverPrm.progressReportInterval = 1.0; % Unless...
	solverPrm.progressReportInterval = 0.0; % Overwrite?
	%
	solverPrm.bestFVarCoeffA = 2.0;
	solverPrm.bestFVarCoeffB = 2.0;
	%
	solverPrm.maxNumRecords = 20; % Unless...
	%solverPrm.maxNumRecords = sizeX;  % Overwrite?
	solverPrm.useQNJ = false; % Unless...
	%solverPrm.useQNJ = true; % Overwrite?
	%
	solverPrm.qnj_basisDropThresh = sqrt(eps); % Risky! So, consier...
	solverPrm.qnj_basisDropThresh = 0.1; % Overwrite?
	solverPrm.qnj_sMaxInit = 3.0;
	solverPrm.qnj_sMax_limitOnGood = true; % Fix 2023-01-15. Formerly called 'pre-limit'.
	solverPrm.qnj_sMaxBT = 0.1;
	solverPrm.qnj_sMaxFT = 2.0;
	%solverPrm.qnj_sMaxLo = 1.0E-7;
	solverPrm.qnj_sMaxLo = [];
	solverPrm.qnj_sMaxHi = 10.0;
	solverPrm.qnj_dMaxInit = [];
	solverPrm.qnj_dMax_limitOnGood = true; % Fix 2023-01-15. Formerly called 'pre-limit'.
	solverPrm.qnj_dMaxBT = 0.1;
	solverPrm.qnj_dMaxFT = 2.0;
	solverPrm.qnj_dMaxLo = 10.0 * solverPrm.xTol;
	solverPrm.qnj_dMaxHi = [];
	%
	solverPrm.qnj_forceGradAfterBad = true;
	solverPrm.qnj_useCtsFromHarvest = true;
return;
endfunction
