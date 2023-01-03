clear;
mydefs;
if ( stopsignalpresent() )
	msg( __FILE__, __LINE__, "ERROR: Stop signal already present." );
	retCode = RETCODE__IMPOSED_STOP;
	return;
endif
%
setprngstates(0);
sizeX = 50
sizeL = min([ sizeX, round(sqrt(sqrt(1E6*sizeX))) ])
vecXCrit = randn(sizeX,1);
%randn(sizeX,1); % Burn to match test_mysgdmom.
fCrit = 1.0;
%cVals = [ 1.0, 1.0E-2, 0.0, 0.0 ]
%cVals = [ 1.0, 1.0E-2, 1.0E-4, 1.0E-6 ]
%cVals = [ 1.0, 1.0E-1, 1.0E-2, 1.0E-3 ]
cVals = [ 0.0, 1.0, 1.0E-2, 1.0E-2 ]
%noisePrm = [ 0.0, 0.0; 0.0, 0.0; 0.0, 0.0 ]
noisePrm = [ 1.0E-12, 1.0E-2; 1.0e-4, 1.0e-4; 1.0e-4, 1.0e-4 ]
%noisePrm = [ 1.0E-8, 0.0; 0.0, 0.0; 0.0, 0.0 ];
%
tic();
msgnnl( __FILE__, __LINE__, "Generating function... " );
matAS = cVals(1)*sparse(eye(sizeX,sizeX)) + cVals(2)*sparse(diag(randn(sizeX,1))) + cVals(3)*sprandn(sizeX,sizeX,sizeL*1.0/sizeX);
matAW = cVals(4)*randn(sizeL,sizeX);
funchFG = @(x) funcQuad1230( x, vecXCrit, fCrit, matAS, matAW, noisePrm );
funchFG_noiseless = @(x) funcQuad1230( x, vecXCrit, fCrit, matAS, matAW, zeros(3,2) );
vecX0 = zeros(sizeX,1);
toc();
%condH = cond( matAS'*matAS + matAW'*matAW )
%
tic();
msgnnl( __FILE__, __LINE__, "Analyzing function... " );
[ f0, vecG0 ] = funchFG_noiseless( vecX0 );
numNLS = 100; % "Noise Level Samples"
[ rvecFNLS, matGNLS ] = funchFG( vecXCrit + zeros(sizeX,numNLS) );
toc();
fNoiseLevel = sqrt( sum( (rvecFNLS-fCrit).^2 ) / numNLS )
gNoiseLevel = sqrt( sum(sum( matGNLS.^2, 1 )) / numNLS )
%
prm = [];
%prm.funch_vecGCrit = @(x)( matH0*(x-vecXCrit) + matA0*(matA0'*(x-vecXCrit)) );
prm.funchFGCrit = funchFG;
prm.funchFGNoiselessCrit = funchFG_noiseless;
%prm.matHCrit = matH0 + matA0*(matA0');
prm.learningRate = 0.1;
prm.momentumFactor = 0.9;
%
prm.xTol = eps^0.9 * ( norm(vecX0) + norm(vecXCrit) );
prm.gTol = 3.0*gNoiseLevel + eps^0.7 * norm(vecG0);
prm.fTol = 3.0*fNoiseLevel + eps^0.7 * f0;
prm.timeLimit = -1;
prm.iterLimit = -1;
prm.progressReportInterval = 0.0;
%
%[ vecXFin, retCode, datOut ] = mysgdmom( funchFG, vecX0, prm );
prm.ledgerLimit = 20;
prm.numFevalPerSuperPt = 20;
[ vecXFin, retCode, datOut ] = sgsolve( funchFG, vecX0, prm );
if ( stopsignalpresent() )
	msg( __FILE__, __LINE__, "Received stop signal." );
	return;
endif
[ fFin, vecGFin ] = funchFG_noiseless( vecXFin );
xRes = norm(vecXFin-vecXCrit)
xTol = prm.xTol
gRes = norm(vecGFin)
gTol = prm.gTol
fRes = fFin - fCrit
fTol = prm.fTol
