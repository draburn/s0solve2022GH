clear;
mydefs;
%
setprngstates(0);
sizeX = 50
vecXCrit = randn(sizeX,1);
fCrit = 1.0;
randn(sizeX,1); % Burn.
matAS = diag(sqrt(abs(randn(sizeX,1))));
matAW = [];
%condHS = cond(matAS'*matAS)
%
%noisePrm = [ 0.0, 0.0; 1.0e-4, 1.0e-4; 1.0e-4, 1.0e-4 ];
noisePrm = [ 1.0E-8, 0.0; 0.0, 0.0; 0.0, 0.0 ];
funchFG = @(x) funcQuad1230( x, vecXCrit, fCrit, matAS, matAW, noisePrm );
funchFG_noiseless = @(x) funcQuad1230( x, vecXCrit, fCrit, matAS, matAW, zeros(3,2) );
%
vecX0 = zeros(sizeX,1);
%
[ f0, vecG0 ] = funchFG_noiseless( vecX0 );
numNLS = 100; % "Noise Level Samples"
[ rvecFNLS, matGNLS ] = funchFG( vecXCrit + zeros(sizeX,numNLS) );
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
prm.xTol = eps^0.7 * ( norm(vecX0) + norm(vecXCrit) );
prm.gTol = 3.0*gNoiseLevel + eps^0.7 * norm(vecG0);
prm.fTol = 3.0*fNoiseLevel + eps^0.7 * f0;
%
%[ vecXFin, retCode, datOut ] = mysgdmom( funchFG, vecX0, prm );
%prm.ledgerLimit = 10;
%prm.numFevalPerSuperPt = 50;
[ vecXFin, retCode, datOut ] = sgsolve( funchFG, vecX0, prm );
[ fFin, vecGFin ] = funchFG_noiseless( vecXFin );
xRes = norm(vecXFin-vecXCrit)
xTol = prm.xTol
gRes = norm(vecGFin)
gTol = prm.gTol
fRes = fFin - fCrit
fTol = prm.fTol
