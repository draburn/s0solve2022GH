error( "This code is currently unsupported." );
clear;
mydefs;
%
setprngstates(0);
sizeX = 50
densityFactor = 0.0;
expVarCoeff = 0.0;
%
%vecXCrit = randn(sizeX,1);
vecXCrit = randn(sizeX,1) .* exp(expVarCoeff*abs(randn(sizeX,1)));
fCrit = abs( randn() * exp(expVarCoeff*abs(randn())) );
%fCrit = 1.0;

%matH0 = mtm( sprandn(sizeX,sizeX,densityFactor) );
%matH0 = mtm( randn(sizeX,sizeX) )
matH0 = full(diag(abs(randn(sizeX,1))));
condHSecret = cond(matH0)

%matA0 = randn(sizeX,1+floor(densityFactor*sizeX));
matA0 = zeros(sizeX,1);

switch (1)
case (0)
	noiseDat = zeros(5,2);
case (1)
	noiseDat = zeros(5,2);
	noiseDat(1,2) = 1.0e-8;
case (10)
	noiseDat = zeros(5,2);
	noiseDat(2,1) = 1.0e-8;
	noiseDat(3,1) = 1.0e-8;
case (100)
	noiseDat = [ ...
	  1.0e-2, 1.0e-8; ...
	  1.0e-2, 1.0e-8; ...
	  1.0e-2, 1.0e-8; ...
	  1.0e-6, 1.0e-12; ...
	  1.0e-6, 1.0e-12 ];
otherwise
	error( "Invalid case." );
endswitch
noiseDensityMax = densityFactor;
funchFG = @(x) funcQuadCustom1229( x, vecXCrit, fCrit, matH0, matA0, noiseDat, noiseDensityMax );
funchFG_noiseless = @(x) funcQuadCustom1229( x, vecXCrit, fCrit, matH0, matA0, zeros(5,2), 0.0 );
%funchFG = @(x) funcSimpleQuad( x, vecXCrit, fCrit, matH0 );
%funchFG_noiseless = @(x) funcSimpleQuad( x, vecXCrit, fCrit, matH0, zeros(7,7) );
%
vecX0 = zeros(sizeX,1);
[ f0, vecG0 ] = funchFG_noiseless( vecX0 );
%
numNLS = 10000; % "Noise Level Samples"
[ rvecFNLS, matGNLS ] = funchFG( vecXCrit + zeros(sizeX,numNLS) );
fNoiseLevel = sqrt( sum( (rvecFNLS-fCrit).^2 ) / numNLS )
gNoiseLevel = sqrt( sum(sum( matGNLS.^2, 1 )) / numNLS )
if (1)
	% Reasonable.
	prm.xTol = eps^0.5 * ( norm(vecX0) + norm(vecXCrit) );
	prm.gTol = 3.0*gNoiseLevel + eps^0.7 * norm(vecG0);
	prm.fTol = 10.0*fNoiseLevel + eps^0.7 * f0;
else
	% Only makes sense if fCrit = 0.0; also: aggressive.
	prm.xTol = eps * ( norm(vecX0) + norm(vecXCrit) );
	prm.gTol = 0.1*gNoiseLevel + eps * norm(vecG0);
	prm.fTol = 3.0*fNoiseLevel + eps * f0;
endif
%
prm.funch_vecGCrit = @(x)( matH0*(x-vecXCrit) + matA0*(matA0'*(x-vecXCrit)) );
prm.funchFGCrit = funchFG;
prm.funchFGNoiselessCrit = funchFG_noiseless;
prm.matHCrit = matH0 + matA0*(matA0');
prm.learningRate = 0.1;
prm.momentumFactor = 0.9;
%
[ vecXFin, retCode, datOut ] = mysgdmom( funchFG, vecX0, prm );
%prm.ledgerLimit = 10;
%prm.numFevalPerSuperPt = 50;
%[ vecXFin, retCode, datOut ] = sgsolve( funchFG, vecX0, prm );
[ fFin, vecGFin ] = funchFG_noiseless( vecXFin );
xRes = norm(vecXFin-vecXCrit)
xTol = prm.xTol
gRes = norm(vecGFin)
gTol = prm.gTol
fRes = fFin - fCrit
fTol = prm.fTol
