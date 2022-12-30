clear;
mydefs;
%
%setprngstates(0);
%setprngstates(1);
%setprngstates(22842624);
%setprngstates(71281280);
%%%setprngstates(53119872); % Rather challenging?
%%%sizeX = 1000
setprngstates(0);
sizeX = 50
expVarCoeff = 0.0
%noiseDat = [ 0.0, 0.0, 0.0; 0.0, 0.0, 0.0 ]
noiseDat = [ 1.0e-8, 0.0, 0.0; 0.0, 0.0, 0.0 ];
%noiseDat = [ 1.0e-4, 0.0, 0.0; 0.0, 0.0, 0.0 ];
%noiseDat = [ 1.0e-4, 1.0e-6, 1.0e-6/sizeX; 1.0e-6, 1.0e-8, 1.0e-8/sizeX ];
prm = [];
%prm.iterLimit = 500;
%
vecXSecret = randn(sizeX,1) .* exp(expVarCoeff*abs(randn(sizeX,1)));
%%%vecXSecret = [1:sizeX]';
%fSecret = max([ exp(expVarCoeff*randn()), 0.0 ]);
fSecret = 1.0;
%matSF = diag(exp(expVarCoeff*randn(sizeX,1)));
%matSX = diag(exp(expVarCoeff*randn(sizeX,1)));
%matA = randn(sizeX,sizeX) .* exp(expVarCoeff*randn(sizeX,sizeX));
%matB = matSF * matA * matSX;
%matHSecret = matB' * matB;
matHSecret = full(diag(abs(randn(sizeX,1))));
%%%matHSecret = diag([1:sizeX]);
condHSecret = cond(matHSecret)
%
funchFG = @(x) funcSimpleQuad( x, vecXSecret, fSecret, matHSecret, noiseDat );
funchFG_noiseless = @(x) funcSimpleQuad( x, vecXSecret, fSecret, matHSecret, zeros(2,3) );
vecX0 = zeros(sizeX,1);
[ f0, vecG0 ] = funchFG_noiseless( vecX0 );
%
if (0)
	%matX = full(eye(sizeX,sizeX+1));
	%matX = randn(sizeX,sizeX+1);
	%matX = [ randn(sizeX,sizeX), zeros(sizeX,1) ];
	matX = randn(sizeX,sizeX+1) + [ 1000.0; zeros(sizeX-1,1) ];
	[ rvecF, matG ] = funchFG(matX);
	[ fFit, vecGFit, matHFit ] = hessfit( matX, rvecF, matG );
	matHFit
	return;
endif
numNLS = 10000; % "Noise Level Samples"
[ rvecFNLS, matGNLS ] = funchFG( vecXSecret + zeros(sizeX,numNLS) );
fNoiseLevel = sqrt( sum( (rvecFNLS-fSecret).^2 ) / numNLS )
gNoiseLevel = sqrt( sum(sum( matGNLS.^2, 1 )) / numNLS )
if (1)
	% Reasonable.
	prm.xTol = eps^0.5 * ( norm(vecX0) + norm(vecXSecret) );
	prm.gTol = 3.0*gNoiseLevel + eps^0.7 * norm(vecG0);
	prm.fTol = 10.0*fNoiseLevel + eps^0.7 * f0;
else
	% Only makes sense if fCrit = 0.0; also: aggressive.
	prm.xTol = eps * ( norm(vecX0) + norm(vecXSecret) );
	prm.gTol = 0.1*gNoiseLevel + eps * norm(vecG0);
	prm.fTol = 3.0*fNoiseLevel + eps * f0;
endif
%
prm.funch_vecGSecret = @(x)( matHSecret*(x-vecXSecret) );
prm.funchFGSecret = funchFG;
prm.funchFGNoiselessSecret = funchFG_noiseless;
prm.matHSecret = matHSecret;
prm.learningRate = 0.1;
prm.momentumFactor = 0.9;
%
%[ vecXFin, retCode, datOut ] = mysgdmom( funchFG, vecX0, prm );
[ vecXFin, retCode, datOut ] = sgsolve( funchFG, vecX0, prm );
[ fFin, vecGFin ] = funchFG_noiseless( vecXFin );
xRes = norm(vecXFin-vecXSecret)
xTol = prm.xTol
gRes = norm(vecGFin)
gTol = prm.gTol
fRes = fFin - fSecret
fTol = prm.fTol
