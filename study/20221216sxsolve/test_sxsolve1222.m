clear;
mydefs;
%
%setprngstates(0);
%setprngstates(1);
%setprngstates(22842624);
setprngstates(71281280);
sizeX = 20
expVarCoeff = 0.0
%noiseDat = [ 0.0, 0.0, 0.0; 0.0, 0.0, 0.0 ]
noiseDat = [ 1.0e-6, 0.0, 0.0; 1.0e-8, 0.0, 0.0 ];
%noiseDat = [ 1.0e-2, 0.0, 0.0; 1.0e-4, 0.0, 0.0 ];
prm = [];
%prm.iterLimit = 500;
%
vecXSecret = randn(sizeX,1) .* exp(expVarCoeff*abs(randn(sizeX,1)));
%fSecret = max([ exp(expVarCoeff*randn()), 0.0 ]);
fSecret = 1.0;
%matSF = diag(exp(expVarCoeff*randn(sizeX,1)));
%matSX = diag(exp(expVarCoeff*randn(sizeX,1)));
%matA = randn(sizeX,sizeX) .* exp(expVarCoeff*randn(sizeX,sizeX));
%matB = matSF * matA * matSX;
%matHSecret = matB' * matB;
matHSecret = full(diag(abs(randn(sizeX,1))));
condHSecret = cond(matHSecret)
%
funchFG = @(x) funcSimpleQuad( x, vecXSecret, fSecret, matHSecret, noiseDat );
funchFG_noiseless = @(x) funcSimpleQuad( x, vecXSecret, fSecret, matHSecret, zeros(2,3) );
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
prm.fTol = 30.0*fNoiseLevel;
prm.gTol = 10.0*gNoiseLevel;
prm.fevalLimit = 1000;
%
vecX0 = zeros(sizeX,1);
prm.funch_vecGSecret = @(x)( matHSecret*(x-vecXSecret) );
prm.funchFGSecret = funchFG;
prm.funchFGNoiselessSecret = funchFG_noiseless;
prm.matHSecret = matHSecret;
%echo__prm = prm
[ vecXCalc, retCode, datOut ] = sxsolve1222( funchFG, vecX0, prm );
assert( reldiff(vecXCalc,vecXSecret) < 0.01 );
fTol = prm.fTol
gTol = prm.gTol
fNoiseLevel
gNoiseLevel
