clear;
mydefs;
%
%setprngstates(0);
setprngstates(1);
%setprngstates(22842624);
sizeX = 5
expVarCoeff = 0.0
noiseDat = [ 0.0, 0.0, 0.0; 0.0, 0.0, 0.0 ]
%noiseDat = [ 1.0e-2, 0.0, 0.0; 1.0e-4, 0.0, 0.0 ];
prm = [];
%
vecXSecret = randn(sizeX,1) .* exp(expVarCoeff*abs(randn(sizeX,1)));
%fSecret = max([ exp(expVarCoeff*randn()), 0.0 ]);
fSecret = 0.0;
%matSF = diag(exp(expVarCoeff*randn(sizeX,1)));
%matSX = diag(exp(expVarCoeff*randn(sizeX,1)));
%matA = randn(sizeX,sizeX) .* exp(expVarCoeff*randn(sizeX,sizeX));
%matB = matSF * matA * matSX;
%matHSecret = matB' * matB;
matHSecret = full(diag(abs(randn(sizeX,1))));
condHSecret = cond(matHSecret)
%
funchFG = @(x) funcSimpleQuad( x, vecXSecret, fSecret, matHSecret, noiseDat );
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
%
vecX0 = zeros(sizeX,1);
prm.funch_vecGSecret = @(x)( matHSecret*(x-vecXSecret) );
prm.matHSecret = matHSecret;
%echo__prm = prm
sxsolve1208( funchFG, vecX0, prm );
