clear;
numFigs = 0;
setprngstates(0);
tic();
%
sizeX = 10;
sizeF = sizeX;
numElemPerCol = 0;
numAddtlElemPerRow = 5;
cEye = 5.0;
c0 = 0.0;
csx = 0.0;
csf = 0.0;
%
matJ0 = cEye*eye(sizeF,sizeX);
for n=1:sizeX
for t=1:numElemPerCol
	m = ceil( sqrt(eps) + (sizeF-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = randn();
endfor
endfor
for m=1:sizeF
for t=1:numAddtlElemPerRow
	n = ceil( sqrt(eps) + (sizeX-2.0*sqrt(eps))*rand() );
	matJ0(m,n) = randn();
endfor
endfor
matJ0 += c0 * randn(sizeF,sizeX);
assert( isrealarray(matJ0,[sizeF,sizeX]) );

matSX = diag(exp(csx*randn(sizeX,1)));
matSF = diag(exp(csf*randn(sizeF,1)));
matJ = matSF*matJ0/matSX;
%

	funchW = @(v)( matJ*v );
	%
	vecX_secret = randn(sizeX,1);
	vecF = matJ*vecX_secret;
	%
	linsolf_prm = [];
	linsolf_prm.tol = 0.01;
	[ vecX, datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
	assert( reldiff(matJ*vecX,-vecF) < 1.1*linsolf_prm.tol );
	%
	%
	vecX_secret = randn(sizeX,1);
	vecF = matJ*vecX_secret;
	%
	sss_prm = [];
	sss_prm.matVR = datOut.matV;
	sss_prm.matWR = datOut.matW;
	sss_prm.tol = 0.01;
	[ vecX, datOut ] = splitspacesolf( funchW, -vecF, sizeX, sss_prm );
	assert( reldiff(matJ*vecX,-vecF) < 1.1*sss_prm.tol );
