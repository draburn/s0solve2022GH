% See Notes 2022-05-30-2100.
clear;
numFigs = 0;
setprngstates(0);
tic();
%
sizeX = 100;
sizeF = sizeX;
numElemPerCol = 0;
numAddtlElemPerRow = 0;
cEye = 1.0;
c0 = 0.1;
csx = 1.0;
csf = 1.0;
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


funchW = @(v)( matJ*v );
%
numRuns = 100;
matJA_compound = eye(sizeF,sizeX);
fevalVals_compound = zeros(1,numRuns);
matVTot_compound = zeros(sizeX,1);
kTotVals_compound = zeros(1,numRuns);
fevalVals_splitspace = zeros(1,numRuns);
kTotVals_splitspace = zeros(1,numRuns);
for n=1:100
	vecX_secret = randn(sizeX,1);
	vecF = matJ*vecX_secret;
	%
	linsolf_prm = [];
	if ( rcond(matJA_compound) < sqrt(eps) )
		msg( __FILE__, __LINE__, sprintf( "rcond(matJA_compound) = %g.", rcond(matJA_compound) ) );
		break;
	endif
	linsolf_prm.matP = inv(matJA_compound);
	linsolf_prm.tol = 0.1;
	[ vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
	assert( reldiff(matJ*vecX,-vecF) < 1.1*linsolf_prm.tol );
	fevalVals_compound(n) = linsolf_datOut.fevalCount;
	matJA_compound += ( linsolf_datOut.matW - (matJA_compound*linsolf_datOut.matV) ) * ( linsolf_datOut.matV' );
	matVTot_compound = utorthdrop( [matVTot_compound, linsolf_datOut.matV] );
	kTotVals_compound(n) = size(matVTot_compound,2);
	%
	if ( 1==n )
		fevalVals_splitspace(n) = fevalVals_compound(n);
		kTotVals_splitspace(n) = kTotVals_compound(n);
		matVR = linsolf_datOut.matV;
		matWR = linsolf_datOut.matW;
	else
		sss_prm = [];
		sss_prm.matVR = matVR;
		sss_prm.matWR = matWR;
		sss_prm.tol = linsolf_prm.tol;
		[ vecX, sss_datOut ] = splitspacesolf( funchW, -vecF, sizeX, sss_prm );
		assert( reldiff(matJ*vecX,-vecF) < 1.1*sss_prm.tol );
		fevalVals_splitspace(n) = sss_datOut.fevalCount;
		matVR = sss_datOut.matVR;
		matWR = sss_datOut.matWR;
		kTotVals_splitspace(n) = size(matVR,2);
	endif
endfor
%
numFigs++; figure(numFigs);
plot( ...
  fevalVals_compound, "o-", "linewidth", 3, "markersize", 15, ...
  fevalVals_splitspace, "x-", "linewidth", 2, "markersize", 10 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  [ 0, cumsum(fevalVals_compound) ], "o-", "linewidth", 3, "markersize", 15, ...
  [ 0, cumsum(fevalVals_splitspace) ], "x-", "linewidth", 2, "markersize", 10 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  kTotVals_compound, "o-", "linewidth", 3, "markersize", 15, ...
  kTotVals_splitspace, "x-", "linewidth", 2, "markersize", 10 );
grid on;

toc();
