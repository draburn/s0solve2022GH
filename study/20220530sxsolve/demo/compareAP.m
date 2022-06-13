clear;
numFigs = 0;
setprngstates(0);
tic();
%
%
% Init calculation stuff.
sizeX = 100;
sizeF = sizeX;
cI = 10.0;
cD = 0.1;
matIX = eye(sizeF,sizeX);
matJ0 = diag(1.0+abs(randn(sizeX)));
matR = randn(sizeF,sizeX);
tol = 0.1;
numRuns = 20;
%
%
% First run: either full space or not...
matJ = cI*matJ0 + matR;
vecX_secret = randn(sizeX,1);
vecF = matJ*vecX_secret;
funchW = @(v)( matJ*v );
%
linsolf_prm = [];
linsolf_prm.tol = tol;
[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
fevalIncr = linsolf_datOut.fevalCount;
matV = linsolf_datOut.matV;
matW = linsolf_datOut.matW;
%
%
% Init APs...
n = 0;
%
condVals(1) = cond(matJ);
%
noAP_fevalCountVals(1) = fevalIncr;
%
compound_fevalCountVals(1) = fevalIncr;
compound_matJA = matIX + (matW-matV)*(matV');
%
compoundFull0_fevalCountVals(1) = sizeX;
compoundFull0_matJA = matJ;
%
%reorth_fevalCountVals(1) = fevalIncr;
%reorth_matV = matV;
%reorth_matW = matW;
%
%reorthFull0_fevalCountVals(1) = sizeX;
%reorthFull0_matV = matIX;
%reorthFull0h_matW = matJ;
%
splitspace_fevalCountVals(1) = fevalIncr;
splitspace_matVR = matV;
splitspace_matWR = matW;
%
splitspaceFull0_fevalCountVals(1) = sizeX;
splitspaceFull0_matVR = matIX;
splitspaceFull0_matWR = matJ;
%
for n=2:numRuns
	matR = ( matR + cD*randn(sizeF,sizeX) ) / sqrt( 1.0 + cD^2 );
	matJ = cI*matJ0 + matR;
	vecX_secret = randn(sizeX,1);
	vecF = matJ*vecX_secret;
	funchW = @(v)( matJ*v );
	%
	condVals(n) = cond(matJ);
	%
	if (1)
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		noAP_fevalCountVals(n) = linsolf_datOut.fevalCount;
	endif
	%
	if ( rcond(compound_matJA) > sqrt(eps) )
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(compound_matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		compound_fevalCountVals(n) = linsolf_datOut.fevalCount;
		compound_matJA += ( linsolf_datOut.matW - (compound_matJA*linsolf_datOut.matV) ) * ( linsolf_datOut.matV' );
	endif
	%
	if ( rcond(compoundFull0_matJA) > sqrt(eps) )
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(compoundFull0_matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		compoundFull0_fevalCountVals(n) = linsolf_datOut.fevalCount;
		compoundFull0_matJA += ( linsolf_datOut.matW - (compoundFull0_matJA*linsolf_datOut.matV) ) * ( linsolf_datOut.matV' );
	endif
	%
	if (1)
		sss_prm = [];
		sss_prm.matVR = splitspace_matVR;
		sss_prm.matWR = splitspace_matWR;
		sss_prm.tol = tol;
		[ sss_vecX, sss_datOut ] = splitspacesolf( funchW, -vecF, sizeX, sss_prm );
		splitspace_fevalCountVals(n) = sss_datOut.fevalCount;
		splitspace_matVR = sss_datOut.matVR;
		splitspace_matWR = sss_datOut.matWR;
	endif
	%
	if (1)
		sss_prm = [];
		sss_prm.matVR = splitspaceFull0_matVR;
		sss_prm.matWR = splitspaceFull0_matWR;
		sss_prm.tol = tol;
		[ sss_vecX, sss_datOut ] = splitspacesolf( funchW, -vecF, sizeX, sss_prm );
		splitspaceFull0_fevalCountVals(n) = sss_datOut.fevalCount;
		splitspaceFull0_matVR = sss_datOut.matVR;
		splitspaceFull0_matWR = sss_datOut.matWR;
	endif
endfor
%
%
%
numFigs++; figure(numFigs);
plot( ...
  noAP_fevalCountVals, "o-", "linewidth", 2, "markersize", 40, ...
  compound_fevalCountVals, "x-", "linewidth", 2, "markersize", 32, ...
  compoundFull0_fevalCountVals, "s-", "linewidth", 2, "markersize", 24, ...
  splitspace_fevalCountVals, "^-", "linewidth", 2, "markersize", 16, ...
  splitspaceFull0_fevalCountVals, "v-", "linewidth", 2, "markersize", 8 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  [ 0.0, cumsum(noAP_fevalCountVals) ], "o-", "linewidth", 2, "markersize", 40, ...
  [ 0.0, cumsum(compound_fevalCountVals) ], "x-", "linewidth", 2, "markersize", 32, ...
  [ 0.0, cumsum(compoundFull0_fevalCountVals) ], "s-", "linewidth", 2, "markersize", 24, ...
  [ 0.0, cumsum(splitspace_fevalCountVals) ], "^-", "linewidth", 2, "markersize", 16, ...
  [ 0.0, cumsum(splitspaceFull0_fevalCountVals) ], "v-", "linewidth", 2, "markersize", 8 );
grid on;
%
numFigs++; figure(numFigs);
semilogy( condVals, "*-", "linewidth", 2, "markersize", 10 );
grid on;
%
%
msg( __FILE__, __LINE__, "Goodbye." );
toc();
return;
