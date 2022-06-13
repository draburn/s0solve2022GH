clear;
numFigs = 0;
setprngstates(0);
tic();
%
%
% Init calculation stuff.
sizeX = 500;
sizeF = sizeX;
cI = 10.0;
cD = 0.01;
matIX = eye(sizeF,sizeX);
matR = randn(sizeF,sizeX);
tol = 0.01;
numRuns = 20;
%
%
% First run: either full space or not...
matJ = cI*matIX + matR;
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
	matJ = cI*matIX + matR;
	vecX_secret = randn(sizeX,1);
	vecF = matJ*vecX_secret;
	funchW = @(v)( matJ*v );
	%
	condVals(n) = cond(matJ);
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
  compound_fevalCountVals, "o-", "linewidth", 2, "markersize", 32, ...
  compoundFull0_fevalCountVals, "s-", "linewidth", 2, "markersize", 24, ...
  splitspace_fevalCountVals, "^-", "linewidth", 2, "markersize", 16, ...
  splitspaceFull0_fevalCountVals, "v-", "linewidth", 2, "markersize", 8 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  [ 0.0, cumsum(compound_fevalCountVals) ], "o-", "linewidth", 2, "markersize", 32, ...
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



funchW = @(v)( matJ*v );
%
numRuns = 100;
matJA_compound = matIX;
fevalVals_compound = zeros(1,numRuns);
matVTot_compound = zeros(sizeX,1);
kTotVals_compound = zeros(1,numRuns);
fevalVals_splitspace = zeros(1,numRuns);
kTotVals_splitspace = zeros(1,numRuns);
for n=1:numRuns
	%%%matJ += 0.01*randn(sizeF,sizeX);
	vecX_secret = randn(sizeX,1);
	vecF = matJ*vecX_secret;
	%
	linsolf_prm = [];
	if ( rcond(matJA_compound) < sqrt(eps) )
		msg( __FILE__, __LINE__, sprintf( "rcond(matJA_compound) = %g.", rcond(matJA_compound) ) );
		break;
	endif
	linsolf_prm.matP = inv(matJA_compound);
	linsolf_prm.tol = 0.01;
	[ vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
	%assert( reldiff(matJ*vecX,-vecF) < 1.1*linsolf_prm.tol );
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
		%assert( reldiff(matJ*vecX,-vecF) < 1.1*sss_prm.tol );
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
  fevalVals_splitspace, "x-", "linewidth", 2, "markersize", 10, ...
  [ sizeX, ones(1,numRuns-1)], "^-", "markersize", 5 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  [ 0, cumsum(fevalVals_compound) ], "o-", "linewidth", 3, "markersize", 15, ...
  [ 0, cumsum(fevalVals_splitspace) ], "x-", "linewidth", 2, "markersize", 10, ...
  [ 0, sizeX+(0:numRuns-1) ], "^-", "markersize", 5 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  kTotVals_compound, "o-", "linewidth", 3, "markersize", 15, ...
  kTotVals_splitspace, "x-", "linewidth", 2, "markersize", 10, ...
  sizeX*ones(1,numRuns), "^-", "markersize", 5 );
grid on;

toc();
