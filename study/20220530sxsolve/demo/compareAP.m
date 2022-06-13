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
%matJ0 = matIX;
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
reorth_fevalCountVals(1) = fevalIncr;
reorth_matVPool = matV;
reorth_matWPool = matW;
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
	if ( 1 )
		matJA = compound_matJA;
		%
		assert( rcond(matJA) > 100.0*eps );
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		%
		compound_fevalCountVals(n) = linsolf_datOut.fevalCount;
		compound_matJA += ( linsolf_datOut.matW - (matJA*linsolf_datOut.matV) ) * ( linsolf_datOut.matV' );
	endif
	%
	if ( 1 )
		matJA = compoundFull0_matJA;
		%
		assert( rcond(matJA) > 100.0*eps );
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		%
		compoundFull0_fevalCountVals(n) = linsolf_datOut.fevalCount;
		compoundFull0_matJA += ( linsolf_datOut.matW - (matJA*linsolf_datOut.matV) ) * ( linsolf_datOut.matV' );
	endif
	%
	%
	if ( 1 )
		matVPool = reorth_matVPool;
		matWPool = reorth_matWPool;
		%
		matJA = matIX + (matWPool-matVPool)*(matVPool');
		assert( rcond(matJA) > 100.0*eps );
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		matVNew = linsolf_datOut.matV;
		matWNew = linsolf_datOut.matW;
		%
		if (0)
			matVNew = [ 1.0; 0.0; 0.0; 0.0 ]
			matVNew'*matVNew
			matVPool = [ 1.0/sqrt(2.0), 1.0/sqrt(3.0); -1.0/sqrt(2.0), 1.0/sqrt(3.0); 0.0, 1.0/sqrt(3.0); 0.0, 0.0 ]
			matVPool'*matVPool
			matVPool0 = [ matVNew, matVPool ]
			matVPool0'*matVPool0
			
			[ matVPool1, rvecDrop ] = utorthdrop( matVPool0, 1.0e-4 )
			% VP1 = VP0 * R. But, VP0'*VP0 is not I.
			% So, R = (VP1'*VP0)^-1.
			matT = pinv(matVPool1'*matVPool0);
			matVPool1
			matVPool0*matT
			
			matWNew = [ 1.0; 2.0; 3.0; 0.0 ]
			matWPool = [ 1.0, 0.0; 0.0, 1.0; 0.0, 0.0; 0.0, 0.0 ]
			matWPool0 = [ matWNew, matWPool ]
			matWPool1 = matWPool0*matR
			%
			return
		endif
		if (0)
			matV = [ 1.0, 1.0/sqrt(2.0); 0.0, 1.0/sqrt(2.0) ]
			matW = [ 1.0, 10.0; 100.0, 1000.0 ]
			%
			matV1 = utorthdrop( matV )
			matW1 = matW/(matV1*matV)
			
			[ matV1, matW1 ] = utorthdrop_pair( matV, matW )
			return
		endif
		%
		%matVPool0 = [ matVNew, matVPool ]
		%[ matVPool1, rvecDrop ] = utorthdrop( matVPool0, 1.0e-4 )
		%matWPool1 = [ matWNew, matWPool ]/(matVPool1'*matVPool0)
		[ matVPool1, matWPool1, rvecDrop ] = utorthdrop_pair( [ matVNew, matVPool ], [ matWNew, matWPool ], 1.0e-4 );
		%return
		
		%
		reorth_fevalCountVals(n) = linsolf_datOut.fevalCount;
		reorth_matVPool = matVPool1;
		reorth_matWPool = matWPool1;
		%size(reorth_matVPool,2)
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
  noAP_fevalCountVals, "o-", "linewidth", 2, "markersize", 15, ...
  compound_fevalCountVals, "x-", "linewidth", 2, "markersize", 15, ...
  compoundFull0_fevalCountVals, "s-", "linewidth", 2, "markersize", 15, ...
  reorth_fevalCountVals, "p-", "linewidth", 2, "markersize", 15, ...
  splitspace_fevalCountVals, "^-", "linewidth", 2, "markersize", 15, ...
  splitspaceFull0_fevalCountVals, "v-", "linewidth", 2, "markersize", 15 );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  [ 0.0, cumsum(noAP_fevalCountVals) ], "o-", "linewidth", 2, "markersize", 15, ...
  [ 0.0, cumsum(compound_fevalCountVals) ], "x-", "linewidth", 2, "markersize", 15, ...
  [ 0.0, cumsum(compoundFull0_fevalCountVals) ], "s-", "linewidth", 2, "markersize", 15, ...
  [ 0.0, cumsum(reorth_fevalCountVals) ], "p-", "linewidth", 2, "markersize", 15, ...
  [ 0.0, cumsum(splitspace_fevalCountVals) ], "^-", "linewidth", 2, "markersize", 15, ...
  [ 0.0, cumsum(splitspaceFull0_fevalCountVals) ], "v-", "linewidth", 2, "markersize", 15 );
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
