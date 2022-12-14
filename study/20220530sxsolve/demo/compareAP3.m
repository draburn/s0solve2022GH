clear;
numFigs = 0;
prngstates = setprngstates(0);
tic();
%
% Init calculation stuff.
sizeX = 50
sizeF = sizeX

cI = 2.0
%cI = 100.0

cR = 0.1
cD = 0.1
matIX = eye(sizeF,sizeX);
%matJ0 = diag(1.0+cR*abs(randn(sizeX,1)));
matJ0 = diag( 1.0+cR*randn(sizeX,1) );
%matJ0 = randn(sizeF,sizeX);
matR = randn(sizeF,sizeX);
matSF = diag(exp(1.0*randn(sizeF,1)));
matSX = diag(exp(1.0*randn(sizeX,1)));
tol = 0.1
numRuns = 20;
strRunName = sprintf( "cI%g_cR%g_cD%g_prng%d_x%d_tol%g", cI, cR, cD, prngstates, sizeX, tol );
%
%matJ = cI*matJ0 + matR;
matJ = matSF*( cI*matJ0 + matR ) /matSX;
vecX_secret = randn(sizeX,1);
vecF = matJ*vecX_secret;
funchW = @(v)( matJ*v );
%
%
% First run: either full space or not...
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
compound_vecX = -matV*(matW\vecF);
%
compoundFull0_fevalCountVals(1) = sizeX;
compoundFull0_matJA = matJ;
%
reorth_fevalCountVals(1) = fevalIncr;
reorth_matVPool = matV;
reorth_matWPool = matW;
%
reorthFull0_fevalCountVals(1) = sizeX;
reorthFull0_matVPool = matIX;
reorthFull0_matWPool = matJ;
%
splitspace_fevalCountVals(1) = fevalIncr;
splitspace_matVR = matV;
splitspace_matWR = matW;
%
splitspaceFull0_fevalCountVals(1) = sizeX;
splitspaceFull0_matVR = matIX;
splitspaceFull0_matWR = matJ;
%
semisplitspace_fevalCountVals(1) = fevalIncr;
semisplitspace_matVR = matV;
semisplitspace_matWR = matW;
%
semisplitspaceFull0_fevalCountVals(1) = sizeX;
semisplitspaceFull0_matVR = matIX;
semisplitspaceFull0_matWR = matJ;
%
semicombo_fevalCountVals(1) = fevalIncr;
semicombo_matVR = matV;
semicombo_matWR = matW;
%
for n=2:numRuns
	if ( stopsignalpresent() )
		msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
		retCode = RETCODE__IMPOSED_STOP;
		break;
	endif
	if (0)
		matR = ( matR + cD*randn(sizeF,sizeX) ) / sqrt( 1.0 + cD^2 );
		%matJ = cI*matJ0 + matR;
		matJ = matSF*( cI*matJ0 + matR ) /matSX;
		vecX_secret = randn(sizeX,1);
		vecF = matJ*vecX_secret;
	else
		matR = ( matR + cD*randn(sizeF,sizeX) ) / sqrt( 1.0 + cD^2 );
		%matJ = cI*matJ0 + matR;
		matJ = matSF*( cI*matJ0 + matR ) /matSX;
		vecF = matJ*(vecX_secret+compound_vecX);
		norm(vecF)
	endif
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
	%
	if ( 1 )
		matJA = compound_matJA;
		%
		if ( rcond(matJA) > 100.0*eps );
		assert( rcond(matJA) > 100.0*eps );
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		%
		compound_fevalCountVals(n) = linsolf_datOut.fevalCount;
		compound_matJA += ( linsolf_datOut.matW - (matJA*linsolf_datOut.matV) ) * ( linsolf_datOut.matV' );
		compound_vecX -= linsolf_datOut.matV*(linsolf_datOut.matW\vecF);
		endif
	endif
	%
	if ( 1 )
		matJA = compoundFull0_matJA;
		%
		if ( rcond(matJA) > 100.0*eps );
		assert( rcond(matJA) > 100.0*eps );
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		%
		compoundFull0_fevalCountVals(n) = linsolf_datOut.fevalCount;
		compoundFull0_matJA += ( linsolf_datOut.matW - (matJA*linsolf_datOut.matV) ) * ( linsolf_datOut.matV' );
		endif
	endif
	%
	%
	if ( 1 )
		matVPool = reorth_matVPool;
		matWPool = reorth_matWPool;
		%
		matJA = matIX + (matWPool-matVPool)*(matVPool');
		if ( rcond(matJA) > 100.0*eps );
		assert( rcond(matJA) > 100.0*eps );
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		matVNew = linsolf_datOut.matV;
		matWNew = linsolf_datOut.matW;
		%
		[ matVPool1, matWPool1, rvecDrop ] = utorthdrop_pair( [ matVNew, matVPool ], [ matWNew, matWPool ], 1.0e-4 );
		reorth_fevalCountVals(n) = linsolf_datOut.fevalCount;
		reorth_matVPool = matVPool1;
		reorth_matWPool = matWPool1;
		endif
	endif
	%
	if ( 1 )
		matVPool = reorthFull0_matVPool;
		matWPool = reorthFull0_matWPool;
		%
		matJA = matIX + (matWPool-matVPool)*(matVPool');
		if ( rcond(matJA) > 100.0*eps );
		assert( rcond(matJA) > 100.0*eps );
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		matVNew = linsolf_datOut.matV;
		matWNew = linsolf_datOut.matW;
		%
		[ matVPool1, matWPool1, rvecDrop ] = utorthdrop_pair( [ matVNew, matVPool ], [ matWNew, matWPool ], 1.0e-4 );
		reorthFull0_fevalCountVals(n) = linsolf_datOut.fevalCount;
		reorthFull0_matVPool = matVPool1;
		reorthFull0_matWPool = matWPool1;
		endif
	endif
	%
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
	%
	%
	if (1)
		semi_prm = [];
		semi_prm.matVR = semisplitspace_matVR;
		semi_prm.matWR = semisplitspace_matWR;
		semi_prm.tol = tol;
		[ semi_vecX, semi_datOut ] = semisplitspacesolf( funchW, -vecF, sizeX, semi_prm );
		semisplitspace_fevalCountVals(n) = semi_datOut.fevalCount;
		semisplitspace_matVR = semi_datOut.matVR;
		semisplitspace_matWR = semi_datOut.matWR;
	endif
	%
	if (1)
		semi_prm = [];
		semi_prm.matVR = semisplitspaceFull0_matVR;
		semi_prm.matWR = semisplitspaceFull0_matWR;
		semi_prm.tol = tol;
		[ semi_vecX, semi_datOut ] = semisplitspacesolf( funchW, -vecF, sizeX, semi_prm );
		semisplitspaceFull0_fevalCountVals(n) = semi_datOut.fevalCount;
		semisplitspaceFull0_matVR = semi_datOut.matVR;
		semisplitspaceFull0_matWR = semi_datOut.matWR;
	endif
	%
	%
	if (1)
		% This doesn't seem to behave correctly.
		matVPool = semicombo_matVR;
		matWPool = semicombo_matWR;
		%rcond(matJ)
		%
		
		%vecFScale_est = sumsq(matWPool,2)/size(matWPool,2);
		%matSF_est = diag(sqrt( vecFScale_est + eps*norm(vecFScale_est) ));
		%matJA = matSF_est + (matWPool-matVPool)*(matVPool'); WRONG!
		%matJA = matSF_est*(matIX-matVPool*(matVPool')) + matWPool*(matVPool'); % RIGHT!
		
		%foo1 = sum(matWPool.*matVPool,2);
		%foo2 = sum(matVPool.*matVPool,2);
		%%%matJA = diag(foo1./foo2) + (matWPool-matVPool)*(matVPool'); WRONG!
		%matJA = diag(foo1./foo2)*(matIX-matVPool*(matVPool')) + matWPool*(matVPool'); % RIGHT?
		%rcond(matJA\matJ)
		%rcond(pinv(matJA)*matJ)
		%rcond(matJ/matJA)
		%rcond(matJ*pinv(matJA))
		
		%s = sqrt( sum(sumsq(matWPool))/sum(sumsq(matVPool)) );
		%s = 0.7
		s = 1.0;
		%%%matJA = s*matIX + (matWPool-matVPool)*(matVPool'); WRONG!
		matJA = s*(matIX-matVPool*(matVPool')) + matWPool*(matVPool'); % RIGHT!
		%rcond(matJA\matJ)
		%rcond(pinv(matJA)*matJ)
		%rcond(matJ/matJA)
		%rcond(matJ*pinv(matJA))
		
		%matJA = matIX + (matWPool-matVPool)*(matVPool');
		%rcond(matJA\matJ)
		%rcond(pinv(matJA)*matJ)
		%rcond(matJ/matJA)
		%rcond(matJ*pinv(matJA))
		%return
		
		
		%assert( rcond(matJA) > 100.0*eps );
		%if ( rcond(matJA) > 100.0*eps );
		if (1)
		linsolf_prm = [];
		linsolf_prm.tol = tol;
		linsolf_prm.matP = pinv(matJA);
		[ linsolf_vecX, linsolf_datOut ] = linsolf( funchW, -vecF, zeros(sizeX,1), linsolf_prm );
		matVNew = linsolf_datOut.matV;
		matWNew = linsolf_datOut.matW;
		%
		semicombo_fevalCountVals(n) = linsolf_datOut.fevalCount;
		
		switch (3)
		case 1
			% This is semicombo concept
			matVPlus = utorthdrop( [ matVNew, matVPool ], 1.0e-10 );
			[ matVTau, matWTau ] = utorthdrop_pair( [ matVNew, matVPool ], [ matWNew, matWPool ], 0.5 );
			matWTemp = matWPool*(matVPool'*matVPlus);
			matZ = matVPlus'*matVTau;
			matWPlus = matWTemp + ( matWTau - matWTemp * matZ ) * (matZ');
		case 2
			% Misc hacks
			matVPlus = utorthdrop( [ matVNew, matVPool ], 1.0e-10 );
			%%%[ matVTau, matWTau ] = utorthdrop_pair( [ matVNew, matVPool ], [ matWNew, matWPool ], 0.9 );
			matVTau = matVNew;
			matWTau = matWNew;
			matWTemp = matWPool*(matVPool'*matVPlus);
			matZ = matVPlus'*matVTau;
			matWPlus = matWTemp + ( matWTau - matWTemp * matZ ) * (matZ');
			
			%%%[ matVPlus, matWPlus ] = utorthdrop_pair( [ matVNew, matVPool], [ matWNew, matWPool ], 1.0e-10 );
			%[ matVPlus, matWPlus ] = utorthdrop_pair( [ matVNew, matVPool], [ matWNew, matWPool ], 0.7 );
			%matZ = matVPlus'*matVNew;
			%matWPlus = matWPlus + ( matWNew - matWPlus*matZ)*(matZ');
		case 3
			% Not semicombo concept, but testing something related ("oblique").
			matVPlus = utorthdrop( [ matVPool, matVNew ] );
			matWTemp = zeros(size(matVPlus));
			matWTemp(:,1:size(matWPool,2)) = matWPool;
			matZ = matVPlus'*matVNew;
			matWPlus = matWTemp + ( matWNew - matWTemp*matZ ) * ( matZ' );
		endswitch
		
		semicombo_matVR = matVPlus;
		semicombo_matWR = matWPlus;
		endif
	endif
endfor
%
%
%
numFigs++; figure(numFigs);
plot( ...
  semicombo_fevalCountVals, "<-", "linewidth", 2, "markersize", 29, ...
  semisplitspace_fevalCountVals, "d-", "linewidth", 2, "markersize", 27, ...
  semisplitspaceFull0_fevalCountVals, "h-", "linewidth", 2, "markersize", 25, ...
  compound_fevalCountVals, "x-", "linewidth", 2, "markersize", 23, ...
  compoundFull0_fevalCountVals, "s-", "linewidth", 2, "markersize", 21, ...
  splitspace_fevalCountVals, "^-", "linewidth", 2, "markersize", 19, ...
  splitspaceFull0_fevalCountVals, "v-", "linewidth", 2, "markersize", 17, ...
  reorth_fevalCountVals, "p-", "linewidth", 2, "markersize", 15, ...
  reorthFull0_fevalCountVals, "*-", "linewidth", 2, "markersize", 13, ...
  noAP_fevalCountVals, "o-", "linewidth", 2, "markersize", 11 );
set( legend( ...
  "semicombo/hack", ...
  "semisplit", ...
  "semisplit full0 ", ...
  "comp", ...
  "comp full0", ...
  "split", ...
  "split full0", ...
  "reorth", ...
  "reorth full0", ...
  "no AP" ), "Interpreter", "none" );
set( title( strRunName ), "Interpreter", "none" );
set( xlabel( "Jacobian index" ), "Interpreter", "none" );
set( ylabel( "feval count" ), "Interpreter", "none" );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  [ 0.0, cumsum(semicombo_fevalCountVals) ], "<-", "linewidth", 2, "markersize", 29, ...
  [ 0.0, cumsum(semisplitspace_fevalCountVals) ], "d-", "linewidth", 2, "markersize", 27, ...
  [ 0.0, cumsum(semisplitspaceFull0_fevalCountVals) ], "h-", "linewidth", 2, "markersize", 25, ...
  [ 0.0, cumsum(compound_fevalCountVals) ], "x-", "linewidth", 2, "markersize", 23, ...
  [ 0.0, cumsum(compoundFull0_fevalCountVals) ], "s-", "linewidth", 2, "markersize", 21, ...
  [ 0.0, cumsum(splitspace_fevalCountVals) ], "^-", "linewidth", 2, "markersize", 19, ...
  [ 0.0, cumsum(splitspaceFull0_fevalCountVals) ], "v-", "linewidth", 2, "markersize", 17, ...
  [ 0.0, cumsum(reorth_fevalCountVals) ], "p-", "linewidth", 2, "markersize", 15, ...
  [ 0.0, cumsum(reorthFull0_fevalCountVals) ], "*-", "linewidth", 2, "markersize", 13, ...
  [ 0.0, cumsum(noAP_fevalCountVals) ], "o-", "linewidth", 2, "markersize", 11 );
set( title( strRunName ), "Interpreter", "none" );
set( xlabel( "Jacobian index" ), "Interpreter", "none" );
set( ylabel( "cumulative feval count" ), "Interpreter", "none" );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  [ cumsum(semicombo_fevalCountVals) ], "<-", "linewidth", 2, "markersize", 29, ...
  [ cumsum(semisplitspace_fevalCountVals) ], "d-", "linewidth", 2, "markersize", 27, ...
  [ cumsum(semisplitspaceFull0_fevalCountVals) ], "h-", "linewidth", 2, "markersize", 25, ...
  [ cumsum(compound_fevalCountVals) ], "x-", "linewidth", 2, "markersize", 23, ...
  [ cumsum(compoundFull0_fevalCountVals) ], "s-", "linewidth", 2, "markersize", 21, ...
  [ cumsum(splitspace_fevalCountVals) ], "^-", "linewidth", 2, "markersize", 19, ...
  [ cumsum(splitspaceFull0_fevalCountVals) ], "v-", "linewidth", 2, "markersize", 17, ...
  [ cumsum(reorth_fevalCountVals) ], "p-", "linewidth", 2, "markersize", 15, ...
  [ cumsum(reorthFull0_fevalCountVals) ], "*-", "linewidth", 2, "markersize", 13 );
set( title( strRunName ), "Interpreter", "none" );
set( xlabel( "Jacobian index" ), "Interpreter", "none" );
set( ylabel( "cumulative feval count" ), "Interpreter", "none" );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  cumsum(semicombo_fevalCountVals)-cumsum(compoundFull0_fevalCountVals), "<-", "linewidth", 2, "markersize", 29, ...
  cumsum(semisplitspace_fevalCountVals)-cumsum(compoundFull0_fevalCountVals), "d-", "linewidth", 2, "markersize", 27, ...
  cumsum(semisplitspaceFull0_fevalCountVals)-cumsum(compoundFull0_fevalCountVals), "h-", "linewidth", 2, "markersize", 25, ...
  cumsum(compound_fevalCountVals)-cumsum(compoundFull0_fevalCountVals), "x-", "linewidth", 2, "markersize", 23, ...
  cumsum(compoundFull0_fevalCountVals)-cumsum(compoundFull0_fevalCountVals), "s-", "linewidth", 2, "markersize", 21, ...
  cumsum(splitspace_fevalCountVals)-cumsum(compoundFull0_fevalCountVals), "^-", "linewidth", 2, "markersize", 19, ...
  cumsum(splitspaceFull0_fevalCountVals)-cumsum(compoundFull0_fevalCountVals), "v-", "linewidth", 2, "markersize", 17 );
set( title( strRunName ), "Interpreter", "none" );
set( xlabel( "Jacobian index" ), "Interpreter", "none" );
set( ylabel( "cumulative feval count - 'cf0'" ), "Interpreter", "none" );
grid on;
%
numFigs++; figure(numFigs);
semilogy( condVals, "*-", "linewidth", 2, "markersize", 10 );
set( title( strRunName ), "Interpreter", "none" );
set( xlabel( "Jacobian index" ), "Interpreter", "none" );
set( ylabel( "condition number" ), "Interpreter", "none" );
grid on;
%
%
msg( __FILE__, __LINE__, "Goodbye." );
toc();
return;
