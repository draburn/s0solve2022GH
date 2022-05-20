% Function...

function [ vecY, datOut ] = findLevPt_0519( vecG, matH, bTrgt=[], matB=[], prmIn=[], datIn=[] )
	mydefs;
	[ matB, prm, dat ] = __init( vecG, matH, bTrgt, matB, prmIn, datIn );
	%
	t1 = 1.0;
	[ b1, bPrime1, vecY1, funchYPrime1 ] = __calc( t1, vecG, matH, matB, prm, dat );
	if ( isempty(bTrgt) )
		msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Full step was requested." );
		vecY = vecY1;
		datOut.retCode = RETCODE__SUCCESS;
		datOut.t = t1;
		datOut.b = b1;
		datOut.bPrime = bPrime1;
		datOut.vecYPrime = funchYPrime1();
		datOut.iterCount = 0;
		return;
	elseif ( b1 <= bTrgt + bTrgt*prm.bRelTol )
		msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Full step is within tolerance." );
		vecY = vecY1;
		datOut.retCode = RETCODE__SUCCESS;
		datOut.t = t1;
		datOut.b = b1;
		datOut.bPrime = bPrime1;
		datOut.vecYPrime = funchYPrime1();
		datOut.iterCount = 0;
		return;
	endif
	%
	%
	t0 = 0.0;
	[ b0, bPrime0, vecY0, funchYPrime0 ] = __calc( t0, vecG, matH, matB, prm, dat );
	if ( b0 >= bTrgt - bTrgt*prm.bRelTol )
		msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Zero step is within tolerance." );
		vecY = vecY0;
		datOut.retCode = RETCODE__SUCCESS;
		datOut.t = t0;
		datOut.b = b0;
		datOut.bPrime = bPrime0;
		datOut.vecYPrime = funchYPrime0();
		datOut.iterCount = 0;
		return;
	endif
	%
	%
	[ vecY, datOut ] = __find( t0, b0, bPrime0, t1, b1, bPrime1, vecG, matH, bTrgt, matB, prm, dat );
	return;
endfunction


function [ matB, prm, dat ] = __init( vecG, matH, bTrgt=[], matB=[], prmIn=[], datIn=[] )
	mydefs;
	%
	if ( isempty(matB) )
		matB = eye(size(matH));
	endif
	%
	prm = struct();
	prm.verbLev = VERBLEV__FLAGGED; prm.valdLev = VALDLEV__ZERO; % Production.
	prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__MEDIUM; % Integration.
	prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Dev.
	prm.bRelTol = 1.0E-4;
	prm.cholRelTol = sqrt(eps);
	prm.extrapThresh0 = 100.0*eps;
	prm.extrapThresh1 = 1.0 - 100.0*eps;
	prm.iterMax = 100;
	prm = overwritefields( prm, prmIn );
	%
	matC = mygetfield( datIn, "matC", [] );
	if (isempty(matC))
		matC = matB'*matB;
	endif
	dat.matC = matC;
	dat.matE0 = sqrt(eps)*matH;
	dat.matEX = sqrt(eps)*eye(size(matH));
	dat.matE1 = sqrt(eps)*matC;
	dat = overwritefields( dat, prmIn );
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( isrealscalar(prm.bRelTol) );
		assert( 0.0 < prm.bRelTol );
		assert( prm.bRelTol < 1.0 );
		assert( isrealscalar(prm.cholRelTol) );
		assert( 0.0 <= prm.cholRelTol )
		assert( prm.cholRelTol <= 1.0 );
		assert( isrealscalar(prm.bRelTol) );
		assert( isrealscalar(prm.extrapThresh0) );
		assert( isrealscalar(prm.extrapThresh1) );
		assert( 0.0 <= prm.extrapThresh0 );
		assert( prm.extrapThresh0 <= prm.extrapThresh1 );
		assert( prm.extrapThresh1 <= 1.0 );
		%
		sz = size( vecG, 1 );
		assert( isrealarray(vecG,[sz,1]) );
		assert( isrealarray(matH,[sz,sz]) );
		assert( issymmetric(matH) );
		if (~isempty(bTrgt))
			assert( isscalar(bTrgt) );
			assert( 0.0 <= bTrgt );			
		endif
		szb = size(matB,1);
		assert( isrealarray(matB,[szb,sz]) );
		assert( isrealarray(dat.matC,[sz,sz]) );
		assert( issymmetric(dat.matC) );
		assert( isrealarray(dat.matE0,[sz,sz]) );
		assert( issymmetric(dat.matE0) );
		assert( isrealarray(dat.matE1,[sz,sz]) );
		assert( issymmetric(dat.matE1) );
		assert( isrealarray(dat.matEX,[sz,sz]) );
		assert( issymmetric(dat.matEX) );
	endif
	if ( prm.valdLev >= VALDLEV__MEDIUM );
		assert( reldiff(dat.matC,matB'*matB) < sqrt(eps) );
	endif
	if ( prm.valdLev >= VALDLEV__VERY_HIGH )
		eigH = eig(matH);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matH): %g ~ %g", min(eigH), max(eigH) ) );
		eigC = eig(dat.matC);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matC): %g ~ %g", min(eigC), max(eigC) ) );
		eigE0 = eig(dat.matE0);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matE0): %g ~ %g", min(eigE0), max(eigE0) ) );
		eigEX = eig(dat.matEX);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matEX): %g ~ %g", min(eigEX), max(eigEX) ) );
		eigE1 = eig(dat.matE1);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matE1): %g ~ %g", min(eigE1), max(eigE1) ) );
		%
		if ( min(eigH) < -sqrt(eps)*max(abs(eigH)) )
			error( "Hessian matrix has a clearly negative eigenvalue." );
		elseif ( min(eigC) < -eps*max(abs(eigC)) )
			error( "Constraint matrix appears to have a non-positive eigenvalue." );
		elseif ( min(eigE0) < -eps*max(abs(eigE0)) )
			error( "Regularization matrix E0 appears to have a non-positive eigenvalue." );
		elseif ( min(eigEX) < -eps*max(abs(eigEX)) )
			error( "Regularization matrix EX appears to have a non-positive eigenvalue." );
		elseif ( min(eigE1) < -eps*max(abs(eigE1)) )
			error( "Regularization matrix E1 appears to have a non-positive eigenvalue." );
		endif
	endif
	return;
endfunction


% If the initial Cholesky factorization fails, try again with extrapolation or checking tolerance.
% Since d/dt( y ) is typically not needed, return merely a function handle allowing calculation.
function [ b, bPrime, vecY, funchYPrime ] = __calc( t, vecG, matH, matB, prm, dat )
	mydefs;
	[ matR, cholFlag ] = chol( t*matH + (1.0-t)*dat.matC );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > prm.cholRelTol * max(abs(diag(matR))) )
		[ b, bPrime, vecY, vecRho ] = __calcFromChol( t, vecG, matR, matB, dat.matC );
		funchYPrime = @()( matR \ vecRho );
		return;
	endif
	endif
	%
	msgif( prm.msgInfo, __FILE__, __LINE__, sprintf( "Extrapolation is necessary for t = %g.", t ) );
	if ( t <= prm.extrapThresh0 )
		[ matR1, cholFlag1 ] = chol( t*matH + (1.0-t)*dat.matC + dat.matE0 );
		[ matR2, cholFlag2 ] = chol( t*matH + (1.0-t)*dat.matC + 2.0*dat.matE0 );
	elseif ( t < prm.extrapThresh1 )
		[ matR1, cholFlag1 ] = chol( t*matH + (1.0-t)*dat.matC + dat.matEX );
		[ matR2, cholFlag2 ] = chol( t*matH + (1.0-t)*dat.matC + 2.0*dat.matEX );
	else
		[ matR1, cholFlag1 ] = chol( t*matH + (1.0-t)*dat.matC + dat.matE1 );
		[ matR2, cholFlag2 ] = chol( t*matH + (1.0-t)*dat.matC + 2.0*dat.matE1 );
	endif
	if ( 0 ~= cholFlag1 )
		error( "Cholesky factorization failed even with regularization matrix." );
	elseif ( 0 ~= cholFlag2 )
		error( "Cholesky factorization (somehow) failed with regularization matrix second time." );
	endif
	%
	[ b1, bPrime1, vecY1, vecRho1 ] = __calcFromChol( t, vecG, matR1, matB, dat.matC );
	[ b2, bPrime2, vecY2, vecRho2 ] = __calcFromChol( t, vecG, matR2, matB, dat.matC );
	b = (2.0*b1) - b2;
	bPrime = (2.0*bPrime1) - bPrime2;
	vecY = (2.0*vecY1) - vecY2;
	funchYPrime = @()(  (2.0*( matR1 \ vecRho1 )) - ( matR2 \ vecRho2 )  );
	return;
endfunction


% Here's how the math works...
%  C  ==  B^T * B
%  M  ==  t*H + (1-t)*C  =  R^T * R
%  d/dt( M^-1 )  =  -M^-1 * ( H - C ) * M^-1
%  gamma  ==  M^-1 * (-g)
%  y  =  t * gamma
%  b  =  || B * y ||  =  t * sqrt( gamma^T * C * gamma )
%  d/dt( y )  =  M^-1 * C * gamma
%  d/dt( b )  =  (C * gamma)^T * M^-1 * (C * gamma)
% Note: "C" used to be called "S".
%
% Using Cholesky factorization of M, and being efficient...
%  M = R^T * R
%  gamma = R \ ( R^T \ (-g) )
%  b0 = || B * gamma ||
%  rho = R^T \ ( C * gamma )
%  y = t * gamma;
%  b = t * b0;
%  d/dt( b ) = || rho ||^2 / b0
%  d/dt( y ) = R \ rho
% But, since d/dt(y) is typically not needed and requires an extra backsub, we return R and rho.
function [ b, bPrime, vecY, vecRho ] = __calcFromChol( t, vecG, matR, matB, matC )
	vecGamma = matR \ ( matR' \ (-vecG) );
	b0 = norm( matB * vecGamma );
	vecRho = matR' \ ( matC * vecGamma );
	b = t * b0;
	bPrime = sumsq( vecRho ) / b0;
	vecY = t * vecGamma;
	return;
endfunction


function [ vecY, datOut ] = __find( tLo, bLo, bPrimeLo, tHi, bHi, bPrimeHi, vecG, matH, bTrgt, matB, prm, dat )
	[ vecY, datOut ] = __find_bisection( tLo, bLo, bPrimeLo, tHi, bHi, bPrimeHi, vecG, matH, bTrgt, matB, prm, dat );
	return;
endfunction


% Placeholder.
function [ vecY, datOut ] = __find_bisection( tLo, bLo, bPrimeLo, tHi, bHi, bPrimeHi, vecG, matH, bTrgt, matB, prm, dat )
	mydefs;
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( tLo < tHi );
		assert( bLo < bTrgt );
		assert( bHi > bTrgt );
		assert( bPrimeLo > 0.0 );
		assert( bPrimeHi > 0.0 );
	endif
	%
	for iterCount = 1 : prm.iterMax
		t = (tLo+tHi)/2.0;
		[ b, bPrime, vecY, funchYPrime ] = __calc( t, vecG, matH, matB, prm, dat );
		if ( abs(b-bTrgt) < bTrgt*prm.bRelTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
			datOut.retCode = RETCODE__SUCCESS;
			datOut.t = t;
			datOut.b = b;
			datOut.bPrime = bPrime;
			datOut.vecYPrime = funchYPrime();
			datOut.iterCount = iterCount;
			return;
		endif
		if ( t <= tLo || t >= tHi )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "NUMERICAL ISSUE: Went out of bounds." );
			datOut.retCode = RETCODE__NUMERICAL_ISSUE;
			datOut.t = t;
			datOut.b = b;
			datOut.bPrime = bPrime;
			datOut.vecYPrime = funchYPrime();
			datOut.iterCount = iterCount;
			return;
		endif
		if ( b < bTrgt )
			tLo = t;
			bLo = b;
			bPrimeLo = bPrime;
		else
			tHi = t;
			bHi = b;
			bPrimeHi = bPrime;
		endif
	endfor
	%
	msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: Reached iterMax." );
	datOut.retCode = RETCODE__IMPOSED_STOP;
	datOut.t = t;
	datOut.b = b;
	datOut.bPrime = bPrime;
	datOut.vecYPrime = funchYPrime();
	datOut.iterCount = iterCount;
	return;
endfunction


%!test
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 5;
%!	%
%!	sizeF = sizeX;
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX);
%!	matB = randn(sizeF,sizeX);
%!	%
%!	matH = matJ'*matJ;
%!	vecG = matJ'*vecF;
%!	bTrgt = abs(randn())*norm(matB*(matH\vecG));
%!	%
%!	prm = [];
%!	dat = [];
%!	[ vecY, datOut ] = findLevPt_0519( vecG, matH, bTrgt, matB, prm, dat );
%!	bActual = norm(matB*vecY);
%!	msg( __FILE__, __LINE__, sprintf( "bDesired = %g,  bActual = %g,  residual = %g.", bTrgt, bActual, bActual - bTrgt ) );


%!test
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 500;
%!	%
%!	sizeF = sizeX;
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX);
%!	matB = randn(sizeF,sizeX);
%!	%
%!	matH = matJ'*matJ;
%!	vecG = matJ'*vecF;
%!	bTrgt = abs(randn())*norm(matB*(matH\vecG));
%!	%
%!	prm = [];
%!	dat = [];
%!	[ vecY, datOut ] = findLevPt_0519( vecG, matH, bTrgt, matB, prm, dat );
%!	bActual = norm(matB*vecY);
%!	msg( __FILE__, __LINE__, sprintf( "bDesired = %g,  bActual = %g,  residual = %g.", bTrgt, bActual, bActual - bTrgt ) );
