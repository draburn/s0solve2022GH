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
	prm.verbLev = VERBLEV__FLAGGED; prm.valdLev = VALDLEV__ZERO; % Production.
	prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__MEDIUM; % Integration.
	%prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Dev.
	prm.bRelTol = 1.0E-4;
	prm.cholRelTol = sqrt(eps);
	prm.extrapThresh0 = 100.0*eps;
	prm.extrapThresh1 = 1.0 - 100.0*eps;
	prm.iterMax = 100;
	prm.minStepCoeff = 0.1;
	prm.minStepFallThresh = 0.1;
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
		assert( isrealscalar(prm.iterMax) );
		assert( abs(prm.iterMax-round(prm.iterMax)) < sqrt(eps) );
		assert( 1 <= prm.iterMax );
		assert( isrealscalar(prm.minStepCoeff) );
		assert( 0.0 <= prm.minStepCoeff )
		assert( prm.minStepCoeff <= 0.5 );
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
		elseif ( min(eigC) < -sqrt(eps)*max(abs(eigC)) )
			error( "Constraint matrix has a clearly negative eigenvalue." );
		elseif ( min(eigE0) < -eps*max(abs(eigE0)) && min(eigC) < -eps*max(abs(eigC)) )
			error( "Constraint and E0 regularization matrices appear to have non-positive eigenvalues." );
		elseif ( min(eigEX) < -eps*max(abs(eigEX)) && min(eigH) < -eps*max(abs(eigH)) && min(eigC) < -eps*max(abs(eigC)) )
			error( "Hessian, constraint, and EX regularization matrices appear to have non-positive eigenvalues." );
		elseif ( min(eigE1) < -eps*max(abs(eigE1)) && min(eigH) < -eps*max(abs(eigH)) )
			error( "Hessian and E1 regularization matrices appear to have non-positive eigenvalues." );
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
	msgif( prm.verbLev >= VERBLEV__INFO, __FILE__, __LINE__, sprintf( "Extrapolation is necessary for t = %g.", t ) );
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
	mydefs;
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( tLo < tHi );
		assert( bLo < bTrgt );
		assert( bHi > bTrgt );
		assert( bPrimeLo > 0.0 );
		assert( bPrimeHi > 0.0 );
	endif
	%
	iterCount = 0;
	applyMinStepConstraint = false;
	while (1)
		if ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "IMPOSED STOP: Reached iterMax (%d).", prm.iterMax ) );
			datOut.retCode = RETCODE__IMPOSED_STOP;
			break;
		endif
		iterCount++;
		%
		if ( applyMinStepConstraint )
			tMin = tLo + prm.minStepCoeff*(tHi-tLo);
			tMax = tHi - prm.minStepCoeff*(tHi-tLo);
		else
			tMin = tLo;
			tMax = tHi;
		endif
		[ t, bModel ] = __cubicRoot( tLo, bLo, bPrimeLo, tHi, bHi, bPrimeHi, bTrgt, tMin, tMax );
		if ( ~isempty(t) )
			stepTypeStr = "c";
		elseif ( abs(bLo-bTrgt) < abs(bHi-bTrgt) )
			t = median([ tMin, tLo + (bTrgt-bLo)/bPrimeLo, (tHi+tLo)/2.0 ]);
			bModel = bLo + bPrimeLo*(t-tLo);
			stepTypeStr = "l";
		else
			t = median([ (tHi+tLo)/2.0, tHi - (bHi-bTrgt)/bPrimeHi, tMax ] );
			bModel = bHi + bPrimeHi*(t-tHi);
			stepTypeStr = "h";
		endif
		%
		if ( t <= tLo || t >= tHi )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( ...
			  "NUMERICAL ISSUE: Guess is out of bounds ( %g ~ %g ~ %g).", tLo, t, tHi )  );
			datOut.retCode = RETCODE__NUMERICAL_ISSUE;
			break;
		endif
		[ b, bPrime, vecY, funchYPrime ] = __calc( t, vecG, matH, matB, prm, dat );
		msgif( prm.verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( ...
		  "  %2d: %s %d;  t: %g ~ %g ~ %g,  1-t: %g ~ %g ~ %g,  b: %g ~ %g (e %g, w %g) ~ %g;  b-bT: %g (e %g).", ...
		  iterCount, stepTypeStr, applyMinStepConstraint, ...
		  tLo, t, tHi, 1.0-tLo, 1.0-t, 1.0-tHi, ...
		  bLo, b, bModel, bTrgt, bHi, b-bTrgt, bModel-bTrgt ) );
		if ( abs(b-bTrgt) < bTrgt*prm.bRelTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
			datOut.retCode = RETCODE__SUCCESS;
			break;
		endif
		if ( b <= bLo || b >= bHi || bPrime <= 0.0 )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( ...
			  "NUMERICAL ISSUE: Function was non-monotonic ( %g ~ %g ~ %g, %g ~ %g ~ %g).", ...
			  bLo, b, bHi, bPrimeLo, bPrime, bPrimeHi )  );
			datOut.retCode = RETCODE__NUMERICAL_ISSUE;
			break;
		endif
		%
		if ( b < bTrgt )
			if ( ~applyMinStepConstraint && abs(b-bModel) > prm.minStepFallThresh*abs(bHi-bLo) )
				applyMinStepConstraint = true;
			else
				applyMinStepConstraint = false;
			endif
			tLo = t;
			bLo = b;
			bPrimeLo = bPrime;
		else
			if ( ~applyMinStepConstraint && abs(b-bModel) > prm.minStepFallThresh*abs(bHi-bLo) )
				applyMinStepConstraint = true;
			else
				applyMinStepConstraint = false;
			endif
			tHi = t;
			bHi = b;
			bPrimeHi = bPrime;
		endif
	endwhile
	%
	datOut.t = t;
	datOut.b = b;
	datOut.bPrime = bPrime;
	datOut.vecYPrime = funchYPrime();
	datOut.iterCount = iterCount;
	return;
endfunction


function [ t, bModel ] = __cubicRoot( tLo, bLo, bPrimeLo, tHi, bHi, bPrimeHi, bTrgt, tMin, tMax )
	[ x, fModel ] = __cubicRootScaled( ...
	  bLo - bTrgt, bPrimeLo*(tHi-tLo), bHi - bTrgt, bPrimeHi*(tHi-tLo), (tMin-tLo)/(tHi-tLo), (tMax-tLo)/(tHi-tLo) );
	if ( isempty(x) )
		t = [];
		bModel = [];
	else
		t = tLo + (tHi-tLo)*x;
		bModel = bTrgt + fModel;
	endif
	return;
endfunction
% f = a*(x^3) + b*(x^2) + fp0*x + f0.
% fp = 3*a*(x^2) + 2*b*x + fp0.
% fpp = 6*a + 2*b
% THIS IS HACKY!
function [ x, fModel ] = __cubicRootScaled( f0, fp0, f1, fp1, xMin, xMax )
	a = (fp1-fp0) - 2.0*(f1-f0-fp0);
	b = 3.0*(f1-f0-fp0) - (fp1-fp0);
	% since fp0, fp1 > 0, f' > 0 where f'' = 0 => no quad roots.
	xe = -b/(3.0*a);
	if ( 3.0*a*(xe^2) + 2.0*b*xe + fp0 <= 0.0 )
		x = [];
		fModel = [];
		return;
	endif
	%
	r = roots([ a, b, fp0, f0 ]);
	x = [];
	for n=1:length(r)
	if ( isreal(r(n)) && 0.0 < r(n) && r(n) < 1.0 )
		if ( isempty(x) )
			x = r(n);
		else
			msg( __FILE__, __LINE__, "Data dump..." );
			[ f0, fp0, f1, fp1 ]
			[ a, b, fp0, f0 ]
			r
			x
			error( "Inconsistency error calculating cubic root." );
		endif
	endif
	endfor
	if ( isempty(x) )
		msg( __FILE__, __LINE__, "Data dump..." );
		[ f0, fp0, f1, fp1 ]
		[ a, b, fp0, f0 ]
		r
		x
		error( "Inconsistency error calculating cubic root." );
	endif
	x = median([ xMin, x, xMax ]);
	fModel = a*(x^3) + b*(x^2) + fp0*x + f0;
	return;
endfunction


%!test
%!	msg( __FILE__, __LINE__, "SKIPPING TEST." ); return;
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
%!	b1 = norm(matB*((matH+eye(size(matH))*max(diag(matH))*sqrt(eps))\vecG));
%!	bTrgt = min([ 0.5*abs(randn())*b1, b1 ]);
%!	%
%!	prm = [];
%!	dat = [];
%!	[ vecY, datOut ] = findLevPt_0519( vecG, matH, bTrgt, matB, prm, dat );
%!	bActual = norm(matB*vecY);
%!	msg( __FILE__, __LINE__, sprintf( "bDesired = %g,  bActual = %g,  residual = %g.", bTrgt, bActual, bActual - bTrgt ) );


%!test
%!	%msg( __FILE__, __LINE__, "SKIPPING TEST." ); return;
%!	numFigs = 0;
%!	setprngstates(47290048);
%!	sizeX = 1000;
%!	%
%!	sizeF = 900;
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX);
%!	matB = randn(sizeX,sizeX);
%!	%
%!	matH = matJ'*matJ;
%!	vecG = matJ'*vecF;
%!	b1 = norm(matB*((matH+eye(size(matH))*max(diag(matH))*sqrt(eps))\vecG));
%!	bTrgt = min([ 0.5*abs(randn())*b1, b1 ]);
%!	%
%!	prm = [];
%!	%prm.bRelTol = 0.1;
%!	dat = [];
%!	%dat.matC = matB'*matB;
%!	tic();
%!	[ vecY, datOut ] = findLevPt_0519( vecG, matH, bTrgt, matB, prm, dat );
%!	toc();
%!	bActual = norm(matB*vecY);
%!	msg( __FILE__, __LINE__, sprintf( "bDesired = %g,  bActual = %g,  residual = %g.", bTrgt, bActual, bActual - bTrgt ) );
