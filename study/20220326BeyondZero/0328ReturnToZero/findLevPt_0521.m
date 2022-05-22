% Function...

function [ vecY, datOut ] = findLevPt_0521( vecG, matH, bTrgt=[], matB=[], prmIn=[], datIn=[] )
	mydefs;
	[ matB, prm, dat ] = __init( vecG, matH, bTrgt, matB, prmIn, datIn );
	%
	if ( bTrgt <= bTrgt*prm.bRelTol )
		[ bPrimeZero, funchYPrimeZero ] = __calcZero( vecG, matH, matB, prm, dat );
		vecY = zeros(size(vecG));
		datOut.mu = +Inf;
		datOut.b = 0.0;
		datOut.bPrime = bPrimeZero;
		datOut.vecYPrime = funchYPrimeZero();
		datOut.iterCount = 0;
		msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Zero step was requested." );
		return;
	endif
	%
	[ bNewt, bPrimeNewt, vecYNewt, funchYPrimeNewt ] = __calcNewt( vecG, matH, matB, prm, dat );
	if ( isempty(bTrgt) || abs(bTrgt-bNewt) <= bTrgt*prm.bRelTol )
		vecY = vecYNewt;
		datOut.mu = 0.0;
		datOut.b = bNewt;
		datOut.bPrime = bPrimeNewt;
		datOut.vecYPrime = funchYPrimeNewt();
		datOut.iterCount = 0;
		if ( isempty(bTrgt) )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Full step was requested." );
		else
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Full step is within tolerance." );
		endif
		return;
	endif
	%
	[ vecY, datOut ] = __find( bNewt, bPrimeNewt, vecG, matH, bTrgt, matB, prm, dat );
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
	prm.verbLev = VERBLEV__DETAILS; prm.valdLev = VALDLEV__HIGH; % Performance testing.
	prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Dev.
	prm.bRelTol = sqrt(eps);
	prm.cholRelTol = sqrt(eps);
	prm.epsReguRel = sqrt(eps);
	prm.iterMax = 100;
	prm = overwritefields( prm, prmIn );
	%
	matC = mygetfield( datIn, "matC", [] );
	if (isempty(matC))
		matC = matB'*matB;
	endif
	dat.matC = matC;
	dat.hScale = norm(diag(matH));
	if ( 0 == dat.hScale )
		error( "Diagonal of Hessian matrix is all zeros." );
	endif
	dat.cScale = norm(diag(matC));
	if ( 0 == dat.cScale )
		error( "Constraint matrix is singular." );
	endif
	dat.epsReguScaled = prm.epsReguRel * dat.hScale / dat.cScale;
	dat = overwritefields( dat, prmIn );
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		sz = size( vecG, 1 );
		assert( isrealarray(vecG,[sz,1]) );
		assert( isrealarray(matH,[sz,sz]) );
		assert( issymmetric(matH) );
		assert( min(diag(matH)) >= 0.0 );
		assert( max(diag(matH)) >= max(max(matH)) );
		if (~isempty(bTrgt))
			assert( isscalar(bTrgt) );
			assert( 0.0 < bTrgt );
		endif
		szb = size(matB,1);
		assert( isrealarray(matB,[szb,sz]) );
		if ( szb > sz )
			msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, ...
			  "NOTE: Leading dimesion of boundary matrix is greater than problem size." );
			msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, ...
			  "  Replacing the input matB with chol(matB'*matB) might provide a speed-up." );
		endif
		%
		assert( isrealscalar(prm.verbLev) );
		assert( isrealscalar(prm.valdLev) );
		assert( isrealscalar(prm.bRelTol) );
		assert( 0.0 < prm.bRelTol );
		assert( prm.bRelTol < 1.0 );
		assert( isrealscalar(prm.cholRelTol) );
		assert( 0.0 <= prm.cholRelTol );
		assert( prm.cholRelTol <= 1.0 );
		assert( 0.0 < prm.epsReguRel );
		assert( prm.epsReguRel <= 1.0 );
		assert( isrealscalar(prm.iterMax) );
		assert( abs(prm.iterMax-round(prm.iterMax)) < sqrt(eps) );
		assert( 1 <= prm.iterMax );
		%
		assert( isrealarray(dat.matC,[sz,sz]) );
		assert( issymmetric(dat.matC) );
		assert( min(diag(dat.matC)) >= 0.0 );
		assert( max(diag(dat.matC)) >= max(max(dat.matC)) );
		assert( isrealscalar(dat.epsReguScaled) );
		assert( 0.0 < dat.epsReguScaled );
	endif
	if ( prm.valdLev >= VALDLEV__MEDIUM );
		assert( reldiff(dat.matC,matB'*matB) < sqrt(eps) );
	endif
	if ( prm.valdLev >= VALDLEV__UNLIMITED )
		eigH = eig(matH);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matH): %0.3e ~ %0.3e", min(eigH), max(eigH) ) );
		eigC = eig(dat.matC);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matC): %0.3e ~ %0.3e", min(eigC), max(eigC) ) );
		%
		if ( min(eigH) < -sqrt(eps)*max(abs(eigH)) )
			error( "Hessian matrix has a clearly negative eigenvalue." );
		elseif ( min(eigC) < -sqrt(eps)*max(abs(eigC)) )
			error( "Constraint matrix has a clearly negative eigenvalue." );
		endif
		if ( min(eigC) <= 0.0 )
			msgif ( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, ...
			  "WARNING: Constraint matrix appears to be non-positive-definite." );
		endif
	endif
	return;
endfunction


% If chol() failes, use linear extrapolation.
% Since calculation of dy/dmu would require an additional backsub, return merely a function to allow its calculation.
function [ b, bPrime, vecY, funchYPrime ] = __calcNewt( vecG, matH, matB, prm, dat )
	mydefs;
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > prm.cholRelTol * max(abs(diag(matR))) )
		[ b, bPrime, vecY, vecRho ] = __calcFromChol( vecG, matR, matB, dat.matC );
		funchYPrime = @()( matR \ vecRho );
		return;
	endif
	endif
	clear matR;
	%
	msgif( prm.verbLev >= VERBLEV__INFO, __FILE__, __LINE__, "Using extrapolation for Newton step." );
	% Note that epsReguScaled should be scaled for matC.
	[ matR1, cholFlag1 ] = chol( matH + dat.epsReguScaled * dat.matC );
	if ( 0 ~= cholFlag1 )
		error( "Cholesky factorization failed even with regularization matrix." );
	endif
	[ matR2, cholFlag2 ] = chol( matH + (2.0*dat.epsReguScaled) * dat.matC );
	if ( 0 ~= cholFlag2 )
		error( "Cholesky factorization (somehow) failed with regularization matrix second time." );
	endif
	%
	[ b1, bPrime1, vecY1, vecRho1 ] = __calcFromChol( vecG, matR1, matB, dat.matC );
	[ b2, bPrime2, vecY2, vecRho2 ] = __calcFromChol( vecG, matR2, matB, dat.matC );
	b = (2.0*b1) - b2;
	bPrime = (2.0*bPrime1) - bPrime2;
	vecY = (2.0*vecY1) - vecY2;
	funchYPrime = @()(  (2.0*( matR1 \ vecRho1 )) - ( matR2 \ vecRho2 )  );
	return;
endfunction


% Math...
%  C = B^T * B
%  M = H + mu*C
%  y = M^-1 * (-g)
%  b = || B *y ||
%  dM/dmu = - M^-1 * C * M^-1
%  dy/dmu = - M^-1 * C * M^-1 * (-g) = - M^-1 * C * y
%  d(b^2)/dmu = -2 * y^T * C * dy/dmu = -2 * || M^-(1/2) * C * y ||^2
%  db/dmu = (d(b^2)/dmu) / ( 2 * b ) = || M^-(1/2) * C * y ||^2 / b
% Use rho = M^(-1/2) * C * y. (I forget why it's called "rho".)
% Use M = R^T * R.
% Note that dy/dmu is not calculated here, but can be calculated from rho and R.
function [ b, bPrime, vecY, vecRho ] = __calcFromChol( vecG, matR, matB, matC )
	vecY = matR \ ( matR' \ (-vecG) );
	b = norm( matB * vecY );
	vecRho = matR' \ ( matC * vecY );
	bPrime = -sumsq( vecRho ) / b;
	return;
endfunction


function [ vecY, datOut ] = __find( bNewt, bPrimeNewt, vecG, matH, bTrgt, matB, prm, dat )
	mydefs;
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( 0.0 < bTrgt );
		assert( bTrgt < bNewt );
		assert( bPrimeNewt < 0.0 );
	endif
	%
	mu0 = 0.0;
	b0 = bNewt;
	bPrime0 = bPrimeNewt;
	haveFinite1 = false;
	mu1 = +Inf;
	b1 = 0.0;
	bPrime1 = 0.0;
	iterCount = 0;
	while (1)
		if ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "IMPOSED STOP: Reached iterMax (%d).", prm.iterMax ) );
			datOut.retCode = RETCODE__IMPOSED_STOP;
			break;
		endif
		iterCount++;
		%
		if (haveFinite1)
			msg( __FILE__, __LINE__, "haveFinite1 is unsupported." );
			break;
		else
			mu = mu0 - ( (b0/bTrgt) - 1.0 ) * b0 / bPrime0;
		endif
		%
		[ mu0, mu, mu1 ]
		if ( mu <= mu0 || mu >= mu1 )
			error( "NUMERICAL ISSUE: Guess went out of bounds." );
		endif
		[ b, bPrime, vecY, funchYPrime ] = __calc( mu, vecG, matH, matB, prm, dat );
		if ( abs(b-bTrgt) < bTrgt*prm.bRelTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
			datOut.retCode = RETCODE__SUCCESS;
			break;
		endif
		if ( b >= b0 || b <= b1 || bPrime >= 0.0 )
			error( "NUMERICAL ISSUE: Function became non-monotonic." );
		endif
		%
		if ( b < bTrgt )
			haveFinite1 = true;
			mu1 = mu;
			b1 = b;
			bPrime1 = bPrime;
		else
			mu0 = mu;
			b0 = b;
			bPrime0 = bPrime;
		endif
	endwhile
	%
	datOut.mu = mu;
	datOut.b = b;
	datOut.bPrime = bPrime;
	datOut.vecYPrime = funchYPrime();
	datOut.iterCount = iterCount;
	return;
endfunction



function [ b, bPrime, vecY, funchYPrime ] = __calc( mu, vecG, matH, matB, prm, dat )
	matR = chol( matH + mu*dat.matC );
	vecY = matR \ ( matR' \ (-vecG) );
	b = norm( matB * vecY );
	vecRho = matR' \ ( dat.matC * vecY );
	bPrime = -sumsq( vecRho ) / b;
	funchYPrime = @()( matR \ vecRho );
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
