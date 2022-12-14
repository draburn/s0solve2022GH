% Function...

function [ vecY, vecYPrime, b, bPrime, iterCount ] = findLevPt_basic( vecG, matH, bMaxTrgt=[], matB=[], prm=[] )
	% PARSE INPUT.
	sz = size( vecG, 1 );
	assert( isrealarray(vecG,[sz,1]) );
	assert( isrealarray(matH,[sz,sz]) );
	assert( issymmetric(matH) );
	if ( ~isempty(bMaxTrgt) )
		assert( isrealscalar(bMaxTrgt) );
		assert( 0.0 < bMaxTrgt );
	endif
	if ( isempty(matB) )
		hScale = norm(diag(matH));
		matB = diag(sqrt( diag(matH) + eps*hScale ));
	endif
	assert( isrealarray(matB) ); % First dim of matB is actually free.
	assert( size(matB,2) == sz );
	assert( size(size(matB),2) == 2 ); % True even is matB is a scalar or vector.
	%
	bRelTol = mygetfield( prm, "bRelTol", 0.01 );
	assert( isrealscalar(bRelTol) );
	assert( 0.0 < bRelTol );
	assert( bRelTol < 1.0 );
	bTol = bMaxTrgt * bRelTol;
	%
	iterMax = mygetfield( prm, "iterMax", 100 );
	assert( isrealscalar(iterMax) );
	assert( abs( iterMax - round(iterMax) ) < sqrt(eps) );
	assert( 1 <= iterMax );
	%
	doMsg = mygetfield( prm, "doMsg", true );
	assert( islogical(doMsg) );
	assert( isscalar(doMsg) );
	%
	%
	% DO PREP.
	iterCount = 0;
	matC = matB'*matB;
	%
	tHi = 1.0;
	[ b, bPrime, vecY, vecYPrime ] = __calc( tHi, vecG, matH, matB, matC, prm );
	if ( isempty(bMaxTrgt) )
		msgif( doMsg, __FILE__, __LINE__, "SUCCESS: Full step was requested." );
		return;
	elseif ( b <= bMaxTrgt + bTol )
		msgif( doMsg, __FILE__, __LINE__, "SUCCESS: Full step is within tolerance." );
		return;
	endif
	bHi = b;
	bPrimeHi = bPrime;
	%
	tLo = 0.0;
	[ b, bPrime, vecY, vecYPrime ] = __calc( tLo, vecG, matH, matB, matC, prm );
	bLo = b; % Should be zero.
	bPrimeLo = bPrime;
	%
	%
	% DO WORK.
	for iterCount = 1 : iterMax
		% Find next guess according to either linear model, but cap to midpoint.
		if ( bPrimeLo * ( tHi - tLo ) <= 2.0 * ( bMaxTrgt - bLo ) )
			deltaTLo = ( tHi - tLo ) / 2.0;
		else
			deltaTLo = ( bMaxTrgt - bLo ) / bPrimeLo;
			if ( deltaTLo < sqrt(eps) )
				deltaTLo = sqrt(eps);
			endif
		endif
		if ( bPrimeHi * ( tHi - tLo ) <= 2.0 * ( bHi - bMaxTrgt ) )
			deltaTHi = ( tHi - tLo ) / 2.0;
		else
			deltaTHi = ( bHi - bMaxTrgt ) / bPrimeHi;
			if ( deltaTHi < sqrt(eps) )
				deltaTHi = sqrt(eps);
			endif
		endif
		deltaTLo = cap( deltaTLo, 0.0, (tHi-tLo)/2.0 );
		deltaTHi = cap( deltaTHi, 0.0, (tHi-tLo)/2.0 );
		%
		% Take whichever guess one is closer to the generating point.. ish.
		if ( deltaTHi < deltaTLo )
			t = tHi - deltaTHi;
		else
			t = tLo + deltaTLo;
		endif
		if ( 0 == mod(iterCount,3) )
			t = (tLo+tHi)/2.0;
			%t = (tLo+tHi+deltaTLo-deltaTHi)/2.0;
		endif
		assert( t > tLo ); % A few checks, just to be safe.
		assert( t < tHi );
		[ b, bPrime, vecY, vecYPrme ] = __calc( t, vecG, matH, matB, matC, prm );
		if ( abs( b - bMaxTrgt ) <= bTol )
			msgif( doMsg, __FILE__, __LINE__, sprintf( "SUCCESS: Converged after %d iterations.", iterCount ) );
			return;
		endif
		msgif( doMsg, __FILE__, __LINE__, sprintf( "  iter %d:  t: %g ~ %g ~ %g;  (1-t: %g ~ %g ~ %g);  b: %g ~ %g ~ %g (%g).", ...
		  iterCount, tLo, t, tHi, 1.0-tLo, 1.0-t, 1.0-tHi, bLo, b, bHi, bMaxTrgt ) );
		if ( b <= bLo || b >= bHi )
			error( "ALGORITHM BREAKDOWN: Iteration went out of bounds; matrices may be poorly scaled." );
		endif
		%
		if ( b < bMaxTrgt )
			assert( t > tLo );
			tLo = t;
			bLo = b;
			bPrimeLo = bPrime;
		else
			assert( t < tHi );
			tHi = t;
			bHi = b;
			bPrimeHi = bPrime;
		endif
	endfor
	error( "IMPOSED STOP: Reached iteration limit."  );
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
function [ b, bPrime, vecY, vecYPrime ] = __calc( t, vecG, matH, matB, matC, prm )
	[ matR, cholFlag ] = chol( t*matH + (1.0-t)*matC );
	if ( 0 ~= cholFlag )
		error( "ALGORITHM BREAKDOWN: Cholesky factorization failed; please ensure matrices are positive definite." );
	endif
	vecGamma = matR \ ( matR' \ (-vecG) );
	b0 = norm( matB * vecGamma );
	vecRho = matR' \ ( matC * vecGamma );
	vecY = t * vecGamma;
	b = t * b0;
	bPrime = sumsq( vecRho ) / b0;
	vecYPrime = matR \ vecRho; % Not used internally, but user may want this.
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
%!	%	
%!	bMaxTrgt = 0.1;
%!	prm = [];
%!	prm.doMsg = true;
%!	[ vecY, vecYPrime, b, bPrime ] = findLevPt_basic( vecG, matH, bMaxTrgt, matB, prm )
%!	bActual = norm(matB*vecY);
%!	msg( __FILE__, __LINE__, sprintf( "bDesired = %g,  bReported = %g,  bActual = %g,  residual = %g.", ...
%!	  bMaxTrgt, b, bActual, bActual - bMaxTrgt ) );


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
%!	%	
%!	bMaxTrgt = 0.1;
%!	prm = [];
%!	prm.doMsg = true;
%!	[ vecY, vecYPrime, b, bPrime ] = findLevPt_basic( vecG, matH, bMaxTrgt, matB, prm );
%!	bActual = norm(matB*vecY);
%!	msg( __FILE__, __LINE__, sprintf( "bDesired = %g,  bReported = %g,  bActual = %g,  residual = %g.", ...
%!	  bMaxTrgt, b, bActual, bActual - bMaxTrgt ) );
