% Function...

function [ vecY, vecYPrime, b, bPrime ] = findLevPt_basic( vecG, matH, bMax=[], matB=[], prm=[] )
	% PARSE INPUT.
	sz = size( vecG, 1 );
	assert( isrealarray(vecG,[sz,1]) );
	assert( isrealarray(matH,[sz,sz]) );
	assert( issymmetric(matH) );
	if ( ~isempty(bMax) )
		assert( isrealscalar(bMax) );
		assert( 0.0 < bMax );
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
	bTol = bMax * bRelTol;
	%
	iterMax = mygetfield( prm, "iterMax", 100 );
	assert( isrealscalar(iterMax) );
	assert( abs( iterMax - round(iterMax) ) < sqrt(eps) );
	assert( 1 <= iterMax );
	%
	doMsg = mygetfield( prm, "doMsg", false );
	assert( islogical(doMsg) );
	assert( isscalar(doMsg) );
	%
	%
	% START WORK.
	matC = matB'*matB;
	%
	tHi = 1.0;
	[ b, bPrime, vecY, vecYPrime ] = __calc( tHi, vecG, matH, matB, matC, prm );
	if ( isempty(bMax) )
		msgif( doMsg, __FILE__, __LINE__, "SUCCESS: Full step was requested." );
		return;
	elseif ( b <= bMax + bTol )
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
	for iterCount = 1 : iterMax
		% Find next t according to either model, but cap to midpoint.
		if ( bPrimeLo * ( tHi - tLo ) <= 2.0 * ( bMax - bLo ) )
			deltaTLo = ( tHi - tLo ) / 2.0;
		else
			deltaTLo = ( bMax - bLo ) / bPrimeLo;
		endif
		if ( bPrimeHi * ( tHi - tLo ) <= 2.0 * ( bHi - bMax ) )
			deltaTHi = ( tHi - tLo ) / 2.0;
		else
			deltaTHi = ( bHi - bMax ) / bPrimeHi;
		endif
		%
		% Take which one is closer to the current point.
		if ( deltaTHi < deltaTLo )
			t = tHi - deltaTHi;
		else
			t = tLo + deltaTLo;
		endif
		assert( t > tLo );
		assert( t < tHi );
		[ b, bPrime, vecY, vecYPrme ] = __calc( t, vecG, matH, matB, matC, prm );
		if ( abs( b - bMax ) <= bTol )
			msgif( doMsg, __FILE__, __LINE__, sprintf( "SUCCESS: Converged after %d iterations.", iterCount ) );
			return;
		endif
		if ( b <= bLo || b >= bHi )
			error( "Iteration went out of bounds; the system may be poorly scaled." );
		endif
		%
		if ( b < bMax )
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
	error( "Reached iteration limit."  );
return;
endfunction


% Here's how the math works...
%  S  ==  B^T * B
%  M  ==  t*H + (1-t)*S  =  R^T * R
%  d/dt( M^-1 )  =  -M^-1 * ( H - S ) * M^-1
%  gamma  ==  M^-1 * (-g)
%  y  =  t * gamma
%  b  =  || B * y ||  =  t * sqrt( gamma^T * S * gamma )
%  d/dt( y )  =  M^-1 * S * gamma
%  d/dt( b )  =  (S * gamma)^T * M^-1 * (S * gamma)
%
% Using Cholesky factorization of M, and being efficient...
%  M = R^T * R
%  gamma = R \ ( R^T \ (-g) )
%  b0 = || B * gamma ||
%  rho = R^T \ ( S * gamma )
%  y = t * gamma;
%  b = t * b0;
%  d/dt( b ) = || rho ||^2 / b0
%  d/dt( y ) = R \ rho
% But, since d/dt(y) is typically not needed and requires an extra backsub, we return R and rho.
function [ b, bPrime, vecY, vecYPrime ] = __calc( t, vecG, matH, matB, matC, prm )
	[ matR, cholFlag ] = chol( t*matH + (1.0-t)*matC );
	if ( 0 ~= cholFlag )
		error( "Cholesky factorization failed; please ensure matrices are positive definite." );
	endif
	vecGamma = matR \ ( matR' \ (-vecG) );
	b0 = norm( matB * vecGamma );
	vecRho = matR' \ ( matC * vecGamma );
	vecY = t * vecGamma;
	b = t * b0;
	bPrime = sumsq( vecRho ) / b0;
	vecYPrime = matR \ vecRho;
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
%!	bMax = 0.1;
%!	prm = [];
%!	prm.doMsg = true;
%!	[ vecY, vecYPrime, b, bPrime ] = findLevPt_basic( vecG, matH, bMax, matB, prm )
%!	bActual = norm(matB*vecY);
%!	msg( __FILE__, __LINE__, sprintf( "bDesired = %g,  bReported = %g,  bActual = %g,  residual = %g.", ...
%!	  bMax, b, bActual, bActual - bMax ) );


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
%!	bMax = 0.1;
%!	prm = [];
%!	prm.doMsg = true;
%!	[ vecY, vecYPrime, b, bPrime ] = findLevPt_basic( vecG, matH, bMax, matB, prm );
%!	bActual = norm(matB*vecY);
%!	msg( __FILE__, __LINE__, sprintf( "bDesired = %g,  bReported = %g,  bActual = %g,  residual = %g.", ...
%!	  bMax, b, bActual, bActual - bMax ) );
