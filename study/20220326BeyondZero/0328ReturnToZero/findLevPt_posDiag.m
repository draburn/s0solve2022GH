% Function...
%  findLev for positive (positive definite) diagonal matrices.
%  Had hoped this would provide a good initial guess for the general case,
%   but, that doesn't seem to be happening.

function [ vecY, vecYPrime, b, bPrime, t ] = findLevPt_posDiag( vecG, matDH, bMax, matB, matDC, prm=[] )
	% PARSE INPUT.
	sz = size( vecG, 1 );
	assert( isrealarray(vecG,[sz,1]) );
	assert( isrealarray(matDH,[sz,sz]) );
	assert( isdiag(matDH) );
	assert( min(diag(matDH)) > 0.0 );
	assert( isrealscalar(bMax) );
	assert( 0.0 < bMax );
	assert( isrealarray(matB) ); % First dim of matB is actually free.
	assert( size(matB,2) == sz );
	assert( size(size(matB),2) == 2 ); % True even is matB is a scalar or vector.
	assert( isrealarray(matDC,[sz,sz]) );
	assert( isdiag(matDC) );
	assert( min(diag(matDC)) > 0.0 );
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
	% DO PREP.
	t = 1.0;
	[ b, bPrime, vecY, vecYPrime ] = __calc( t, vecG, matDH, matB, matDC, prm );
	if ( isempty(bMax) )
		msgif( doMsg, __FILE__, __LINE__, "SUCCESS: Full step was requested." );
		return;
	elseif ( b <= bMax + bTol )
		msgif( doMsg, __FILE__, __LINE__, "SUCCESS: Full step is within tolerance." );
		return;
	endif
	tHi = t;
	bHi = b;
	bPrimeHi = bPrime;
	%
	t = 0.0;
	[ b, bPrime, vecY, vecYPrime ] = __calc( t, vecG, matDH, matB, matDC, prm );
	tLo = t;
	bLo = b; % Should be zero.
	bPrimeLo = bPrime;
	%
	%
	% DO WORK.
	for iterCount = 1 : iterMax
		% Find next guess according to either linear model, but cap to midpoint.
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
		% Take whichever guess one is closer to the generating point.
		if ( deltaTHi < deltaTLo )
			t = tHi - deltaTHi;
		else
			t = tLo + deltaTLo;
		endif
		assert( t > tLo ); % A few checks, just to be safe.
		assert( t < tHi );
		[ b, bPrime, vecY, vecYPrme ] = __calc( t, vecG, matDH, matB, matDC, prm );
		if ( abs( b - bMax ) <= bTol )
			msgif( doMsg, __FILE__, __LINE__, sprintf( "SUCCESS: Converged after %d iterations.", iterCount ) );
			return;
		endif
		if ( b <= bLo || b >= bHi )
			error( "ALGORITHM BREAKDOWN: Iteration went out of bounds; matrices may be poorly scaled." );
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
	error( "IMPOSED STOP: Reached iteration limit."  );
return;
endfunction


% Here's how the math works for the general pos-def symmetric case.
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
function [ b, bPrime, vecY, vecYPrime ] = __calc( t, vecG, matDH, matB, matDC, prm )
	matM = t*matDH + (1.0-t)*matDC;
	vecGamma = matM \ (-vecG);
	b0 = norm( matB * vecGamma );
	vecY = t * vecGamma;
	vecRho = sqrt(matM) \ ( matDC * vecGamma );
	vecY = t * vecGamma;
	b = t * b0;
	bPrime = sumsq( vecRho ) / b0;
	vecYPrime = matM \ ( matDC * vecGamma ); % Not used internally, but user may want this.
return;
endfunction


%!test
%!	%msg( __FILE__, __LINE__, "Skipping this test." ); return;
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
%!	matC = matB'*matB;
%!	matDH = diag(diag(matH));
%!	matDC = diag(diag(matC));
%!	hScale = norm( matB*(matDH\vecG) ) / norm( matB*(matH\vecG) );
%!	cScale = norm( matB*(matDC\vecG) ) / norm( matB*(matC\vecG) );
%!	matDH *= hScale;
%!	matDC *= cScale;
%!	%[ norm(matB*(matC\vecG)), norm(matB*(matDC\vecG)) ]
%!	%[ norm(matB*(matH\vecG)), norm(matB*(matDH\vecG)) ]
%!	%
%!	bMax = 0.1;
%!	prm = [];
%!	prm.doMsg = true;
%!	[ vecY, vecYPrime, b, bPrime, t ] = findLevPt_posDiag( vecG, matDH, bMax, matB, matDC, prm );
%!	bActual = norm(matB*vecY);
%!	vecYShouldBe = ( t*matDH + (1.0-t)*matDC ) \ (-t*vecG);
%!	vecYWouldBe = ( t*matH + (1.0-t)*matC ) \ (-t*vecG);
%!	bWouldBe = norm( matB*vecYWouldBe );
%!	bShouldBe = norm( matB*vecYShouldBe );
%!	if ( sizeX <= 5 )
%!		g_ddh_ddc_y_ywb_yp  = [ vecG, diag(matDH), diag(matDC), vecY, vecYShouldBe, vecYWouldBe, vecYPrime ]
%!	endif
%!	msg( __FILE__, __LINE__, sprintf( "bPrime = %g;  t = %g.", bPrime, t ) );
%!	msg( __FILE__, __LINE__, sprintf( "bWouldBe = %f,  bShouldBe = %f,  bReported = %g,  bDesired = %g,  bActual = %g;  residual = %g.", ...
%!	  bWouldBe, bShouldBe, b, bMax, bActual, bActual - bMax ) );
%!	printf( "\n\n" );



%!test
%!	%msg( __FILE__, __LINE__, "Skipping this test." ); return;
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 5;
%!	%
%!	sizeF = sizeX;
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX) .* exp(3.0*randn(sizeF,sizeX)) * exp(3.0*randn());
%!	matB = randn(sizeF,sizeX) .* exp(3.0*randn(sizeF,sizeX)) * exp(3.0*randn());
%!	%
%!	matH = matJ'*matJ;
%!	vecG = matJ'*vecF;
%!	matC = matB'*matB;
%!	matDH = diag(diag(matH));
%!	matDC = diag(diag(matC));
%!	hScale = norm( matB*(matDH\vecG) ) / norm( matB*(matH\vecG) );
%!	cScale = norm( matB*(matDC\vecG) ) / norm( matB*(matC\vecG) );
%!	matDH *= hScale;
%!	matDC *= cScale;
%!	%
%!	bMax = 0.1;
%!	prm = [];
%!	prm.doMsg = true;
%!	[ vecY, vecYPrime, b, bPrime, t ] = findLevPt_posDiag( vecG, matDH, bMax, matB, matDC, prm );
%!	bActual = norm(matB*vecY);
%!	vecYShouldBe = ( t*matDH + (1.0-t)*matDC ) \ (-t*vecG);
%!	vecYWouldBe = ( t*matH + (1.0-t)*matC ) \ (-t*vecG);
%!	bWouldBe = norm( matB*vecYWouldBe );
%!	bShouldBe = norm( matB*vecYShouldBe );
%!	if ( sizeX <= 5 )
%!		g_ddh_ddc_y_ywb_yp  = [ vecG, diag(matDH), diag(matDC), vecY, vecYShouldBe, vecYWouldBe, vecYPrime ]
%!	endif
%!	msg( __FILE__, __LINE__, sprintf( "bPrime = %g;  t = %g.", bPrime, t ) );
%!	msg( __FILE__, __LINE__, sprintf( "bWouldBe = %f,  bShouldBe = %f,  bReported = %g,  bDesired = %g,  bActual = %g;  residual = %g.", ...
%!	  bWouldBe, bShouldBe, b, bMax, bActual, bActual - bMax ) );
%!	printf( "\n\n" );



%!test
%!	%msg( __FILE__, __LINE__, "Skipping this test." ); return;
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
%!	matC = matB'*matB;
%!	matDH = diag(diag(matH));
%!	matDC = diag(diag(matC));
%!	hScale = norm( matB*(matDH\vecG) ) / norm( matB*(matH\vecG) );
%!	cScale = norm( matB*(matDC\vecG) ) / norm( matB*(matC\vecG) );
%!	matDH *= hScale;
%!	matDC *= cScale;
%!	%
%!	bMax = 0.1;
%!	prm = [];
%!	prm.doMsg = true;
%!	[ vecY, vecYPrime, b, bPrime, t ] = findLevPt_posDiag( vecG, matDH, bMax, matB, matDC, prm );
%!	bActual = norm(matB*vecY);
%!	vecYShouldBe = ( t*matDH + (1.0-t)*matDC ) \ (-t*vecG);
%!	vecYWouldBe = ( t*matH + (1.0-t)*matC ) \ (-t*vecG);
%!	bWouldBe = norm( matB*vecYWouldBe );
%!	bShouldBe = norm( matB*vecYShouldBe );
%!	if ( sizeX <= 5 )
%!		g_ddh_ddc_y_ywb_yp  = [ vecG, diag(matDH), diag(matDC), vecY, vecYShouldBe, vecYWouldBe, vecYPrime ]
%!	endif
%!	msg( __FILE__, __LINE__, sprintf( "bPrime = %g;  t = %g.", bPrime, t ) );
%!	msg( __FILE__, __LINE__, sprintf( "bWouldBe = %f,  bShouldBe = %f,  bReported = %g,  bDesired = %g,  bActual = %g;  residual = %g.", ...
%!	  bWouldBe, bShouldBe, b, bMax, bActual, bActual - bMax ) );
%!	printf( "\n\n" );



%!test
%!	%msg( __FILE__, __LINE__, "Skipping this test." ); return;
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 500;
%!	%
%!	sizeF = sizeX;
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX).*exp(3.0*randn(sizeF,sizeX));
%!	matB = randn(sizeF,sizeX).*exp(3.0*randn(sizeF,sizeX))*exp(3.0*randn());
%!	%
%!	matH = matJ'*matJ;
%!	vecG = matJ'*vecF;
%!	matC = matB'*matB;
%!	matDH = diag(diag(matH));
%!	matDC = diag(diag(matC));
%!	hScale = norm( matB*(matDH\vecG) ) / norm( matB*(matH\vecG) );
%!	cScale = norm( matB*(matDC\vecG) ) / norm( matB*(matC\vecG) );
%!	matDH *= hScale;
%!	matDC *= cScale;
%!	%
%!	bMax = 0.1;
%!	prm = [];
%!	prm.doMsg = true;
%!	[ vecY, vecYPrime, b, bPrime, t ] = findLevPt_posDiag( vecG, matDH, bMax, matB, matDC, prm );
%!	bActual = norm(matB*vecY);
%!	vecYShouldBe = ( t*matDH + (1.0-t)*matDC ) \ (-t*vecG);
%!	vecYWouldBe = ( t*matH + (1.0-t)*matC ) \ (-t*vecG);
%!	bWouldBe = norm( matB*vecYWouldBe );
%!	bShouldBe = norm( matB*vecYShouldBe );
%!	if ( sizeX <= 5 )
%!		g_ddh_ddc_y_ywb_yp  = [ vecG, diag(matDH), diag(matDC), vecY, vecYShouldBe, vecYWouldBe, vecYPrime ]
%!	endif
%!	msg( __FILE__, __LINE__, sprintf( "bPrime = %g;  t = %g.", bPrime, t ) );
%!	msg( __FILE__, __LINE__, sprintf( "bWouldBe = %f,  bShouldBe = %f,  bReported = %g,  bDesired = %g,  bActual = %g;  residual = %g.", ...
%!	  bWouldBe, bShouldBe, b, bMax, bActual, bActual - bMax ) );
%!	printf( "\n\n" );



%!test
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 500;
%!	%
%!	sizeF = sizeX;
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX).*exp(3.0*randn(sizeF,sizeX));
%!	matB = randn(sizeF,sizeX).*exp(3.0*randn(sizeF,sizeX))*exp(3.0*randn());
%!	matJ = diag(diag(matJ));
%!	matB = diag(diag(matB));
%!	%
%!	matH = matJ'*matJ;
%!	vecG = matJ'*vecF;
%!	matC = matB'*matB;
%!	matDH = diag(diag(matH));
%!	matDC = diag(diag(matC));
%!	hScale = norm( matB*(matDH\vecG) ) / norm( matB*(matH\vecG) );
%!	cScale = norm( matB*(matDC\vecG) ) / norm( matB*(matC\vecG) );
%!	matDH *= hScale;
%!	matDC *= cScale;
%!	assert( reldiff(matH,matDH) < sqrt(eps) );
%!	assert( reldiff(matC,matDC) < sqrt(eps) );
%!	%	
%!	bMax = 0.1;
%!	prm = [];
%!	prm.doMsg = true;
%!	[ vecY, vecYPrime, b, bPrime, t ] = findLevPt_posDiag( vecG, matDH, bMax, matB, matDC, prm );
%!	bActual = norm(matB*vecY);
%!	vecYShouldBe = ( t*matDH + (1.0-t)*matDC ) \ (-t*vecG);
%!	vecYWouldBe = ( t*matH + (1.0-t)*matC ) \ (-t*vecG);
%!	bWouldBe = norm( matB*vecYWouldBe );
%!	bShouldBe = norm( matB*vecYShouldBe );
%!	if ( sizeX <= 5 )
%!		g_ddh_ddc_y_ywb_yp  = [ vecG, diag(matDH), diag(matDC), vecY, vecYShouldBe, vecYWouldBe, vecYPrime ]
%!	endif
%!	msg( __FILE__, __LINE__, sprintf( "bPrime = %g;  t = %g.", bPrime, t ) );
%!	msg( __FILE__, __LINE__, sprintf( "bWouldBe = %f,  bShouldBe = %f,  bReported = %g,  bDesired = %g,  bActual = %g;  residual = %g.", ...
%!	  bWouldBe, bShouldBe, b, bMax, bActual, bActual - bMax ) );
%!	printf( "\n\n" );
