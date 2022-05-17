% Function...

function [ vecY, vecYPrime, b, bPrime ] = findLevPt( vecG, matH, bMax=[], matB=[], prm=[] )
	sz = size( vecG );
	tLo = 0.0;
	tHi = 1.0;
	%
	if ( isempty(matB) )
		matB = eye(sz,sz);
	endif
	matBTB = mygetfield( prm, "matBTB", [] );
	if ( isempty(matBTB) )
		matBTB = matB' * matB;
	endif
	matRegu = mygetfield( prm, "matRegu", [] );
	if ( isempty(matRegu) )
		matRegu = (eps^0.4)*diag(abs(diag(matBTB))) + (eps^0.7)*eye(sz,sz);
	endif
	%
	%
	if ( mygetfield( prm, "useDogLeg", false ) )
		error( "Not implemented!" );
	endif
	%
	%
	cholRelThresh = mygetfield( prm, "cholRelThresh", sqrt(eps) );
	bTol = mygetfield( prm, "bTol", 0.1 );
	%
	[ b, bPrime, vecY, funchYPrime ] = __calc( tHi, vecG, matH, matB, matBTB, matRegu, cholRelThresh );
	if ( isempty(bMax) )
		vecYPrime = funchYPrime();
		return;
	elseif ( b <= bMax + bTol )
		vecYPrime = funchYPrime();
		return;
	endif
	bHi = b;
	bPrimeHi = bPrime;
	%
	[ b, bPrime, vecY, funchYPrime ] = __calc( tLo, vecG, matH, matB, matBTB, matRegu, cholRelThresh );
	if ( b >= bMax - bTol )
		if ( b > bMax + bTol )
			msg( __FILE__, __LINE__, "WARNING: tLo was too large for bMax; returning point for tLo." );
		endif
		vecYPrime = funchYPrime();
		return;
	endif
	bLo = b;
	bPrimeLo = bPrime;
	%
	iterMax = mygetfield( prm, "iterMax", 100 );
	for n = 1 : iterMax
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
		if ( deltaTHi < deltaTLo )
			t = tHi - deltaTHi;
		else
			t = tLo + deltaTLo;
		endif
		assert( t > tLo );
		assert( t < tHi );
		%tVals = [ tLo, tLo + deltaTLo, t, tHi - deltaTHi, tHi ]
		[ b, bPrime, vecY, funchYPrime ] = __calc( t, vecG, matH, matB, matBTB, matRegu, cholRelThresh );
		%bVals = [ bLo, b, bHi, bMax ]
		if ( abs( b - bMax ) <= bTol )
			vecYPrime = funchYPrime();
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


% If the initial Cholesky factorization fails (or is skipped),
%  try it again while adding 1 or 2 times E (the regularization matrix),
%  then extrapolate to 0 times E.
% Since d/dt( y ) is typically not needed,
%  return only a function handle allowing it to be calculated. (This works, right?)
function [ b, bPrime, vecY, funchYPrime ] = __calc( t, vecG, matH, matB, matS, matE, cholRelThresh )
	if ( 0.0 < cholRelThresh && cholRelThresh < 1.0 )
		[ matR, cholFlag ] = chol( t*matH + (1.0-t)*matS );
		if ( 0 == cholFlag )
		if ( min(diag(matR)) > cholRelThresh * max(abs(diag(matR))) )
			[ vecY, b, bPrime, vecRho ] = __calcGivenR( t, vecG, matR, matB, matS );
			funchYPrime = @()( matR \ vecRho );
			return;
		endif
		endif
		clear matR;
	endif
	msg( __FILE__, __LINE__, "Extrapolating!" );
	%
	[ matR1, cholFlag ] = chol( t*matH + (1.0-t)*matS + matE );
	if ( 0 ~= cholFlag )
		error( "Cholesky factorization failed even with regularization matrix." );
	endif
	[ matR2, cholFlag ] = chol( t*matH + (1.0-t)*matS + (2.0*matE) );
	if ( 0 ~= cholFlag )
		error( "Cholesky factorization failed with regularization matrix second time." );
	endif
	[ vecY1, b1, bPrime1, vecRho1 ] = __calcGivenR( t, vecG, matR1, matB, matS );
	[ vecY2, b2, bPrime2, vecRho2 ] = __calcGivenR( t, vecG, matR2, matB, matS );
	vecY = (2.0*vecY1) - vecY2;
	b = (2.0*b1) - b2;
	bPrime = (2.0*bPrime1) - bPrime2;
	funchYPrime = @()(  (2.0*( matR1 \ vecRho1 )) - ( matR2 \ vecRho2 )  );
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
function [ vecY, b, bPrime, vecRho ] = __calcGivenR( t, vecG, matR, matB, matS )
	vecGamma = matR \ ( matR' \ (-vecG) );
	b0 = norm( matB * vecGamma );
	vecRho = matR' \ ( matS * vecGamma );
	vecY = t * vecGamma;
	b = t * b0;
	bPrime = sumsq( vecRho ) / b0;
return;
endfunction


%!test
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 200;
%!	sizeF = 100;
%!	%
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX);
%!	if (1)
%!		vecPhi = randn(sizeX,1);
%!		vecPhi /= norm(vecPhi);
%!		matJ -= (matJ*vecPhi)*(vecPhi');
%!	endif
%!	matB_gen = randn(sizeF,sizeX);
%!	matB_unscaled = matB_gen'*matB_gen;
%!	%matB_unscaled = diag(diag(matJ'*matJ));
%!	bMax_unscaled = sqrt(2.0);
%!	%
%!	%
%!	matH = matJ'*matJ;
%!	vecG = matJ'*vecF;
%!	%
%!	% Regularize things?
%!	matB_unscaled += sqrt(eps)*max(diag(matB_unscaled))*eye(sizeX,sizeX);
%!	matH += sqrt(eps)*max(diag(matH))*eye(sizeX,sizeX);
%!	%
%!	% Scale B?
%!	hobScale =  sqrt(max(diag(matH))) / max(diag(matB_unscaled));
%!	matB_scaled = matB_unscaled * hobScale;
%!	bMax_scaled = bMax_unscaled * hobScale;
%!	%
%!	%
%!	prm = [];
%!	prm.bTol = 0.001;
%!	[ vecY, vecYPrime, b, bPrime ] = findLevPt( vecG, matH, bMax_scaled, matB_scaled, prm );
%!	[ norm(matB_scaled*vecY), bMax_scaled, norm(matB_scaled*vecY) - bMax_scaled ]
%!	[ norm(matB_unscaled*vecY), bMax_unscaled, norm(matB_unscaled*vecY) - bMax_unscaled ]
