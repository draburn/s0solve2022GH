% Function...

function [ vecY, vecYPrime, b, bPrime, n ] = findLevPt( vecG, matH, bMax=[], matB=[], prm=[] )
	sz = size( vecG );
	tLo = 0.0;
	tHi = 1.0;
	n = 0;
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
	
	%[ norm( matH\vecG ), vecG'*matH*vecG, norm( matH*vecG), norm( matH*randn(sz,1) ), norm( matH*randn(sz,1) ), norm( matH*randn(sz,1) ) ]
	%[ norm( matBTB\vecG ), vecG'*matBTB*vecG, norm( matBTB*vecG), norm( matBTB*randn(sz,1) ), norm( matBTB*randn(sz,1) ), norm( matBTB*randn(sz,1) ) ]
	%foo_b1 = norm( matB * (matH\vecG) )
	%foo_b0prime = norm( matB * (matBTB\vecG) )
	%msg( __FILE__, __LINE__, "Goodbye!" ); return;
	
	%
	%
	cholRelThresh = mygetfield( prm, "cholRelThresh", sqrt(eps) );
	bTol = mygetfield( prm, "bTol", 0.01*bMax );
	%
	[ b, bPrime, vecY, funchYPrime ] = __calc( tHi, vecG, matH, matB, matBTB, matRegu, cholRelThresh );
	if ( isempty(bMax) )
		vecYPrime = funchYPrime();
		msg( __FILE__, __LINE__, "Hey!" );
		return;
	elseif ( b <= bMax + bTol )
		vecYPrime = funchYPrime();
		msg( __FILE__, __LINE__, "Hey!" );
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
		b
		bMax
		bTol
		vecYPrime = funchYPrime();
		msg( __FILE__, __LINE__, "Hey!" );
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
		deltaTLo = cap( deltaTLo, 0.0, (tHi-tLo)/2.0 );
		deltaTHi = cap( deltaTHi, 0.0, (tHi-tLo)/2.0 );
		%
		%if ( 1 == n )
		%	% Force an over-step on first iteration.
		%	deltaTHi = max([ 2.0*deltaTHi, (tHi-tLo)*0.0001 ]);
		%	deltaTLo = max([ 2.0*deltaTLo, (tHi-tLo)*0.0001 ]);
		%	deltaTHi = min([ deltaTHi, (tHi-tLo)/2.0 ]);
		%	deltaTLo = min([ deltaTLo, (tHi-tLo)/2.0 ]);
		%endif
		if ( deltaTHi < deltaTLo )
			t = tHi - deltaTHi;
		else
			t = tLo + deltaTLo;
		endif
		%
		if (0)
		if ( 0 == mod(n,3) )
			% ... or, periodically force a bisection?
			t = ( tLo + tHi ) / 2.0;
		endif
		endif
		%
		if (0)
		if ( 1 == n )
			matC = matBTB + sqrt(eps)*diag(diag(matBTB));
			matDH = diag( diag(matH) + eps*norm(diag(matH)) );
			matDC = diag( diag(matC) + eps*norm(diag(matC)) );
			matDH *= norm( matB*(matDH\vecG) ) / norm( matB*(matH\vecG) );
			matDC *= norm( matB*(matDC\vecG) ) / norm( matB*(matC\vecG) );
			posDiag_prm = [];
			posDiag_prm.doMsg = true;
			posDiag_prm = mygetfield( prm, "posDiag_prm", posDiag_prm );
			[ posDiag_vecY, posDiag_vecYPrime, posDiag_b, posDiag_bPrime, t ] = findLevPt_posDiag( vecG, matDH, bMax, matB, matDC, posDiag_prm );
			[ bMax, norm( t*matB*( (t*matDH + (1.0-t)*matDC)\vecG) ), norm( t*matB*( (t*matH + (1.0-t)*matC)\vecG) ) ]
			t
			msg( __FILE__, __LINE__, "HALT!" ); return;
		endif
		endif
		%
		if (0)
		if ( 1 == n )
			matC = matBTB + sqrt(eps)*diag(diag(matBTB));
			matDH = diag( diag(matH) + eps*norm(diag(matH)) );
			matDC = diag( diag(matC) + eps*norm(diag(matC)) );
			matDH *= norm( matB*(matDH\vecG) ) / norm( matB*(matH\vecG) );
			bHalfzies = norm(matB*( (matH+matC)\vecG ));
			posDiag_prm = [];
			posDiag_prm.doMsg = true;
			posDiag_prm = mygetfield( prm, "posDiag_prm", posDiag_prm );
			[ posDiag_vecY, posDiag_vecYPrime, posDiag_b, posDiag_bPrime, t ] = findLevPt_posDiag( vecG, matDH, bHalfzies, matB, matDC, posDiag_prm );
			cScl = (1.0-t)/t
			matDC *= cScl;
			posDiag_prm = [];
			posDiag_prm.doMsg = true;
			posDiag_prm = mygetfield( prm, "posDiag_prm", posDiag_prm );
			[ posDiag_vecY, posDiag_vecYPrime, posDiag_b, posDiag_bPrime, t ] = findLevPt_posDiag( vecG, matDH, bMax, matB, matDC, posDiag_prm );
			[ bMax, norm( t*matB*( (t*matDH + (1.0-t)*matDC)\vecG) ), norm( t*matB*( (t*matH + (1.0-t)*matC)\vecG) ) ]
			t
			%msg( __FILE__, __LINE__, "HALT!" ); return;
		endif
		endif
		%
		if (1)
		if ( 1 == n )
			mu = norm( matB * (matBTB\vecG) ) / bMax;
			t = 1.0/(1.0+mu)
		endif
		endif
		%
		assert( t > tLo );
		assert( t < tHi );
		%tVals = [ tLo, tLo + deltaTLo, t, tHi - deltaTHi, tHi ]
		%t
		[ b, bPrime, vecY, funchYPrime ] = __calc( t, vecG, matH, matB, matBTB, matRegu, cholRelThresh );
		%bVals = [ bLo, b, bHi, bMax ]
		if ( abs( b - bMax ) <= bTol )
			vecYPrime = funchYPrime();
			t
			omt = 1.0 - t
			n
			msg( __FILE__, __LINE__, "Hey!" );
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
	%myquat = b0/norm(vecRho)
	%error( "HALT!" );
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
%!	%hobScale =  sqrt(max(diag(matH))) / max(diag(matB_unscaled));
%!	%hobScale =  1.0e-4*sqrt(max(diag(matH))) / max(diag(matB_unscaled))
%!	%bPrime0 = norm(matB_unscaled*( (matB_unscaled'*matB_unscaled) \ vecG ))
%!	%b1 = norm(matB_unscaled*( matH \ vecG ))
%!	%hobScale =  sqrt(bPrime0/b1)
%!	%hobScale = sqrt( norm(matH\vecG) / norm(matH\(matB_unscaled'*(matB_unscaled*(matH\vecG)))) )
%!	vecHUG = matH\vecG;
%!	vecSHUG = matB_unscaled' * ( matB_unscaled*vecHUG );
%!	if (0)
%!	vecHUSHUG = matH \ vecSHUG;
%!	[ norm( vecSHUG ), norm(vecHUSHUG) ]
%!	shugthushug = vecSHUG' * vecHUSHUG
%!	hobScale2 = norm( matB_unscaled*vecHUG ) / sqrt( shugthushug )
%!	a = sumsq( matB_unscaled * vecHUSHUG )
%!	b = 2.0*( vecSHUG' * vecHUSHUG )
%!	c = sumsq( matB_unscaled * vecHUG ) * 3.0/4.0
%!	discrim = b^2 - 4 * a * c
%!	a*(hobScale2^4)
%!	b*(hobScale2^2)
%!	return;
%!	hobScale = hobScale2
%!	else
%!		matS_foo = matB_unscaled'*matB_unscaled;
%!		matS_foo *= sqrt(eps) * max(diag(matH)) / max(diag(matS_foo));
%!		matR1 = chol( matH + matS_foo );		
%!		matR2 = chol( matH + 2.0*matS_foo );
%!		vecRTUSHUG = (2.0*( matR1' \ vecSHUG )) - ( matR2' \ vecSHUG );
%!		hobScale = norm( matB_unscaled*vecHUG ) / norm( vecRTUSHUG )
%!		%hobScale = 1.0;
%!		%matR = chol( matH );
%!		%vecRTUSHUG = matR' \ vecSHUG;
%!		%hobScale = norm( matB_unscaled*vecHUG ) / norm( vecRTUSHUG )
%!		%return;
%!	endif
%!	matB_scaled = matB_unscaled * hobScale;
%!	bMax_scaled = bMax_unscaled * hobScale;
%!	%rcond( matH + matB_scaled'*matB_scaled )
%!	%
%!	%
%!	prm = [];
%!	prm.bTol = 0.001*bMax_scaled;
%!	[ vecY, vecYPrime, b, bPrime ] = findLevPt( vecG, matH, bMax_scaled, matB_scaled, prm );
%!	[ norm(matB_scaled*vecY), bMax_scaled, norm(matB_scaled*vecY) - bMax_scaled ]
%!	[ norm(matB_unscaled*vecY), bMax_unscaled, norm(matB_unscaled*vecY) - bMax_unscaled ]
