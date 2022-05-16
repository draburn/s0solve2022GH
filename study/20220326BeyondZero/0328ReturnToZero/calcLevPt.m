% Function...

function [ vecY, vecYPrime, s, sPrime ] = calcLevPt( ...
  vecG, matH, t = 1.0, matS = [], matRegu = [], cholTol = sqrt(eps), doChecks = false )
	% Note:
	%  prime = d/dt; may be unreliable if matrix is singular.
	%  s = sqrt( y^T * S * y ).
	%  deltaOmega = g^T * y + 0.5 * y^T * H * y.
	% Maybe inform if singular?
	%
	sz = size(vecG,1);
	if ( doChecks )
		assert( isscalar(doChecks) );
		assert( islogical(doChecks) );
		assert( isrealarray(vecG,[sz,1]) );
		assert( isrealarray(matH,[sz,sz]) );
		assert( issymmetric(matH) );
		assert( isrealscalar(t) );
		assert( 0.0 <= t );
		assert( t <= 1.0 );
	endif
	if ( isempty(matS) )
		matS = max(diag(matH)) * eye(sz,sz);
	endif
	if (isempty(matRegu) )
		matRegu = sqrt(eps) * ( matS + sqrt(eps)*max(diag(matS))*eye(sz,sz) );
	endif
	if ( doChecks )
		assert( isrealarray(matS,[sz,sz]) );
		assert( issymmetric(matS) );
		% Also, t*matH + (1.0-t)*matS must be positive definite or positive semi-definite.
		assert( isrealarray(matRegu,[sz,sz]) );
		assert( issymmetric(matRegu) );
		% Also, matS must be positive definite.
	endif
	%
	matM = t*matH + (1.0-t)*matS;
	if ( cholTol > 0.0 )
		% Try without regularization.
		[ matR, cholFlag ] = chol( matM );
		if ( 0 == cholFlag && min(diag(matR)) > cholTol*max(abs(diag(matR))) )
			if ( 1 == nargout )
				vecY = matR \ ( matR' \ (-t*vecG) );
			else
				[ vecY, vecYPrime, s, sPrime ] = __fromChol( vecG, matH, t, matS, matR );
			endif
			return;
		endif
	endif
	%
	[ matR, cholFlag ] = chol( matM + matRegu );
	if ( 0 ~= cholFlag )
		error( "Cholesky factorization failed even with regularization matrix." );
	endif
	if ( 1 == nargout )
		vecY1 = matR \ ( matR' \ (-t*vecG) );
	else
		[ vecY1, vecYPrime1, s1, sPrime1 ] = __fromChol( vecG, matH, t, matS, matR );
	endif
	%
	[ matR, cholFlag ] = chol( matM + 2.0*matRegu );
	if ( 0 ~= cholFlag )
		error( "Cholesky factorization failed with regularization matrix second time." );
	endif
	switch (nargout)
	case 1
		vecY2 = matR \ ( matR' \ (-t*vecG) );
		vecY = 2.0*vecY1 - vecY2;
	otherwise
		[ vecY2, vecYPrime2, s2, sPrime2 ] = __fromChol( vecG, matH, t, matS, matR );
		vecY = 2.0*vecY1 - vecY2;
		vecYPrime = 2.0*vecYPrime1 - vecYPrime2;
		s = 2.0*s1 - s2;
		sPrime = 2.0*sPrime1 - sPrime2;
	endswitch
	return;
endfunction


function [ vecY, vecYPrime, s, sPrime ] = __fromChol( vecG, matH, t, matS, matR )
	vecMING = matR \ ( matR' \ (-vecG) ); % Called "gamma" in notes 2022-05-15-2400.
	vecY = t * vecMING;
	vecSMING = matS * vecMING;
	vecYPrime = matR \ ( matR' \ vecSMING );
	s0 = sqrt(max([ 0.0, vecMING' * vecSMING ]));
	s = t * s0;
	sPrime = max([ 0.0, ( vecYPrime' * vecSMING ) / s0 ]);
	% Note: If we didn't bother calculating vecYPrime, we could get sPrime using:
	%sPrime_alt = sumsq( matR' \ ( matS * vecMING ) ) / s0;
	%assert( reldiff( sPrime, sPrime_alt ) < sqrt(eps) );
return;
endfunction


%!test
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 5;
%!	sizeF = 5;
%!	%
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX);
%!	if (0)
%!		vecPhi = randn(sizeX,1);
%!		vecPhi /= norm(vecPhi);
%!		matJ -= (matJ*vecPhi)*(vecPhi');
%!	endif
%!	matH = matJ'*matJ;
%!	vecG = matJ'*vecF;
%!	matS = diag(diag(matH));
%!	matE = sqrt(eps)*matS + (eps*max(diag(matH)))*eye(sizeX,sizeX);
%!	matS_sing = matS; matS_sing(1:2:end,1:2:end) = 0.0;
%!	%
%!	vecY = calcLevPt( vecG, matH, 0.0 );
%!	vecY = calcLevPt( vecG, matH, 1.0 );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 0.0 );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 1.0 );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 0.0, matS );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 1.0, matS );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 0.0, matS, matS );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 1.0, matS, matS );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 0.0, matS, matE, 0.0, true );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 1.0, matS, matE, 0.0, true );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 0.0, matS_sing, [], sqrt(eps), true );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 1.0, matS_sing, [], sqrt(eps), true );
%!	[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, 0.5, matS_sing, [], sqrt(eps), true );
%!	%
%!	hadError = false;
%!	try
%!		% If t*H + (1-t)*S is npd and E is singular, we should get an error.
%!		vecY = calcLevPt( vecG, matH, 0.0, matS_sing, matS_sing );
%!	catch
%!		hadError = true;
%!		msg( __FILE__, __LINE__, "calcLevPt correctly threw an error." );
%!	end
%!	if (~hadError)
%!		error( "calcLevPt failed to throw an error when it should have." );
%!	endif
%!	return;
%!	%
%!	numVals = 101;
%!	foo = linspace( 1.0, 0.0, numVals );
%!	tVals = ( 1.0 - (foo.^2) ).^2;
%!	%
%!	for n=1:numVals
%!		[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, tVals(n) );
%!		vecYVals(:,n) = vecY;
%!		vecYPrimeVals(:,n) = vecYPrime;
%!		sVals(n) = s;
%!		sPrimeVals(n) = sPrime;
%!	endfor
%!	n = 1+round( (numVals-1)*0.2 ); sModelVals1 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	n = 1+round( (numVals-1)*0.5 ); sModelVals2 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	n = 1+round( (numVals-1)*0.8 ); sModelVals3 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  tVals, sVals, 'o-', ...
%!	  tVals, cap( sModelVals1, 0.0, max(sVals) ), 'x-', ...
%!	  tVals, cap( sModelVals2, 0.0, max(sVals) ), '^-', ...
%!	  tVals, cap( sModelVals3, 0.0, max(sVals) ), 'v-' );
%!	grid on;
%!	numFigs++; figure(numFigs);
%!	plot( tVals, sPrimeVals, 'o-' );
%!	grid on;
%!	%
%!	for n=1:numVals
%!		[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, tVals(n), matS );
%!		vecYVals(:,n) = vecY;
%!		vecYPrimeVals(:,n) = vecYPrime;
%!		sVals(n) = s;
%!		sPrimeVals(n) = sPrime;
%!	endfor
%!	n = 1+round( (numVals-1)*0.2 ); sModelVals1 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	n = 1+round( (numVals-1)*0.5 ); sModelVals2 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	n = 1+round( (numVals-1)*0.8 ); sModelVals3 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  tVals, sVals, 'o-', ...
%!	  tVals, cap( sModelVals1, 0.0, max(sVals) ), 'x-', ...
%!	  tVals, cap( sModelVals2, 0.0, max(sVals) ), '^-', ...
%!	  tVals, cap( sModelVals3, 0.0, max(sVals) ), 'v-' );
%!	grid on;
%!	numFigs++; figure(numFigs);
%!	plot( tVals, sPrimeVals, 'o-' );
%!	grid on;
%!	%
%!	for n=1:numVals
%!		[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, tVals(n), matS, matE );
%!		vecYVals(:,n) = vecY;
%!		vecYPrimeVals(:,n) = vecYPrime;
%!		sVals(n) = s;
%!		sPrimeVals(n) = sPrime;
%!	endfor
%!	n = 1+round( (numVals-1)*0.2 ); sModelVals1 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	n = 1+round( (numVals-1)*0.5 ); sModelVals2 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	n = 1+round( (numVals-1)*0.8 ); sModelVals3 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  tVals, sVals, 'o-', ...
%!	  tVals, cap( sModelVals1, 0.0, max(sVals) ), 'x-', ...
%!	  tVals, cap( sModelVals2, 0.0, max(sVals) ), '^-', ...
%!	  tVals, cap( sModelVals3, 0.0, max(sVals) ), 'v-' );
%!	grid on;
%!	numFigs++; figure(numFigs);
%!	plot( tVals, sPrimeVals, 'o-' );
%!	grid on;
%!	%
%!	for n=1:numVals
%!		[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, tVals(n), matS_sing );
%!		vecYVals(:,n) = vecY;
%!		vecYPrimeVals(:,n) = vecYPrime;
%!		sVals(n) = s;
%!		sPrimeVals(n) = sPrime;
%!	endfor
%!	n = 1+round( (numVals-1)*0.2 ); sModelVals1 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	n = 1+round( (numVals-1)*0.5 ); sModelVals2 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	n = 1+round( (numVals-1)*0.8 ); sModelVals3 = sVals(n) + sPrimeVals(n)*( tVals - tVals(n) );
%!	numFigs++; figure(numFigs);
%!	plot( ...
%!	  tVals, sVals, 'o-', ...
%!	  tVals, cap( sModelVals1, 0.0, max(sVals) ), 'x-', ...
%!	  tVals, cap( sModelVals2, 0.0, max(sVals) ), '^-', ...
%!	  tVals, cap( sModelVals3, 0.0, max(sVals) ), 'v-' );
%!	grid on;
%!	numFigs++; figure(numFigs);
%!	plot( tVals, sPrimeVals, 'o-' );
%!	grid on;
