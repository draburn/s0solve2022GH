% Function...

function [ vecY, vecYPrime, s, sPrime, deltaOmega, deltaOmegaPrime ] = calcLevPt( ...
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
		matS = eye(sz,sz);
	endif
	if (isempty(matRegu) )
		matRegu = 0.1*sqrt(eps)*eye(sz,sz);
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
			switch (nargout)
			case 1
				vecY = __fromChol_1( vecG, matH, t, matS, matR );
			case 2
				[ vecY, vecYPrime ] = __fromChol_2( vecG, matH, t, matS, matR );
			case { 3, 4 }
				[ vecY, vecYPrime, s, sPrime ] = __fromChol_4( vecG, matH, t, matS, matR );
			otherwise
				[ vecY, vecYPrime, s, sPrime, deltaOmega, deltaOmegaPrime ] = __fromChol_6( vecG, matH, t, matS, matR );
			endswitch
			return;
		endif
	endif
	%
	[ matR, cholFlag ] = chol( matM + matRegu );
	if ( 0 ~= cholFlag )
		error( "Cholesky factorization failed even with regularization matrix." );
	endif
	switch (nargout)
	case 1
		vecY1 = __fromChol_1( vecG, matH, t, matS, matR );
	case 2
		[ vecY1, vecYPrime1 ] = __fromChol_2( vecG, matH, t, matS, matR );
	case { 3, 4 }
		[ vecY1, vecYPrime1, s1, sPrime1 ] = __fromChol_4( vecG, matH, t, matS, matR );
	otherwise
		[ vecY1, vecYPrime1, s1, sPrime1, deltaOmega1, deltaOmegaPrime1 ] = __fromChol_6( vecG, matH, t, matS, matR );
	endswitch
	%
	[ matR, cholFlag ] = chol( matM + 2.0*matRegu );
	if ( 0 ~= cholFlag )
		error( "Cholesky factorization failed with regularization matrix second time." );
	endif
	switch (nargout)
	case 1
		vecY2 = __fromChol_1( vecG, matH, t, matS, matR );
		vecY = vecY1 - (2.0*vecY2);
	case 2
		[ vecY2, vecYPrime2 ] = __fromChol_2( vecG, matH, t, matS, matR );
		vecY = vecY1 - (2.0*vecY2);
		vecYPrime = vecYPrime1 - (2.0*vecYPrime2);
	case { 3, 4 }
		[ vecY2, vecYPrime2, s2, sPrime2 ] = __fromChol_4( vecG, matH, t, matS, matR );
		vecY = vecY1 - (2.0*vecY2);
		vecYPrime = vecYPrime1 - (2.0*vecYPrime2);
		s = s1 - (2.0*s2);
		sPrime = sPrime1 - (2.0*sPrime2);
	otherwise
		[ vecY2, vecYPrime2, s2, sPrime2, deltaOmega2, deltaOmegaPrime2 ] = __fromChol_6( vecG, matH, t, matS, matR );
		vecY = vecY1 - (2.0*vecY2);
		vecYPrime = vecYPrime1 - (2.0*vecYPrime2);
		s = s1 - (2.0*s2);
		sPrime = sPrime1 - (2.0*sPrime2);
		deltaOmega = deltaOmega1 - (2.0*deltaOmega2);
		deltaOmegaPrime = deltaOmegaPrime1 - (2.0*deltaOmegaPrime2);
	endswitch
	return;
endfunction


function vecY = __fromChol_1( vecG, matH, t, matS, matR )
	vecMING = matR \ ( matR' \ (-vecG) );
	vecY = t*vecMING;
return;
endfunction


function [ vecY, vecYPrime ] = __fromChol_2( vecG, matH, t, matS, matR )
	vecMING = matR \ ( matR' \ (-vecG) );
	vecY = t * vecMING;
	vecSMING = matS * vecMING;
	vecYPrime = matR \ ( matR' \ vecSMING );
return;
endfunction


function [ vecY, vecYPrime, s, sPrime ] = __fromChol_4( vecG, matH, t, matS, matR )
	vecMING = matR \ ( matR' \ (-vecG) );
	vecY = t * vecMING;
	vecSMING = matS * vecMING;
	vecYPrime = matR \ ( matR' \ vecSMING );
	s = todo
return;
endfunction


%!test
%!	setprngstates(0);
%!	sizeX = 5;
%!	sizeF = 5;
%!	%
%!	vecF = randn(sizeF,1);
%!	matJ = randn(sizeF,sizeX);
%!	if (1)
%!		vecPhi = randn(sizeX,1);
%!		vecPhi /= norm(vecPhi);
%!		matJ -= (matJ*vecPhi)*(vecPhi');
%!	endif
%!	matH = matJ'*matJ;
%!	vecG = matJ'*vecF;
%!	%
%!	[ vecY, vecYPrime ] = calcLevPt( vecG, matH, 0.5 )
