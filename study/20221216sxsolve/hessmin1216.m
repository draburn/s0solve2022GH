function [ vecDelta, datOut ] = hessmin1216( f0, vecG, matH, matB=[], bMax=[], prm=[] )
	datOut = [];
	sizeX = size(vecG,1);
	assert( 1 <= sizeX );
	assert( isrealscalar(f0) );
	assert( f0 >= 0.0 );
	assert( isrealarray(vecG,[sizeX,1]) );
	assert( isrealarray(matH,[sizeX,sizeX]) );
	assert( issymmetric(matH) );
	if ( isempty(matB) )
		matB = eye(sizeX,sizeX);
	endif
	assert( isrealarray(matB) );
	assert( issymmetric(matB) );
	if ( ~isdiag(matB) )
		msg( __FILE__, __LINE__, "WARNING: Boundary matrix is not diagonal. Performance may be very slow." );
	endif
	%
	if ( 0.0 == f0 )
		msg( __FILE__, __LINE__, "WARNING: Initial objective function value is zero." );
		vecX = zeros(sizeX,1);
		return;
	endif
	if ( 0.0 == norm(vecG) )
		msg( __FILE__, __LINE__, "WARNING: Initial gradient is zero." );
		vecX = zeros(sizeX,1);
		return;
	endif
	matRB = chol(matB);
	vecG = matRB \ ( matRB' \ vecG );
	matH = ( (matRB\(matRB'\matH)) / matRB ) / (matRB');
	matI = eye(sizeX,sizeX);
	%
	h = sqrt(sum(sum(matH.^2,1))); % Or not?
	if ( 0.0 == h )
		msg( __FILE__, __LINE__, "WARNING: Hessian matrix is zero." );
		h = 1.0;
	endif
	c0 = f0/h;
	vecB = (-1.0/h)*vecG;
	matA = (matH/h) - matI;
	%
	[ sR, vecDeltaR, cR, vecDeltaPrimeR, cPrimeR ] = getFarPt( c0, vecB, matA, prm )
	vecDelta = matRB \ ( matRB' \ vecDeltaR );
	return;
	%
	error( "END OF VALID CODE." );

return;
endfunction

% After scaling...
%  c0 = f0/h
%  vecB = -vecG/h
%  matA = matH/h - I
function [ vecDelta, c, vecDeltaPrime, cPrime ] = calcDelta( s, c0, vecB, matA )
	[ matR, cholFlag ] = chol( eye(size(matA)) + s * matA );
	if ( 0 ~= cholFlag )
		vecDelta = [];
		c = [];
		vecDeltaPrime = [];
		cPrime = [];
		return;
	endif
	vecMInvB = matR \ ( matR' \ vecB );
	vecDelta = s * vecMInvB;
	if ( 1 >= nargout )
		return;
	endif
	vecIPADelta = ( eye(size(matA)) + matA ) * vecDelta;
	c = c0 - (vecB'*vecDelta) + 0.5*(vecDelta'*vecIPADelta);
	if ( 2 >= nargout )
		return;
	endif
	%
	vecMInv2B = matR \ ( matR' \ vecMInvB );
	vecDeltaPrime = vecMInv2B;
	cPrime = -(vecB'*vecDeltaPrime) + (vecDeltaPrime'*vecIPADelta);
	%
	%delta = norm(vecDelta);
	%deltaPrime = (vecMInv2B'*vecMInvB)/norm(vecMInvB);
return;
endfunction;
function [ vecDelta, c, vecDeltaPrime, cPrime ] = calcDelta_withExtrap( s, c0, vecB, matA, epsS=sqrt(eps), rdTol=0.1 )
	[ vecDelta1, c1, vecDeltaPrime1, cPrime1 ] = calcDelta( s-epsS, c0, vecB, matA );
	if ( isempty(vecDelta1) )
		vecDelta = [];
		c = [];
		vecDeltaPrime = [];
		cPrime = [];
	endif
	[ vecDelta2, c2, vecDeltaPrime2, cPrime2 ] = calcDelta( s-2.0*epsS, c0, vecB, matA );
	if ( isempty(vecDelta2) )
		vecDelta = [];
		c = [];
		vecDeltaPrime = [];
		cPrime = [];
	endif
	vecDelta = 2.0*vecDelta1 - vecDelta2;
	vecDeltaPrime = 2.0*vecDeltaPrime1 - vecDeltaPrime;
	if ( reldiff(vecDelta,vecDelta1) > rdTol || reldiff(vecDeltaPrime,vecDeltaPrime1) > rdTol )
		vecDelta = [];
		c = [];
		vecDeltaPrime = [];
		cPrime = [];
	endif
	vecIPADelta = ( eye(size(matA)) + matA ) * vecDelta;
	c = c0 - (vecB'*vecDelta) + 0.5*(vecDelta'*vecIPADelta);
	cPrime = -(vecB'*vecDeltaP) + (vecDeltaP'*vecIPADelta);
return;
endfunction;

function [ s, vecDelta, c, vecDeltaPrime, cPrime ] = getFarPt( c0, vecB, matA, prm )
	[ vecDelta, c, vecDeltaPrime, cPrime ] = calcDelta( 1.0, c0, vecB, matA );
	if ( ~isempty(vecDelta) )
		s = 1.0;
		return;
	endif
	[ vecDelta, c, vecDeltaPrime, cPrime ] = calcDelta_withExtrap( 1.0, c0, vecB, matA )
	if ( ~isempty(vecDelta) )
		s = 1.0;
		return;
	endif
	eigvecA = eig(matA);
	assert( min(eigvecA) < -0.999 );
	s = min([ -1.0/min(eigvecA), 1.0-sqrt(eps) ]);
	[ vecDelta, c, vecDeltaPrime, cPrime ] = calcDelta( s, c0, vecB, matA );
	assert( ~isempty(vecDelta) );
return;
endfunction
