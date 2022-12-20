function [ vecDelta, datOut ] = hessmin( f0, vecG, matH, matB=[], bMax=[], prm=[] )
	error( "This code is abandoned, though it may still work." );
	% Do initialization.
	datOut = [];
	sizeX = size(vecG,1);
	assert( 1 <= sizeX );
	matI = eye(sizeX,sizeX);
	assert( isrealscalar(f0) );
	fMinAllowed = mygetfield( prm, "fMinAllowed", -1.0*f0 );
	if ( ~isempty(fMinAllowed) )
		assert( isrealscalar(fMinAllowed) );
		assert( fMinAllowed < f0 );
	endif
	assert( isrealarray(vecG,[sizeX,1]) );
	if ( 0.0 == norm(vecG) )
		msg( __FILE__, __LINE__, "WARNING: Initial gradient is zero." );
		vecDelta = zeros(sizeX,1);
		return;
	endif
	assert( isrealarray(matH,[sizeX,sizeX]) );
	assert( issymmetric(matH) );
	if ( isempty(matB) )
		matB = matI;
	endif
	assert( isrealarray(matB) );
	assert( issymmetric(matB) );
	if ( ~isdiag(matB) )
		msg( __FILE__, __LINE__, "WARNING: Boundary matrix is not diagonal. Performance may be very slow." );
	endif
	if ( ~isempty(bMax) )
		assert( isrealscalar(bMax) );
		assert( 0.0 <= bMax );
		if ( bMax == 0.0 )
			msg( __FILE__, __LINE__, "WARNING: A zero step was requested." );
			vecDelta = zeros(sizX,1);
			return;
		endif
	endif
	%
	% Apply scaling.
	matRB = mygetfield( prm, "matRB", [] );
	if ( isempty(matRB) )
		matRB = chol(matB);
	else
		rd = reldiff( matB, matRB'*matRB );
		if ( rd > sqrt(eps) )
			msg( __FILE__, __LINE__, "WARNING: Provided Cholesky factorization of barrier matrix may be incorrect." );
		endif
	
	endif
	vecG = matRB \ ( matRB' \ vecG );
	matH = ( (matRB\(matRB'\matH)) / matRB ) / (matRB');
	%
	% Do main work.
	[ vecDelta, datOut ] = __solve( f0, vecG, matH, bMax, fMinAllowed, prm );
	%
	% Properly scale result.
	vecDelta = matRB \ ( matRB' \ vecDelta );
return;
endfunction


function [ vecDelta, f, vecDeltaPrime, fPrime, matR ] = __getDelta( s, f0, vecG, matH, h, cholTol=sqrt(eps) )
	if ( 1.0 == s )
		[ matR, cholFlag ] = chol( matH );
	else
		[ matR, cholFlag ] = chol( s*matH + ((h*(s-1.0))*eye(size(matH))) );
	endif
	if ( 0 ~= cholFlag || min(diag(matR)) < cholTol * max(abs(diag(matR))) )
		vecDelta = [];
		f = [];
		vecDeltaPrime = [];
		fPrime = [];
		return;
	endif
	error( "Check me." );
	vecMInvG = matR \ ( matR' \ vecG );
	vecDelta = (-s)*matMinvG;
	vecHDelta = matH*vecDelta;
	f = f0 + vecDelta'*vecG + (vecDelta'*vecHDelta)/2.0;
	vecMInv2G = matR \ ( matR' \ vecMInvG ):
	vecDeltaPrime = vecMInv2G;
	fPrime = vecDeltaPrime'*( vecG + vecHDelta );
return;
endfunction


function [ vecDelta, datOut ] = __solve( f0, vecG, matH, bMax, fMinAllowed, prm )
	datOut = [];
	hScl = sqrt(sum(sum(matH.^2)));
	error( "END OF VALID CODE." );
	if ( 0.0 == hScl )
		msg( __FILE__, __LINE__, "WARNING: Hessian matrix is zero." );
		if ( isempty(bMax) && isempty(fMinAllowed) )
			error( "Hessian is zero, gradient is not, and no constraints have been specified. There is no solution." );
		elseif ( isempty(bMax) )
			error( "DO WE WANT TO GO TO F = 0 IN THIS CASE???" );
			% So, fMinAllowed functions more like a tolerance?
			% bTol?
		elseif ( isempty(fMinAllowed) )
		else
		endif
	endif
	
	%
	
	%
	cholTol = mygetfield( prm, "cholTol", sqrt(eps) );
	assert( isrealscalar(cholTol) );
	
	%
	%if ( ~isempty(f0) )
	%	fMinAllowed = mygetfield( prm, "fMinAllowed", -0.1*f0 );
	%	if ( ~isempty(fMinAllowed) )
	%		assert( isrealscalar(fMinAllowed) )
	%		assert( fMinAllowed < f0 );
	%	endif
	%endif
	%
	hScl = sqrt(sum(sum(matH.^2)));
	if ( 0.0 == hScl )
		if ( isempty(bMax) )
		endif
	endif
	%
	s1 = 1.0;
	[ vecDelta1, f1, vecDeltaPrime1, fPrime1, matR1 ] = __getDelta( s1, f0, vecG, matH, hScl, prm.cholTol );
	if ( isempty(vecDelta1) )
		vecLambda = eig(matH);
		muCrit = max([ 0.0, -vecLambda ]);
		muSafe = muCrit + sqrt(eps)*max(abs(vecLambda));
		s1 = 1.0 / ( h*muSafe + 1.0 );
		[ vecDelta1, f1, vecDeltaPrime1, fPrime1, matR1 ] = __getDelta( s1, f0, vecG, matH, hScl, prm.cholTol );
	else
		muCrit = 0.0;
		muSafe = 0.0;
	endif
	%
	% Even if bMax is empty, we still may need to BT for fMinAllowed.
	% Note: If bMax is empty, we might still consider extrapolting to s1 = 1.0;
	error( "What next?" );
return;
endfunction

% After scaling...
%  c0 = f0/h
%  vecB = -vecG/h
%  matA = matH/h - I
function [ vecDelta, c, vecDeltaPrime, cPrime ] = __calcDelta( s, c0, vecB, matA, prm )
	[ matR, cholFlag ] = chol( eye(size(matA)) + s * matA );
	if ( 0 ~= cholFlag )
		vecDelta = [];
		c = [];
		vecDeltaPrime = [];
		cPrime = [];
		return;
	elseif ( min(diag(matR)) < prm.cholTol*max(abs(diag(matR))) )
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
function [ vecDelta, c, vecDeltaPrime, cPrime ] = __calcDelta_withExtrap( s, c0, vecB, matA, prm )
	[ vecDelta1, c1, vecDeltaPrime1, cPrime1 ] = __calcDelta( s-prm.epsS, c0, vecB, matA, prm );
	if ( isempty(vecDelta1) )
		vecDelta = [];
		c = [];
		vecDeltaPrime = [];
		cPrime = [];
		return;
	endif
	[ vecDelta2, c2, vecDeltaPrime2, cPrime2 ] = __calcDelta( s-2.0*prm.epsS, c0, vecB, matA, prm );
	if ( isempty(vecDelta2) )
		vecDelta = [];
		c = [];
		vecDeltaPrime = [];
		cPrime = [];
		return;
	endif
	vecDelta = 2.0*vecDelta1 - vecDelta2;
	vecDeltaPrime = 2.0*vecDeltaPrime1 - vecDeltaPrime2;
	[ vecDelta, vecDelta1, vecDelta2 ]
	[ vecDeltaPrime, vecDeltaPrime1, vecDeltaPrime2 ]
	reldiff(vecDelta,vecDelta1)
	reldiff(vecDeltaPrime,vecDeltaPrime1)
	if ( reldiff(vecDelta,vecDelta1) > prm.extrap_rdTol )
		vecDelta = [];
		c = [];
		vecDeltaPrime = [];
		cPrime = [];
		return;
	endif
	vecIPADelta = ( eye(size(matA)) + matA ) * vecDelta;
	c = c0 - (vecB'*vecDelta) + 0.5*(vecDelta'*vecIPADelta);
	cPrime = -(vecB'*vecDeltaPrime) + (vecDeltaPrime'*vecIPADelta);
return;
endfunction;

function [ s, vecDelta, c, vecDeltaPrime, cPrime ] = __getFarPt( c0, vecB, matA, prm )
	[ vecDelta, c, vecDeltaPrime, cPrime ] = __calcDelta( 1.0, c0, vecB, matA, prm );
	if ( ~isempty(vecDelta) )
		s = 1.0;
		return;
	endif
	[ vecDelta, c, vecDeltaPrime, cPrime ] = __calcDelta_withExtrap( 1.0, c0, vecB, matA, prm );
	if ( ~isempty(vecDelta) )
		s = 1.0;
		return;
	endif
	eigvecA = eig(matA)
	assert( min(eigvecA) < -0.999 );
	s = min([ -1.0/min(eigvecA), 1.0 ]) * ( 1.0 - prm.epsS );
	[ vecDelta, c, vecDeltaPrime, cPrime ] = __calcDelta( s, c0, vecB, matA, prm );
	assert( ~isempty(vecDelta) );
return;
endfunction

function [ res, c, resPrime, cPrime ] = __funcFC_forMyFZero( s, c0, vecB, matA, bMax, prm )
	matR = chol( eye(size(matA)) + s * matA );
	vecMInvB = matR \ ( matR' \ vecB );
	vecDelta = s * vecMInvB;
	res = norm(vecDelta) - bMax;
	vecIPADelta = ( eye(size(matA)) + matA ) * vecDelta;
	c = c0 - (vecB'*vecDelta) + 0.5*(vecDelta'*vecIPADelta);
	%
	vecMInv2B = matR \ ( matR' \ vecMInvB );
	vecDeltaPrime = vecMInv2B;
	resPrime = (vecMInv2B'*vecMInvB)/norm(vecMInvB);
	cPrime = -(vecB'*vecDeltaPrime) + (vecDeltaPrime'*vecIPADelta);
return;
endfunction;
