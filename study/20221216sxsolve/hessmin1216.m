function [ vecDelta, datOut ] = hessmin1216( f0, vecG, matH, matB=[], bMax=[], prm=[] )
	msg( __FILE__, __LINE__, "ERROR: THIS CODE HAS BEEN ABANDONED AND PROBABLY DOES NOT WORK." );
	msg( __FILE__, __LINE__, "  12-16-2330: We should simply modify matH first;" );
	msg( __FILE__, __LINE__, "  This works both to ensure matH is pos-def and to satisfy the fMinAllowed constraint." );
	error( "This code has been abandoned." )
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
	[ vecDelta, datOut ] = __solve( c0, vecB, matA, bMax, prm );
	vecDelta = matRB \ ( matRB' \ vecDelta );
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

function [ vecDelta, datOut ] = __solve( c0, vecB, matA, bMax, prm )
	datOut = [];
	%
	prm.cholTol = mygetfield( prm, "cholTol", sqrt(eps) );
	prm.epsS = mygetfield( prm, "epsS", sqrt(eps) );
	prm.extrap_rdTol = mygetfield( prm, "extrap_rdTol", 0.1 );
	prm.cTolCoeff = mygetfield( prm, "cTolCoeff", 0.1 );
	prm.bTol = mygetfield( prm, "bTol", 1.0e-4 );
	%
	[ sR, vecDeltaR, cR, vecDeltaPrimeR, cPrimeR ] = __getFarPt( c0, vecB, matA, prm );
	bTol = prm.bTol;
	if ( isempty(prm.cTolCoeff) )
		cTol = [];
	else
		cTol = c0 * prm.cTolCoeff;
	endif
	if ( isempty(bMax) || norm(vecDeltaR) <= bMax + bTol )
	if ( isempty(cTol) || cR >= -cTol )
		msg( __FILE__, __LINE__, "SUCCESS: Initial far point is satisfactory." );
		vecDelta = vecDeltaR;
		return;
	endif
	endif
	%
	sL = 0.0;
	if ( isempty(bMax) )
		% myfzero is going to require some value of bMax.
		% We must be here because the constraint isn't satisfied,
		%  likely because the Hessian has a negative eigenvalue.
		assert( ~isempty(cTol) );
		assert( cR < -cTol );
		bMax = -1.0;
	endif
	mfzPrm = [];
	mfzPrm.fTol = bTol;
	mfzPrm.cTol = cTol;
	%mfzPrm.fR = norm(vecDeltaR) - bMax;
	%mfzPrm.cR = cR;
	%mfzPrm.fpR = deltaPrimeR OR zero.
	%
	funchFC = @(s)( __funcFC_forMyFZero( s, c0, vecB, matA, bMax, prm ) );
	[ s, mfzDat ] = myfzero1216( funchFC, [sL,sR], mfzPrm )
	vecDelta = __calcDelta( s, c0, vecB, matA, prm )
return;
endfunction
