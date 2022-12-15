function [ vecX, datOut ] = myhessmin_gperp( f, vecG, matH, matD=[], dMax=[], prm=[] )
	datOut = [];
	sz = size(vecG,1);
	if ( 1 )
		assert( 2 <= sz );
		assert( isrealscalar(f) );
		assert( 0.0 < f );
		assert( isrealarray(vecG,[sz,1]) );
		assert( isrealarray(matH,[sz,sz]) );
		assert( issymmetric(matH,1.0e-4) );
		if ( ~isempty(matD) )
			assert( isdiag(matD) );
			assert( 0.0 < min(diag(matD)) );
		endif
		if ( ~isempty(dMax) )
			assert( isrealscalar(dMax) );
			assert( 0.0 < dMax );
		endif
	endif
	if ( isempty(matD) )
		[ vecX, datOut ] = __solve( f, vecG, matH, dMax, prm );
	else
		matDInv = diag(1.0./diag(matD));
		[ vecXScaled, datOut ] = __solve( f, matDInv*vecG, matDInv*matH*matDInv, dMax, prm );
		vecX = matDInv * vecXScaled;
	endif
return;
endfunction


function [ vecX, datOut ] = __solve( f, vecG, matH, xMax, prm )
	datOut = [];
	sz = size(vecG,1);
	debugMode = mygetfield( prm, "debugMode", false );
	%
	fMin = mygetfield( prm, "fMin", [] );
	if ( ~isempty(fMin) )
		error( "ERROR: fMin is not (yet) supported." );
		fMinTrigger = mygetfield( prm, "fMinTrigger", -0.1*f );
		if ( 1 )
			assert( isrealscalar(fMin) );
			assert( fMin < f );
			assert( isrealscalar(fMinTrigger) );
			assert( fMin <= fMinTrigger );
			assert( fMinTrigger <= f );
		endif
	endif
	%
	gSq = sumsq(vecG);
	if ( 0.0 == gSq )
		msg( __FILE__, __LINE__, "WARNING: vecG is zero." );
		vecX = zeros(sz,1);
		return;
	endif
	%
	hScl = max(max(abs(matH)));
	if ( 0.0 == hScl )
		msg( __FILE__, __LINE__, "WARNING: matH is zero." );
		if ( ~isempty(xMax) )
			vecX = vecG * ( -xMax / g );
		else
			error( "ERROR: matH is zero and neither xMax nor fMin is specified." );
		endif
	endif
	%
	%
	% We'll start by making sure the Hessian is pos-semi-def and calculating the Newton step.
	hDiagMin = min(diag(matH));
	gthg = vecG'*matH*vecG;
	vecXNewton = [];
	if ( hDiagMin > 0.0  &&  gthg > 0.0 )
		% Attempt pos-def handling.
		cholTol = mygetfield( prm, "cholTol", 1.0e-6 );
		[ matR, cholFlag ] = chol( matH );
		if ( 0 == cholFlag )
		if ( min(diag(matR)) > cholTol*sqrt(hScl) )
			vecXNewton = matR \ ( matR' \ (-vecG) );
			msgif( debugMode, __FILE__, __LINE__, "Accepted pos-def form." );
		endif
		endif
	endif
	%
	if ( isempty(vecXNewton) && hDiagMin >= 0.0 && gthg >= 0.0 )
		% Try extrapolation.
		epsExtrap = mygetfield( prm, "epsExtrap", 1.0e-10 );
		[ matR1, cholFlag ] = chol( matH + epsExtrap*hScl*eye(sz,sz) );
		if ( 0 == cholFlag )
			matR2 = chol( matH + 2.0*epsExtrap*hScl*eye(sz,sz) );
			matR3 = chol( matH + 3.0*epsExtrap*hScl*eye(sz,sz) );
			vecXN1 = matR1 \ ( matR1' \ (-vecG) );
			vecXN2 = matR2 \ ( matR2' \ (-vecG) );
			vecXN3 = matR3 \ ( matR3' \ (-vecG) );
			vecXNewton = (3.0*(vecXN1-vecXN2)) + vecXN3;
			vecXNAlt = 2.0*vecXN1-vecXN2;
			if ( reldiff(vecXNAlt,vecXNewton) > 0.1 )
				vecXNewton = [];
				msgif( debugMode, __FILE__, __LINE__, "Rejected extrapolation form." );
			else
				msgif( debugMode, __FILE__, __LINE__, "Accepted extrapolation form." );
			endif
		endif
	endif
	%
	if ( isempty(vecXNewton) )
		% No choice but to perturb Hessian.
		epsPerturb = mygetfield( prm, "epsPerturb", 1.0e-8 );
		vecPosDefDiagMin = sum(abs(matH),2) - diag(abs(matH)) + epsPerturb*sz*hScl; % Scalar autobroadcast.
		vecHDiag = diag(matH);
		vecHDiagMod = vecHDiag;
		vecHDiagMod( vecHDiag < vecPosDefDiagMin ) = vecPosDefDiagMin( vecHDiag < vecPosDefDiagMin );
		matH += diag( vecHDiagMod - vecHDiag );
		matR = chol( matH );
		vecXNewton = matR \ ( matR' \ (-vecG) );
		msgif( debugMode,  __FILE__, __LINE__, "Accepted perturbation form." );
		gthg = vecG'*matH*vecG;
	endif
	%
	assert( gthg > 0.0 );
	vecXCauchy = vecG * (-gSq/gthg);
	%
	% Now we get to the "_gperp" aspect:
	%  the last component is considered "special" compared to the rest of the components.
	% In particular, the last component may be the local gradient direction perpendicular to the rest of the space,
	%  the corresponding elements of the Hessian having a larger uncertainty and the last diagonal element
	%  being just a guess.
	% We need to make sure things behave well whether that guess is too large or too small.
	% We consider two "steepest descent" directions:
	%  one in the actual (scaled) gradient descent direction, which includes the perpendicular component,
	%  and one which is fully in the rest of the space.
	% I'll use "wayPt" to refer to the min in the basis { perpendicular gradient direction, steepest-descent in rest of space }.
	matV = [ vecG(1:sz-1), zeros(sz-1,1); 0.0, vecG(sz) ];
	assert( norm(matV(:,1)) > 0.0 );
	matV(:,1) /= norm(matV(:,1));
	assert( norm(matV(:,2)) > 0.0 );
	matV(:,2) /= norm(matV(:,2));
	% Vectors are already clearly orthogonal.
	minivecG = matV'*vecG;
	minimatH = matV'*matH*matV;
	minimatR = chol(minimatH); % Validate pos-def.
	minivecXN = minimatR \ ( minimatR' \ (-minivecG) );
	vecXWayPt = matV * minivecXN;
	% And, we have 3 segments: 0 -> Cauchy -> wayPt -> Newton.
	% I'm not totally sure they are well ordered, but, let's check...
	if ( 0.0 >= norm(vecXCauchy) )
		msg( __FILE__, __LINE__, "WARNING: 0.0 >= norm(vecXCauchy)." );
	endif
	if ( norm(vecXCauchy)*(1.0+sqrt(eps)) >= norm(vecXWayPt)*(1.0-sqrt(eps)) )
		msgif( debugMode, __FILE__, __LINE__, "WARNING: norm(vecXCauchy) >= norm(vecXWayPt)." );
	endif
	if ( norm(vecXWayPt)*(1.0+sqrt(eps)) >= norm(vecXNewton)*(1.0-sqrt(eps)) )
		msgif( debugMode, __FILE__, __LINE__, "WARNING: norm(vecXWayPt) >= norm(vecXNewton)." );
	endif
	assert( norm(vecXCauchy)*(1.0-sqrt(eps)) <= norm(vecXNewton)*(1.0+sqrt(eps)) ); % Require this one.
	%
	if ( isempty(xMax) )
		msgif( debugMode, __FILE__, __LINE__, "Using full Newton step because xMax is empty." );
		vecX = vecXNewton;
	elseif ( xMax > norm(vecXNewton) )
		msgif( debugMode, __FILE__, __LINE__, "Using full Newton step because it is closer than xMax." );
		vecX = vecXNewton;
	elseif ( xMax < norm(vecXCauchy) )
		msgif( debugMode, __FILE__, __LINE__, "Using first leg (0 -> Cauchy)." );
		vecX = vecXCauchy * xMax / norm(vecXCauchy);
	elseif (  norm(vecXCauchy) < norm(vecXWayPt)  &&  norm(vecXWayPt) < norm(vecXNewton)  )
		% This is my assumption for the typical case.
		if ( xMax < norm(vecXWayPt) )
			msgif( debugMode, __FILE__, __LINE__, "Using second leg (Cauchy -> way-point) ." );
			vecXStart = vecXCauchy;
			vecDXLeg = vecXWayPt - vecXCauchy;
			% a*s^2 + b*s + c = 0
			a = sumsq( vecDXLeg );
			b = 2.0 * ( vecXStart' * vecDXLeg );
			c = sumsq( vecXStart ) - xMax^2;
			discrim = (b^2) - (4.0*a*c);
			assert( discrim >= 0.0 );
			assert( a > 0.0 );
			assert( c <= 0.0 );
			s = (-b+sqrt(discrim))/(2.0*a);
			vecX = vecXStart + s*vecDXLeg;
		else
			msgif( debugMode, __FILE__, __LINE__, "Using third leg (way-point -> Newton)." );
			vecXStart = vecXWayPt;
			vecDXLeg = vecXNewton - vecXWayPt;
			% a*s^2 + b*s + c = 0
			a = sumsq( vecDXLeg );
			b = 2.0 * ( vecXStart' * vecDXLeg );
			c = sumsq( vecXStart ) - xMax^2;
			discrim = (b^2) - (4.0*a*c);
			assert( discrim >= 0.0 );
			assert( a > 0.0 );
			assert( c <= 0.0 );
			s = (-b+sqrt(discrim))/(2.0*a);
			vecX = vecXStart + s*vecDXLeg;
		endif
	else
		% The expected ordering does not hold.
		% I'm not sure if this scenario should be possible.
		% Ignore the way pt.
		msgif( debugMode, __FILE__, __LINE__, "Using original leg (Cauchy -> Newton)." );
		vecXStart = vecXCauchy;
		vecDXLeg = vecXNewton - vecXCauchy;
		% a*s^2 + b*s + c = 0
		a = sumsq( vecDXLeg );
		b = 2.0 * ( vecXStart' * vecDXLeg );
		c = sumsq( vecXStart ) - xMax^2;
		discrim = (b^2) - (4.0*a*c);
		assert( discrim >= 0.0 );
		assert( a > 0.0 );
		assert( c <= 0.0 );
		s = (-b+sqrt(discrim))/(2.0*a);
		vecX = vecXStart + s*vecDXLeg;
	endif
return;
endfunction

%!test
%!	setprngstates(0);
%!	sizeX = 2 + ceil(20*rand());
%!	size1 = sizeX-1;
%!	size2 = 1;
%!	fMin = 0.0;
%!	vecXMin = randn(sizeX,1);
%!	foo = randn(size1,sizeX); matH1 = foo'*foo; clear foo;
%!	foo = randn(size2,sizeX); matH2 = foo'*foo; clear foo;
%!	%
%!	vecX = zeros(sizeX,1);
%!	matD = eye(sizeX,sizeX);
%!	%
%!	%
%!	%
%!	matH = matH1 + matH2;
%!	vecG = matH * ( vecX - vecXMin );
%!	fMin = 0.0;
%!	f = fMin + ((vecX-vecXMin)'*matH*(vecX-vecXMin))/2.0
%!	%
%!	%
%!	dMax = 1.0
%!	vecDelta = myhessmin( f, vecG, matH, matD, dMax );
%!	d = norm(matD*vecDelta)
%!	assert( d*(1.0-sqrt(eps)) <= dMax*(1.0+sqrt(eps)) );
%!	vecXNew = vecX + vecDelta;
%!	fNew = fMin + ((vecXNew-vecXMin)'*matH*(vecXNew-vecXMin))/2.0
%!	assert( fNew < f );
%!	%
%!	dMax = 2.0
%!	vecDelta = myhessmin( f, vecG, matH, matD, dMax );
%!	d = norm(matD*vecDelta)
%!	assert( d*(1.0-sqrt(eps)) <= dMax*(1.0+sqrt(eps)) );
%!	vecXNew = vecX + vecDelta;
%!	fNew = fMin + ((vecXNew-vecXMin)'*matH*(vecXNew-vecXMin))/2.0
%!	assert( fNew < f );
%!	%
%!	dMax = []
%!	vecDelta = myhessmin( f, vecG, matH, matD, dMax );
%!	d = norm(matD*vecDelta)
%!	vecXNew = vecX + vecDelta;
%!	fNew = fMin + ((vecXNew-vecXMin)'*matH*(vecXNew-vecXMin))/2.0
%!	assert( fNew < f );
%!	%
%!	%
%!	%
%!	matH = matH1;
%!	vecG = matH * (vecX-vecXMin);
%!	fMin = 0.0;
%!	f = fMin + ((vecX-vecXMin)'*matH*(vecX-vecXMin))/2.0
%!	%
%!	%
%!	dMax = 0.1
%!	vecDelta = myhessmin( f, vecG, matH, matD, dMax );
%!	d = norm(matD*vecDelta)
%!	assert( d*(1.0-sqrt(eps)) <= dMax*(1.0+sqrt(eps)) );
%!	vecXNew = vecX + vecDelta;
%!	fNew = fMin + ((vecXNew-vecXMin)'*matH*(vecXNew-vecXMin))/2.0
%!	assert( fNew < f );
%!	%
%!	dMax = 2.0
%!	vecDelta = myhessmin( f, vecG, matH, matD, dMax );
%!	d = norm(matD*vecDelta)
%!	assert( d*(1.0-sqrt(eps)) <= dMax*(1.0+sqrt(eps)) );
%!	vecXNew = vecX + vecDelta;
%!	fNew = fMin + ((vecXNew-vecXMin)'*matH*(vecXNew-vecXMin))/2.0
%!	assert( fNew < f );
%!	%
%!	dMax = []
%!	vecDelta = myhessmin( f, vecG, matH, matD, dMax );
%!	d = norm(matD*vecDelta)
%!	vecXNew = vecX + vecDelta;
%!	fNew = fMin + ((vecXNew-vecXMin)'*matH*(vecXNew-vecXMin))/2.0
%!	assert( fNew < f );
%!	%
%!	%
%!	%
%!	matH = matH1 - matH2;
%!	vecG = matH * ( vecX - vecXMin );
%!	fMin = 100.0;
%!	f = fMin + ((vecX-vecXMin)'*matH*(vecX-vecXMin))/2.0
%!	%
%!	%
%!	dMax = 0.01
%!	vecDelta = myhessmin( f, vecG, matH, matD, dMax );
%!	d = norm(matD*vecDelta)
%!	assert( d*(1.0-sqrt(eps)) <= dMax*(1.0+sqrt(eps)) );
%!	vecXNew = vecX + vecDelta;
%!	fNew = fMin + ((vecXNew-vecXMin)'*matH*(vecXNew-vecXMin))/2.0
%!	assert( fNew < f );
%!	%
%!	dMax = 0.07
%!	vecDelta = myhessmin( f, vecG, matH, matD, dMax );
%!	d = norm(matD*vecDelta)
%!	assert( d*(1.0-sqrt(eps)) <= dMax*(1.0+sqrt(eps)) );
%!	vecXNew = vecX + vecDelta;
%!	fNew = fMin + ((vecXNew-vecXMin)'*matH*(vecXNew-vecXMin))/2.0
%!	assert( fNew < f );
%!	%
%!	dMax = []
%!	vecDelta = myhessmin( f, vecG, matH, matD, dMax );
%!	d = norm(matD*vecDelta)
%!	vecXNew = vecX + vecDelta;
%!	fNew = fMin + ((vecXNew-vecXMin)'*matH*(vecXNew-vecXMin))/2.0
%!	assert( fNew < f );

