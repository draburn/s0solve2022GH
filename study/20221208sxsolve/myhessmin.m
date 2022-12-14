function [ vecX, datOut ] = myhessmin( f, vecG, matH, matD=[], dMax=[], prm=[] )
	datOut = [];
	sz = size(vecG,1);
	if ( 1 )
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
			msg( __FILE__, __LINE__, "Accepted pos-def form." );
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
				msg( __FILE__, __LINE__, "Rejected extrapolation form." );
			else
				msg( __FILE__, __LINE__, "Accepted extrapolation form." );
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
		msg( __FILE__, __LINE__, "Accepted perturbation form." );
		gthg = vecG'*matH*vecG;
		clear matH;
	endif
	%
	assert( gthg > 0.0 );
	vecXCauchy = vecG * (-gSq/gthg);
	vecDXLeg = vecXNewton - vecXCauchy;
	assert( 0.0 < norm(vecXCauchy) );
	assert( norm(vecXCauchy)*(1.0-sqrt(eps)) <= norm(vecXNewton)*(1.0+sqrt(eps)) );
	assert( (vecDXLeg'*vecXCauchy) >= -sqrt(eps)*( sumsq(vecDXLeg) + sumsq(vecXCauchy) ) );
	%
	if ( isempty(xMax) )
		vecX = vecXNewton;
		return;
	elseif ( xMax < norm(vecXCauchy) )
		vecX = vecXCauchy * xMax / norm(vecXCauchy);
		return;
	elseif ( xMax > norm(vecXNewton) )
		vecX = vecXNewton;
		return;
	endif
	%
	msg( __FILE__, __LINE__, "Using second leg." );
	% a*s^2 + b*s + c = 0
	a = sumsq( vecDXLeg );
	b = 2.0 * ( vecXCauchy'*vecDXLeg);
	c = sumsq( vecXCauchy ) - xMax^2;
	discrim = (b^2) - (4.0*a*c);
	assert( discrim >= 0.0 );
	assert( a > 0.0 );
	assert( c <= 0.0 );
	s = (-b+sqrt(discrim))/(2.0*a);
	vecX = vecXCauchy + s*vecDXLeg;
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

