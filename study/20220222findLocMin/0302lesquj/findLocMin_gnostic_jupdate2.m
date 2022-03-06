% Function...
%  A 'gnostic for studying how various forms of Jacobian updating impact convergence.

function [ vecX, datOut ] = findLocMin_gnostic_jupdate2( vecX0, funchF, prm=[] )
	commondefs;
	findLocMin_gnostic_jupdate2__defs;
	time0 = time();
	fevalCount = 0;
	jevalCount = 0;
	collected_vecXVals = [];
	collected_vecFVals = [];
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__UNLIMITED );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(verbLev) );
		assert( isrealscalar(valdLev) );
	endif
	%
	sizeX = size(vecX0,1);
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isrealarray(vecX0,[sizeX,1]) );
	endif
	%
	vecF0 = funchF( vecX0 );
	fevalCount++;
	sizeF = size(vecF0,1);
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isrealarray(vecF0,[sizeF,1]) );
	endif
	collected_vecXVals = [ collected_vecXVals, vecX0 ];
	collected_vecFVals = [ collected_vecFVals, vecF0 ];
	collected_indexCurrent = size(collected_vecXVals,2);
	%
	matJ0 = jacobs( vecX0, funchF ); jevalCount++;
	if ( valdLev >= VALDLEV__MEDIUM )
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
	endif
	%
	omega0 = sumsq(vecF0)/2.0;
	vecG0 = matJ0'*vecF0;
	gNorm0 = norm(vecG0);
	matH0 = matJ0'*matJ0;
	hNorm0 = max(abs(diag(matH0)));
	assert( 0.0 < omega0 );
	assert( 0.0 < gNorm0 );
	fNorm0 = norm(vecF0);
	jNorm0 = sqrt(hNorm0);
	assert( 0.0 < jNorm0 );
	%
	%
	%
	omegaMin = mygetfield( prm, "omegaMin", omega0*(eps^2) );
	gNormMin = mygetfield( prm, "gNormMin", gNorm0*eps );
	iterMax = mygetfield( prm, "iterMax", 300 );
	omegaMax = mygetfield( prm, "omegaMax", omega0/eps );
	distMax = mygetfield( prm, "distMax", fNorm0/(eps*jNorm0) );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(distMax) );
		assert( isrealscalar(omegaMin) );
		assert( isrealscalar(omegaMax) );
		assert( isrealscalar(gNormMin) );
		assert( isrealscalar(iterMax) );
	endif
	%
	stepType = mygetfield( prm, "stepType", STEP_TYPE__BLIND_NEWTON );
	jupdateType = mygetfield( prm, "jupdateType", JUPDATE_TYPE__RECALC );
	cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	muRegu = mygetfield( prm, "muRegu", sqrt(eps) );
	sMin = mygetfield( prm, "sMin", 0.0001 );
	allowUphillSteps = mygetfield( prm, "allowUphillSteps", false );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(stepType) );
		assert( isrealscalar(jupdateType) );
		assert( isrealscalar(cholSafeTol) );
		assert( 0.0 < cholSafeTol );
		assert( isrealscalar(muRegu) );
		assert( 0.0 < muRegu );
		assert( isrealscalar(sMin) );
		assert( 0.0 <= sMin );
		assert( isscalar(allowUphillSteps) );
		assert( isbool(allowUphillSteps) );
	endif
	%
	omegaFallTol = mygetfield( prm, "omegaFallTol", omega0*(eps^2) );
	deltaJNormTol = mygetfield( prm, "deltaJNormTol", jNorm0*eps );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(omegaFallTol) );
		assert( isrealscalar(deltaJNormTol) );
	endif
	%
	%
	%
	vecX_best = vecX0;
	vecF_best = vecF0;
	omega_best = omega0;
	if ( nargout >= 2 )
		datOut = [];
	endif
	%
	%
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	iterCount = 0;
	doMainLoop = true;
	while (1)
		% Set quantities.
		if ( valdLev >= VALDLEV__MEDIUM )
			assert( isrealarray(vecX,[sizeX,1]) );
			assert( isrealarray(vecF,[sizeF,1]) );
			assert( isrealarray(matJ,[sizeF,sizeX]) );
		endif
		omega = sumsq(vecF)/2.0;
		vecG = matJ'*vecF;
		matH = matJ'*matJ;
		gNorm = norm(vecG);
		hNorm = max(abs(diag(matH)));
		%
		%
		% Log progress.
		if ( verbLev >= VERBLEV__PROGRESS )
		if ( abs( iterCount - round(sqrt(iterCount))^2 ) < 0.001 )
		if ( 0 == iterCount )
			msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
			  time()-time0, iterCount, ...
			  fevalCount, jevalCount, ...
			  sumsq(vecF)/2.0, -1.0 ) );
		else
			msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
			  time()-time0, iterCount, ...
			  fevalCount, jevalCount, ...
			  sumsq(vecF)/2.0, (sumsq(vecF_prev)-sumsq(vecF))/2.0 ) );
		endif
		endif
		endif
		%
		%
		% Check pre-iter stop crit.
		if ( omega <= omegaMin )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: omega <= omegaMin." );
			doMainLoop = false;
		endif
		if ( gNorm <= gNormMin )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "WEAK SUCCESS: gNorm <= gNormMin." );
			doMainLoop = false;
		endif
		if ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			doMainLoop = false;
		endif
		if ( omega >= omegaMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: omega >= omegaMax." );
			doMainLoop = false;
		endif
		if ( norm(vecX-vecX0) >= distMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecX-vecX0) >= distMax." );
			doMainLoop = false;
		endif
		if ( ~doMainLoop )
			break;
		endif
		%
		%
		% Do iteration.
		iterCount++;
		%
		%
		% Generate delta.
		switch( stepType )
		case STEP_TYPE__BLIND_NEWTON
			[ matR, cholFlag ] = chol( matH );
			if ( 0~=cholFlag || min(diag(matR)) <= cholSafeTol*max(abs(diag(matR))) )
				assert( 0.0 ~= hNorm );
				[ matR, cholFlag ] = chol( matH + muRegu*hNorm*eye(sizeX,sizeX) );
				if ( 0~=cholFlag || min(diag(matR)) <= cholSafeTol*max(abs(diag(matR))) )
					error( "Cholesky factorization failed after regularization; this should be impossible." );
				endif
			endif
			vecDelta = matR \ ( matR' \ (-vecG) );
			% This can lead to crazy large steps that can cause issues elsewhere.
		case STEP_TYPE__BLIND_GRAD_MIN
			error( "STEP_TYPE__BLIND_GRAD_MIN is not implemented yet." );
		case STEP_TYPE__SCAN_LEV_MIN
			funchDeltaOfS = @(s)( ( s*matH + (1.0-s)*hNorm*eye(sizeX,sizeX) ) \ ( -s*vecG ) );
			funchLevOmegaOfS = @(s)(0.5*sumsq(funchF( vecX + funchDeltaOfS(s) )));
			[ matR, cholFlag ] = chol( matH );
			if ( 0~=cholFlag || min(diag(matR)) <= cholSafeTol*max(abs(diag(matR))) )
				sOfMin = fminbnd( funchLevOmegaOfS, sMin, 0.9999 );
			else
				sOfMin = fminbnd( funchLevOmegaOfS, sMin, 1.0 );
			endif
			vecDelta = funchDeltaOfS(sOfMin);
			% This needs validation.
		otherwise
			error( "Invalid value of stepType." );
		endswitch
		%
		%
		% Look at delta.
		if ( valdLev >= VALDLEV__MEDIUM )
			assert( isrealarray(vecDelta,[sizeX,1]) );
		endif
		vecX_trial = vecX + vecDelta;
		vecDelta = vecX_trial - vecX; % This can be important due to FPfx.
		%deltaNormSq = sumsq(vecDelta);
		deltaNormSq = vecDelta'*vecDelta; % This can be important too???
		deltaNorm = sqrt(deltaNormSq);
		assert( 0.0 < deltaNorm );
		%
		vecF_trial = funchF( vecX_trial );
		fevalCount++;
		if ( valdLev >= VALDLEV__MEDIUM )
			assert( isrealarray(vecF_trial,[sizeX,1]) );
		endif
		omega_trial = sumsq(vecF_trial)/2.0;
		if ( omega_trial < omega_best )
			vecX_best = vecX_trial;
			vecF_best = vecF_trial;
			omega_best = omega_trial;
		endif
		msgif( verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, sprintf( "||vecDelta|| = %10.3e,  omega_trial = %10.3e, fall = %10.3e", deltaNorm, omega_trial, omega-omega_trial ) );
		%
		%
		% Decide whether or not to move to the new point.
		if ( omega_trial < omega )
			% We could consider Wolfe or compare-to-model criteria, but, this is enough for now.
			msgif( verbLev >= VERBLEV__PROGRESS+15, __FILE__, __LINE__, "Accepting delta as next because it decreases omega." );
			acceptDeltaAsNext = true;
		elseif ( allowUphillSteps )
			msgif( verbLev >= VERBLEV__PROGRESS+15, __FILE__, __LINE__, "Accepting delta as next because uphill steps are allowed." );
			acceptDeltaAsNext = true;
		else
			msgif( verbLev >= VERBLEV__PROGRESS+15, __FILE__, __LINE__, "Rejecting delta as next." );
			acceptDeltaAsNext = false;
		endif
		if ( acceptDeltaAsNext )
			vecX_next = vecX_trial;
			vecF_next = vecF_trial;
			omegaFall = omega - omega_trial; % Negative if uphill.
			collected_vecXVals = [ collected_vecXVals, vecX_trial ];
			collected_vecFVals = [ collected_vecFVals, vecF_trial ];
			collected_indexOfNext = size(collected_vecXVals,2);
		else
			vecX_next = vecX;
			vecF_next = vecF;
			omegaFall = 0.0;
			collected_vecXVals = [ collected_vecXVals, vecX_trial ];
			collected_vecFVals = [ collected_vecFVals, vecF_trial ];
		endif
		%
		%
		% Decide whether or not to update our Jacobian.
		switch( jupdateType )
		case JUPDATE_TYPE__NONE
			matJ_next = matJ;
		case JUPDATE_TYPE__BROYDEN
			fooX = vecX_trial - vecX;
			fooY = vecF_trial - ( vecF + matJ*fooX );
			fooXSq = fooX'*fooX;
			assert( 0.0 ~= fooXSq );
			fooJ = (fooY*(fooX'))/fooXSq;
			matJ_next = matJ + fooJ;
			if ( valdLev >= VALDLEV__HIGH )
				assert( reldiff( matJ_next*(vecX_trial-vecX), vecF_trial-vecF ) <= sqrt(eps) );
			endif
		case JUPDATE_TYPE__BROYDEN_ALT
			% This could behave differently due to finite precision effects.
			matJ_next = matJ + ( vecF_trial - (vecF+matJ*vecDelta) )*(vecDelta'/deltaNormSq);
			% Instead dividing the full matrix by the scalar works better in one test case.
			%matJ_next = matJ + ((( vecF_trial - (vecF+matJ*vecDelta) )*(vecDelta'))/deltaNormSq);
			if ( valdLev >= VALDLEV__HIGH )
				assert( reldiff( matJ_next*(vecX_trial-vecX), vecF_trial-vecF ) <= sqrt(eps) );
			endif
		case JUPDATE_TYPE__SECANT_REORTHONORM
			error( "JUPDATE_TYPE__SECANT_REORTHONORM is not implemented yet." );
		case JUPDATE_TYPE__LESQUJ_PRIMAL
			assert( reldiff(vecX_next,collected_vecXVals(:,collected_indexOfNext)) <= eps );
			assert( reldiff(vecF_next,collected_vecFVals(:,collected_indexOfNext)) <= eps );
			lesquj_prm = [];
			lesquj_prm.indexOfPt0 = collected_indexOfNext;
			lesquj_prm.jevalDat(1).vecX = vecX0;
			lesquj_prm.jevalDat(1).vecF = vecF0;
			lesquj_prm.jevalDat(1).matJ = matJ0;
			[ lesquj_vecX0, lesquj_vecF0, lesquj_matJ0, lesquj_datOut ] = calcLesquj_basic( collected_vecXVals, collected_vecFVals, lesquj_prm );
			assert( reldiff(lesquj_vecX0,vecX_next) <= eps );
			assert( reldiff(lesquj_vecF0,vecF_next) <= eps );
			matJ_next = lesquj_matJ0;
			if ( valdLev >= VALDLEV__HIGH )
				assert( isrealarray(matJ_next,[sizeF,sizeX]) );
			endif
		case JUPDATE_TYPE__RECALC
			matJ_next = jacobs( vecX_next, funchF ); jevalCount++;
		otherwise
			error( "Invalid value of jupdateType." );
		endswitch
		deltaJNorm = sqrt(sum(sumsq(matJ_next-matJ)));
		%
		%
		vecX_prev = vecX;
		vecF_prev = vecF;
		matJ_prev = matJ;
		vecX = vecX_next;
		vecF = vecF_next;
		matJ = matJ_next;
		%
		%
		if ( (omegaFall < omegaFallTol) && (deltaJNorm < deltaJNormTol) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: (omegaFall < omegaFallTol) && (deltaJNorm < deltaJNormTol)." );
			doMainLoop = false;
		endif
	endwhile
	%
	if ( verbLev >= VERBLEV__MAIN )
	if ( 0 == iterCount )
		msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
		  time()-time0, iterCount, ...
		  fevalCount, jevalCount, ...
		  sumsq(vecF)/2.0, -1.0 ) );
	else
		msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
		  time()-time0, iterCount, ...
		  fevalCount, jevalCount, ...
		  sumsq(vecF)/2.0, (sumsq(vecF_prev)-sumsq(vecF))/2.0 ) );
	endif
	endif
endfunction

%!function [ vecF, matJ ] = funcFJ_cubyDiagTest( vecX, c )
%!	sizeX = size(vecX,1);
%!	vecXE = (1:sizeX)';
%!	vecF = (vecX-vecXE) + c*(vecX-vecXE).^3;
%!	if ( nargout >= 2 )
%!		matJ = eye(sizeX) + c*3.0*diag((vecX-vecXE).^2);
%!	endif
%!endfunction


%!test
%!	clear;
%!	commondefs;
%!	findLocMin_gnostic_jupdate2__defs;
%!	setprngstates(0);
%!	numFigs = 0;
%!	%
%!	caseNum = 40;
%!	msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
%!	switch (caseNum)
%!	case 0
%!		sizeX = 2;
%!		c_cuby = 0.0;
%!	case 10
%!		sizeX = 15;
%!		c_cuby = 0.0;
%!	case 20
%!		sizeX = 2;
%!		c_cuby = 0.01;
%!	case 30
%!		sizeX = 15;
%!		c_cuby = 0.01;
%!	case 40
%!		sizeX = 15;
%!		c_cuby = 0.1;
%!	case 100
%!		sizeX = 15;
%!		c_cuby = 1.0;
%!	case 200
%!		sizeX = 20;
%!		c_cuby = 1.0;
%!	otherwise
%!		error( "Ivalid caseNum." );
%!	endswitch
%!	funchFJ = @(dummyX)( funcFJ_cubyDiagTest( dummyX, c_cuby ) );
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	msg( __FILE__, __LINE__, "" );
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN + default ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__BROYDEN;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm );
%!	%
%!	msg( __FILE__, __LINE__, "" );
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__BROYDEN;
%!	prm.stepType = STEP_TYPE__SCAN_LEV_MIN;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm );
%!	%
%!	msg( __FILE__, __LINE__, "" );
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__LESQUJ_PRIMAL + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__LESQUJ_PRIMAL;
%!	prm.stepType = STEP_TYPE__SCAN_LEV_MIN;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm );
%!	%
%!	%
%!	msg( __FILE__, __LINE__, "" );
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__RECALC + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__RECALC;
%!	prm.stepType = STEP_TYPE__SCAN_LEV_MIN;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm );
