function [ vecX, retCode, datOut ] = findObjSurfMin_onOff( vecX0, funchSurf, funchOmega, prm=[] )
	commondefs;
	thisFile = "findObjSurfMin_onOff";
	valdLev = mygetfield( prm, "valdLev", VALDLEV__MEDIUM );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__WARN );
	%valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	assert( isrealscalar(valdLev) );
	assert( isrealscalar(verbLev) );
	msg_copious( verbLev, thisFile, __LINE__, "Welcome." );
	dumpData = false;
	%
	% Validate main input.
	sizeX = size(vecX0,1);
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealarray(vecX0,[sizeX,1]) );
	end
	[ vecS0, vecNHat0, vecUHat0, matNablaST0 ] = funchSurf( vecX0 );
	[ omega0, vecNablaOmega0 ] = funchOmega( vecX0 );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealarray(vecS0,[sizeX,1]) );
		assert( isrealarray(vecNHat0,[sizeX,1]) );
		assert( isrealarray(vecUHat0,[sizeX,1]) );
		assert( isrealarray(matNablaST0,[sizeX,sizeX]) );
		assert( abs(norm(vecNHat0)-1.0) < eps075*sizeX );
		assert( abs(norm(vecUHat0)-1.0) < eps075*sizeX );
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecNablaOmega0,[sizeX,1]) );
	end
	%
	tauX = mygetfield( prm, "tauX", 1e-3 );
	epsX = mygetfield( prm, "epsX", 1e-4*tauX );
	h0_within = mygetfield( prm, "h0_within", norm(vecNablaOmega0)/tauX );
	h0_onSurf = mygetfield( prm, "h0_onSurf", 1e1*h0_within );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(tauX) );
		assert( isrealscalar(epsX) );
		assert( isrealscalar(h0_within) );
		assert( isrealscalar(h0_onSurf) );
		assert( 0.0 < tauX );
		assert( 0.0 < epsX );
		assert( 0.0 <= h0_within );
		assert( 0.0 < h0_onSurf );
	end
	funchBigF_within = @(vecX)( funcOmega_withinSurf( vecX, funchSurf, funchOmega, tauX, h0_within ) );
	%%%funchBigF_onSurf = @(vecX)( funcOmega_onSurf( vecX, funchSurf, funchOmega, h0_onSurf ) );
	funchBigF_onSurf = @(vecX)( funcOmega_onSurf_lg( vecX, funchSurf, funchOmega, tauX, h0_within ) );
	%
	iterLimit = mygetfield( prm, "iterLimit", 10 );
	%%%tolMagNablaOmega = mygetfield( prm, "tolMagNablaOmega", eps025*norm(vecNablaOmega0) + eps050*sizeX );
	tolAbsDeltaOmega = mygetfield( prm, "tolAbsDeltaOmega", eps050*abs(omega0) + eps075 );
	tolRelDeltaOmega = mygetfield( prm, "tolRelDeltaOmega", eps075 );
	tolMagNablaOmega = mygetfield( prm, "tolMagNablaOmega", eps025*norm(vecNablaOmega0) + eps050*sizeX ...
	  + sqrt( tolAbsDeltaOmega * norm(vecNablaOmega0) / tauX ) );
	tolAbsMagStep = mygetfield( prm, "tolAbsMagStep", eps050*norm(vecX0-vecS0) + eps075*(norm(vecX0)+norm(vecS0)) + eps*sizeX );
	tolRelMagStep = mygetfield( prm, "tolRelMagStep", eps075 );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(iterLimit) );
		assert( isrealscalar(tolMagNablaOmega) );
		assert( isrealscalar(tolAbsDeltaOmega) );
		assert( isrealscalar(tolRelDeltaOmega) );
	end
	%
	if ( nargout >= 3 )
		datOut = [];
	end
	%
	%
	iterCount = 0;
	vecX = vecX0;
	vecS = vecS0;
	vecNHat = vecNHat0;
	vecUHat = vecUHat0;
	matNablaST = matNablaST0;
	omega = omega0;
	vecNablaOmega = vecNablaOmega0;
	while (1)
		reCalcEveryIter = true;
		if (reCalcEveryIter)
			[ vecS, vecNHat, vecUHat, matNablaST ] = funchSurf( vecX );
			[ omega, vecNablaOmega ] = funchOmega( vecX );
		end
		%
		% Check pre-iter stop crit.
		if ( norm(vecNablaOmega) <= tolMagNablaOmega )
			msg_main( verbLev, thisFile, __LINE__, "Success: tolMagNablaOmega." );
			retCode = RETCODE__SUCCESS;
			return;
		end
		%
		iterCount++;
		if ( iterCount > iterLimit )
			msg_warn( verbLev, thisFile, __LINE__, "Imposed stop: iterLimit." );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		msg_progress( verbLev, thisFile, __LINE__, sprintf( "Starting iteration %d.", iterCount ) );
		%
		%
		%
		vecD = vecX-vecS;
		u = vecUHat'*vecD;
		if (dumpData)
			echo_u = u
			ntNablaOmega = vecNHat'*vecNablaOmega
		end
		if ( u > tauX )
			msg_error( verbLev, thisFile, __LINE__, "Point is well outside surface." );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			return;
		elseif ( u > -epsX && vecNHat'*vecNablaOmega <= 0.0 )
			msg_progress( verbLev, thisFile, __LINE__, "Point is on surface." );
			% Near surface and gradient is inward (descent is outward).
			if (dumpData)
				echo_foo1 = matNablaST'*vecNablaOmega
				echo_foo2 = norm(matNablaST'*vecNablaOmega)
				echo_foo3 = vecNablaOmega - vecNHat*(vecNHat'*vecNablaOmega)
				echo_foo4 = norm(echo_foo3)
				[ f, vecNablaBigF ] = funchBigF_onSurf( vecX )
			end
			if ( norm(matNablaST'*vecNablaOmega) <= tolMagNablaOmega )
				msg_main( verbLev, thisFile, __LINE__, "Success: surface tolMagNablaOmega." );
				retCode = RETCODE__SUCCESS;
				return;
				% This is the primary return condition.
			end
			if (0)
				normNablaOmega = norm(vecNablaOmega)
				ntNablaOmega = vecNHat'*vecNablaOmega
				echo_foo1 = matNablaST'*vecNablaOmega
				echo_foo2 = norm(matNablaST'*vecNablaOmega)
				echo_foo3 = vecNablaOmega - vecNHat*(vecNHat'*vecNablaOmega)
				echo_foo4 = norm(echo_foo3)
				echo__tolNablaOmega = tolMagNablaOmega
				[ f, vecNablaBigF ] = funchBigF_onSurf( vecX )
			end
			%%%assert( 2 >= iterCount );
			msg_flagged( verbLev, thisFile, __LINE__, sprintf( "Iter %d: performing onSurf sovle.", iterCount ) );
			funchBigF = funchBigF_onSurf;
			wasOnSurf = true;
		else
			msg_progress( verbLev, thisFile, __LINE__, "Point is inside surface." );
			%%%assert( 1 == iterCount );
			msg_flagged( verbLev, thisFile, __LINE__, sprintf( "Iter %d: performing within sovle.", iterCount ) );
			funchBigF = funchBigF_within;
			wasOnSurf = false;
		end
		%
		%
		% Do work.
		useProvidedGradients = true;
		if (useProvidedGradients)
			fminunc_opts = optimset( 'GradObj', 'on', 'TolX', 1e-12, 'TolFun', 1e-12 );
			if (1)
			vecX_trial = fminunc( funchBigF, vecX, fminunc_opts );
			else
				[ vecX_trial, fminunc_f, fminunc_info, fminunc_output, fminunc_grad, fminunc_hess ] = fminunc( funchBigF, vecX, fminunc_opts )
			end
		else
			vecX_trial = fminunc( funchBigF, vecX );
		end
		if ( valdLev >= VALDLEV__MEDIUM )
			assert( isrealarray(vecX_trial,[sizeX,1]) );
		end
		%
		%
		% Look at new point.
		[ vecS_trial, vecNHat_trial, vecUHat_trial, matNablaST_trial ] = funchSurf( vecX_trial );
		if ( valdLev >= VALDLEV__MEDIUM )
			assert( isrealarray(vecS_trial,[sizeX,1]) );
			assert( isrealarray(vecNHat_trial,[sizeX,1]) );
			assert( isrealarray(vecUHat_trial,[sizeX,1]) );
			assert( abs(norm(vecNHat_trial)-1.0) < eps075*sizeX );
			assert( abs(norm(vecUHat_trial)-1.0) < eps075*sizeX );
		end
		vecD_trial = vecX_trial - vecS_trial;
		u_trial = vecUHat_trial' * vecD_trial;
		if (0)
			msg( thisFile, __LINE__, "FEVAL original!" );
			echo__vecX = vecX
			echo__omega = omega
			echo__vecX_trial = vecX_trial
			echo__tauX = tauX
			echo__h0_within = h0_within
			echo__h0_onSurf = h0_onSurf
			[ vecS_trial, vecNHat_trial, vecUHat_trial, matNablaST_trial ] = funchSurf( vecX_trial )
			echo__vecD_trial = vecX_trial - vecS_trial
			[ omega_trial, vecNablaOmega_trial ] = funchOmega( vecX_trial )
			[ f, vecNablaBigF ] = funchBigF_within( vecX_trial )
			[ f, vecNablaBigF ] = funchBigF_onSurf( vecX_trial )
			echo__proposed_tol = sqrt( h0_within * tolAbsDeltaOmega )
		end
		%if ( u_trial > epsX && abs(u_trial) > tauX )
		if ( u_trial > 0.0 )
			msg_progress( verbLev, thisFile, __LINE__, "Pulling result to surface." );
			vecX_trial = vecS_trial;
			vecD_trial = vecX_trial - vecS_trial; % Zero.
			u_trial = vecUHat_trial' * vecD_trial; % Zero.
			if (0)
				msg( thisFile, __LINE__, "FEVAL pulled!" );
				echo__vecX = vecX
				echo__omega = omega
				echo__vecX_trial = vecX_trial
				echo__tauX = tauX
				echo__h0_within = h0_within
				echo__h0_onSurf = h0_onSurf
				[ vecS_trial, vecNHat_trial, vecUHat_trial, matNablaST_trial ] = funchSurf( vecX_trial )
				echo__vecD_trial = vecX_trial - vecS_trial
				[ omega_trial, vecNablaOmega_trial ] = funchOmega( vecX_trial )
				[ f, vecNablaBigF ] = funchBigF_within( vecX_trial )
				[ f, vecNablaBigF ] = funchBigF_onSurf( vecX_trial )
				echo__proposed_tol = sqrt( h0_within * tolAbsDeltaOmega )
			end
		end
		%
		[ omega_trial, vecNablaOmega_trial ] = funchOmega( vecX_trial );
		if ( valdLev >= VALDLEV__MEDIUM )
			assert( isrealscalar(omega_trial) );
			assert( isrealarray(vecNablaOmega_trial,[sizeX,1]) );
		end
		if ( omega_trial >= omega )
			msg_warn( verbLev, thisFile, __LINE__, "Failed to decrease omega." );
			if (0)
				echo__vecX = vecX_trial;
				echo__vecS = vecS_trial;
				echo__echo__vecNHat = vecNHat_trial;
				echo__vecUHat = vecUHat_trial;
				echo__matNablaST = matNablaST_trial;
				echo__omega = omega_trial;
				echo__vecNablaOmega = vecNablaOmega_trial;
			end
			assert(0);
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			return
		end
		%
		%
		%
		% Accept step.
		%msg( thisFile, __LINE__, "Accepting step." );
		omega_prev = omega;
		vecX = vecX_trial;
		vecS = vecS_trial;
		vecNHat = vecNHat_trial;
		vecUHat = vecUHat_trial;
		matNablaST = matNablaST_trial;
		omega = omega_trial;
		vecNablaOmega = vecNablaOmega_trial;
		if (0)
		if ( omega >= omega_prev - tolAbsDeltaOmega )
			msg_main( verbLev, thisFile, __LINE__, "Reached tolAbsDeltaOmega" );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		if ( abs(omega-omega_prev) < tolRelDeltaOmega*(abs(omega)+abs(omega_prev)) )
			msg_main( verbLev, thisFile, __LINE__, "Reached tolRelDeltaOmega" );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		end
	end
	error( "Loop broke." );
end


%!test
%!	commondefs;
%!	thisFile = "findObjSurfMin_onOff test 1";
%!	setprngstates(0);
%!	numFigs = 0;
%!	%
%!	sizeX = 5;
%!	%
%!	bigR = 10000.0;
%!	vecC = zeros(sizeX,1);
%!	funchSurf = @(x)( funcSurf_ellip( x, bigR, vecC ) );
%!	%
%!	vecXRoot = (1:sizeX)';
%!	h0 = 1.0;
%!	funchOmega = @(x)( funcOmega_ellip( x, h0, vecXRoot ) );
%!	%	
%!	vecX0 = zeros(sizeX,1);
%!	vecX = findObjSurfMin_onOff( vecX0, funchSurf, funchOmega );
%!	assert( norm(vecX-vecXRoot) <= (eps^0.50)*(norm(vecX)+norm(vecXRoot)) );


%!test
%!	commondefs;
%!	thisFile = "findObjSurfMin_onOff test 2";
%!	setprngstates(10801488); % Nice.
%!	numFigs = 0;
%!	%
%!	sizeX = 2;
%!	%
%!	bigR = 1.0 + abs(randn());
%!	vecC = randn(sizeX,1);
%!	funchSurf = @(x)( funcSurf_ellip( x, bigR, vecC ) );
%!	%
%!	vecXRoot = randn(sizeX,1);
%!	h0 = abs(randn());
%!	funchOmega = @(x)( funcOmega_ellip( x, h0, vecXRoot ) );
%!	%
%!	theta = 2*pi*rand();
%!	vecX0 = vecC;
%!	[ omega0, vecNablaOmega0 ] = funchOmega( vecX0 );
%!	tauX = 1e-2;
%!	h0 = norm(vecNablaOmega0)/tauX;
%!	funchF = @(x)( funcOmega_withinSurf( x, funchSurf, funchOmega, tauX, h0 ) );
%!	%
%!	vecXF = findObjSurfMin_onOff( vecX0, funchSurf, funchOmega );
%!	%
%!	numPts = 101;
%!	thetaVals = linspace(0.0,2.0*pi,numPts);
%!	vecXVals = vecC + 3.0*bigR*[ cos(thetaVals); sin(thetaVals) ];
%!	for n=1:numPts
%!		vecSVals(:,n) = funchSurf( vecXVals(:,n) );
%!	end
%!	funchZ = @(x,y)(log( eps*omega0 + funchF([x;y]) ) );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourfunch( funchZ );
%!	axis equal;
%!	hold on;
%!	plot( vecSVals(1,:), vecSVals(2,:), 'ko-' );
%!	plot( ...
%!	  [ vecX0(1), vecXF(1) ], [ vecX0(2), vecXF(2) ], '-', ...
%!	  vecX0(1), vecX0(2), 's', 'linewidth', 3, 'markersize', 15, ...
%!	  vecXF(1), vecXF(2), 'x', 'linewidth', 3, 'markersize', 15 );
%!	hold off;
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figure(s) look correct. ***" );



%!test
%!	commondefs;
%!	thisFile = "findObjSurfMin_onOff test 2";
%!	setprngstates();
%!	numFigs = 1;
%!	sizeX = 2;
%!	%
%!	sizeF_surf = 2;
%!	bigR_surf = 1.0 + abs(randn());
%!	vecXCent_surf = randn(sizeX,1);
%!	matA_surf = randn(sizeF_surf,sizeX);
%!	funchSurf = @(x)( funcSurf_ellip( x, bigR_surf, vecXCent_surf, matA_surf ) );
%!	%
%!	vecXCent_omega = randn(sizeX,1);
%!	h0_omega = abs(randn());
%!	funchOmega = @(x)( funcOmega_ellip( x, h0_omega, vecXCent_omega ) );
%!	%
%!	theta = 2*pi*rand();
%!	vecX0 = vecXCent_surf;
%!	[ omega0, vecNablaOmega0 ] = funchOmega( vecX0 );
%!	tauX_combo = 1e-2;
%!	h0_combo = norm(vecNablaOmega0)/tauX_combo;
%!	funchF = @(x)( funcOmega_withinSurf( x, funchSurf, funchOmega, tauX_combo, h0_combo ) );
%!	%
%!	vecXF = findObjSurfMin_onOff( vecX0, funchSurf, funchOmega );
%!	%
%!	numPts = 101;
%!	thetaVals = linspace(0.0,2.0*pi,numPts);
%!	vecXVals = vecXCent_surf + 3.0*bigR_surf*[ cos(thetaVals); sin(thetaVals) ];
%!	for n=1:numPts
%!		vecSVals(:,n) = funchSurf( vecXVals(:,n) );
%!	end
%!	funchZ = @(x,y)(log( eps*omega0 + funchF([x;y]) ) );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourfunch( funchZ );
%!	axis equal;
%!	hold on;
%!	plot( vecSVals(1,:), vecSVals(2,:), 'ko-' );
%!	plot( ...
%!	  [ vecX0(1), vecXF(1) ], [ vecX0(2), vecXF(2) ], '-', ...
%!	  vecX0(1), vecX0(2), 's', 'linewidth', 3, 'markersize', 15, ...
%!	  vecXF(1), vecXF(2), 'x', 'linewidth', 3, 'markersize', 15 );
%!	hold off;
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figure(s) look correct. ***" );
