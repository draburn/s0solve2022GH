% Function...


function [ omegaVals, vecNablaOmegaVals ] = funcOmegaWithinSurf( vecXVals, funchOmegaBase, funchSurf, prm=[] )
	if ( 3 > nargin || 4 < nargin )
		msg( __FILE__, __LINE__, "Bad nargin." );
		print_usage();
		return; % Superfluous?
	elseif ( 2 < nargout )
		msg( __FILE__, __LINE__, "Bad nargout." );
		print_usage();
		return; % Superfluous?
	end
	%
	sizeX = size(vecXVals,1);
	numVals = size(vecXVals,2);
	%
	h0 = mygetfield( prm, "h0", 0.01 );
	tau = mygetfield( prm, "tau", 0.01 );
	epsFD = mygetfield( prm, "epsFD", tau/100.0 );
	%
	debugMode = mygetfield( prm, "debugMode", false );
	if (debugMode)
		assert( isrealarray(vecXVals,[sizeX,numVals]) );
		assert( isrealscalar(h0) );
		assert( isrealscalar(tau) );
		assert( isrealscalar(epsFD) );
		assert( 0.0 < epsFD );
	end
	%
	if ( 1 < numVals )
		[ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = funchSurf( vecXVals );
		vecDVals = vecXVals - vecSVals;
		switch (nargout)
		case 1
			% This could probably be improved by limiting the "outside the surface" calculations
			% to those that are actually outside the surface -- rather than just the result --
			% but, meh.
			omegaVals_i = funchOmegaBase( vecXVals );
			[ omegaVals_s, vecNablaOmegaVals_s ] = funchOmegaBase( vecSVals );
			omegaVals_o = omegaVals_s + sum( vecDVals .* vecNablaOmegaVals_s, 1 ) ...
			  + 0.5 * sumsq( vecDVals, 1 ) .* ( h0 + sqrt(sumsq( vecNablaOmegaVals_s, 1 )) / tau );
			outFlagVals = ( sum(vecDVals.*vecNHatVals,1) > 0.0 );
			omegaVals = omegaVals_i;
			omegaVals(outFlagVals) = omegaVals_o(outFlagVals);
		case 2
			[ omegaVals_i, vecNablaOmegaVals_i ] = funchOmegaBase( vecXVals );
			[ omegaVals_s, vecNablaOmegaVals_s ] = funchOmegaBase( vecSVals );
			omegaVals_o = omegaVals_s + sum( vecDVals .* vecNablaOmegaVals_s, 1 ) ...
			  + 0.5 * sumsq( vecDVals, 1 ) .* ( h0 + sqrt(sumsq( vecNablaOmegaVals_s, 1 )) / tau );
			outFlagVals = ( sum(vecDVals.*vecNHatVals,1) > 0.0 );
			omegaVals = omegaVals_i;
			omegaVals(outFlagVals) = omegaVals_o(outFlagVals);
			%
			vecNablaOmegaVals = vecNablaOmegaVals_i;
			for n=1:numVals
				if ( ~outFlagVals(n) )
					continue;
				end
				vecD = vecDVals(:,n);
				vecNablaOmega_s = vecNablaOmegaVals_s(:,n);
				vecNHat = vecNHatVals(:,n);
				vecS = vecSVals(:,n);
				matNablaST = matNablaSTVals(:,:,n);
				%
				vecXi = vecD + ( vecNablaOmega_s * sumsq(vecD) / (2.0*tau*norm(vecNablaOmega_s)) );
				% Use finite-differencing to approximate vecNabla2OmegaBase * vecXi.
				if ( vecNHat'*vecXi > 0.0 ) % Try to make point be inside.
					[ omega_sp, vecNablaOmega_sp ] = funchOmegaBase( vecS - epsFD * vecXi );
					vecNabala2OmegaBaseXi = ( vecNablaOmega_sp - vecNablaOmega_s ) / (-epsFD);
				else
					[ omega_sp, vecNablaOmega_sp ] = funchOmegaBase( vecS + epsFD * vecXi );
					vecNabala2OmegaBaseXi = ( vecNablaOmega_sp - vecNablaOmega_s ) / epsFD;
				end
				%
				vecNablaOmegaVals(:,n) = vecNablaOmega_s + ( matNablaST * vecNabala2OmegaBaseXi ) ...
				  + ( vecD - (matNablaST*vecD) ) * ( (norm(vecNablaOmega_s)/tau) + h0 );
			end
			%
			return;
			%
			vecXiVals = vecDVals + ( vecNablaOmegaVals_s .* sumsq(vecDVals,1) ./ (2.0*tau*sqrt(sumsq(vecNablaOmegaVals_s,1))) );
			epsFDVals = epsFD * ( ones(1,numVals) - 2.0*(sum(vecNHatVals.*vecXiVals,1)>0.0) );
			[ omegaVals_sp, vecNablaOmegaVals_sp ] = funchOmegaBase( vecSVals + epsFDVals.*vecXiVals );
			vecNabala2OmegaBaseXiVals = ( vecNablaOmegaVals_sp - vecNablaOmegaVals_s ) ./ epsFDVals;
			%
			vec0Vals = vecDVals .* ( h0 + sqrt(sumsq(vecNablaOmegaVals_s,1))/tau );
			vecNablaOmegaVals = vecNablaOmegaVals_s + vec0Vals;
			vec1Vals = vecNabala2OmegaBaseXiVals - vec0Vals;
			parfor n=1:numVals
				vecNablaOmegaVals(:,n) += matNablaSTVals(:,:,n) * vec1Vals(:,n);
			end
		otherwise
		error( "Impossible case." );
		end
	return;
	end
	%
	% numVals is 1.
	[ vecS, vecNHat, vecUHat, matNablaST ] = funchSurf( vecXVals );
	vecD = vecXVals - vecS;
	if ( vecD'*vecNHat <= 0.0 )
		% We're inside.
		switch (nargout)
		case 1
			omegaVals = funchOmegaBase( vecXVals );
		case 2
			[ omegaVals, vecNablaOmegaVals ] = funchOmegaBase( vecXVals );
		otherwise
			error( "Impossible case." );
		end
	return;
	end
	%
	% numVals is 1 and we're outside.
	switch (nargout)
	case 1
		[ omega_s, vecNablaOmega_s ] = funchOmegaBase( vecS );
		omegaVals = omega_s + vecD'*vecNablaOmega_s + 0.5 * sumsq(vecD) * ( h0 + norm(vecNablaOmega_s)/tau );
	case 2
		[ omega_s, vecNablaOmega_s ] = funchOmegaBase( vecS );
		omegaVals = omega_s + vecD'*vecNablaOmega_s + 0.5 * sumsq(vecD) * ( h0 + norm(vecNablaOmega_s)/tau );
		%
		vecXi = vecD + ( vecNablaOmega_s * sumsq(vecD) / (2.0*tau*norm(vecNablaOmega_s)) );
		% Use finite-differencing to approximate vecNabla2OmegaBase * vecXi.
		if ( vecNHat'*vecXi > 0.0 ) % Try to make point be inside.
			[ omega_sp, vecNablaOmega_sp ] = funchOmegaBase( vecS - epsFD * vecXi );
			vecNabala2OmegaBaseXi = ( vecNablaOmega_sp - vecNablaOmega_s ) / (-epsFD);
		else
			[ omega_sp, vecNablaOmega_sp ] = funchOmegaBase( vecS + epsFD * vecXi );
			vecNabala2OmegaBaseXi = ( vecNablaOmega_sp - vecNablaOmega_s ) / epsFD;
		end
		%
		vecNablaOmegaVals = vecNablaOmega_s + ( matNablaST * vecNabala2OmegaBaseXi ) ...
		  + ( vecD - (matNablaST*vecD) ) * ( (norm(vecNablaOmega_s)/tau) + h0 );
	otherwise
		error( "Impossible case." );
	end
return;
end


%!test
%!	msg( __FILE__, __LINE__, "Performing basic execution test." );
%!	setprngstates(0);
%!	%
%!	for trialIndex=1:5
%!	%
%!	%
%!	sizeX = 2 + round(2.0*abs(randn()));
%!	numVals = 10 + round(5.0*abs(randn()));
%!	%
%!	vecXCent_surf = randn(sizeX,1);
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent_surf );
%!	vecXCent_base = randn(sizeX,1);
%!	funchOmegaBase = @(dummyX) funcOmegaEllip( dummyX, vecXCent_base );
%!	%
%!	vecXVals = vecXCent_surf + 0.4*randn(sizeX,numVals);
%!	%
%!	%
%!	% Test with minimal input.
%!	for n=1:numVals
%!		omega = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf );
%!		assert( isrealscalar(omega) );
%!		[ omega, vecNablaOmega ] = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf );
%!		assert( isrealscalar(omega) );
%!		assert( isrealarray(vecNablaOmega,[sizeX,1]) );
%!	end
%!	%
%!	%
%!	% Test with maximal input.
%!	bigR = 0.01 + abs(randn);
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA_surf = (0.01*eye(sizeX,sizeX)) + (matA0'*matA0);
%!	debugMode_surf = true;
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent_surf, bigR, matA_surf, debugMode_surf );
%!	%
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA_base = (0.01*eye(sizeX,sizeX)) + (matA0'*matA0);
%!	omega0 = abs(randn);
%!	omega1 = 0.01 + abs(randn);
%!	debugMode_base = true;
%!	funchOmegaBase = @(dummyX) funcOmegaEllip( dummyX, vecXCent_base, matA_base, omega0, omega1, debugMode_base );
%!	%
%!	prm = [];
%!	omegaVals = funcOmegaWithinSurf( vecXVals, funchOmegaBase, funchSurf, prm );
%!	assert( isrealarray(omegaVals,[1,numVals]) );
%!	[ omegaVals, vecNablaOmegaVals ] = funcOmegaWithinSurf( vecXVals, funchOmegaBase, funchSurf, prm );
%!	assert( isrealarray(omegaVals,[1,numVals]) );
%!	assert( isrealarray(vecNablaOmegaVals,[sizeX,numVals]) );
%!	%
%!	%
%!	% Test vectorization.
%!	epsOmega = sqrt(eps*sumsq(reshape(omegaVals,[],1)))/numVals;
%!	epsNablaOmega = sqrt(eps*sumsq(reshape(vecNablaOmegaVals,[],1)))/numVals;
%!	for n=1:numVals
%!		omega = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf, prm );
%!		assert( reldiff(omega,omegaVals(n),epsOmega) < sqrt(eps) );
%!		[ omega, vecNablaOmega ] = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf );
%!		assert( reldiff(omega,omegaVals(n),epsOmega) < sqrt(eps) );
%!		assert( reldiff(vecNablaOmega,vecNablaOmegaVals(:,n),epsNablaOmega) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );


%!test
%!	msg( __FILE__, __LINE__, "Performing finite-differencing test." );
%!	%msg( __FILE__, __LINE__, "SKIPPING TEST." ); return;
%!	setprngstates(0);
%!	%
%!	for trialIndex=1:5
%!	%
%!	%
%!	sizeX = 2 + round(2.0*abs(randn()));
%!	%
%!	vecXCent_surf = randn(sizeX,1);
%!	bigR = 0.1 + 3.0*abs(randn);
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA_surf = (eye(sizeX,sizeX)) + (matA0'*matA0);
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent_surf, bigR, matA_surf );
%!	%
%!	vecXCent_base = randn(sizeX,1);
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA_base = (eye(sizeX,sizeX)) + (matA0'*matA0);
%!	omega0 = abs(randn);
%!	omega1 = 0.1 + abs(randn);
%!	debugMode_base = true;
%!	funchOmegaBase = @(dummyX) funcOmegaEllip( dummyX, vecXCent_base, matA_base, omega0, omega1 );
%!	%
%!	prm = [];
%!	prm.h0 = 0.01;
%!	prm.tau = 0.01;
%!	funchOmega = @(dummyX) funcOmegaWithinSurf( dummyX, funchOmegaBase, funchSurf, prm );
%!	%
%!	numVals = 10 + round(5.0*abs(randn()));
%!	vecXVals = randn(sizeX,numVals);
%!	%
%!	%
%!	[ omegaVals, vecNablaOmegaVals ] = funchOmega( vecXVals );
%!	epsOmega = sqrt(eps*sumsq(reshape(omegaVals,[],1)))/numVals;
%!	epsNablaOmega = sqrt(eps*sumsq(reshape(vecNablaOmegaVals,[],1)))/numVals;
%!	for n=1:numVals
%!		[ omega, vecNablaOmega ] = funchOmega( vecXVals(:,n) );
%!		epsX = 1e-6;
%!		vecNablaOmega_fd = zeros(sizeX,1);
%!		for m=1:sizeX
%!			vecXP = vecXVals(:,n); vecXP(m) += epsX;
%!			vecXM = vecXVals(:,n); vecXM(m) -= epsX;
%!			omegaP = funchOmega( vecXP );
%!			omegaM = funchOmega( vecXM );
%!			vecNablaOmega_fd(m) = ( omegaP - omegaM ) / (2.0*epsX);
%!		end
%!		assert( reldiff(vecNablaOmega_fd,vecNablaOmegaVals(:,n),epsNablaOmega) < 100.0*epsX );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );

%!test
%!	msg( __FILE__, __LINE__, "Generating visualization." );
%!	%msg( __FILE__, __LINE__, "SKIPPING TEST." ); return;
%!	setprngstates();
%!	numFigs0 = 0;
%!	numFigs = numFigs0;
%!	setAxisEqual = true;
%!	%
%!	sizeX = 2;
%!	%
%!	vecXCent_surf = randn(sizeX,1);
%!	bigR = 0.1 + 3.0*abs(randn);
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA_surf = (eye(sizeX,sizeX)) + (matA0'*matA0);
%!	funchSurf = @(dummyX) funcSurfEllip( dummyX, vecXCent_surf, bigR, matA_surf );
%!	%
%!	% Trace surface with "pts".
%!	numPts = 1001;
%!	thetaPts = linspace(0.0,2.0*pi,numPts);
%!	vecXPts = vecXCent_surf + bigR*[ cos(thetaPts); sin(thetaPts) ];
%!	vecSPts = funchSurf( vecXPts );
%!	%
%!	x1gc = sum(vecSPts(1,:))/numPts;
%!	x2gc = sum(vecSPts(2,:))/numPts;
%!	x1gv = max(abs(vecSPts(1,:)-x1gc));
%!	x2gv = max(abs(vecSPts(2,:)-x2gc));
%!	xgv = max([ x1gv, x2gv ]);
%!	ax = [ x1gc, x1gc, x2gc, x2gc ] + 1.5*xgv*[-1,1,-1,1];
%!	%
%!	vecXCent_base = randn(sizeX,1);
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA_base = (eye(sizeX,sizeX)) + (matA0'*matA0);
%!	omega0 = abs(randn);
%!	omega1 = 0.1 + abs(randn);
%!	debugMode_base = true;
%!	funchOmegaBase = @(dummyX) funcOmegaEllip( dummyX, vecXCent_base, matA_base, omega0, omega1 );
%!	%
%!	prm = [];
%!	prm.h0 = 0.01;
%!	prm.tau = 0.01;
%!	funchOmega = @(dummyX) funcOmegaWithinSurf( dummyX, funchOmegaBase, funchSurf, prm );
%!	%
%!	isVectorized = true;
%!	%ax = [ -5.0, 5.0, -5.0, 5.0 ];
%!	numXVals = [ 51, 55 ];
%!	[ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F ] = ...
%!	  genVizGrids( funchOmega, isVectorized, ax, numXVals );
%!	%
%!	%
%!	numFigs++; figure(numFigs);
%!	%gridZ = sqrt(sqrt(gridF)); strZ = "sqrt(sqrt(omega))";
%!	gridZ = log( 1.0 + gridF.^2 ); strZ = "log( 1 + omega^2 )";
%!	contourf( gridX1, gridX2, gridZ );
%!	colormap( 0.3 + 0.7*colormap("default") );
%!	hold on
%!	plot( vecSPts(1,:), vecSPts(2,:), 'ro-', 'markersize', 2 );
%!	hold off;
%!	if (setAxisEqual)
%!		axis equal;
%!		axis equal; % Needed twice b/c of bug in Octave?
%!	end
%!	grid on;
%!	title(sprintf("%s vs (x1,x2); %10.3e ~ %10.3e", strZ, min(min(gridZ)), max(max(gridZ)) ) );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = gridD1F; strZ = "d/dx1 omega";
%!	contourf( gridCX1, gridCX2, gridZ );
%!	colormap( 0.3 + 0.7*colormap("default") );
%!	hold on
%!	plot( vecSPts(1,:), vecSPts(2,:), 'ro-', 'markersize', 2 );
%!	hold off;
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	grid on;
%!	title(sprintf("%s vs (x1,x2); %10.3e ~ %10.3e", strZ, min(min(gridZ)), max(max(gridZ)) ) );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = gridD2F; strZ = "d/dx2 omega";
%!	contourf( gridCX1, gridCX2, gridZ );
%!	colormap( 0.3 + 0.7*colormap("default") );
%!	hold on
%!	plot( vecSPts(1,:), vecSPts(2,:), 'ro-', 'markersize', 2 );
%!	hold off;
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	grid on;
%!	title(sprintf("%s vs (x1,x2); %10.3e ~ %10.3e", strZ, min(min(gridZ)), max(max(gridZ)) ) );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	%
%!	msg( __FILE__, __LINE__, sprintf("Please check figures %d ~ %d for reasonableness.", numFigs0+1, numFigs) );
%!	return
