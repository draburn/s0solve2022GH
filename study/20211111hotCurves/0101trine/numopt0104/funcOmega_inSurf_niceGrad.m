function [ f, vecNablaF ] = funcOmega_inSurf_niceGrad( vecX, funchSurf, funchOmega, deltaR, h0, prm=[] )
	[ vecS, vecU, vecV, matNablaST ] = funchSurf( vecX );
	if ( vecV'*(vecX-vecS) < 0.0 )
		% We're inside or the surface.
		if ( 1 == nargout )
			f = funchOmega( vecX );
		elseif ( 2 == nargout )
			[ f, vecNablaF ] = funchOmega( vecX );
		end
	return;
	end
	%
	% We're on or outside the surface.
	[ omega, vecNablaOmega ] = funchOmega( vecS );
	vecD = vecX - vecS;
	%f = omega + ( vecD' * vecNablaOmega) + ( 0.5 * ( vecD' * vecD ) ) * ( norm(vecNablaOmega)/deltaR + h0 );
	f = omega + ( vecD' * vecNablaOmega) + ( 0.5 * ( vecD' * vecD ) ) * ( norm(vecNablaOmega)/deltaR + h0 );
	assert( nargout == 1 );
return;
end


%!test
%!	thisFile = "funcOmega_inSurf_niceGrad test: runs";
%!	setprngstates(); ax=[ -5, 5, -5, 5 ];
%!	%setprngstates(11208128); ax = [ -1, 0.5, 1.5, 3.0 ] % Forked min. 
%!	%setprngstates(56196480); ax = [ -5, 5, -5, 5 ] % Pretty.
%!	%setprngstates(1283408); ax = [ -5, 5, -5, 5 ] % Local max? contourfunch issue?
%!	%setprngstates(85681808);  ax = [ -0.5, 2.5, -2.5, 0.5 ] % Nice. Clean. Simple.
%!	numFigs = 0;
%!	sizeX = 2;
%!	%
%!	sizeF_surf = 2;
%!	bigR_surf = abs(randn());
%!	vecXCent_surf = randn(sizeX,1);
%!	matA_surf = randn(sizeF_surf,sizeX);
%!	vecS = funcSurf_ellip( [0.0;0.0], bigR_surf, vecXCent_surf, matA_surf );
%!	funchSurf = @(x)( funcSurf_ellip( x, bigR_surf, vecXCent_surf, matA_surf ) );
%!	vecS = funchSurf( [0.0;0.0] );
%!	%
%!	sizeF_omega = 2;
%!	h0_omega = abs(randn());
%!	vecXCent_omega = randn(sizeX,1);
%!	matA_omega = randn(sizeF_omega,sizeX);
%!	omega = funcOmega_ellip( [0.0;0.0], h0_omega, vecXCent_omega, matA_omega );
%!	funchOmega = @(x)( funcOmega_ellip( x, h0_omega, vecXCent_omega, matA_omega ) );
%!	omega = funchOmega( [0.0;0.0] );
%!	%
%!	deltaR = bigR_surf / 10.0;
%!	[ omega0, vecNablaOmega0 ] = funchOmega( vecXCent_surf );
%!	h0 = norm(vecNablaOmega0) / deltaR;
%!	%
%!	numPts = 101;
%!	thetaVals = linspace(0.0,2.0*pi,numPts);
%!	vecXVals = vecXCent_surf + 3.0*bigR_surf*[ cos(thetaVals); sin(thetaVals) ];
%!	for n=1:numPts
%!		vecSVals(:,n) = funcSurf_ellip( vecXVals(:,n), bigR_surf, vecXCent_surf, matA_surf );
%!	end
%!	%
%!	funchF = @(x)( funcOmega_inSurf_niceGrad( x, funchSurf, funchOmega, deltaR, h0 ) );
%!	isVectorized = false;
%!	numXVals = [ 51, 53 ];
%!	[ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F ] = genVizGrids( funchF, isVectorized, ax, numXVals );
%!	%
%!	numCLevs = 31;
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridX1, gridX2, sqrt(gridF), numCLevs );
%!	colormap(mycmap);
%!	hold on;
%!	plot( ...
%!	  vecSVals(1,:), vecSVals(2,:), 'ko-', ...
%!	  vecXCent_surf(1), vecXCent_surf(2), 's', 'linewidth', 3, 'markersize', 15, ...
%!	  vecXCent_omega(1), vecXCent_omega(2), 'x', 'linewidth', 3, 'markersize', 15 );
%!	hold off;
%!	grid on;
%!	title( "sqrt(F) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, gridD1F, numCLevs );
%!	colormap(mycmap);
%!	grid on;
%!	title("dF/dx1 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, gridD2F, numCLevs );
%!	colormap(mycmap);
%!	grid on;
%!	title("dF/dx2 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, sqrt(sqrt( gridD1F.^2 + gridD2F.^2 )), numCLevs );
%!	colormap(mycmap);
%!	hold on;
%!	plot( ...
%!	  vecSVals(1,:), vecSVals(2,:), 'ko-', ...
%!	  vecXCent_surf(1), vecXCent_surf(2), 's', 'linewidth', 3, 'markersize', 15, ...
%!	  vecXCent_omega(1), vecXCent_omega(2), 'x', 'linewidth', 3, 'markersize', 15 );
%!	hold off;
%!	grid on;
%!	title("sqrt(||nablaF||) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( gridCX1', gridCX2', atan2( gridD2F, gridD1F )' );
%!	hold on;
%!	plot( ...
%!	  vecSVals(1,:), vecSVals(2,:), 'ko-', ...
%!	  vecXCent_surf(1), vecXCent_surf(2), 's', 'linewidth', 3, 'markersize', 15, ...
%!	  vecXCent_omega(1), vecXCent_omega(2), 'x', 'linewidth', 3, 'markersize', 15 );
%!	hold off;
%!	set(gca,'ydir','normal');
%!	colormap(hsv);
%!	grid on;
%!	title( "theta(nablaF) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figure(s) look correct. ***" );


%!test
%!	thisFile = "funcOmega_ellip test: numerical";
%!	setprngstates();
%!	numFigs = 3;
%!	%
%!	msg( thisFile, __LINE__, "~~~ Compare analytic and numerical gradients! ~~~" );
