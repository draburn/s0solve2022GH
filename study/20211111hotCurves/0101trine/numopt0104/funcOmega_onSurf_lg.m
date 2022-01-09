function [ f, vecNablaF ] = funcOmega_onSurf_lg( vecX, funchSurf, funchOmega, tauX, h0, prm=[] )
	% "lg" stands for "Local gradient".
	thisFile = "funcOmega_onSurf_lg";
	if ( 1 == nargout )
		
	end
	[ vecS, vecNHat, vecUHat, matNablaST ] = funchSurf( vecX );
	[ omega, vecNablaOmega ] = funchOmega( vecS );
	%
	vecD = vecX-vecS;
	h1 = h0 + (norm(vecNablaOmega)/tauX);
	f = omega + ( 0.5 * h1 * sumsq(vecD) );
	if ( 2 <= nargout )
		vecXi = vecNablaOmega * ( sumsq(vecD) / (2.0*tauX*norm(vecNablaOmega)) );
		epsFD = 1E-4;
		if ( vecNHat'*vecXi > 0.0 )
			epsFD = -epsFD;
		end
		[ omegaP, vecNablaOmegaP ] = funchOmega( vecS + epsFD * vecXi );
		vecNabla2OmegaXi = ( vecNablaOmegaP - vecNablaOmega ) / epsFD;
		vecNablaF = h1 * vecD + ( matNablaST * ( vecNablaOmega + vecNabla2OmegaXi - (h1*vecD) ) );
	end
return;
end


%!test
%!	thisFile = "funcOmega_onSurf_lg test: runs";
%!	setprngstates(0); ax=[ -3, 3, -3, 3 ];
%!	%setprngstates(); ax=[ -5, 5, -5, 5 ];
%!	%setprngstates(93978080); ax=[ -5, 5, -5, 5 ];
%!	%setprngstates(11208128); ax = [ -1, 0.5, 1.5, 3.0 ] % Forked min. 
%!	%setprngstates(56196480); ax = [ -5, 5, -5, 5 ] % Pretty.
%!	%setprngstates(1283408); ax = [ -5, 5, -5, 5 ] % Local max? contourfunch issue?
%!	%setprngstates(85681808);  ax = [ -0.5, 2.5, -2.5, 0.5 ] % Nice. Clean. Simple.
%!	%setprngstates(72144288); ax=[ -5, 5, -5, 5 ]; % Nice and complex.
%!	%setprngstates(3859184); ax=[ -5, 5, -5, 5 ]; % Cute boxy thing
%!	% 
%!	numFigs = 0;
%!	sizeX = 2;
%!	%
%!	switch 2
%!	case 1
%!		bigR_surf = 1.0;
%!		vecXCent_surf = zeros(sizeX,1);
%!		matA_surf = eye(sizeX,sizeX);
%!		funchSurf = @(x)( funcSurf_ellip( x, bigR_surf, vecXCent_surf, matA_surf ) );
%!	case 2
%!		sizeF_surf = 2;
%!		bigR_surf = abs(randn());
%!		vecXCent_surf = randn(sizeX,1);
%!		matA_surf = randn(sizeF_surf,sizeX);
%!		funchSurf = @(x)( funcSurf_ellip( x, bigR_surf, vecXCent_surf, matA_surf ) );
%!	end
%!	%
%!	switch 2
%!	case 1
%!		h0_omega = 1.0;
%!		vecXCent_omega = ones(sizeX,1);
%!		matA_omega = eye(sizeX,sizeX);
%!		funchOmega = @(x)( funcOmega_ellip( x, h0_omega, vecXCent_omega, matA_omega ) );
%!	case 2
%!		sizeF_omega = 2;
%!		h0_omega = abs(randn());
%!		vecXCent_omega = randn(sizeX,1);
%!		matA_omega = randn(sizeF_omega,sizeX);
%!		funchOmega = @(x)( funcOmega_ellip( x, h0_omega, vecXCent_omega, matA_omega ) );
%!	end
%!	omega = funchOmega( [0.0;0.0] );
%!	%
%!	tauX = bigR_surf / 10.0;
%!	[ omega0, vecNablaOmega0 ] = funchOmega( vecXCent_surf );
%!	h0 = norm(vecNablaOmega0) / tauX;
%!	%
%!	numPts = 101;
%!	thetaVals = linspace(0.0,2.0*pi,numPts);
%!	vecXVals = vecXCent_surf + 3.0*bigR_surf*[ cos(thetaVals); sin(thetaVals) ];
%!	for n=1:numPts
%!		vecSVals(:,n) = funcSurf_ellip( vecXVals(:,n), bigR_surf, vecXCent_surf, matA_surf );
%!	end
%!	%
%!	funchF = @(x)( funcOmega_onSurf_lg( x, funchSurf, funchOmega, tauX, h0 ) );
%!	%
%!	%[ omega00, vecNablaOmega ] = funchF( [1.0;0.0] )
%!	%
%!	numCLevs = 31;
%!	%
%!	isVectorized = false;
%!	numXVals = [ 51, 53 ];
%!	[ gridX1, gridX2, gridOmega ] = genVizGrids( funchOmega, isVectorized, ax, numXVals );
%!	numFigs++; figure(numFigs);
%!	contourf( gridX1, gridX2, log(eps+gridOmega), numCLevs );
%!	colormap(mycmap);
%!	grid on;
%!	title("log( eps + omega ) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	[ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F ] = genVizGrids( funchF, isVectorized, ax, numXVals );
%!	[ gridX1, gridX2, gridOmega ] = genVizGrids( funchOmega, isVectorized, ax, numXVals );
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
%!	contourf( gridCX1, gridCX2, atan2( gridD2F, gridD1F ), numCLevs );
%!	colormap(hsv);
%!	hold on;
%!	plot( ...
%!	  vecSVals(1,:), vecSVals(2,:), 'ko-', ...
%!	  vecXCent_surf(1), vecXCent_surf(2), 's', 'linewidth', 3, 'markersize', 15, ...
%!	  vecXCent_omega(1), vecXCent_omega(2), 'x', 'linewidth', 3, 'markersize', 15 );
%!	hold off;
%!	grid on;
%!	title("theta(nablaF) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( gridCX1(:,1), gridCX2(1,:), atan2( gridD2F, gridD1F )' );
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
%!	numFigs++; figure(numFigs);
%!	imagesc( gridCX1(:,1), gridCX2(1,:), xygrids2img( gridD1F, gridD2F ) );
%!	hold on;
%!	contour( gridCX1, gridCX2, sqrt(sqrt( gridD1F.^2 + gridD2F.^2 )), 10 );
%!	colormap(zeros(64,3));
%!	plot( ...
%!	  vecSVals(1,:), vecSVals(2,:), 'ko-', ...
%!	  vecXCent_surf(1), vecXCent_surf(2), 's', 'linewidth', 3, 'markersize', 15, ...
%!	  vecXCent_omega(1), vecXCent_omega(2), 'x', 'linewidth', 3, 'markersize', 15 );
%!	hold off;
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "nablaF vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figure(s) look correct. ***" );
%!	%
%!	numTestVals = 100;
%!	vecXTestVals = vecXCent_surf + randn(sizeX,numTestVals);
%!	epsX = 1e-4;
%!	for n=1:numTestVals
%!		vecX00 = vecXTestVals(:,n);
%!		vecXP0 = vecX00; vecXP0(1) += epsX;
%!		vecXM0 = vecX00; vecXM0(1) -= epsX;
%!		vecX0P = vecX00; vecX0P(2) += epsX;
%!		vecX0M = vecX00; vecX0M(2) -= epsX;
%!		[ omega00, vecNablaOmega ] = funchF( vecX00 );
%!		omegaP0 = funchF( vecXP0 );
%!		omegaM0 = funchF( vecXM0 );
%!		omega0P = funchF( vecX0P );
%!		omega0M = funchF( vecX0M );
%!		vecNablaOmegaFD = [ (omegaP0-omegaM0)/(2.0*epsX); (omega0P-omega0M)/(2.0*epsX) ];
%!		assert( norm(vecNablaOmegaFD-vecNablaOmega) < 1e-6*(norm(vecNablaOmegaFD)+norm(vecNablaOmega)) );
%!	end
%!	%
%!	msg( thisFile, __LINE__, "" );
%!	msg( thisFile, __LINE__, "*** Consider adding Hessian check. ***" );
%!	msg( thisFile, __LINE__, "" );
