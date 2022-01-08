function [ vecX, retCode, datOut ] = findObjSurfMin( vecX0, funchSurf, funchOmega, prm=[] )
	commondefs;
	thisFile = "findObjSurfMin";
	valdLev = mygetfield( prm, "valdLev", VALDLEV__MEDIUM );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__WARN );
	%valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	assert( isrealscalar(valdLev) );
	assert( isrealscalar(verbLev) );
	msg_copious( verbLev, thisFile, __LINE__, "Welcome." );
	%
	% Validate main input.
	sizeX = size(vecX0,1);
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealarray(vecX0,[sizeX,1]) );
	end
	if ( valdLev >= VALDLEV__MEDIUM )
		[ vecS0, vecU0, vecV0, matNablaST0 ] = funchSurf( vecX0 );
		[ omega0, vecNablaOmega0 ] = funchOmega( vecX0 );
		assert( isrealarray(vecS0,[sizeX,1]) );
		assert( isrealarray(vecU0,[sizeX,1]) );
		assert( isrealarray(vecV0,[sizeX,1]) );
		assert( isrealarray(matNablaST0,[sizeX,sizeX]) );
		assert( abs(norm(vecU0)-1.0) < eps075*sizeX );
		assert( abs(norm(vecV0)-1.0) < eps075*sizeX );
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecNablaOmega0,[sizeX,1]) );
	end
	%
	%
	% Set up the objective function that combines the input funchSurf and funchOmega.
	tauX = mygetfield( prm, "tauX", 1e-2 );
	h0 = mygetfield( prm, "h0", norm(vecNablaOmega0)/tauX );
	epsX = mygetfield( prm, "epsX", 1e-5 );
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealscalar(tauX) );
		assert( isrealscalar(h0) );
		assert( isrealscalar(epsX) );
	end
	funchBigF = @(vecX)( funcOmega_withinSurf( vecX, funchSurf, funchOmega, tauX, h0, epsX ) );
	%
	%
	% Do work.
	useProvidedGradients = true;
	if (useProvidedGradients)
		fminunc_opts = optimset( 'GradObj', 'off' );
		vecX = fminunc( funchBigF, vecX0, fminunc_opts );
	else
		vecX = fminunc( funchBigF, vecX0 );
	end
	if ( valdLev >= VALDLEV__LOW )
		assert( isrealarray(vecX,[sizeX,1]) );
	end
return;
end


%!test
%!	commondefs;
%!	thisFile = "findObjSurfMin test 1";
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
%!	vecX = findObjSurfMin( vecX0, funchSurf, funchOmega );
%!	assert( norm(vecX-vecXRoot) <= (eps^0.50)*(norm(vecX)+norm(vecXRoot)) );


%!test
%!	commondefs;
%!	thisFile = "findObjSurfMin test 2";
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
%!	vecXF = findObjSurfMin( vecX0, funchSurf, funchOmega );
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
%!	thisFile = "findObjSurfMin test 2";
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
%!	vecXF = findObjSurfMin( vecX0, funchSurf, funchOmega );
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
