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
	h0 = mygetfield( prm, "h0", 0.01 );
	tau = mygetfield( prm, "tau", 0.01 );
	numVals = size(vecXVals,2);
	%
	if ( 1 < numVals )
		[ vecSVals, vecNHatVals, vecUHatVals, matNablaSTVals ] = funchSurf( vecXVals );
		vecDVals = vecXVals - vecSVals;
		switch (nargout)
		case 1
			omegaVals_i = funchOmegaBase( vecXVals );
			[ omegaVals_s, vecNablaOmegaVals_s ] = funchOmegaBase( vecSVals );
			omegaVals_o = omegaVals_s + sum( vecDVals .* vecNablaOmegaVals_s, 1 ) ...
			  + 0.5 * sumsq( vecDVals, 1 ) .* ( h0 + sqrt(sumsq( vecNablaOmegaVals_s, 1 )) / tau );
			outFlagVals = ( sum(vecDVals.*vecNHatVals,1) > 0.0 );
			omegaVals = omegaVals_i;
			omegaVals(outFlagVals) = omegaVals_o(outFlagVals);
		case 2
			[ omegaVals_i, vecNablaOmegaVals_i, matNabla2OmegaVals_i ] = funchOmegaBase( vecXVals );
			[ omegaVals_s, vecNablaOmegaVals_s, matNabla2OmegaVals_s ] = funchOmegaBase( vecSVals );
			omegaVals_o = omegaVals_s + sum( vecDVals .* vecNablaOmegaVals_s, 1 ) ...
			  + 0.5 * sumsq( vecDVals, 1 ) .* ( h0 + sqrt(sumsq( vecNablaOmegaVals_s, 1 )) / tau );
			outFlagVals = ( sum(vecDVals.*vecNHatVals,1) > 0.0 );
			omegaVals = omegaVals_i;
			omegaVals(outFlagVals) = omegaVals_o(outFlagVals);
			msg( __FILE__, __LINE__, "The rest of this case is not implemented." );
			error ( "To-do." );
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
		[ omega_s, vecNablaOmega_s, matNabla2Omega_s ] = funchOmegaBase( vecS );
		omegaVals = omega_s + svecD'*vecNablaOmega_s + 0.5 * sumsq(vecD) * ( h0 + norm(vecNablaOmega_s)/tau );
		msg( __FILE__, __LINE__, "The rest of this case is not implemented." );
		error ( "To-do." );
	otherwise
		error( "Impossible case." );
	end
return;
end


%!test
%!	msg( __FILE__, __LINE__, "Performing basic execution test." );
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = 2 + round(2.0*abs(randn()));
%!	numVals = 50 + round(10.0*abs(randn()));
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
%!		%[ omega, vecNablaOmega ] = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf );
%!		%assert( isrealscalar(omega) );
%!		%assert( isrealarray(vecNablaOmega,[sizeX,1]) );
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
%!	%[ omegaVals, vecNablaOmegaVals ] = funcOmegaWithinSurf( vecXVals, funchOmegaBase, funchSurf, prm );
%!	%assert( isrealarray(omegaVals,[1,numVals]) );
%!	%assert( isrealarray(vecNablaOmegaVals,[sizeX,numVals]) );
%!	%
%!	%
%!	% Test vectorization.
%!	epsOmega = sqrt(eps*sumsq(reshape(omegaVals,[],1)))/numVals;
%!	%epsNablaOmega = sqrt(eps*sumsq(reshape(vecNablaOmegaVals,[],1)))/numVals;
%!	for n=1:numVals
%!		omega = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf, prm );
%!		assert( reldiff(omega,omegaVals(n),epsOmega) < sqrt(eps) );
%!		%[ omega, vecNablaOmega ] = funcOmegaWithinSurf( vecXVals(:,n), funchOmegaBase, funchSurf );
%!		%assert( reldiff(omega,omegaVals(n),epsOmega) < sqrt(eps) );
%!		%assert( reldiff(vecNablaOmega,vecNablaOmegaVals(:,n),epsNablaOmega) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );


%!test
%!	msg( __FILE__, __LINE__, "Generating visualization." );
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


%!	error("NOT A TEST.");
%!	gridZ11 = (gridD1F(3:end,2:end-1) - gridD1F(1:end-2,2:end-1)) ./ ...
%!	  ( (gridCX1(3:end,2:end-1) - gridCX1(1:end-2,2:end-1)) .* sqrt(1.0+(gridD1F(2:end-1,2:end-1).^2)) );
%!	gridCCX1 = (gridCX1(3:end,2:end-1) + gridCX1(1:end-2,2:end-1))/2.0;
%!	gridCCX2 = (gridCX2(3:end,2:end-1) + gridCX2(1:end-2,2:end-1))/2.0;
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = gridZ11; strZ = "Z11";
%!	contourf( gridCCX1, gridCCX2, gridZ );
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
%!	%
%!	msg( __FILE__, __LINE__, sprintf("Please check figures %d ~ %d for reasonableness.", numFigs0+1, numFigs) );
%!	return


%!	error("NOT A TEST.");
%!	%
%!	%
%!	gridG1 = gridD1F;
%!	gridG2 = gridD2F;
%!	gridD1G1 = (gridG1(3:end,2:end-1) - gridG1(1:end-2,2:end-1))./(gridCX1(3:end,2:end-1) - gridCX1(1:end-2,2:end-1));
%!	gridD1G2 = (gridG2(3:end,2:end-1) - gridG2(1:end-2,2:end-1))./(gridCX1(3:end,2:end-1) - gridCX1(1:end-2,2:end-1));
%!	gridD2G1 = (gridG1(2:end-1,3:end) - gridG1(2:end-1,1:end-2))./(gridCX2(2:end-1,3:end) - gridCX2(2:end-1,1:end-2));
%!	gridD2G2 = (gridG2(2:end-1,3:end) - gridG2(2:end-1,1:end-2))./(gridCX2(2:end-1,3:end) - gridCX2(2:end-1,1:end-2));
%!	gridCCX1 = (gridCX1(3:end,2:end-1) + gridCX1(1:end-2,2:end-1))/2.0;
%!	gridCCX2 = (gridCX2(3:end,2:end-1) + gridCX2(1:end-2,2:end-1))/2.0;
%!	%
%!	numFigs++; figure(numFigs);
%!	%gridZ = gridD1G1.^2 + gridD1G2.^2 + gridD2G1.^2 + gridD2G2.^2; strZ = "||nabla^2 omega||";
%!	gridZ = gridD1G1; strZ = "d2 omega";
%!	contourf( gridCCX1, gridCCX2, gridZ );
%!	colormap( 0.3 + 0.7*colormap("default") );
%!	%hold on
%!	%plot( vecSPts(1,:), vecSPts(2,:), 'ro-', 'markersize', 2 );
%!	%hold off;
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	grid on;
%!	title(sprintf("%s vs (x1,x2); %10.3e ~ %10.3e", strZ, min(min(gridZ)), max(max(gridZ)) ) );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	%
%!	%
%!	return
%!	%
%!	%
%!	numFigs++; figure(numFigs);
%!	n2 = round( (1+numXVals(2))/2.0 );
%!	plot( gridX1(:,n2), gridF(:,n2), 'o-' );
%!	grid on;
%!	%
%!	numFigs++; figure(numFigs);
%!	n2 = round( (1+numXVals(2))/2.0 );
%!	plot( ...
%!	  cent(gridX1(:,n2)), diff(gridF(:,n2))./diff(gridX1(:,n2)), 'o-', ...
%!	  cent(gridX1(:,n2)), diff(log(1.0+gridF(:,n2).^2))./diff(gridX1(:,n2)), 'x-' );
%!	grid on;
%!	%
%!	%
%!	x = gridX1(:,n2);
%!	f = gridF(:,n2);
%!	df = diff(f)./diff(x);
%!	cx = cent(x);
%!	ddf = diff(df)./diff(cx);
%!	ccx = cent(cx);
%!	numFigs++; figure(numFigs);
%!	plot( ccx, ddf, 'o-', ccx, diff(log(1.0+df.^2))./diff(cx) );
%!	grid on;
%!	%
%!	%
%!	grid_foo_F = log( 0.01 + gridD1F.^2 +gridD2F.^2 );
%!	grid_foo_X = gridCX1;
%!	grid_foo_Y = gridCX2;
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = grid_foo_F; strZ = "loggy thing 0 of omega";
%!	contourf( grid_foo_X, grid_foo_Y, grid_foo_F );
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
%!	%
%!	for foon=1:1
%!		grid_foo_DFDX = ...
%!		   ( grid_foo_F(3:end,2:end-1) - grid_foo_F(1:end-2,2:end-1) ) ...
%!		 ./( grid_foo_X(3:end,2:end-1) - grid_foo_X(1:end-2,2:end-1) );
%!		grid_foo_DFDY = ...
%!		   ( grid_foo_F(2:end-1,3:end) - grid_foo_F(2:end-1,1:end-2) ) ...
%!		 ./( grid_foo_Y(2:end-1,3:end) - grid_foo_Y(2:end-1,1:end-2) );
%!		grid_foo_F = log( 0.01 + grid_foo_DFDX.^2 + grid_foo_DFDY.^2 );
%!		grid_foo_X = ( grid_foo_X(3:end,2:end-1) + grid_foo_X(1:end-2,2:end-1) )/2.0;
%!		grid_foo_Y = ( grid_foo_Y(2:end-1,3:end) + grid_foo_Y(2:end-1,1:end-2) )/2.0;
%!	end
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = grid_foo_F; strZ = "loggy thing 1 of omega";
%!	contourf( grid_foo_X, grid_foo_Y, grid_foo_F );
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
%!	return
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = log(gridD1F.^2+gridD2F.^2); strZ = "log(||nablaOmega||)";
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
%!	%
%!	gridD1D1F = (gridD1F(3:end,2:end-1) - gridD1F(1:end-2,2:end-1))./(gridCX1(3:end,2:end-1) - gridCX1(1:end-2,2:end-1));
%!	gridD1D2F = (gridD2F(3:end,2:end-1) - gridD2F(1:end-2,2:end-1))./(gridCX1(3:end,2:end-1) - gridCX1(1:end-2,2:end-1));
%!	gridD2D1F = (gridD1F(2:end-1,3:end) - gridD1F(2:end-1,1:end-2))./(gridCX2(2:end-1,3:end) - gridCX2(2:end-1,1:end-2));
%!	gridD2D2F = (gridD2F(2:end-1,3:end) - gridD2F(2:end-1,1:end-2))./(gridCX2(2:end-1,3:end) - gridCX2(2:end-1,1:end-2));
%!	gridCCX1 = (gridCX1(3:end,2:end-1) + gridCX1(1:end-2,2:end-1))/2.0;
%!	gridCCX2 = (gridCX2(3:end,2:end-1) + gridCX2(1:end-2,2:end-1))/2.0;
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = log(gridD1D1F.^2 + gridD1D2F.^2 + gridD2D1F.^2 + gridD2D2F.^2 ); strZ = "2*log(||nabla2 Omega||)";
%!	contourf( gridCCX1, gridCCX2, gridZ );
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
%!	msg( __FILE__, __LINE__, sprintf("Please check figures %d ~ %d for reasonableness.", numFigs0+1, numFigs) );
