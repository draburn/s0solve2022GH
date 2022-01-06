function [ f, vecNablaF ] = funcOmega_surfGradOmega( vecX, funchSurf, funchOmega, deltaR, h0, prm=[] )
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
%!	thisFile = "funcOmega_surfGradOmega test: runs";
%!	commondefs;
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
%!	funchF = @(x)( funcOmega_surfGradOmega( x, funchSurf, funchOmega, deltaR, h0 ) );
%!	funchZ = @(x,y)( log(eps*omega0+funchF( [x;y] )) );
%!	echo_Z = funchZ(0.0,0.0);
%!	%
%!	numPts = 101;
%!	thetaVals = linspace(0.0,2.0*pi,numPts);
%!	vecXVals = vecXCent_surf + 3.0*bigR_surf*[ cos(thetaVals); sin(thetaVals) ];
%!	for n=1:numPts
%!		vecSVals(:,n) = funcSurf_ellip( vecXVals(:,n), bigR_surf, vecXCent_surf, matA_surf );
%!	end
%!	%
%!	%contourfunch( funchZ );
%!	contourfunch( funchZ, ax );
%!	axis equal;
%!	hold on;
%!	plot( vecSVals(1,:), vecSVals(2,:), 'ko-' );
%!	hold off;
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figure(s) look correct. ***" );
