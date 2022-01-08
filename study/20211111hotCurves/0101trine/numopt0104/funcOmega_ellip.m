function [ omega, vecNablaOmega, matNabla2Omega, aryNabla3Omega ] = funcOmega_ellip( vecX, h0, vecXCent, matA )
	if ( 3 >= nargin )
		vecD = vecX - vecXCent;
		omega = ( vecD' * vecD ) * h0 / 2.0;
		if ( 2 <= nargout )
			vecNablaOmega = vecD;
			if ( 3 <= nargout )
				matNabla2Omega = h0 * eye( length(vecX), length(vecX) );
				if ( 4 <= nargout )
					aryNabla3Omega = zeros( length(vecX), length(vecX), length(vecX) );
				end
			end
		end
	else
		vecB = matA * ( vecX - vecXCent );
		omega = ( vecB' * vecB ) * h0 / 2.0;
		if ( 2 <= nargout )
			vecNablaOmega = h0 * ( matA' * vecB );
			if ( 3 <= nargout )
				matNabla2Omega = h0 * ( matA' * matA );
				if ( 4 <= nargout )
					aryNabla3Omega = zeros( length(vecX), length(vecX), length(vecX) );
				end
			end
		end
	end
return;
end


%!test
%!	thisFile = "funcOmega_ellip: test execution";
%!	h0 = 1.0;
%!	vecXCent = [ 0.0; 0.0 ];
%!	%
%!	vecX = [ 1.0; 2.0 ];
%!	[ f ] = funcOmega_ellip( vecX, h0, vecXCent );
%!	[ f, vecDF ] = funcOmega_ellip( vecX, h0, vecXCent );
%!	[ f, vecDF, matD2F ] = funcOmega_ellip( vecX, h0, vecXCent );
%!	[ f, vecDF, matD2F, aryD3F ] = funcOmega_ellip( vecX, h0, vecXCent );
%!	%
%!	matA = [ 1.0, 0.0; 1.0, 2.0 ];
%!	[ f, vecDF, matD2F, aryD3F ] = funcOmega_ellip( vecX, h0, vecXCent, matA );


%!test
%!	thisFile = "funcOmega_ellip test: viz";
%!	setprngstates();
%!	numFigs = 0;
%!	%
%!	sizeX = 2;
%!	%
%!	switch 3
%!	case 1
%!		sizeF = sizeX;
%!		h0 = 1.0
%!		vecXCent = zeros(sizeX,1)
%!		matA = eye(sizeX,sizeX)
%!	case 2
%!		sizeF = sizeX;
%!		h0 = 2.0
%!		vecXCent = [ 1.0; 1.0 ]
%!		matA = [ 2.0, 1.0; 0.0, 1.0 ]
%!	case 3
%!		sizeF = 2;
%!		h0 = abs(randn());
%!		vecXCent = randn(sizeX,1);
%!		matA = randn(sizeF,sizeX);
%!	end
%!	setAxisEqual = true;
%!	%
%!	numTestVals = 10;
%!	vecXTestVals = vecXCent + randn(sizeX,numTestVals);
%!	epsX = 1e-4;
%!	%for n=1:numTestVals
%!	for n=1:1
%!		vecX00 = vecXTestVals(:,n);
%!		vecXP0 = vecXTestVals(:,n);
%!		vecXM0 = vecXTestVals(:,n);
%!		vecX0P = vecXTestVals(:,n);
%!		vecX0M = vecXTestVals(:,n);
%!		vecXP0(1) += epsX;
%!		vecXM0(1) -= epsX;
%!		vecX0P(2) += epsX;
%!		vecX0M(2) -= epsX;
%!		[ omega00, vecNablaOmega ] = funcOmega_ellip( vecX00, h0, vecXCent, matA );
%!		omegaP0 = funcOmega_ellip( vecXP0, h0, vecXCent, matA );
%!		omegaM0 = funcOmega_ellip( vecXM0, h0, vecXCent, matA );
%!		omega0P = funcOmega_ellip( vecX0P, h0, vecXCent, matA );
%!		omega0M = funcOmega_ellip( vecX0M, h0, vecXCent, matA );
%!		vecNablaOmegaFD = [ (omegaP0-omegaM0)/(2.0*epsX); (omega0P-omega0M)/(2.0*epsX) ];
%!		assert( norm(vecNablaOmegaFD-vecNablaOmega) < 1e-8*(norm(vecNablaOmegaFD)+norm(vecNablaOmega)) );
%!	end
%!	%
%!	msg( thisFile, __LINE__, "" );
%!	msg( thisFile, __LINE__, "*** Consider adding Hessian check. ***" );
%!	msg( thisFile, __LINE__, "" );
%!	%
%!	isVectorized = false;
%!	ax = [ -5.0, 5.0, -5.0, 5.0 ];
%!	numXVals = [ 31, 33 ];
%!	[ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F ] = genVizGrids( ...
%!	  @(x)( funcOmega_ellip( x, h0, vecXCent, matA ) ), ...
%!	  isVectorized, ax, numXVals );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridX1, gridX2, sqrt(gridF) );
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	grid on;
%!	title( "sqrt(F) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%%% Why do I need to call setAxisEqual twice???
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, sqrt( gridD1F.^2 + gridD2F.^2 ) );
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	grid on;
%!	title("||nablaF|| vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, atan2( gridD2F, gridD1F ) );
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	grid on;
%!	title("theta(nablaF) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( gridCX1(:,1), gridCX2(1,:), atan2( gridD2F, gridD1F )' );
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	set(gca,'ydir','normal');
%!	colormap(hsv);
%!	grid on;
%!	title( "theta(nablaF) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( gridCX1(:,1), gridCX2(1,:), xygrids2img( gridD1F, gridD2F ) );
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	set(gca,'ydir','normal');
%!	grid on;
%!	title( "nablaF vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figure(s) look correct. ***" );
