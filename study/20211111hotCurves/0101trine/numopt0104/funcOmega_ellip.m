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
			vecNablaOmega = matA' * vecB;
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
%!	commondefs;
%!	setprngstates();
%!	numFigs = 0;
%!	%
%!	sizeX = 2;
%!	sizeF = 3;
%!	h0 = abs(randn());
%!	vecXCent = randn(sizeX,1);
%!	matA = randn(sizeF,sizeX);
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
%!	grid on;
%!	title( "sqrt(F) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, sqrt( gridD1F.^2 + gridD2F.^2 ) );
%!	grid on;
%!	title("||nablaF|| vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	imagesc( gridCX1', gridCX2', atan2( gridD2F, gridD1F )' );
%!	set(gca,'ydir','normal');
%!	colormap(hsv);
%!	grid on;
%!	title( "theta(nablaF) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figure(s) look correct. ***" );
