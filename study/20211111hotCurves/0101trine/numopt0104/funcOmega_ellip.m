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
%!	thisFile = "funcOmega_ellip test: runs";
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
%!	sizeF = 2;
%!	h0 = abs(randn());
%!	vecXCent = randn(sizeX,1);
%!	matA = randn(sizeF,sizeX);
%!	%
%!	funchZ = @(x,y)( sqrt(funcOmega_ellip( [x;y], h0, vecXCent, matA )) );
%!	numFigs++; figure(numFigs);
%!	contourfunch( funchZ );
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figure(s) look correct. ***" );
