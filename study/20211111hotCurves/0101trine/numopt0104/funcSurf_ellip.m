% Usage examples...
%    vecS = funchSurf_ellip( vecX, bigR, vecXCent );
%    [ vecS, vecNHat, vecUHat, matNablaST ] = funcSurf_ellip( vecX, bigR, vecXCent, matA )
% Calculate the "surface function" for the ellipse specified by vecXCent and matA / bigR,
% with vecS = vecXCent + ( vecX - vecXCent ) * bigR / norm( matA * ( vecX - vecXCent ) ).
% vecNHat is the outward surface normal.
% vecUHat is the outward direction along which vecS is locally constant.

function [ vecS, vecNHat, vecUHat, matNablaST ] = funcSurf_ellip( vecX, bigR, vecXCent, matA )
	if ( nargin <= 2 )
		print_usage( "funcSurf_ellip" );
		return;
	elseif ( nargin > 4 )
		print_usage( "funcSurf_ellip" );
		return;
	end
	%
	thisFile = "funcSurf_ellip";
	%
	if ( 3 >= nargin )
		vecD = vecX - vecXCent;
		d = norm(vecD);
		if ( 0.0 == d )
			vecD(1) = 1.0;
			d = 1.0;
		end
		vecNHat = vecD/d; % Surface normal for a circle.
		vecS = vecXCent + (bigR*vecNHat);
		if ( 3 <= nargout )
			vecUHat = vecNHat; % Also surface normal.
			if ( 4 <= nargout )
				matNablaST = ( (bigR/d) * eye(length(vecX),length(vecX)) ) - ( ((bigR/(d^3))*vecD) * (vecD') );
				matNablaST = ( matNablaST' + matNablaST ) / 2.0;
			end
		end
	else
		vecD = vecX - vecXCent;
		d = norm(vecD);
		if ( 0.0 == d )
			vecD(1) = 1.0;
			d = 1.0;
		end
		%
		vecB = matA * vecD;
		b = norm(vecB);
		if ( 0.0 == b )
			msg( thisFile, __LINE__, "Input matA is singular." );
		end
		vecS = vecXCent + (vecD*(bigR/b));
		%
		if ( 3 <= nargout )
			vecP = matA' * vecB;
			p = norm(vecP);
			if ( 0.0 == p )
				msg( thisFile, __LINE__, "Input matA is singular." );
			end
			vecNHat = vecP/p; % Surface normal.
			vecUHat = vecD/d; % vecS is constant if vecX moves in this direction.
			if ( 4 <= nargout )
				matNablaST = ((bigR/b)*eye(length(vecX),length(vecX))) - ( ((bigR/(b^3))*vecP) * (vecD') );
			end
		end
	end
end


%!test
%!	thisFile = "funcOmega_ellip: execution test";
%!	%
%!	sizeX = 2;
%!	vecX = [ 1.0; 1.0 ];
%!	bigR = 1.0;
%!	vecC = [ 0.0; 0.0 ];
%!	[ vecS, vecNHat, vecUHat, matNablaST ] = funcSurf_ellip( vecX, bigR, vecC );
%!	assert( isrealarray(vecS,[sizeX,1]) );
%!	assert( isrealarray(vecNHat,[sizeX,1]) );
%!	assert( isrealarray(vecUHat,[sizeX,1]) );
%!	assert( isrealarray(matNablaST,[sizeX,sizeX]) );
%!	%
%!	matB = [ 1.0, 0.0; 0.0, 2.0 ];
%!	[ vecS, vecNHat, vecUHat, matNablaST ] = funcSurf_ellip( vecX, bigR, vecC, matB );
%!	assert( isrealarray(vecS,[sizeX,1]) );
%!	assert( isrealarray(vecNHat,[sizeX,1]) );
%!	assert( isrealarray(vecUHat,[sizeX,1]) );
%!	assert( isrealarray(matNablaST,[sizeX,sizeX]) );


%!test
%!	thisFile = "funcOmega_ellip: manual visual test";
%!	commondefs;
%!	setprngstates();
%!	numFigs = 0;
%!	%
%!	sizeX = 2;
%!	sizeF = 2;
%!	bigR = abs(randn);
%!	vecXCent = randn(sizeX,1);
%!	matA = randn(sizeF,sizeX);
%!	%
%!	%
%!	numPts = 100;
%!	vecXVals = vecXCent + randn(sizeX,numPts);
%!	parfor n=1:numPts
%!		vecSVals(:,n) = funcSurf_ellip( vecXVals(:,n), bigR, vecXCent );
%!	end
%!	%
%!	numFigs++; figure(numFigs);
%!	n = 1;
%!	plot( [ vecXVals(1,:); vecSVals(1,:) ], [ vecXVals(2,:); vecSVals(2,:) ], 'x-', 'markersize', 10 );
%!	axis equal;
%!	hold on;
%!	plot( vecSVals(1,:), vecSVals(2,:), 'ko', 'linewidth', 2 );
%!	hold off;
%!	grid on;
%!	%
%!	%
%!	parfor n=1:numPts
%!		vecSVals(:,n) = funcSurf_ellip( vecXVals(:,n), bigR, vecXCent, matA );
%!	end
%!	%
%!	numFigs++; figure(numFigs);
%!	n = 1;
%!	plot( [ vecXVals(1,:); vecSVals(1,:) ], [ vecXVals(2,:); vecSVals(2,:) ], 'x-', 'markersize', 10 );
%!	axis equal;
%!	hold on;
%!	plot( vecSVals(1,:), vecSVals(2,:), 'ko', 'linewidth', 2 );
%!	hold off;
%!	grid on;
%!	%
%!	%
%!	numPts = 101;
%!	thetaVals = linspace(0.0,2.0*pi,numPts);
%!	vecXVals = vecXCent + 3.0*bigR*[ cos(thetaVals); sin(thetaVals) ];
%!	for n=1:numPts
%!		vecSVals(:,n) = funcSurf_ellip( vecXVals(:,n), bigR, vecXCent, matA );
%!	end
%!	distScale = max([ max(abs( vecSVals(1,:) - vecXCent(1) )), max(abs( vecSVals(2,:) - vecXCent(2) )) ]);
%!	%
%!	if ( 0 == n )
%!		n = 1;
%!	end
%!	n = ceil( numPts * rand() );
%!	vecX = vecXVals(:,n);
%!	[ vecS, vecNHat, vecUHat, matNablaST ] = funcSurf_ellip( vecX, bigR, vecXCent, matA );
%!	%
%!	% Look at surface tangent.
%!	vecZ = 1e-4*(2.0*rand(sizeX,1)-1.0);
%!	[ vecSAtPlusZ ] = funcSurf_ellip( vecX + vecZ, bigR, vecXCent, matA );
%!	[ vecSAtMinusZ ] = funcSurf_ellip( vecX - vecZ, bigR, vecXCent, matA );
%!	vecTM = (vecSAtPlusZ+vecSAtMinusZ)/2.0 - distScale*(vecSAtPlusZ-vecSAtMinusZ)/norm(vecSAtPlusZ-vecSAtMinusZ);
%!	vecTP = (vecSAtPlusZ+vecSAtMinusZ)/2.0 + distScale*(vecSAtPlusZ-vecSAtMinusZ)/norm(vecSAtPlusZ-vecSAtMinusZ);
%!	%
%!	numFigs++; figure(numFigs);
%!	plot( vecSVals(1,:), vecSVals(2,:), 'ko-' );
%!	axis equal;
%!	hold on;
%!	plot( ...
%!	  vecXCent(1), vecXCent(2), 'k+', 'linewidth', 5, 'markersize', 10, ...
%!	  [ vecX(1); vecS(1) ], [ vecX(2); vecS(2) ], 's-', 'linewidth', 8, 'markersize', 10, ...
%!	  [ vecX(1); vecXCent(1) ], [ vecX(2); vecXCent(2) ], 'x-', 'linewidth', 4, 'markersize', 10 );
%!	plot( [ vecS(1); vecS(1)+distScale*vecUHat(1) ], [ vecS(2); vecS(2)+distScale*vecUHat(2) ], '^-', 'linewidth', 3, 'markersize', 8 );
%!	plot( ...
%!	  [ vecS(1); vecS(1)+distScale*vecNHat(1) ], [ vecS(2); vecS(2)+distScale*vecNHat(2) ], 'v-', 'linewidth', 3, 'markersize', 8, ...
%!	  [ vecTM(1); vecTP(1) ], [ vecTM(2); vecTP(2) ], '-', 'linewidth', 2 );
%!	hold off;
%!	legend( ...
%!	  'surface', ...
%!	  'center', 'sample pt pulled to surf', 'sample pt pulled to center', ...
%!	  's-no-change direction', ...
%!	  'surface normal direction', 'surface tangent', ...
%!	  'location', 'NorthWestOutside' );
%!	grid on;
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figures look correct. ***" );


%!test
%!	thisFile = "funcOmega_ellip test: numerical test";
%!	commondefs;
%!	setprngstates();
%!	numFigs = 0;
%!	%
%!	sizeX = 2;
%!	sizeF = 2;
%!	numPts = 100;
%!	bigR = abs(randn);
%!	vecXCent = randn(sizeX,1);
%!	matA = randn(sizeF,sizeX);
%!	%
%!	vecX = randn(sizeX,1);
%!	%
%!	%
%!	[ vecS, vecNHat, vecUHat, matNablaST ] = funcSurf_ellip( vecX, bigR, vecXCent, matA );
%!	%
%!	[ vecS_atS, vecNHat_atS, vecUHat_atS ] = funcSurf_ellip( vecS, bigR, vecXCent, matA );
%!	assert( norm(vecS-vecS_atS) < 1e-8*(norm(vecS)+norm(vecS_atS)) );
%!	assert( norm(vecNHat-vecNHat_atS) < 1e-8*(norm(vecNHat)+norm(vecNHat_atS)) );
%!	assert( norm(vecUHat-vecUHat_atS) < 1e-8*(norm(vecUHat)+norm(vecUHat_atS)) );
%!	%
%!	[ vecSAtPlusU ] = funcSurf_ellip( vecX + vecUHat, bigR, vecXCent, matA );
%!	assert( norm(vecS-vecSAtPlusU) < (eps^0.75)*(norm(vecS)+norm(vecSAtPlusU)) );
%!	vecZ = 1e-4*(2.0*rand(sizeX,1)-1.0);
%!	[ vecSAtPlusZ ] = funcSurf_ellip( vecX + vecZ, bigR, vecXCent, matA );
%!	[ vecSAtMinusZ ] = funcSurf_ellip( vecX - vecZ, bigR, vecXCent, matA );
%!	assert( abs(vecNHat'*(vecSAtPlusZ-vecSAtMinusZ)) < (1e-8)*( norm(vecNHat) + norm(vecSAtPlusZ) + norm(vecSAtMinusZ) ) );
%!	%
%!	[ vecS, vecNHat, vecUHat, matNablaST ] = funcSurf_ellip( vecX, bigR, vecXCent, matA );
%!	[ vecSAtPlusU ] = funcSurf_ellip( vecX + vecUHat, bigR, vecXCent, matA );
%!	assert( norm(vecS-vecSAtPlusU) < (eps^0.75)*(norm(vecS)+norm(vecSAtPlusU)) );
%!	vecZ = 1e-6*(2.0*rand(sizeX,1)-1.0);
%!	[ vecSAtPlusZ ] = funcSurf_ellip( vecX + vecZ, bigR, vecXCent, matA );
%!	[ vecSAtMinusZ ] = funcSurf_ellip( vecX - vecZ, bigR, vecXCent, matA );
%!	%
%!	if ( abs(vecNHat'*(vecSAtPlusZ-vecSAtMinusZ)) > (1e-8)*( norm(vecNHat) + norm(vecSAtPlusZ-vecSAtMinusZ) ) )
%!		msg( thisFile, __LINE__, "" );
%!		msg( thisFile, __LINE__, "Failing surface normal check..." );
%!		msg( thisFile, __LINE__, "  This could happen by chance, becuase the points used are determined sloppily." );
%!		msg( thisFile, __LINE__, "  But, if this happens frequently, there is likely a bug." );
%!		echo__vecNHat = vecNHat
%!		echo__vecSAtPlusZ = vecSAtPlusZ
%!		echo__vecSAtMinusZ = vecSAtMinusZ
%!		echo__absDotProd = abs(vecNHat'*(vecSAtPlusZ-vecSAtMinusZ))
%!		echo__normRHS = (1e-8)*( norm(vecNHat) + norm(vecSAtPlusZ-vecSAtMinusZ) )
%!	end
%!	assert( abs(vecNHat'*(vecSAtPlusZ-vecSAtMinusZ)) < (1e-8)*( norm(vecNHat) + norm(vecSAtPlusZ-vecSAtMinusZ) ) );
