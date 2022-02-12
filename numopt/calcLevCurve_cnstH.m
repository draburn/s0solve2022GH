% Function ...
%  [ vecXPts, datOut ] = calcLevCurve_cnstH( vecX0, omega0, vecG0, matH, prm=[] )
% Calculates points along the Levenberg curve for the specified quantities.

function [ vecXPts, datOut ] = calcLevCurve_cnstH( vecX0, omega0, vecG0, matH, prm=[] )
	%
	% DRaburn 2022.02.05:
	%  This code is inelegant and hack-ish, though less so than calcGradCurve_cnstH.
	%  Oh well.
	%
	sizeX = size(vecX0,1);
	debugMode = mygetfield( prm, "debugMode", false );
	matS = mygetfield( prm, "matS", eye(sizeX,sizeX) );
	omegaMin = mygetfield( prm, "omegaMin", 0.0 );
	numPts = mygetfield( prm, "numPts", 101 );
	if (debugMode)
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecG0,[sizeX,1]) );
		assert( isrealarray(matH,[sizeX,sizeX]) );
		assert( issymmetric(matH) );
		assert( isrealarray(matS,[sizeX,sizeX]) );
		assert( isrealscalar(omegaMin) );
		assert( isrealscalar(numPts) );
	end
	%
	datOut = [];
	%
	matD = matS'*matS;
	eps075 = eps^0.75;
	%
	% Do main loop.
	vecXPts(:,1) = vecX0;
	for n=2:numPts
		s = (n-1.0)/(numPts-1.0);
		matM = (s*matH) + ((1.0-s)*matD);
		[ matR, cholFlag ] = chol( matM );
		if (0~=cholFlag)
			if (debugMode)
				msg( __FILE__, __LINE__, sprintf( "chol() failed for s = %f.", s ) );
			end
			break;
		end
		vecDelta = matR \ ( matR' \ ( -s*vecG0) );
		%
		omega = omega0 + vecG0'*vecDelta + 0.5*vecDelta'*matH*vecDelta;
		if ( omega < omegaMin )
			if (debugMode)
				msg( __FILE__, __LINE__, sprintf( "omega = %f.", omega ) );
			end
			break;
		end
		%
		vecXPts(:,n) = vecX0 + vecDelta;
	end
	%
	if ( 1+numPts == n )
		% All points were accepted.
		return;
	elseif ( 2 == n )
		% This is bad.
		msg( __FILE__, __LINE__, "Only one point was accepted." );
		return;
	end
	vecXPts = vecXPts(:,1:n-1);
	%
	%
	% DRaburn 2022.02.06:
	%  Extrapolate to omegaMin using last two points.
	vecXA = vecXPts(:,end-1);
	vecXB = vecXPts(:,end);
	vecZ = vecXB-vecXA; % We'll take a step in this direction...
	%
	vecDB = vecXB - vecX0;
	omegaB = omega0 + (vecG0'*vecDB) + 0.5*(vecDB'*matH*vecDB);
	vecGB = vecG0 + matH*vecDB;
	%
	%
	% DRaburn 2022.01.25:
	%  I'm not 100%, given the above, that the s here will always hit omegaMin;
	%  but, when it doesn't, this should still be a meaningful point... Right?
	s = calcLinishRootOfQuad( 0.5*(vecZ'*matH*vecZ), vecGB'*vecZ, omegaB-omegaMin );
	assert( isrealscalar(s) );
	if ( s < 0.0 )
		msg( __FILE__, __LINE__, "Extrapolation goes in wrong direction!" );
	elseif ( s < sqrt(eps) )
		msg( __FILE__, __LINE__, "Extrapolation would be too small!" );
	elseif ( s > 1.0/sqrt(eps) )
		msg( __FILE__, __LINE__, "Extrapolation would be too large!" );
	end
	vecXPts(:,end+1) = vecXB + s*vecZ;
return;
end



%!test
%!	sizeX = 2;
%!	sizeF = 2;
%!	switch (1)
%!	case 0
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = ones(sizeF,1);
%!		matJ0 = eye(sizeF,sizeX);
%!	case 1
%!		setprngstates();
%!		vecX0 = randn(sizeX,1);
%!		vecF0 = randn(sizeF,1);
%!		matJ0 = randn(sizeF,sizeX);
%!	end
%!	omega0 = (vecF0'*vecF0)/2.0;
%!	vecG0 = matJ0'*vecF0;
%!	matH0 = matJ0'*matJ0;
%!	funchF = @(dummyX)( vecF0 + matJ0*(dummyX-vecX0) );
%!	funchOmega = @(dummyX)( sumsq(funchF(dummyX),1)/2.0 );
%!	%
%!	%
%!	prm = [];
%!	prm.debugMode = true;
%!	prm.matS = eye(sizeX,sizeX);
%!	vecXPts = calcLevCurve_cnstH( vecX0, omega0, vecG0, matH0, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	%
%!	x1Lo = min(vecXPts(1,:));
%!	x1Hi = max(vecXPts(1,:));
%!	x2Lo = min(vecXPts(2,:));
%!	x2Hi = max(vecXPts(2,:));
%!	x1Diff = max([ 0.1, x1Hi-x1Lo ]);
%!	x2Diff = max([ 0.1, x2Hi-x2Lo ]);
%!	x1Lo = x1Lo - 0.3*x1Diff;
%!	x1Hi = x1Hi + 0.3*x1Diff;
%!	x2Lo = x2Lo - 0.3*x2Diff;
%!	x2Hi = x2Hi + 0.3*x2Diff;
%!	numX1Vals = 101;
%!	numX2Vals = 81;
%!	%
%!	x1Vals = linspace( x1Lo, x1Hi, numX1Vals );
%!	x2Vals = linspace( x2Lo, x2Hi, numX2Vals );
%!	[ x1Mesh, x2Mesh ] = meshgrid( x1Vals, x2Vals );
%!	vecXVals = [ reshape(x1Mesh,1,[]); reshape(x2Mesh,1,[]) ];
%!	omegaVals = funchOmega( vecXVals );
%!	omegaMesh = reshape( omegaVals, numX2Vals, numX1Vals );
%!	%
%!	contourf( x1Mesh, x2Mesh, sqrt(omegaMesh) );
%!	colormap( mycmap(256) );
%!	hold on;
%!	plot( ...
%!	  vecX0(1), vecX0(2), 'p', 'linewidth', 4, 'markersize', 30, ...
%!	  vecXPts(1,:), vecXPts(2,:), 'o-', 'linewidth', 2, 'markersize', 10, ...
%!	  vecXPts(1,end), vecXPts(2,end), 'x-', 'linewidth', 4, 'markersize', 25 );
%!	hold off;
%!	grid on;
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	title( "sqrt(omega) vs x1, x2" );
