% Function ...
%  [ vecXPts, datOut ] = calcGradCurve_cnstH( vecX0, omega0, vecG0, matH, prm=[] )
% Calculates points along the gradient-descent curve for the specified quantities.

function [ vecXPts, datOut ] = calcGradCurve_cnstH( vecX0, omega0, vecG0, matH, prm=[] )
	%
	% DRaburn 2022.02.05:
	%  This code is inelegant and hack-ish.
	%  Oh well.
	%
	sizeX = size(vecX0,1);
	debugMode = mygetfield( prm, "debugMode", false );
	omegaMin = mygetfield( prm, "omegaMin", 0.0 );
	numPts = mygetfield( prm, "numPts", 401 );
	if (debugMode)
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecG0,[sizeX,1]) );
		assert( isrealarray(matH,[sizeX,sizeX]) );
		assert( issymmetric(matH) );
		assert( isrealscalar(omegaMin) );
		assert( isrealscalar(numPts) );
	end
	%
	matS = mygetfield( prm, "matS", [] );
	if ( ~isempty(matS) )
		% Draburn 2022.02.06:
		%  This recursive call is clunky, but works.
		%msg( __FILE__, __LINE__, "Calling with scaled values." );
		prmMod = prm;
		prmMod.matS = [];
		if (debugMode)
			assert( isrealarray(matS,[sizeX,sizeX]) );
		end
		matSInv = inv(matS);
		vecGMod = matSInv*vecG0;
		matHMod = matSInv'*matH*matSInv;
		vecXPts = calcGradCurve_cnstH( vecX0, omega0, vecGMod, matHMod, prmMod );
		vecXPts = vecX0 + matSInv*(vecXPts-vecX0);
		return;
	end
	%
	datOut = [];
	%
	% Do main calculations.
	% DRaburn 2022.02.06:
	%  The current code involves Inf and/or NaN calculations for non-positive eigenvalues.
	%  Also, some of the results might not be well-scaled.
	%  But, improving this is not a priority.
	[ matPsi, matLambda ] = eig( matH );
	vecLambda = diag(matLambda);
	%
	vecY0 = matPsi'*vecX0;
	vecGamma = (matPsi'*vecG0) - (matLambda*vecY0);
	%
	lambdaAbsMin = min(abs(vecLambda))+sqrt(eps)*max(abs(vecLambda));
	sPts = (linspace( 1.0, 0.0, numPts )).^(1.0/lambdaAbsMin);
	%
	sizeY = sizeX;
	vecGOL = vecGamma./vecLambda;
	vecYPts = zeros(sizeY,numPts);
	for n=1:sizeY
	if ( 0.0 == vecLambda(n) )
		vecYPts(n,:) = vecY0(n) + vecGamma(n)*log(sPts);
	else
		vecYPts(n,:) = ( vecY0(n) + vecGOL(n) ) * (sPts.^vecLambda(n)) - vecGOL(n);
	end
	end
	vecXPts = matPsi*vecYPts;
	vecDPts = vecXPts - vecX0;
	omegaPts = omega0 + (vecG0'*vecDPts) + 0.5 * sum( vecDPts.*(matH*vecDPts), 1 );
	%
	% DRaburn 2022.02.06:
	%  Now get rid of the points that were Inf or NaN.
	%  The sumsq is particularly tacky, but, whatever.
	ySumSqPts = sumsq( vecYPts, 1 );
	nanPts = ( (isnan(ySumSqPts)) | isinf(ySumSqPts) );
	skipPts = nanPts | (omegaPts<omegaMin);
	if ( 0 == sum(double(skipPts)) )
		return;
	end
	clear vecDPts;
	clear omegaPts;
	%
	%msg( __FILE__, __LINE__, "Discarding invalid points." );
	vecXPts = vecXPts(:,~skipPts);
	%
	numPts = size(vecXPts,2);
	if ( 1 == numPts )
		msg( __FILE__, __LINE__, "Only one point was accepted." );
		return;
	end
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
%!	prm = [];
%!	prm.debugMode = true;
%!	prm.matS = eye(sizeX,sizeX);
%!	vecXPts = calcGradCurve_cnstH( vecX0, omega0, vecG0, matH0, prm );
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
