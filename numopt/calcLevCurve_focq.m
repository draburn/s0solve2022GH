% Function ...
%  [ vecXPts, datOut ] = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm=[] )
% Calculates points along the "F One Component Quadratic" Levenberg curve for the specified quantities.
% The model function F is given by...
%  vecF = vecF + matJ0*vecDelta + vecEta*(vecPhiHat'*vecDelta)^2
% where vecDelta = vecX - vecX0.
% For values of "p" for which there are more than one solution,
%  the returned vecX is selected per prm.curveSelector.

function [ vecXPts, datOut ] = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm=[] )
	%
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	debugMode = mygetfield( prm, "debugMode", false );
	curveSelector = mygetfield( prm, "curveSelector", 0 );
	numPts = mygetfield( prm, "numPts", 101 );
	sz = mygetfield( prm, "sz", 1.0 );
	matSY = mygetfield( prm, "matSY", eye(sizeX-1,sizeX-1) );
	if ( debugMode );
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isscalar(debugMode) )
		assert( isbool(debugMode) );
		assert( 2 <= sizeX );
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		assert( isrealarray(vecPhiHat,[sizeX,1]) );
		assert( abs(norm(vecPhiHat)-1.0) < (eps^0.75)*sizeX );
		assert( isrealarray(vecEta,[sizeF,1]) );
		%
		assert( isrealscalar(curveSelector) );
		assert( isrealscalar(sz) );
		assert( isrealarray(matSY,[sizeX-1,sizeX-1]) );
	end
	%
	vecXPts = [];
	datOut = [];
	%
	%
	% DRaburn 2022.01.29:
	%  Contrary to notes, take matPsi to span all of sizeX except vecPhiHat.
	%  Or, rather, assume matJ*matPsi is non-singular.
	matIX = eye(sizeX,sizeX);
	matPsi = orth( matIX - (vecPhiHat*(vecPhiHat')) );
	if ( debugMode )
		assert( isrealarray(matPsi,[sizeX,sizeX-1]) );
		matShy = [ vecPhiHat, matPsi ];
		assert( sum(sum(abs(matShy*(matShy')-matIX))) < (eps^0.50)*sizeX );
		assert( sum(sum(abs((matShy')*matShy-matIX))) < (eps^0.50)*sizeX );
		clear matShy;
	end
	vecLambda = matJ0 * vecPhiHat;
	matW = matJ0 * matPsi;
	szsq = sz*sz;
	matD = matSY'*matSY;
	matWTW = matW'*matW;
	%
	%
	vecXPts = zeros(sizeX,numPts);
	vecXPts(:,1) = vecX0;
	for n=2:numPts
		p = (n-1.0)/(numPts-1.0);
		matM = ((1.0-p)*matD) + (p*matWTW);
		matA = matIX - (p*(matW*(matM\(matW'))));
		%
		c0 = p * ( vecF0' * matA * vecLambda );
		c1 = (1.0-p)*szsq + p*( vecLambda' * matA * vecLambda ) + 2.0*p*( vecF0' * matA * vecEta );
		c2 = 3.0*p*( vecLambda' * matA * vecEta );
		c3 = 2.0*p*( vecEta' * matA * vecEta );
		%
		zRoots = calcCubicRoots( c0, c1, c2, c3 );
		numRoots = size(zRoots,2);
		if (debugMode)
			assert( isrealarray(zRoots,[1,numRoots]) );
		end
		vecZetaRoots = vecF0 + (vecLambda*zRoots) + (vecEta*(zRoots.^2));
		vecYRoots = -p*( matM\(matW'*vecZetaRoots) );
		vecXRoots = vecX0 + (vecPhiHat*zRoots) + (matPsi*vecYRoots);
		switch (numRoots)
		case 1
			vecXPts(:,n) = vecXRoots;
		case 2
			msg( __FILE__, __LINE__, "Case with two roots is not fully supported." );
			vecXPts(:,n) = vecXRoots(:,1);
		case 3
			if ( curveSelector < -0.5 )
				vecXPts(:,n) = vecXRoots(:,1);
			elseif ( curveSelector < 0.5 )
				vecXPts(:,n) = vecXRoots(:,2);
			else
				vecXPts(:,n) = vecXRoots(:,3);
			end
		otherwise
			error("Invalid case.");
		end
	end
	%
	% CONSIDER EXTRAPOLATION HERE,
	%  AS PER OTHER CURVES!
return;
end


%!test
%!	sizeX = 2;
%!	sizeF = 2;
%!	switch (1000)
%!	case 0
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = ones(sizeF,1);
%!		matJ0 = eye(sizeF,sizeX);
%!	case 1
%!		setprngstates();
%!		vecX0 = randn(sizeX,1);
%!		vecF0 = randn(sizeF,1);
%!		matJ0 = randn(sizeF,sizeX);
%!	case 1000
%!		setprngstates(22027120);
%!		vecX0 = randn(sizeX,1);
%!		vecF0 = randn(sizeF,1);
%!		matJ0 = randn(sizeF,sizeX);
%!	end
%!	omega0 = (vecF0'*vecF0)/2.0;
%!	vecG0 = matJ0'*vecF0;
%!	matH0 = matJ0'*matJ0;
%!	[ matPsi0, matLambda0 ] = eig(matH0);
%!	[ lambda0AbsMin0, nOfAbsMin0 ] = min(abs(diag(matLambda0)));
%!	vecPhiHat = matPsi0(:,nOfAbsMin0);
%!	vecEta = randn(sizeF,1);
%!	%
%!	funchF = @(dummyX)( vecF0 + matJ0*(dummyX-vecX0) + vecEta*sumsq(vecPhiHat'*(dummyX-vecX0),1) );
%!	funchOmega = @(dummyX)( sumsq(funchF(dummyX),1)/2.0 );
%!	%
%!	prm = [];
%!	prm.debugMode = true;
%!	prm.matS = eye(sizeX,sizeX);
%!	%
%!	%
%!	prm.curveSelector = -1;
%!	vecXPts= calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	vecXPts_l = vecXPts;
%!	%
%!	prm.curveSelector = 0;
%!	vecXPts= calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	vecXPts_c = vecXPts;
%!	%
%!	prm.curveSelector = +1;
%!	vecXPts= calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	vecXPts_r = vecXPts;
%!	%
%!	%
%!	%
%!	x1Lo = min([ min(vecXPts_l(1,:)), min(vecXPts_c(1,:)), min(vecXPts_r(1,:)) ])
%!	x1Hi = max([ max(vecXPts_l(1,:)), max(vecXPts_c(1,:)), max(vecXPts_r(1,:)) ])
%!	x2Lo = min([ min(vecXPts_l(2,:)), min(vecXPts_c(2,:)), min(vecXPts_r(2,:)) ])
%!	x2Hi = max([ max(vecXPts_l(2,:)), max(vecXPts_c(2,:)), max(vecXPts_r(2,:)) ])
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
%!	  vecX0(1), vecX0(2), 'kp', 'linewidth', 5, 'markersize', 40, ...
%!	  vecXPts_l(1,:), vecXPts_l(2,:), 'rv-', 'linewidth', 2, 'markersize', 10, ...
%!	  vecXPts_l(1,end), vecXPts_l(2,end), 'rv-', 'linewidth', 4, 'markersize', 25, ...
%!	  vecXPts_c(1,:), vecXPts_c(2,:), 'gs-', 'linewidth', 2, 'markersize', 15, ...
%!	  vecXPts_c(1,end), vecXPts_c(2,end), 'gs-', 'linewidth', 4, 'markersize', 30, ...
%!	  vecXPts_r(1,:), vecXPts_r(2,:), 'b^-', 'linewidth', 2, 'markersize', 20, ...
%!	  vecXPts_r(1,end), vecXPts_r(2,end), 'b^-', 'linewidth', 4, 'markersize', 35 );
%!	hold off;
%!	grid on;
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	title( "sqrt(omega) vs x1, x2" );
