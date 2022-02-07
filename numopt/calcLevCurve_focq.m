% Function ...
%  [ vecXPts, datOut ] = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm=[] )
% Calculates points along the "F One Component Quadratic" Levenberg curve for the specified quantities.
% The model function F is given by...
%  vecF = vecF0 + matJ0*vecDelta + vecEta*(vecPhiHat'*vecDelta)^2
% where vecDelta = vecX - vecX0. (Note the absence of "0.5*" in front of vecEta.)
% For values of "p" for which there are more than one solution,
%  the returned vecX is selected per prm.curveSelector.

function [ vecXPts, datOut ] = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm=[] )
	%
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	debugMode = mygetfield( prm, "debugMode", false );
	sz = mygetfield( prm, "sz", 1.0 );
	matSY = mygetfield( prm, "matSY", eye(sizeX-1,sizeX-1) );
	curveSelector = mygetfield( prm, "curveSelector", 2 );
	pPts = mygetfield( prm, "pPts", [] );
	numPts = mygetfield( prm, "numPts", 101 );
	if (isempty(pPts))
		pPts = linspace( 0.0, 1.0, numPts );
	else
		numPts = size(pPts,2);
	end
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
		assert( isrealscalar(sz) );
		assert( isrealarray(matSY,[sizeX-1,sizeX-1]) );
		%
		assert( isrealscalar(curveSelector) );
		assert( isrealarray(pPts,[1,numPts]) );
	end
	%
	vecXPts = [];
	if ( nargout >= 2 )
		datOut = [];
		datOut.multiRootCounter = 0;
	end
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
	for n=1:numPts
		p = pPts(n);
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
			switch (curveSelector)
			case -2
				vecDeltaRoots = vecXRoots-vecX0;
				vecFRoots = vecF0 + matJ0*vecDeltaRoots + vecEta*((vecPhiHat'*vecDeltaRoots).^2);
				xiRoots = p*sumsq(vecFRoots,1) + (1.0-p)*( (sz*zRoots).^2 + sumsq(matSY*vecYRoots,1) );
				if ( xiRoots(1) < xiRoots(3) )
					vecXPts(:,n) = vecXRoots(:,1);
				else
					vecXPts(:,n) = vecXRoots(:,3);
				end
				% If we're at the roots (1==p), then both solutions should have 0 == xi;
				% but, finite precision issues may mess things up.
				% So, for the last point, force us to stay with the current root,
				% as per 2 == curveSelector.
				if ( abs(p-1.0)<sqrt(eps) )
					vecXPrev = vecXPts(:,n-1);
					if ( norm(vecXRoots(:,3)-vecXPrev) < norm(vecXRoots(:,1)-vecXPrev) )
						vecXPts(:,n) = vecXRoots(:,3);
					else
						vecXPts(:,n) = vecXRoots(:,1);
					end
				end
			case -1
				vecXPts(:,n) = vecXRoots(:,1);
			case 0
				vecXPts(:,n) = vecXRoots(:,2);
			case 1
				vecXPts(:,n) = vecXRoots(:,3);
			case 2
				vecXPrev = vecXPts(:,n-1);
				if ( norm(vecXRoots(:,3)-vecXPrev) < norm(vecXRoots(:,1)-vecXPrev) )
					vecXPts(:,n) = vecXRoots(:,3);
				else
					vecXPts(:,n) = vecXRoots(:,1);
				end
			otherwise
				error("Invalid value of curveSelector.");
			end
			if ( nargout >= 2 )
				datOut.multiRootCounter++;
				datOut.multiRoot(datOut.multiRootCounter).n = n;
				datOut.multiRoot(datOut.multiRootCounter).p = p;
				datOut.multiRoot(datOut.multiRootCounter).vecXRoots = vecXRoots;
			end
		otherwise
			error("Number of roots is not 1, 2, or 3. This should be impossible.");
		end
	end
	%
	% CONSIDER EXTRAPOLATION HERE,
	%  AS PER OTHER CURVES!
return;
end


%!test
%!	numFigs = 0;
%!	sizeX = 2;
%!	sizeF = 2;
%!	switch (10011)
%!	case 0
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = ones(sizeF,1);
%!		matJ0 = eye(sizeF,sizeX);
%!		vecEta = ones(sizeX,1);
%!	case 1
%!		setprngstates();
%!		vecX0 = randn(sizeX,1);
%!		vecF0 = randn(sizeF,1);
%!		matJ0 = randn(sizeF,sizeX);
%!		vecEta = randn(sizeF,1);
%!	case 10010
%!		% "1->3->1 detour".
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = [ -2.0; 1.0 ];
%!		matJ0 = [ 1.0, 0.0; 0.0, 0.0 ];
%!		vecEta = [ 2.0; 1.0 ];
%!	case 10011
%!		% "1->3->1 disconnected loop".
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = [ -2.0; 1.0 ];
%!		matJ0 = [ 1.0, 0.0; 0.0, 0.02 ];
%!		vecEta = [ 2.0; 1.0 ];
%!	case 10020
%!		% "1->3->1 cts -> 3"
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = [ -100.0; 1.0 ];
%!		matJ0 = [ 3.1, 0.0; 0.0, -3.0 ];
%!		vecEta = [ 0.1; 1.0 ];
%!	case 10050
%!		setprngstates(22027120);
%!		vecX0 = randn(sizeX,1);
%!		vecF0 = randn(sizeF,1);
%!		matJ0 = randn(sizeF,sizeX);
%!		vecEta = randn(sizeF,1);
%!	end
%!	omega0 = (vecF0'*vecF0)/2.0;
%!	vecG0 = matJ0'*vecF0;
%!	matH0 = matJ0'*matJ0;
%!	% We'll take vecPhiHat to be the eigenvec for the absmin eigval of matH0.
%!	[ matPsi0, matLambda0 ] = eig(matH0);
%!	[ lambda0AbsMin0, nOfAbsMin0 ] = min(abs(diag(matLambda0)));
%!	vecPhiHat = matPsi0(:,nOfAbsMin0);
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
%!	vecXPts = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	vecXPts_l = vecXPts;
%!	%
%!	prm.curveSelector = 0;
%!	vecXPts = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	vecXPts_c = vecXPts;
%!	%
%!	prm.curveSelector = +1;
%!	vecXPts = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	vecXPts_r = vecXPts;
%!	%
%!	prm.curveSelector = +2;
%!	vecXPts = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	vecXPts_cts = vecXPts;
%!	%
%!	%
%!	%
%!	x1Lo = min([ min(vecXPts_l(1,:)), min(vecXPts_c(1,:)), min(vecXPts_r(1,:)) ]);
%!	x1Hi = max([ max(vecXPts_l(1,:)), max(vecXPts_c(1,:)), max(vecXPts_r(1,:)) ]);
%!	x2Lo = min([ min(vecXPts_l(2,:)), min(vecXPts_c(2,:)), min(vecXPts_r(2,:)) ]);
%!	x2Hi = max([ max(vecXPts_l(2,:)), max(vecXPts_c(2,:)), max(vecXPts_r(2,:)) ]);
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
%!	numFigs++; figure(numFigs);
%!	contourf( x1Mesh, x2Mesh, sqrt(omegaMesh) );
%!	colormap( mycmap(256) );
%!	hold on;
%!	plot( ...
%!	  vecX0(1), vecX0(2), 'kp', 'linewidth', 5, 'markersize', 60, ...
%!	  vecXPts_l(1,:), vecXPts_l(2,:), 'rv-', 'linewidth', 2, 'markersize', 20, ...
%!	  vecXPts_l(1,end), vecXPts_l(2,end), 'rv-', 'linewidth', 4, 'markersize', 52, ...
%!	  vecXPts_c(1,:), vecXPts_c(2,:), 'gs-', 'linewidth', 2, 'markersize', 15, ...
%!	  vecXPts_c(1,end), vecXPts_c(2,end), 'gs-', 'linewidth', 4, 'markersize', 44, ...
%!	  vecXPts_r(1,:), vecXPts_r(2,:), 'b^-', 'linewidth', 2, 'markersize', 10, ...
%!	  vecXPts_r(1,end), vecXPts_r(2,end), 'b^-', 'linewidth', 4, 'markersize', 36, ...
%!	  vecXPts_cts(1,:), vecXPts_cts(2,:), 'ko-', 'linewidth', 2, 'markersize', 4, ...
%!	  vecXPts_cts(1,end), vecXPts_cts(2,end), 'ko-', 'linewidth', 4, 'markersize', 28 );
%!	hold off;
%!	grid on;
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	title( "sqrt(omega) vs x1, x2 -- all 4 curves" );


%!test
%!	numFigs = 1;
%!	sizeX = 2;
%!	sizeF = 2;
%!	switch (10020)
%!	case 0
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = ones(sizeF,1);
%!		matJ0 = eye(sizeF,sizeX);
%!		vecEta = ones(sizeX,1);
%!	case 1
%!		setprngstates();
%!		vecX0 = randn(sizeX,1);
%!		vecF0 = randn(sizeF,1);
%!		matJ0 = randn(sizeF,sizeX);
%!		vecEta = randn(sizeF,1);
%!	case 10010
%!		% "1->3->1 detour".
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = [ -2.0; 1.0 ];
%!		matJ0 = [ 1.0, 0.0; 0.0, 0.0 ];
%!		vecEta = [ 2.0; 1.0 ];
%!	case 10011
%!		% "1->3->1 disconnected loop".
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = [ -2.0; 1.0 ];
%!		matJ0 = [ 1.0, 0.0; 0.0, 0.02 ];
%!		vecEta = [ 2.0; 1.0 ];
%!	case 10020
%!		% "1->3->1 cts -> 3"
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = [ -100.0; 1.0 ];
%!		matJ0 = [ 3.1, 0.0; 0.0, -3.0 ];
%!		vecEta = [ 0.1; 1.0 ];
%!	case 10050
%!		setprngstates(22027120);
%!		vecX0 = randn(sizeX,1);
%!		vecF0 = randn(sizeF,1);
%!		matJ0 = randn(sizeF,sizeX);
%!		vecEta = randn(sizeF,1);
%!	end
%!	omega0 = (vecF0'*vecF0)/2.0;
%!	vecG0 = matJ0'*vecF0;
%!	matH0 = matJ0'*matJ0;
%!	% We'll take vecPhiHat to be the eigenvec for the absmin eigval of matH0.
%!	[ matPsi0, matLambda0 ] = eig(matH0);
%!	[ lambda0AbsMin0, nOfAbsMin0 ] = min(abs(diag(matLambda0)));
%!	vecPhiHat = matPsi0(:,nOfAbsMin0);
%!	%
%!	funchF = @(dummyX)( vecF0 + matJ0*(dummyX-vecX0) + vecEta*sumsq(vecPhiHat'*(dummyX-vecX0),1) );
%!	funchOmega = @(dummyX)( sumsq(funchF(dummyX),1)/2.0 );
%!	%
%!	prm = [];
%!	%prm.debugMode = true;
%!	prm.pPts = linspace(0.0,1.0,201).^2;
%!	%
%!	[ vecXPts, datOut ] = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	%
%!	vecXPts_l = zeros(sizeX,datOut.multiRootCounter);
%!	vecXPts_c = zeros(sizeX,datOut.multiRootCounter);
%!	vecXPts_r = zeros(sizeX,datOut.multiRootCounter);
%!	for n=1:datOut.multiRootCounter
%!		vecXPts_l(:,n) = datOut.multiRoot(n).vecXRoots(:,1);
%!		vecXPts_c(:,n) = datOut.multiRoot(n).vecXRoots(:,2);
%!		vecXPts_r(:,n) = datOut.multiRoot(n).vecXRoots(:,3);
%!	end
%!	%
%!	%
%!	x1Lo = min([ min(vecXPts_l(1,:)), min(vecXPts_c(1,:)), min(vecXPts_r(1,:)), min(vecXPts(1,:)) ]);
%!	x1Hi = max([ max(vecXPts_l(1,:)), max(vecXPts_c(1,:)), max(vecXPts_r(1,:)), max(vecXPts(1,:)) ]);
%!	x2Lo = min([ min(vecXPts_l(2,:)), min(vecXPts_c(2,:)), min(vecXPts_r(2,:)), min(vecXPts(2,:)) ]);
%!	x2Hi = max([ max(vecXPts_l(2,:)), max(vecXPts_c(2,:)), max(vecXPts_r(2,:)), max(vecXPts(2,:)) ]);
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
%!	numFigs++; figure(numFigs);
%!	contourf( x1Mesh, x2Mesh, sqrt(omegaMesh) );
%!	colormap( mycmap(256) );
%!	hold on;
%!	plot( ...
%!	  vecX0(1), vecX0(2), 'kp', 'linewidth', 5, 'markersize', 50, ...
%!	  vecXPts(1,:), vecXPts(2,:), 'ko-', 'linewidth', 2, 'markersize', 20, ...
%!	  vecXPts(1,end), vecXPts(2,end), 'ko-', 'linewidth', 4, 'markersize', 30, ...
%!	  vecXPts_l(1,:), vecXPts_l(2,:), 'rv', 'linewidth', 2, 'markersize', 15, ...
%!	  vecXPts_c(1,:), vecXPts_c(2,:), 'gs', 'linewidth', 2, 'markersize', 10, ...
%!	  vecXPts_r(1,:), vecXPts_r(2,:), 'b^', 'linewidth', 2, 'markersize', 5 );
%!	hold off;
%!	grid on;
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	title( "sqrt(omega) vs x1, x2 -- cts curve + pts" );



%!test
%!	numFigs = 2;
%!	sizeX = 2;
%!	sizeF = 2;
%!	switch (10020)
%!	case 0
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = ones(sizeF,1);
%!		matJ0 = eye(sizeF,sizeX);
%!		vecEta = ones(sizeX,1);
%!	case 1
%!		setprngstates();
%!		vecX0 = randn(sizeX,1);
%!		vecF0 = randn(sizeF,1);
%!		matJ0 = randn(sizeF,sizeX);
%!		vecEta = randn(sizeF,1);
%!	case 10010
%!		% "1->3->1 detour".
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = [ -2.0; 1.0 ];
%!		matJ0 = [ 1.0, 0.0; 0.0, 0.0 ];
%!		vecEta = [ 2.0; 1.0 ];
%!	case 10011
%!		% "1->3->1 disconnected loop".
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = [ -2.0; 1.0 ];
%!		matJ0 = [ 1.0, 0.0; 0.0, 0.02 ];
%!		vecEta = [ 2.0; 1.0 ];
%!	case 10020
%!		% "1->3->1 cts -> 3"
%!		vecX0 = zeros(sizeX,1);
%!		vecF0 = [ -100.0; 1.0 ];
%!		matJ0 = [ 3.1, 0.0; 0.0, -3.0 ];
%!		vecEta = [ 0.1; 1.0 ];
%!	case 10050
%!		setprngstates(22027120);
%!		vecX0 = randn(sizeX,1);
%!		vecF0 = randn(sizeF,1);
%!		matJ0 = randn(sizeF,sizeX);
%!		vecEta = randn(sizeF,1);
%!	end
%!	omega0 = (vecF0'*vecF0)/2.0;
%!	vecG0 = matJ0'*vecF0;
%!	matH0 = matJ0'*matJ0;
%!	% We'll take vecPhiHat to be the eigenvec for the absmin eigval of matH0.
%!	[ matPsi0, matLambda0 ] = eig(matH0);
%!	[ lambda0AbsMin0, nOfAbsMin0 ] = min(abs(diag(matLambda0)));
%!	vecPhiHat = matPsi0(:,nOfAbsMin0);
%!	%
%!	funchF = @(dummyX)( vecF0 + matJ0*(dummyX-vecX0) + vecEta*sumsq(vecPhiHat'*(dummyX-vecX0),1) );
%!	funchOmega = @(dummyX)( sumsq(funchF(dummyX),1)/2.0 );
%!	%
%!	prm = [];
%!	%prm.debugMode = true;
%!	prm.pPts = linspace(0.0,1.0,201).^2;
%!	prm.curveSelector = -2;
%!	%
%!	[ vecXPts, datOut ] = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, prm );
%!	numPts = size(vecXPts,2);
%!	assert( isrealarray(vecXPts,[sizeX,numPts]) );
%!	%
%!	vecXPts_l = zeros(sizeX,datOut.multiRootCounter);
%!	vecXPts_c = zeros(sizeX,datOut.multiRootCounter);
%!	vecXPts_r = zeros(sizeX,datOut.multiRootCounter);
%!	for n=1:datOut.multiRootCounter
%!		vecXPts_l(:,n) = datOut.multiRoot(n).vecXRoots(:,1);
%!		vecXPts_c(:,n) = datOut.multiRoot(n).vecXRoots(:,2);
%!		vecXPts_r(:,n) = datOut.multiRoot(n).vecXRoots(:,3);
%!	end
%!	%
%!	%
%!	x1Lo = min([ min(vecXPts_l(1,:)), min(vecXPts_c(1,:)), min(vecXPts_r(1,:)), min(vecXPts(1,:)) ]);
%!	x1Hi = max([ max(vecXPts_l(1,:)), max(vecXPts_c(1,:)), max(vecXPts_r(1,:)), max(vecXPts(1,:)) ]);
%!	x2Lo = min([ min(vecXPts_l(2,:)), min(vecXPts_c(2,:)), min(vecXPts_r(2,:)), min(vecXPts(2,:)) ]);
%!	x2Hi = max([ max(vecXPts_l(2,:)), max(vecXPts_c(2,:)), max(vecXPts_r(2,:)), max(vecXPts(2,:)) ]);
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
%!	numFigs++; figure(numFigs);
%!	contourf( x1Mesh, x2Mesh, sqrt(omegaMesh) );
%!	colormap( mycmap(256) );
%!	hold on;
%!	plot( ...
%!	  vecX0(1), vecX0(2), 'kp', 'linewidth', 5, 'markersize', 50, ...
%!	  vecXPts(1,:), vecXPts(2,:), 'ko-', 'linewidth', 2, 'markersize', 20, ...
%!	  vecXPts(1,end), vecXPts(2,end), 'ko-', 'linewidth', 4, 'markersize', 30, ...
%!	  vecXPts_l(1,:), vecXPts_l(2,:), 'rv', 'linewidth', 2, 'markersize', 15, ...
%!	  vecXPts_c(1,:), vecXPts_c(2,:), 'gs', 'linewidth', 2, 'markersize', 10, ...
%!	  vecXPts_r(1,:), vecXPts_r(2,:), 'b^', 'linewidth', 2, 'markersize', 5 );
%!	hold off;
%!	grid on;
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	title( "sqrt(omega) vs x1, x2 -- xi min curve + pts" );
