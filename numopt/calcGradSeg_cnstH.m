% Function ...
%  [ vecXPts, datOut ] = calcGradSeg_cnstH( vecX0, omega0, vecG0, matH, prm=[] )
% Calculates points along the gradient-descent line segment for the specified quantities.

function [ vecXPts, datOut ] = calcGradSeg_cnstH( vecX0, omega0, vecG0, matH, prm=[] )
	%
	sizeX = size(vecX0,1);
	debugMode = mygetfield( prm, "debugMode", false );
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
		assert( isrealscalar(omegaMin) );
		assert( isrealscalar(numPts) );
	end
	%
	matS = mygetfield( prm, "matS", [] );
	if ( ~isempty(matS) )
		% Draburn 2022.02.11:
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
		vecXPts = calcGradSeg_cnstH( vecX0, omega0, vecGMod, matHMod, prmMod );
		vecXPts = vecX0 + matSInv*(vecXPts-vecX0);
		return;
	end
	%
	datOut = [];
	%
	% Do main calculations.
	% Model:
	%  omega = omega0 + vecG0'*vecDelta + 0.5*vecDelta'*matH*vecDelta,
	% with:
	%  vecDelta = -p*vecG0.
	% So:
	%  omega = omega0 + (-vecG0'*vecG0)*p + (0.5*vecG0'*matH*vecG0)*(p^2).
	p = calcLinishRootOfQuad( 0.5*(vecG0'*matH*vecG0), -(vecG0'*vecG0), omega0 );
	vecXPts = vecX0 - vecG0*linspace(0.0,p,numPts);
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
%!	vecXPts = calcGradSeg_cnstH( vecX0, omega0, vecG0, matH0, prm );
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
