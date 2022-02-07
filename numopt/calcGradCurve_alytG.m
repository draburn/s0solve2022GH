% Function...
%  [ vecXPts, datOut ] = calcGradCurve_alytG( vecX0, funchOmegaG, prm=[] )
% Calculates points along the gradient-descent curve for funcOmega using LSODE().
% funchOmegaG must support the following interface:
%  [ omega, vecG ] = funchOmegaG( vecX )

function [ vecXPts, datOut ] = calcGradCurve_alytG( vecX0, funchOmegaG, prm=[] )
	%
	debugMode = mygetfield( prm, "debugMode", false );
	sizeX = size(vecX0,1);
	if (debugMode)
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealarray(vecX0,[sizeX,1]) );
	end
	%
	[ omega0, vecG0 ] = funchOmegaG( vecX0 );
	if (debugMode)
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecG0,[sizeX,1]) );
	end
	normG0 = norm(vecG0);
	assert( 0.0 ~= normG0 );
	%
	targetStepSize = mygetfield( prm, "targetStepSize", 0.1 );
	numPtsPerChunk = mygetfield( prm, "numPtsPerChunk", 101 );
	maxNumChunks = mygetfield( prm, "maxNumChunks", 100 );
	deltaTPerChunk = mygetfield( prm, "deltaTPerChunk", targetStepSize*(numPtsPerChunk-1) );
	deltaNormTol = mygetfield( prm, "deltaNormTol", targetStepSize/100.0 );
	gScale = mygetfield( prm, "gScale", (eps^0.25)*normG0 );
	matS = mygetfield( prm, "matS", [] );
	if (debugMode)
		if (~isempty(matS))
			assert( isrealarray(matS,[sizeX,sizeX]) );
		end
	end
	%
	datOut = [];
	%
	% We only actually care about funchXDot_lsode(), not funchOmegaG() per se.
	function vecG = funchG( dummyX )
		[ omega, vecG ] = funchOmegaG(dummyX);
	end
	funchXDotOfG = @(g)( -g/sqrt( gScale^2 + (g'*g) ) );
	if ( isempty(matS) )
		funchXDot_lsode = @(x,t)( funchXDotOfG(funchG(x)) );
	else
		matD = matS'*matS;
		funchXDot_lsode = @(x,t)( funchXDotOfG(matD\funchG(x)) );
	end
	vecT = deltaTPerChunk*linspace(0.0,1.0,numPtsPerChunk);
	%
	% Main loop.
	vecX = vecX0;
	vecXPts(:,1) = vecX;
	numPts = 1;
	for n=1:maxNumChunks
		vecXPts_new = lsode( funchXDot_lsode, vecX, vecT )';
		numPts_new = size(vecXPts_new,2);
		assert( isrealarray(vecXPts_new,[sizeX,numPts_new]) );
		vecXPts(:,numPts+1:numPts+numPts_new-1) = vecXPts_new(:,2:end);
		numPts += numPts_new-1;
		vecX = vecXPts_new(:,end);
		if ( norm(vecXPts_new(:,end)-vecXPts_new(:,end-1)) <= deltaNormTol )
			break;
		end
	end
	%
	% Drop some points from last chunk.
	vecXF = vecXPts(:,end);
	keepPts = ( sumsq(vecXPts-vecXF,1) >= (deltaNormTol^2) );
	vecXPts = vecXPts(:,keepPts);
	vecXPts(:,end+1) = vecXF;
return;
end


%!function [ omega, vecG ] = funcOmegaG( vecX, vecX0, omega0, vecG0, matH0 )
%!	vecY = vecX-vecX0;
%!	omega = omega0 + vecG0'*vecY + 0.5*vecY'*matH0*vecY;
%!	if ( 2 >= nargout )
%!		vecG = vecG0 + matH0*vecY;
%!	end
%!endfunction


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
%!	%
%!	funchF = @(dummyX)( vecF0 + matJ0*(dummyX-vecX0) );
%!	funchOmega = @(dummyX)( sumsq(funchF(dummyX),1)/2.0 );
%!	funchOmegaG = @(dummyX) funcOmegaG( dummyX, vecX0, omega0, vecG0, matH0 );
%!	%
%!	prm = [];
%!	prm.debugMode = true;
%!	prm.matS = eye(sizeX,sizeX);
%!	vecXPts = calcGradCurve_alytG( vecX0, funchOmegaG, prm );
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
