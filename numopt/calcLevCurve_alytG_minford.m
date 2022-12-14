% Function...
%  [ vecXPts, datOut ] = calcLevCurve_alytG_minford( vecX0, funchOmegaG, prm=[] )
% Calculates points along the Minford (minimum for given distance) Levenberg curve for funcOmega using FMINUNC().
% funchOmegaG must support the following interface:
%  [ omega, vecG ] = funchOmegaG( vecX )
% Note: That this function uses funcOmegaWithinSurf; alternative approaches to "minford" are possible.

function [ matX, datOut ] = calcLevCurve_alytG_minford( vecX0, funchOmegaG, prm=[] )
	%
	sizeX = size(vecX0,1);
	debugMode = mygetfield( prm, "debugMode", false );
	maxNumPts = mygetfield( prm, "maxNumPts", 10000 );
	targetStepSize = mygetfield( prm, "targetStepSize", 0.1 );
	matS = mygetfield( prm, "matS", [] );
	if (debugMode)
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert(isrealarray(vecX0,[sizeX,1]));
		if (~isempty(matS))
			assert(isrealarray(matS,[sizeX,sizeX]));
		end
	end
	%
	[ omega0, vecG0 ] = funchOmegaG( vecX0 );
	if (debugMode)
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecG0,[sizeX,1]) );
	end
	h0 = norm(vecG0)/targetStepSize;
	%
	vecXC = vecX0;
	datOut = [];
	%
	bigR = 0.0;
	vecX = vecX0;
	matX(:,1) = vecX;
	for n=2:maxNumPts
		bigR += targetStepSize;
		if ( isempty(matS) )
			funchSurf = @(x)( funcSurfEllip( x, vecXC, bigR ) );
		else
			funchSurf = @(x)( funcSurfEllip( x, vecXC, bigR, matS ) );
		end
		%
		% Do work.
		prm_fowis = [];
		prm_fowis.h0 = h0;
		prm_fowis.tau = targetStepSize/10.0;
		prm_fowis.epsX = targetStepSize/100.0;
		fcn = @(dummyX)( funcOmegaWithinSurf( dummyX, funchOmegaG, funchSurf, prm_fowis ) );
		fminunc_opts = optimset( 'GradObj', 'on' );
		vecX = fminunc( fcn, vecX, fminunc_opts );
		if ( debugMode )
			assert( isrealarray(vecX,[sizeX,1]) );
		end
		%
		if ( isempty(matS) )
			s_next = norm(vecX-vecXC);
		else
			s_next = norm(matS*(vecX-vecXC));
		end
		if ( s_next < bigR - (0.5*targetStepSize) )
			%msg( __FILE__, __LINE__, "Reached local min." );
			matX(:,n) = vecX;
			return;
		end
		if ( s_next > bigR )
			if ( debugMode )
			if ( s_next > bigR + targetStepSize )
				msg( __FILE__, __LINE__, "WARNING: Generated vecX is well outside surface. This should never happen." );
			end
			end
			% Pull point back to surface.
			vecX = vecXC + bigR*(vecX-vecXC)/s_next;
		end
		matX(:,n) = vecX;
	end
	%
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
%!	vecXPts = calcLevCurve_alytG_minford( vecX0, funchOmegaG, prm );
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
