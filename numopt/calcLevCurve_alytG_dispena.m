% Function...
%  [ vecXPts, datOut ] = calcLevCurve_alytG_dispena( vecX0, funchOmegaG, prm=[] )
% Calculates points along the Dispena (distance penalty) Levenberg curve for funcOmega using FMINUNC().
% funchOmegaG must support the following interface:
%  [ omega, vecG ] = funchOmegaG( vecX )

function [ vecXPts, datOut ] = calcLevCurve_alytG_dispena( vecX0, funchOmegaG, prm=[] )
	%
	sizeX = size(vecX0,1);
	debugMode = mygetfield( prm, "debugMode", false );
	numPts =  mygetfield( prm, "numPts", 101 );
	matS = mygetfield( prm, "matS", [] );
	if (debugMode)
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealarray(vecX0,[sizeX,1]) );
		[ omegaTest, vecGTest ] = funchOmegaG( vecX0 );
		assert( isrealscalar(omegaTest) );
		assert( isrealarray(vecGTest,[sizeX,1]) );
		if (~isempty(matS) )
			assert( isrealarray(matS,[sizeX,sizeX]) );
		end
	end
	pPts = linspace( 0.0, 1.0, numPts );
	% Force p values to be concentrated near start and end.
	pPts = 1.0 - ((1.0-(pPts.^4)).^4);
	%
	if (isempty(matS))
		matD = [];
	else
		matD = matS'*matS;
	end
	%
	function [ omegaPts, vecNablaOmegaPts ] = funchOmegaLev_sansS( vecXPts, funchOmegaGBase, p, vecXCent )
		if ( 1 == nargout )
			omegaBasePts = funchOmegaGBase( vecXPts );
			omegaPts = p*omegaBasePts + 0.5*(1.0-p)*sumsq(vecXPts-vecXCent,1);
		else
			[ omegaBasePts, vecNablaOmegaBasePts ] = funchOmegaGBase( vecXPts );
			omegaPts = p*omegaBasePts + 0.5*(1.0-p)*sumsq(vecXPts-vecXCent,1);
			vecNablaOmegaPts = p*vecNablaOmegaBasePts + (1.0-p)*(vecXPts-vecX0);
		end
	return;
	end
	%
	function [ omegaPts, vecNablaOmegaPts ] = funchOmegaLev_withS( vecXPts, funchOmegaGBase, p, vecXCent, matD )
		vecDPts = vecXPts - vecXCent;
		if ( 1 == nargout )
			omegaBasePts = funchOmegaGBase( vecXPts );
			omegaPts = p*omegaBasePts + 0.5*(1.0-p)*sum( vecDPts'.*matD*vecDPts, 1 );
		else
			[ omegaBasePts, vecNablaOmegaBasePts ] = funchOmegaGBase( vecXPts );
			omegaPts = p*omegaBasePts + 0.5*(1.0-p)*sum( vecDPts'.*matD*vecDPts, 1 );
			vecNablaOmegaPts = p*vecNablaOmegaBasePts + (1.0-p)*matD*vecDPts;
		end
	return;
	end
	%
	vecX = vecX0;
	n = 1;
	p = pPts(n);
	vecXPts(:,n) = vecX0;
	%
	for n=2:numPts
		p = pPts(n);
		if ( isempty(matD) )
			fcn = @(dummyX)( funchOmegaLev_sansS( dummyX, funchOmegaG, p, vecX0 ) );
		else
			fcn = @(dummyX)( funchOmegaLev_withS( dummyX, funchOmegaG, p, vecX0, matD ) );
		end
		opts = optimset( 'GradObj', 'on' );
		vecX = fminunc( fcn, vecX, opts );
		vecXPts(:,n) = vecX;
	end
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
%!	vecXPts = calcLevCurve_alytG_dispena( vecX0, funchOmegaG, prm );
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
