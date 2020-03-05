function [ sizeX, sizeY, gridX, gridY, matV, xp, xq, yp, yq ] = spanspace( ...
  vecP, vecQ, prm=[] )
	sizeP = size(vecP,1);
	assert( isrealarray(vecP,[sizeP,1]) );
	assert( isrealarray(vecQ,[sizeP,1]) );
	%
	perpTol = 10.0*eps*(sizeP^2);
	xRangeScaleLo = 1.5;
	xRangeScaleHi = 1.5;
	yRangeScaleLo = 1.5;
	yRangeScaleHi = 1.5;
	sizeX = 50;
	sizeY = 51;
	%
	normP = sqrt(vecP'*vecP);
	assert( 0.0 < normP );
	normQ = sqrt(vecQ'*vecQ);
	assert( 0.0 < normQ );
	vecX = vecP / normP;
	vecT = vecQ - (vecX*(vecX'*vecQ));
	normT = sqrt(vecT'*vecT);
	assert( normT > perpTol*normQ );
	vecY = vecT / normT;
	%
	xp = vecX'*vecP;
	xq = vecX'*vecQ;
	yp = vecY'*vecP; % Should be zero.
	yq = vecY'*vecQ;
	%
	xLo = min([ 0.0, xp, xq ]);
	xHi = max([ 0.0, xp, xq ]);
	yLo = min([ 0.0, yp, yq ]);
	yHi = max([ 0.0, yp, yq ]);
	%
	xVar = (xHi-xLo)/2.0;
	xMid = (xHi+xLo)/2.0;
	yVar = (yHi-yLo)/2.0;
	yMid = (yHi+yLo)/2.0;
	%
	xMin = xMid - (xRangeScaleLo*xVar);
	xMax = xMid + (xRangeScaleLo*xVar);
	yMin = yMid - (yRangeScaleLo*yVar);
	yMax = yMid + (yRangeScaleLo*yVar);
	xVals = linspace( xMin, xMax, sizeX );
	yVals = linspace( yMin, yMax, sizeY );
	%
	[ gridX, gridY ] = meshgrid( xVals, yVals );
	rvecX = reshape( gridX, 1, [] );
	rvecY = reshape( gridY, 1, [] );
	matV = (vecX*rvecX) + (vecY*rvecY);
return;
end

%!test
%!	sizeX = 10;
%!	vecU1 = randn(sizeX,1);
%!	vecU2 = 3*randn(sizeX,1);
%!	[ numS1Vals, numS2Vals, gridS1, gridS2, matX, xp, xq, yp, yq ] = spanspace( vecU1, vecU2 );
%!	numPts = size(matX,2)
%!	for n=1:numPts
%!		rvecZ(n) = sqrt(min([ ...
%!		  10*sum(matX(:,n).^2), ...
%!		  sum((matX(:,n)-vecU1).^2), ...
%!		  3*sum((matX(:,n)-vecU2).^2) ]));
%!	end
%!	gridZ = reshape( rvecZ, numS2Vals, numS1Vals );
%!	contour( gridS1, gridS2, gridZ );
%!	hold on;
%!	plot( [0,xp], [0,yp], 'ko-' );
%!	plot( [0,xq], [0,yq], 'kx-' );
%!	grid on;
%!	%axis equal;
%!	hold off;
%!	xlabel( "x" );
%!	ylabel( "y" );
