thisFile = "xpost";

msg( thisFile, __LINE__, "Performing post-processing..." );
postTime0 = time();

%ax = 2*[-1,1,-1,1];
contourPlot_numX1Vals = 51;
contourPlot_numX2Vals = 61;
contourPlot_numLevels = 30;
contourPlot_numColors = 64;
contourPlot_axEq = false;
autoAx_borderCoeff = 0.1;

curveColors = 0.6*jet(numCurves);
curveSymbols = [ 'o', 'x', '^', 'v', '*', 's', 'p' ];

autoAx_x1Lo = vecX0(1);
autoAx_x1Hi = vecX0(1);
autoAx_x2Lo = vecX0(2);
autoAx_x2Hi = vecX0(2);

for n=1:numCurves
	matX = curveDat(n).matX;
	%
	assert( sizeX == size(matX,1) );
	matF = funchF(matX);
	%
	numPts = size(matX,2);
	msg( thisFile, __LINE__, sprintf( "Curve %d has %d points.", n, numPts ) );
	curveDat(n).numPts = numPts;
	curveDat(n).matBigDelta = matX - repmat( vecX0, 1, numPts );
	curveDat(n).matF = matF;
	curveDat(n).vecOmega = 0.5*sum( matF.^2, 1 );
	curveDat(n).vecDist = sqrt(sum( curveDat(n).matBigDelta.^2, 1 ));
	%
	curveDat(n).matBigDelta_fromEnd = matX - repmat( matX(:,end), 1, numPts );
	curveDat(n).vecDist_fromEnd = sqrt(sum( curveDat(n).matBigDelta_fromEnd.^2, 1 ));
	%
	curveDat(n).matSmlDelta_cent = matX(:,2:end) - matX(:,1:end-1);
	curveDat(n).matX_cent = ( matX(:,2:end) + matX(:,1:end-1) ) / 2.0;
	curveDat(n).matF_cent = ( matF(:,2:end) + matF(:,1:end-1) ) / 2.0;
	curveDat(n).vecOmega_cent = ( curveDat(n).vecOmega(2:end) + curveDat(n).vecOmega(1:end-1) ) / 2.0;
	curveDat(n).vecStepSize_cent = sqrt(sum( curveDat(n).matSmlDelta_cent.^2, 1 ));
	%
	curveDat(n).plot_color = curveColors(n,:);
	curveDat(n).plot_markerSize = 10;
	curveDat(n).plot_markerStyle = curveSymbols(1+mod(n-1,max(size(curveSymbols))));
	curveDat(n).plot_lineWidth = 2;
	curveDat(n).plot_lineStyle = '-';
	curveDat(n).plot_bigMarkerSize = 30;
	curveDat(n).plot_bigLineWidth = 4;
	%
	autoAx_x1Lo = min([ autoAx_x1Lo, min( matX(1,:) ) ]);
	autoAx_x1Hi = max([ autoAx_x1Hi, max( matX(1,:) ) ]);
	autoAx_x2Lo = min([ autoAx_x2Lo, min( matX(2,:) ) ]);
	autoAx_x2Hi = max([ autoAx_x2Hi, max( matX(2,:) ) ]);
end

if (isempty(ax))
	autoAx_border = autoAx_borderCoeff*max([ autoAx_x1Hi-autoAx_x1Lo, autoAx_x2Hi-autoAx_x2Lo ]);
	ax = [ ...
	  autoAx_x1Lo - autoAx_border, ...
	  autoAx_x1Hi + autoAx_border, ...
	  autoAx_x2Lo - autoAx_border, ...
	  autoAx_x2Hi + autoAx_border ];
end

x1Vals = linspace(ax(1),ax(2),contourPlot_numX1Vals);
x2Vals = linspace(ax(3),ax(4),contourPlot_numX2Vals);
[ x1Mesh, x2Mesh ] = meshgrid( x1Vals, x2Vals );
matX = [ reshape(x1Mesh,1,[]); reshape(x2Mesh,1,[]) ];
matF = funchF(matX,testFuncPrm);
f1Mesh = reshape(matF(1,:),contourPlot_numX2Vals,contourPlot_numX1Vals);
f2Mesh = reshape(matF(2,:),contourPlot_numX2Vals,contourPlot_numX1Vals);
omegaMesh = 0.5*( f1Mesh.^2 + f2Mesh.^2 );
%
msg( thisFile, __LINE__, sprintf( "F1 scale: %g to %g.", min(min((f1Mesh))), max(max((f1Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "F2 scale: %g to %g.", min(min((f2Mesh))), max(max((f2Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "omega scale: %g to %g.", min(min((omegaMesh))), max(max((omegaMesh))) ) );

msg( thisFile, __LINE__, sprintf( "Post processing took %0.3fs.", time()-postTime0 ) );

thisFile = "RETURN FROM xpost";
return;
