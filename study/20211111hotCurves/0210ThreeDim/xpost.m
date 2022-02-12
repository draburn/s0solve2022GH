%msg( __FILE__, __LINE__, "Performing post-processing..." );
postTime0 = time();
%
%
curveSymbols = [ 'o', '^', 'v', 's', 'p' ];
curveColors = 0.7*hsv(numCurves);
vecZ1 = vecXE - vecX0;
vecZ2 = -vecG0;
autoAx_borderCoeff = 0.2;
contourPlot_numZ1Vals = 51;
contourPlot_numZ2Vals = 61;
contourPlot_numLevels = 30;
contourPlot_numColors = 64;
%
%
%
assert( 0.0 ~= norm(vecZ1) );
vecZ1Hat = vecZ1/norm(vecZ1);
vecZ2 = vecZ2 - vecZ1Hat*(vecZ1Hat'*vecZ2);
assert( 0.0 ~= norm(vecZ2) );
vecZ2Hat = vecZ2/norm(vecZ2);
%
autoAx_z1Lo = vecX0(1);
autoAx_z1Hi = vecX0(1);
autoAx_z2Lo = vecX0(2);
autoAx_z2Hi = vecX0(2);
for n=1:numCurves
	curveDat(n).vecFVals = funchF( curveDat(n).vecXVals );
	curveDat(n).omegaVals = 0.5*sumsq(curveDat(n).vecFVals,1);
	curveDat(n).z1Vals = vecZ1Hat'*( curveDat(n).vecXVals - vecX0 );
	curveDat(n).z2Vals = vecZ2Hat'*( curveDat(n).vecXVals - vecX0 );
	curveDat(n).vecZRVals = curveDat(n).vecXVals - vecZ1Hat*curveDat(n).z1Vals - vecZ2Hat*curveDat(n).z2Vals;
	curveDat(n).zdVals = sqrt(sumsq( curveDat(n).vecZRVals, 1));
	%
	%
	curveDat(n).plot_color = curveColors(n,:);
	curveDat(n).plot_markerSize = 8+3*(numCurves-n);
	curveDat(n).plot_markerStyle = curveSymbols(1+mod(n-1,max(size(curveSymbols))));
	curveDat(n).plot_lineWidth = 2;
	curveDat(n).plot_lineStyle = '-';
	curveDat(n).plot_bigMarkerSize = 10+4*(2*numCurves-n);
	curveDat(n).plot_bigLineWidth = 4;
	%
	%
	autoAx_z1Lo = min([ autoAx_z1Lo, min( curveDat(n).z1Vals ) ]);
	autoAx_z1Hi = max([ autoAx_z1Hi, max( curveDat(n).z1Vals ) ]);
	autoAx_z2Lo = min([ autoAx_z2Lo, min( curveDat(n).z2Vals ) ]);
	autoAx_z2Hi = max([ autoAx_z2Hi, max( curveDat(n).z2Vals ) ]);
end
%
%
%
if (isempty(ax))
	autoAx_border = autoAx_borderCoeff*max([ autoAx_z1Hi-autoAx_z1Lo, autoAx_z2Hi-autoAx_z2Lo ]);
	ax = [ ...
	  autoAx_z1Lo - autoAx_border, ...
	  autoAx_z1Hi + autoAx_border, ...
	  autoAx_z2Lo - autoAx_border, ...
	  autoAx_z2Hi + autoAx_border ];
end
%
contourPlot_z1Vals = linspace( ax(1), ax(2), contourPlot_numZ1Vals );
contourPlot_z2Vals = linspace( ax(3), ax(4), contourPlot_numZ2Vals );
[ contourPlot_z1Mesh, contourPlot_z2Mesh ] = meshgrid( contourPlot_z1Vals, contourPlot_z2Vals );
contourPlot_vecXVals = vecX0 + vecZ1Hat*reshape(contourPlot_z1Mesh,1,[]) + vecZ2Hat*reshape(contourPlot_z2Mesh,1,[]);
contourPlot_vecFVals = funchF( contourPlot_vecXVals );
contourPlot_omegaVals = 0.5*sumsq(contourPlot_vecFVals,1);
contourPlot_omegaMesh = reshape( contourPlot_omegaVals, contourPlot_numZ2Vals, contourPlot_numZ1Vals );
%
%
%
msg( __FILE__, __LINE__, sprintf( "Post processing took %0.3fs.", time()-postTime0 ) );
