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
vecZ3Hat = zeros(sizeX,1); % Unless...
if ( 3 <= sizeX )
for n=1:3
	vecZ3 = zeros(sizeX,1);
	vecZ3(n) = 1.0;
	vecZ3 = vecZ3 - vecZ1Hat*(vecZ1Hat'*vecZ3);
	vecZ3 = vecZ3 - vecZ2Hat*(vecZ2Hat'*vecZ3);
	if ( 0.0 == norm(vecZ3) )
		continue;
	end
	vecZ3Hat = vecZ3/norm(vecZ3);
	vecZ3 = vecZ3 - vecZ1Hat*(vecZ1Hat'*vecZ3);
	vecZ3 = vecZ3 - vecZ2Hat*(vecZ2Hat'*vecZ3);
	if ( 0.0 == norm(vecZ3) )
		continue;
	end
	vecZ3Hat = vecZ3/norm(vecZ3);
end
end
%
autoAx_z1Lo = 0.0;
autoAx_z1Hi = 0.0;
autoAx_z2Lo = 0.0;
autoAx_z2Hi = 0.0;
for n=1:numCurves
	curveDat(n).vecFVals = funchF( curveDat(n).vecXVals );
	curveDat(n).omegaVals = 0.5*sumsq(curveDat(n).vecFVals,1);
	curveDat(n).vecBigDeltaVals = curveDat(n).vecXVals - vecX0;
	curveDat(n).z1Vals = vecZ1Hat'*(curveDat(n).vecBigDeltaVals);
	curveDat(n).z2Vals = vecZ2Hat'*(curveDat(n).vecBigDeltaVals);
	curveDat(n).z3Vals = vecZ3Hat'*(curveDat(n).vecBigDeltaVals);
	curveDat(n).vecZRVals = curveDat(n).vecBigDeltaVals - vecZ1Hat*curveDat(n).z1Vals - vecZ2Hat*curveDat(n).z2Vals;
	curveDat(n).zdVals = sqrt(sumsq( curveDat(n).vecZRVals, 1));
	curveDat(n).vecFVals_cnstJ = vecF0 + matJ0*curveDat(n).vecBigDeltaVals;
	curveDat(n).modelInaccuracyVals_cnstJ = sqrt(sumsq( curveDat(n).vecFVals_cnstJ - curveDat(n).vecFVals, 1 ));
	curveDat(n).omegaVals_cnstJ = 0.5*sumsq(curveDat(n).vecFVals_cnstJ,1);
	curveDat(n).bigDeltaNormVals = sqrt(sumsq(curveDat(n).vecBigDeltaVals,1));
	curveDat(n).dacVals = [ 0.0, cumsum(sqrt(sumsq(diff(curveDat(n).vecXVals,1,2),1))) ];
	%
	%
	curveDat(n).plot_color = curveColors(n,:);
	curveDat(n).plot_markerSize = 4+3*(numCurves-n);
	curveDat(n).plot_markerStyle = curveSymbols(1+mod(n-1,max(size(curveSymbols))));
	curveDat(n).plot_lineWidth = 2;
	curveDat(n).plot_lineStyle = '-';
	curveDat(n).plot_bigMarkerSize = 4+3*(2*numCurves-n);
	curveDat(n).plot_bigLineWidth = 2;
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
contourPlot_vecFVals_cnstJ = vecF0 + matJ0*contourPlot_vecXVals;
contourPlot_omegaVals_cnstJ = 0.5*sumsq(contourPlot_vecFVals_cnstJ,1);
contourPlot_omegaMesh_cnstJ = reshape( contourPlot_omegaVals_cnstJ, contourPlot_numZ2Vals, contourPlot_numZ1Vals );
%
%
%
msg( __FILE__, __LINE__, sprintf( "Post processing took %0.3fs.", time()-postTime0 ) );
