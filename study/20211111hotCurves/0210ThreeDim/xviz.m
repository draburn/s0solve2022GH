numFigs = 0;
vizTime0 = time();
%msg( __FILE__, __LINE__, "Drawing graphs..." );
%
useAxEq = false;
plotLegend = true;
plotCurveZ1Z2 = true;
plotCurveZ1Z2_cnstJ = false;
plotCurveZ1ZD = true;
plotOmegaVsZ1 = true;
plotOmegaVsZ1_cnstJ = true;
plotOmegaVsN = true;
plotZ1VsDAC = true;
%
%
%
if ( plotLegend )
	numFigs++; figure(numFigs);
	n = 1;
	cellAry_legend(1,n) = { curveDat(n).strName };
	plot( ...
	  curveDat(n).z1Vals, curveDat(n).z2Vals, ...
	  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
	  'linewidth', curveDat(n).plot_lineWidth, ...
	  'color', curveDat(n).plot_color, ...
	  'markersize', curveDat(n).plot_markerSize );
	hold on;
	for n=2:numCurves
		cellAry_legend(1,n) = { curveDat(n).strName };
		plot( ...
		  curveDat(n).z1Vals, curveDat(n).z2Vals, ...
		  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
		  'linewidth', curveDat(n).plot_lineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_markerSize );
	end
	legObj = legend( cellAry_legend, "location", "northwestoutside" );
	set( legObj, 'Interpreter', 'none' );
	for n=1:numCurves
		plot( ...
		  curveDat(n).z1Vals([1,end]), curveDat(n).z2Vals([1,end]), ...
		  curveDat(n).plot_markerStyle, ...
		  'linewidth', curveDat(n).plot_bigLineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_bigMarkerSize );
	end
	plot( [0.0,norm(vecXE-vecX0)], [0.0,0.0], 'k-' );
	hold off;
	title( "Legend" );
	grid on;
end
%
%
%
if ( plotCurveZ1Z2 )
	numFigs++; figure(numFigs);
	funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
	contourf( contourPlot_z1Mesh, contourPlot_z2Mesh, sqrt(contourPlot_omegaMesh), contourPlot_numLevels );
	colormap( mycmap(contourPlot_numColors) );
	xlabel( "z1" );
	ylabel( "z2" );
	title( "sqrt(omega) vs ( z1, z2 )" );
	hold on;
	plot( [0.0,norm(vecXE-vecX0)], [0.0,0.0], 'k-' );
	if ( useAxEq )
		axis equal;
	end
	for n=1:numCurves
		plot( ...
		  curveDat(n).z1Vals, curveDat(n).z2Vals, ...
		  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
		  'linewidth', curveDat(n).plot_lineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_markerSize );
		plot( ...
		  curveDat(n).z1Vals([1,end]), curveDat(n).z2Vals([1,end]), ...
		  curveDat(n).plot_markerStyle, ...
		  'linewidth', curveDat(n).plot_bigLineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_bigMarkerSize );
	end
	hold off;
	grid on;
end
if ( plotCurveZ1Z2_cnstJ )
	numFigs++; figure(numFigs);
	funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
	contourf( contourPlot_z1Mesh, contourPlot_z2Mesh, sqrt(contourPlot_omegaMesh_cnstJ), contourPlot_numLevels );
	colormap( mycmap(contourPlot_numColors) );
	xlabel( "z1" );
	ylabel( "z2" );
	title( "sqrt(omega cnstJ) vs ( z1, z2 )" );
	hold on;
	plot( [0.0,norm(vecXE-vecX0)], [0.0,0.0], 'k-' );
	if ( useAxEq )
		axis equal;
	end
	for n=1:numCurves
		plot( ...
		  curveDat(n).z1Vals, curveDat(n).z2Vals, ...
		  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
		  'linewidth', curveDat(n).plot_lineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_markerSize );
		plot( ...
		  curveDat(n).z1Vals([1,end]), curveDat(n).z2Vals([1,end]), ...
		  curveDat(n).plot_markerStyle, ...
		  'linewidth', curveDat(n).plot_bigLineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_bigMarkerSize );
	end
	hold off;
	grid on;
end
%
%
%
if ( plotCurveZ1ZD )
	numFigs++; figure(numFigs);
	funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
	plot( [0.0,norm(vecXE-vecX0)], [0.0,0.0], 'k-' );
	xlabel( "z1" );
	ylabel( "zd" );
	title( "curve in (z1,zd) space" );
	if ( useAxEq )
		axis equal;
	end
	hold on;
	for n=1:numCurves
		plot( ...
		  curveDat(n).z1Vals, curveDat(n).zdVals, ...
		  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
		  'linewidth', curveDat(n).plot_lineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_markerSize );
		plot( ...
		  curveDat(n).z1Vals([1,end]), curveDat(n).zdVals([1,end]), ...
		  curveDat(n).plot_markerStyle, ...
		  'linewidth', curveDat(n).plot_bigLineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_bigMarkerSize );
	end
	hold off;
	grid on;
end
%
%
%
if ( plotOmegaVsZ1 )
	numFigs++; figure(numFigs);
	funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
	plot( [0.0,norm(vecXE-vecX0)], [0.0,0.0], 'k-' );
	xlabel( "z1" );
	ylabel( "sqrt(omega)" );
	title( "sqrt(omega) vs z1" );
	if ( useAxEq )
		axis equal;
	end
	hold on;
	for n=1:numCurves
		plot( ...
		  curveDat(n).z1Vals, sqrt(curveDat(n).omegaVals), ...
		  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
		  'linewidth', curveDat(n).plot_lineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_markerSize );
		plot( ...
		  curveDat(n).z1Vals([1,end]), sqrt(curveDat(n).omegaVals([1,end])), ...
		  curveDat(n).plot_markerStyle, ...
		  'linewidth', curveDat(n).plot_bigLineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_bigMarkerSize );
	end
	hold off;
	grid on;
end
if ( plotOmegaVsZ1_cnstJ )
	numFigs++; figure(numFigs);
	funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
	plot( [0.0,norm(vecXE-vecX0)], [0.0,0.0], 'k-' );
	xlabel( "z1" );
	ylabel( "sqrt(omega cnstJ)" );
	title( "sqrt(omega cnstJ) vs z1" );
	if ( useAxEq )
		axis equal;
	end
	hold on;
	for n=1:numCurves
		plot( ...
		  curveDat(n).z1Vals, sqrt(curveDat(n).omegaVals_cnstJ), ...
		  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
		  'linewidth', curveDat(n).plot_lineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_markerSize );
		plot( ...
		  curveDat(n).z1Vals([1,end]), sqrt(curveDat(n).omegaVals_cnstJ([1,end])), ...
		  curveDat(n).plot_markerStyle, ...
		  'linewidth', curveDat(n).plot_bigLineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_bigMarkerSize );
	end
	hold off;
	grid on;
end
%
%
%
if ( plotOmegaVsN )
	numFigs++; figure(numFigs);
	funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
	plot( [0.0,norm(vecXE-vecX0)], [0.0,0.0], 'k-' );
	xlabel( "point index" );
	ylabel( "omega" );
	title( "omega vs point index" );
	hold on;
	for n=1:numCurves
		plot( ...
		  [0.0:curveDat(n).numPts-1.0], sqrt(curveDat(n).omegaVals), ...
		  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
		  'linewidth', curveDat(n).plot_lineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_markerSize );
		plot( ...
		  [0.0,curveDat(n).numPts-1.0], sqrt(curveDat(n).omegaVals([1,end])), ...
		  curveDat(n).plot_markerStyle, ...
		  'linewidth', curveDat(n).plot_bigLineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_bigMarkerSize );
	end
	hold off;
	grid on;
end
%
%
%
if ( plotZ1VsDAC )
	numFigs++; figure(numFigs);
	funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
	plot( [0.0,norm(vecXE-vecX0)], [0.0,0.0], 'k-' );
	xlabel( "DAC" );
	ylabel( "z1" );
	title( "z1 vs DAC" );
	if ( useAxEq )
		axis equal;
	end
	hold on;
	for n=1:numCurves
		plot( ...
		  curveDat(n).dacVals, curveDat(n).z1Vals, ...
		  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
		  'linewidth', curveDat(n).plot_lineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_markerSize );
		plot( ...
		  curveDat(n).dacVals([1,end]), curveDat(n).z1Vals([1,end]), ...
		  curveDat(n).plot_markerStyle, ...
		  'linewidth', curveDat(n).plot_bigLineWidth, ...
		  'color', curveDat(n).plot_color, ...
		  'markersize', curveDat(n).plot_bigMarkerSize );
	end
	hold off;
	grid on;
end
%
%
%
msg( __FILE__, __LINE__, sprintf( "Drawing graphs took %0.3fs.", time()-vizTime0 ) );
