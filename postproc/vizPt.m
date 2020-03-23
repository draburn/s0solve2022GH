function [ figIndex, retCode, datOut ] = vizPt( studyPtDat, plotList, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commoninit;
	vizPt_init;
	thisFile = "vizPt";
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	%
	curveDat = studyPtDat.curveDat;
	numCurves = size(curveDat,2);
	numPlots = max(size(plotList));
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for p=1:numPlots
	figIndex++; figure(figIndex);
	clf;
	strLegend = [];
	hold on
	switch (plotList(p))
	case {VIZPT_PLOT__OMEGA_VS_DELTANORM}
		for n=1:numCurves
			plot( ...
			  curveDat(n).rvecDeltaNorm, curveDat(n).rvecOmega, ...
			  "o-", "color", 0.8*curveDat(n).col );
			strLegend = [ strLegend; curveDat(n).curveName ];
		end
		legend( strLegend, "location", "northeastoutside" );
		plot( 0.0, curveDat(1).rvecOmega(1), "k+", "linewidth", 2, "markersize", 20 );
		for n=1:numCurves
			m = curveDat(n).indexOfMin;
			plot( ...
			  curveDat(n).rvecDeltaNorm(end), curveDat(n).rvecOmega(end), ...
			  "x", "color", 0.4*curveDat(n).col, "linewidth", 2, "markersize", 10, ...
			  curveDat(n).rvecDeltaNorm(m), curveDat(n).rvecOmega(m), ...
			  "s", "color", 0.4*curveDat(n).col, "linewidth", 2, "markersize", 10 );
		end
		xlabel( "deltaNorm" );
		ylabel( "omea" );
		title( "omega vs deltaNorm" );
		%
		%
	case {VIZPT_PLOT__Y1_Y2}
		k1 = 1;
		k2 = 2;
		for n=1:numCurves
			stepType = curveDat(n).stepType;
			m = curveDat(n).indexOfMin;
			plot( ...
			  curveDat(n).matY(k1,:), curveDat(n).matY(k2,:), ...
			  "o-", "color", 0.8*curveDat(n).col );
			strLegend = [ strLegend; curveDat(n).curveName ];
		end
		legend( strLegend, "location", "northeastoutside" );
		plot( 0.0, 0.0, "k+", "linewidth", 2, "markersize", 20 );
		for n=1:numCurves
			m = curveDat(n).indexOfMin;
			plot( ...
			  curveDat(n).matY(k1,end), curveDat(n).matY(k2,end), ...
			  "x", "color", 0.4*curveDat(n).col, "linewidth", 2, "markersize", 10, ...
			  curveDat(n).matY(k1,m), curveDat(n).matY(k2,m), ...
			  "s", "color", 0.4*curveDat(n).col, "linewidth", 2, "markersize", 10 );
		end
		xlabel(sprintf( "y(%d)", k1 ));
		ylabel(sprintf( "y(%d)", k2 ));
		title(sprintf( "curve components y(%d) and y(%d)", k1, k2 ));
		%
		%
	otherwise
		error(sprintf("Invalid value of plotList(%d) (%g).", p, plotList(p) ));
	end
	grid on;
	hold off;
	end
	%
	%
retCode = RETCODE__SUCCESS;
return;
end

%!test
%!	test_studyPt;
