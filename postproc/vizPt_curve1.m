function vizPt_curve1( studyPtDat, indexB, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commoninit;
	thisFile = "vizPt_curve1";
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% INIT.
	%
	if ( indexB < 0 )
		m = studyPtDat.curveDat(abs(indexB)).indexOfMin;
		vecBU = studyPtDat.curveDat(abs(indexB)).matDelta(:,m);
	else
		vecBU = studyPtDat.curveDat(indexB).matDelta(:,end);
	end
	%
	numCurves = max(size(studyPtDat.curveDat));
	vecBV = myorth(vecBU);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CALC.
	%
	for n=1:numCurves
		matDelta = studyPtDat.curveDat(n).matDelta;
		rvecVizX = vecBV' * matDelta;
		matRes = matDelta - (vecBV * rvecVizX);
		rvecVizY = sqrt(sum((matRes.^2),1));
		%
		vizDat(n).rvecX = rvecVizX;
		vizDat(n).rvecY = rvecVizY;
		vizDat(n).col = studyPtDat.curveDat(n).col;
		vizDat(n).curveName = studyPtDat.curveDat(n).curveName;
		vizDat(n).indexOfMin = studyPtDat.curveDat(n).indexOfMin;
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PLOT.
	%
	figIndex++; figure(figIndex);
	clf();
	strLegend = [];
	hold on;
	%
	for n=1:numCurves
		plot( ...
		  vizDat(n).rvecX, vizDat(n).rvecY, "o-", ...
		  "color", 0.8*vizDat(n).col, "markersize", 4+2*(numCurves-n) );
		strLegend = [ strLegend; vizDat(n).curveName ];
	end
	legend( strLegend, "location", "northeastoutside" );
	%
	plot( ...
	  0.0, 0.0, ...
	  "k+", "linewidth", 2, "markersize", 30 );
	for n=1:numCurves
		m = vizDat(n).indexOfMin;
		plot( ...
		  vizDat(n).rvecX(end), vizDat(n).rvecY(end), "s", ...
		  "color", 0.6*vizDat(n).col, "linewidth", 2, "markersize", 10+2*(numCurves-n), ...
		  vizDat(n).rvecX(m), vizDat(n).rvecY(m), "x", ...
		  "color", 0.6*vizDat(n).col, "linewidth", 2, "markersize", 20+2*(numCurves-n) );
	end
	%
	if ( 0 < indexB )
		strXCoord = sprintf( "full %s step", ...
		  studyPtDat.curveDat(indexB).curveName );
	else
		strXCoord = sprintf( "omegaMin %s step", ...
		  studyPtDat.curveDat(abs(indexB)).curveName );
	end
	%
	xlabel(sprintf( "dist along %s", strXCoord ));
	ylabel(sprintf( "ortho dist"));
	title(sprintf("1D Curve Plot: %s",strXCoord));
	%
	grid on;
	hold off;
	%
	%
return;
end
