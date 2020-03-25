function vizPt_curve1d( studyPtDat, indexB, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commoninit;
	thisFile = "vizPt_curve1d";
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% INIT.
	%
	numCurves = max(size(studyPtDat.curveDat));
	if ( indexB < 0 )
		m = studyPtDat.curveDat(abs(indexB)).indexOfMin;
		vecBU = studyPtDat.curveDat(abs(indexB)).matDelta(:,m);
	else
		vecBU = studyPtDat.curveDat(indexB).matDelta(:,end);
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CALC.
	%
	vecBV = myorth(vecBU);
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
	clf();
	strLegend = [];
	hold on;
	%
	for n=1:numCurves
		if (abs(indexB)==n)
			tempLineWidth = 1;
			tempFMT = "^-";
		else
			tempLineWidth = 1;
			tempFMT = "o-";
		end
		plot( ...
		  vizDat(n).rvecX, vizDat(n).rvecY, tempFMT, ...
		  "color", 0.8*vizDat(n).col, ...
		  "markersize", 4+2*(numCurves-n), ...
		  "linewidth", tempLineWidth );
		strLegend = [ strLegend; vizDat(n).curveName ];
	end
	legend( strLegend, "location", "northeastoutside" );
	%
	plot( ...
	  0.0, 0.0, ...
	  "k+", "linewidth", 2, "markersize", 30 );
	for n=1:numCurves
		m = vizDat(n).indexOfMin;
		if (abs(indexB)==n)
			tempLineWidth = 2;
		else
			tempLineWidth = 2;
		end
		plot( ...
		  vizDat(n).rvecX(end), vizDat(n).rvecY(end), "s", ...
		  "color", 0.6*vizDat(n).col, ...
		  "linewidth", tempLineWidth, ...
		  "markersize", 10+2*(numCurves-n), ...
		  vizDat(n).rvecX(m), vizDat(n).rvecY(m), "x", ...
		  "color", 0.6*vizDat(n).col, ...
		  "linewidth", tempLineWidth, ...
		  "markersize", 20+2*(numCurves-n) );
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
	ylabel(sprintf( "ortho dist" ));
	title(sprintf( "1D Curve Plot: %s", strXCoord ));
	%
	grid on;
	hold off;
	%
	%
return;
end

%!test
%!	test_studyPt;
%!	vizPt_curve1d( studyPtDat, 1 );
