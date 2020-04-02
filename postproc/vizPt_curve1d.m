function vizPt_curve1d( studyPtDat, matIndexB, prm=[], datIn=[] )
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
	sizeB = size(matIndexB,2);
	for n=1:sizeB
		indexB = matIndexB(:,n);
		if ( indexB > 0 )
			m = studyPtDat.curveDat(abs(indexB)).indexOfMin;
			vecBU = studyPtDat.curveDat(abs(indexB)).matDelta(:,m);
		else
			vecBU = studyPtDat.curveDat(abs(indexB)).matDelta(:,end);
		end
		matBU(:,n) = vecBU;
		clear indexB;
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CALC.
	%
	matBV = myorth(matBU);
	%
	for n=1:numCurves
		matDelta = studyPtDat.curveDat(n).matDelta;
		matVizX = matBV' * matDelta;
		rvecVizX = sqrt(sum(matVizX.^2,1));
		matRes = matDelta - (matBV * matVizX);
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
		%if (abs(indexB)==n)
		if (0)
			tempLineWidth = 1;
			tempFMT = "v-";
		else
			tempLineWidth = 1;
			tempFMT = "o-";
		end
		plot( ...
		  vizDat(n).rvecX, vizDat(n).rvecY, tempFMT, ...
		  "color", 0.7*vizDat(n).col, ...
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
		%if (abs(indexB)==n)
		if (0)
			tempLineWidth = 2;
		else
			tempLineWidth = 2;
		end
		plot( ...
		  vizDat(n).rvecX(end), vizDat(n).rvecY(end), "s", ...
		  "color", 0.5*vizDat(n).col, ...
		  "linewidth", tempLineWidth, ...
		  "markersize", 10+2*(numCurves-n), ...
		  vizDat(n).rvecX(m), vizDat(n).rvecY(m), "x", ...
		  "color", 0.5*vizDat(n).col, ...
		  "linewidth", tempLineWidth, ...
		  "markersize", 20+2*(numCurves-n) );
	end
	%
	if (0)
	if ( indexB > 0 )
		strXCoord = sprintf( "omegaMin %s step", ...
		  studyPtDat.curveDat(abs(indexB)).curveName );
	else
		strXCoord = sprintf( "full %s step", ...
		  studyPtDat.curveDat(abs(indexB)).curveName );
	end
	%
	xlabel(sprintf( "dist along %s", strXCoord ));
	ylabel(sprintf( "ortho dist" ));
	title(sprintf( "1D Curve Plot: %s", strXCoord ));
	end
	%
	xlabel(sprintf( "dist in space" ));
	ylabel(sprintf( "ortho dist" ));
	title(sprintf( "1D Curve Plot" ));
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
