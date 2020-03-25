function vizPt_curve2d( studyPtDat, indexB1, indexB2, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commoninit;
	thisFile = "vizPt_curve2d";
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% INIT.
	%
	numCurves = max(size(studyPtDat.curveDat));
	if ( indexB1 < 0 )
		m = studyPtDat.curveDat(abs(indexB1)).indexOfMin;
		vecBU1 = studyPtDat.curveDat(abs(indexB1)).matDelta(:,m);
	else
		vecBU1 = studyPtDat.curveDat(indexB1).matDelta(:,end);
	end
	if ( indexB2 < 0 )
		m = studyPtDat.curveDat(abs(indexB2)).indexOfMin;
		vecBU2 = studyPtDat.curveDat(abs(indexB2)).matDelta(:,m);
	else
		vecBU2 = studyPtDat.curveDat(indexB2).matDelta(:,end);
	end
	matBU = [ vecBU1, vecBU2 ];
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CALC.
	%
	matBV = myorth(matBU);
	vecBV1 = matBV(:,1);
	vecBV2 = matBV(:,2);
	%
	for n=1:numCurves
		matDelta = studyPtDat.curveDat(n).matDelta;
		rvecVizX = vecBV1' * matDelta;
		rvecVizY = vecBV2' * matDelta;
		matRes = matDelta - (vecBV1 * rvecVizX) - (vecBV2 * rvecVizY);
		rvecVizOD = sqrt(sum((matRes.^2),1));
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
		if ( abs(indexB1)==n && abs(indexB2)==n )
			% Only possible if comparing a curve's end with omegaMin.
			tempLineWidth = 1;
			tempFMT = "^-";
		elseif (abs(indexB1)==n)
			tempLineWidth = 1;
			tempFMT = "^-";
		elseif (abs(indexB2)==n)
			tempLineWidth = 1;
			tempFMT = "v-";
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
	plot( 0.0, 0.0, "k+", "linewidth", 2, "markersize", 30 );
	for n=1:numCurves
		m = vizDat(n).indexOfMin;
		if (abs(indexB1)==n)
			tempLineWidth = 2;
		elseif (abs(indexB2)==n)
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
	if ( 0 < indexB1 )
		strXCoord = sprintf( "full %s step", ...
		  studyPtDat.curveDat(indexB1).curveName );
	else
		strXCoord = sprintf( "omegaMin %s step", ...
		  studyPtDat.curveDat(abs(indexB1)).curveName );
	end
	if ( 0 < indexB2 )
		strYCoord = sprintf( "full %s step", ...
		  studyPtDat.curveDat(indexB2).curveName );
	else
		strYCoord = sprintf( "omegaMin %s step", ...
		  studyPtDat.curveDat(abs(indexB2)).curveName );
	end
	%
	xlabel(sprintf( "dist along %s", strXCoord ));
	ylabel(sprintf( "dist along ortho %s", strYCoord ));
	title(sprintf( "2D Curve Plot: %s, %s", strXCoord, strYCoord ));
	%
	grid on;
	hold off;
	%
	borderzoom(0.1);
	%
	%
return;
end

%!test
%!	test_studyPt;
%!	vizPt_curve2d( studyPtDat, 1, 2 );
