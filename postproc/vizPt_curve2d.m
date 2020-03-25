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
	vecX0 = studyPtDat.vecX0;
	funchF = studyPtDat.funchF;
	funchFSupportsMultiArg = studyPtDat.funchFSupportsMultiArg;
	vecF0 = studyPtDat.vecF0;
	matV = studyPtDat.matV;
	matW = studyPtDat.matW;
	%
	VIZPT_CONTOUR__NONE = 0;
	VIZPT_CONTOUR__F = 1;
	VIZPT_CONTOUR__FLIN = 2;
	vizContour = VIZPT_CONTOUR__FLIN;
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
	%
	%
	s1Max = 0.0;
	s1Min = 0.0;
	s2Max = 0.0;
	s2Min = 0.0;
	for n=1:numCurves
	%if ( abs(indexB1)==n || abs(indexB2)==n )
	if (1)
		s1Max = max([ s1Max, max(vizDat(n).rvecX) ]);
		s1Min = min([ s1Min, min(vizDat(n).rvecX) ]);
		s2Max = max([ s2Max, max(vizDat(n).rvecY) ]);
		s2Min = min([ s2Min, min(vizDat(n).rvecY) ]);
	end
	end
	%
	s1Var = s1Max - s1Min;
	s2Var = s2Max - s2Min;
	s1Max = s1Max + (0.2*s1Var);
	s1Min = s1Min - (0.2*s1Var);
	s2Max = s2Max + (0.2*s2Var);
	s2Min = s2Min - (0.2*s2Var);
	%
	numS1Vals = 51;
	numS2Vals = 51;
	numVals = numS1Vals * numS2Vals;
	s1Vals = linspace( s1Min, s1Max, numS1Vals );
	s2Vals = linspace( s2Min, s2Max, numS2Vals );
	[ gridS1, gridS2 ] = ndgrid( s1Vals, s2Vals );
	%
	rvecGridS1 = reshape( gridS1, 1, numVals );
	rvecGridS2 = reshape( gridS2, 1, numVals );
	%
	matDelta = matBV * [ rvecGridS1; rvecGridS2 ];
	matX = repmat(vecX0,[1,numVals]) + matDelta;
	%
	switch (vizContour)
	case {VIZPT_CONTOUR__NONE}
		% Nothing to do.
	case {VIZPT_CONTOUR__F}
		if ( funchFSupportsMultiArg )
			matF = funchF(matX);
		else
			parfor n=1:numVals
				matF(:,n) = funchF(matX(:,n));
			end
		end
		rvecOmega = 0.5*sum(matF.^2,1);
		gridZ = reshape( rvecOmega, [ numS1Vals, numS2Vals ] ).^0.5;
	case {VIZPT_CONTOUR__FLIN}
		matFLin = repmat(vecF0,[1,numVals]) + (matW*(matV'*matDelta));
		rvecOmegaLin = 0.5*sum(matFLin.^2,1);
		gridZ = reshape( rvecOmegaLin, [ numS1Vals, numS2Vals ] ).^0.5;
	otherwise
		error(sprintf( "Invalid value of vizContour!"  ));
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
			tempFMT = "v-";
		elseif (abs(indexB2)==n)
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
	if ( VIZPT_CONTOUR__NONE ~= vizContour )
		contourf( gridS1, gridS2, gridZ, 30 );
		colormap( 0.7 + 0.3*jet(256) );
	end
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
	%
return;
end

%!test
%!	test_studyPt;
%!	vizPt_curve2d( studyPtDat, 1, 2 );
