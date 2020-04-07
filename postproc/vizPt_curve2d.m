function vizPt_curve2d( studyPtDat, indexB1, indexB2, strContour="omega", prm=[], datIn=[] )
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
	if ( strcmpi("colorbar",strContour) || strcmpi("colorbarf",strContour) )
		vizPt_curve2d_colorbar( studyPtDat, indexB1, indexB2, strContour, prm, datIn );
		return;
	end
	%
	numCurves = max(size(studyPtDat.curveDat));
	if ( indexB1 > 0 )
		m = studyPtDat.curveDat(abs(indexB1)).indexOfMin;
		vecBU1 = studyPtDat.curveDat(abs(indexB1)).matDelta(:,m);
	else
		vecBU1 = studyPtDat.curveDat(abs(indexB1)).matDelta(:,end);
	end
	if ( indexB2 > 0 )
		m = studyPtDat.curveDat(abs(indexB2)).indexOfMin;
		vecBU2 = studyPtDat.curveDat(abs(indexB2)).matDelta(:,m);
	else
		vecBU2 = studyPtDat.curveDat(abs(indexB2)).matDelta(:,end);
	end
	matBU = [ vecBU1, vecBU2 ];
	%
	vecX0 = studyPtDat.vecX0; % X of starting point.
	funchF = studyPtDat.funchF;
	funchFSupportsMultiArg = studyPtDat.funchFSupportsMultiArg;
	vecF0 = studyPtDat.vecF0;
	matV = studyPtDat.matV;
	matW = studyPtDat.matW;
	%
	funchZ = @(omega)( omega.^0.5 );
	numContours = 30;
	vecXP = mygetfield(prm,"vecXP",vecX0); % X of origin of plane.
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CALC.
	%
	matBV = myorth(matBU);
	vecBV1 = matBV(:,1);
	vecBV2 = matBV(:,2);
	%
	omegaLo = studyPtDat.omegaLo;
	omega0 = studyPtDat.omega0;
	omegaVizMin = 0.0;
	omegaVizMax = 2.0*omega0;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get vizDat.
	%
	for n=1:numCurves
		vizDat(n).numPts = size(studyPtDat.curveDat(n).matX,2);
		vizDat(n).matDelta = studyPtDat.curveDat(n).matX ...
		  - repmat(vecXP,[1,vizDat(n).numPts]);
		vizDat(n).rvecS1 = vecBV1' * vizDat(n).matDelta;
		vizDat(n).rvecS2 = vecBV2' * vizDat(n).matDelta;
		vizDat(n).col = studyPtDat.curveDat(n).col;
		vizDat(n).curveName = studyPtDat.curveDat(n).curveName;
		vizDat(n).indexOfMin = studyPtDat.curveDat(n).indexOfMin;
		%
		vizDat(n).matRes = vizDat(n).matDelta ...
		  - (vecBV1*vizDat(n).rvecS1) - (vecBV2*vizDat(n).rvecS2);
		vizDat(n).rvecSX = sqrt(sum((vizDat(n).matRes).^2,1));
		%
		% Distance to nearest neighbor.
		vizDat(n).rvecDTNN = studyPtDat.curveDat(n).rvecDTNN;
		vizDat(n).rvecInPlane = ( vizDat(n).rvecSX <= sqrt(eps) * vizDat(n).rvecDTNN(1) );
		vizDat(n).allInPlane = sum(vizDat(n).rvecInPlane) == vizDat(n).numPts;
		% DRaburn 2020.04.06...
		%  We could get a more detailed picture of when a curve touches the plane
		%  considering interpolation between the points, and finding the closest
		%  point. But, that would take time, and this is merely viz.
	end
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get gridS1, gridS2.
	%
	s1Max = 0.0;
	s1Min = 0.0;
	s2Max = 0.0;
	s2Min = 0.0;
	if (1)
	% Maximal area of interest.
	for n=1:numCurves
		s1Max = max([ s1Max, max(vizDat(n).rvecS1) ]);
		s1Min = min([ s1Min, min(vizDat(n).rvecS1) ]);
		s2Max = max([ s2Max, max(vizDat(n).rvecS2) ]);
		s2Min = min([ s2Min, min(vizDat(n).rvecS2) ]);
	end
	else
	% Minimal area of interest.
	for n=1:numCurves
	if ( abs(indexB1)==n || abs(indexB2)==n )
		m = vizDat(n).indexOfMin;
		s1Max = max([ s1Max, max(vizDat(n).rvecS1(1:m)) ]);
		s1Min = min([ s1Min, min(vizDat(n).rvecS1(1:m)) ]);
		s2Max = max([ s2Max, max(vizDat(n).rvecS2(1:m)) ]);
		s2Min = min([ s2Min, min(vizDat(n).rvecS2(1:m)) ]);
	end
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
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get gridOmega...
	%
	rvecGridS1 = reshape( gridS1, 1, numVals );
	rvecGridS2 = reshape( gridS2, 1, numVals );
	%
	matX = repmat(vecXP,[1,numVals]) + matBV * [ rvecGridS1; rvecGridS2 ];
	%
	if (strcmpi(strContour,"none"))
		gridOmega = [];
	elseif (strcmpi(strContour,"omega"))
		if ( funchFSupportsMultiArg )
			matF = funchF(matX);
		else
			parfor n=1:numVals
				matF(:,n) = funchF(matX(:,n));
			end
		end
		rvecOmega = 0.5*sum(matF.^2,1);
		gridOmega = reshape( rvecOmega, [ numS1Vals, numS2Vals ] );
	elseif (strcmpi(strContour,"omegaLin"))
		matFLin = repmat(vecF0,[1,numVals]) ...
		  + (matW*(matV'*(matX-repmat(vecX0,[1,numVals]))));
		rvecOmegaLin = 0.5*sum(matFLin.^2,1);
		gridOmega = reshape( rvecOmegaLin, [ numS1Vals, numS2Vals ] );
		%
		if ( sum(sum( ((matV*(matV'*matBV)) - matBV).^2 )) > eps )
			msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
			  "Warning: Contour basis space is outside linear model space." ));
		end
	else
		error(sprintf( "Invalid value of strContour!"  ));
	end
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PLOT.
	%
	strLegend = char([]);
	clf();
	hold off;
	%
	hold on;
	for n=1:numCurves
		tempLineWidth = 2;
		tempFMT = "-";
		if (vizDat(n).allInPlane )
			tempLineWidth = 7;
		end
		plot( ...
		  vizDat(n).rvecS1, vizDat(n).rvecS2, tempFMT, ...
		  "color", 0.7*vizDat(n).col, ...
		  "markersize", 8+4*(numCurves-n), ...
		  "linewidth", tempLineWidth );
		strLegend = [ strLegend; vizDat(n).curveName ];
	end
	legend( strLegend, "location", "northeastoutside" );
	%
	if ( ~strcmpi(strContour,"none") )
		vizPt_curve2d_contour( ...
		  gridS1, ...
		  gridS2, ...
		  gridOmega, ...
		  numContours, ...
		  omegaVizMin, ...
		  omegaLo, ...
		  omega0, ...
		  omegaVizMax, ...
		  @(omega_dummy)( (omega_dummy).^0.5 ), ...
		  "contourf" );
		%image( s1Vals, s2Vals, gridZV' );
		%set(get(gcf,"children"),"ydir","normal");
	end
	%
	plot( 0.0, 0.0, "k+", "linewidth", 2, "markersize", 30 );
	for n=1:numCurves
		if (~vizDat(n).allInPlane )
			plot( ...
			  vizDat(n).rvecS1(vizDat(n).rvecInPlane), ...
			  vizDat(n).rvecS2(vizDat(n).rvecInPlane), "o", ...
			  "color", 0.7*vizDat(n).col, ...
			  "linewidth", 2, ...
			  "markersize", 10+5*(numCurves-n) );
		end
		%
		m = vizDat(n).indexOfMin;
		plot( ...
		  vizDat(n).rvecS1(m), vizDat(n).rvecS2(m), "s", ...
		  "color", 0.7*vizDat(n).col, ...
		  "linewidth", 2, ...
		  "markersize", 10+5*(numCurves-n), ...
		  vizDat(n).rvecS1(m), vizDat(n).rvecS2(m), "x", ...
		  "color", 0.7*vizDat(n).col, ...
		  "linewidth", 2, ...
		  "markersize", 10+5*(numCurves-n) );
	end
	%
	if ( indexB1 > 0 )
		strXCoord = sprintf( "omegaMin %s step", ...
		  studyPtDat.curveDat(abs(indexB1)).curveName );
	else
		strXCoord = sprintf( "full %s step", ...
		  studyPtDat.curveDat(abs(indexB1)).curveName );
	end
	if ( indexB2 > 0 )
		strYCoord = sprintf( "omegaMin %s step", ...
		  studyPtDat.curveDat(abs(indexB2)).curveName );
	else
		strYCoord = sprintf( "full %s step", ...
		  studyPtDat.curveDat(abs(indexB2)).curveName );
	end
	%
	xlabel(sprintf( "dist along %s", strXCoord ));
	ylabel(sprintf( "dist along ortho %s", strYCoord ));
	title(sprintf( "2D Curve Plot (%s): %s, %s", strContour, strXCoord, strYCoord ));
	%
	axis square;
	grid on;
	hold off;
	%
	%
return;
end

%!test
%!	test_studyPt;
%!	vizPt_curve2d( studyPtDat, 1, 2 );
