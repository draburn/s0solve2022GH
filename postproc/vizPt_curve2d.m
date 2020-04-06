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
	vecX0 = studyPtDat.vecX0;
	funchF = studyPtDat.funchF;
	funchFSupportsMultiArg = studyPtDat.funchFSupportsMultiArg;
	vecF0 = studyPtDat.vecF0;
	matV = studyPtDat.matV;
	matW = studyPtDat.matW;
	%
	funchZ = @(omega)( omega.^0.5 );
	numCol_fullScale = 51;
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% CALC.
	%
	matBV = myorth(matBU);
	vecBV1 = matBV(:,1);
	vecBV2 = matBV(:,2);
	%matBV' * matBU
	%
	omegaLo = studyPtDat.curveDat(1).rvecOmega(1);
	omega0 = studyPtDat.curveDat(1).rvecOmega(1);
	dMax = 0.0;
	for n=1:numCurves
		matDelta = studyPtDat.curveDat(n).matDelta;
		rvecVizX = vecBV1' * matDelta;
		rvecVizY = vecBV2' * matDelta;
		matRes = matDelta - (vecBV1 * rvecVizX) - (vecBV2 * rvecVizY);
		rvecVizR = sqrt(sum((matRes.^2),1));
		%
		vizDat(n).rvecX = rvecVizX;
		vizDat(n).rvecY = rvecVizY;
		vizDat(n).col = studyPtDat.curveDat(n).col;
		vizDat(n).curveName = studyPtDat.curveDat(n).curveName;
		vizDat(n).indexOfMin = studyPtDat.curveDat(n).indexOfMin;
		vizDat(n).rvecR = rvecVizR;
		%
		omegaLo = min([ omegaLo, studyPtDat.curveDat(n).rvecOmega ]);
		omegaLo = min([ omegaLo, studyPtDat.curveDat(n).rvecOmegaLin ]);
		dMax = max([ dMax, max(max(matDelta)) ]);
	end
	zLo = funchZ(omegaLo);
	z0 = funchZ(omega0);
	if (zLo>z0/numCol_fullScale)
		zCapMin = zLo * ((zLo/z0)^0.0);
	else
		zCapMin = 0.0;
	end
	%
	cmap_fullScale = 0.5 + ( 0.5 * jet(numCol_fullScale) );
	indexC_z0 = 1 + round( (numCol_fullScale-1.0) * 0.63 );
	zCapMax = zCapMin + ((z0-zCapMin)*(numCol_fullScale-1.0)/(indexC_z0-1.0));
	%
	%
	if (1)
		cmap_fullScale(1,:) *= 0.4;
		cmap_fullScale(end,:) *= 0.4;
		cmap_fullScale(indexC_z0,:) *= 0.4;
		cmap_fullScale(1,:) += 0.4;
		cmap_fullScale(end,:) += 0.4;
		cmap_fullScale(indexC_z0,:) += 0.4;
	else
	% DRaburn 2020.04.02...
	% This gets at the idea of "flagged" values, but
	% should be made *exact* using explicit contour levels.
	n0 = ceil(numCol_fullScale/50);
	for n=1:n0
		cmap_fullScale(n,:) *= 0.4 + 0.5*(n-1.0)/(n0-1.0);
		cmap_fullScale(end+1-n,:) *= 0.4 + 0.5*(n-1.0)/(n0-1.0);
		cmap_fullScale(indexC_z0+n,:) *= 0.4 + 0.5*(n-1.0)/(n0-1.0);
		cmap_fullScale(indexC_z0-n,:) *= 0.4 + 0.5*(n-1.0)/(n0-1.0);
	end
	cmap_fullScale(indexC_z0,:) *= 0.3;
	end
	%
	%
	% DRaburn 2020.04.02...
	% This zUnderMin and zOverMax are just guesses.
	zUnderMin = zCapMin - (zCapMax-zCapMin)*0.1/(numCol_fullScale);
	zOverMax = zCapMax + (zCapMax-zCapMin)*0.1/(numCol_fullScale);
	clev_fullScale = linspace(zUnderMin,zOverMax,numCol_fullScale+1);
	%
	%
	if (strcmpi(strContour,"colorbar"))
		xVals = (0:1);
		yVals = linspace(zCapMin,zCapMax,numCol_fullScale);
		[ gridX, gridY ] = ndgrid(xVals,yVals);
		gridZ = gridY;
		contourf( gridX, gridY, gridZ, 50 );
		colormap( cmap_fullScale );
		set(get(gcf,"children"),"xticklabel",[]);
		grid on;
		ylabel( "value" );
		title( "colorbar" );
		return;
	elseif (strcmpi(strContour,"fullColorMap"))
		xVals = (0:1);
		yVals = linspace(zCapMin,zCapMax,numCol_fullScale);
		[ gridX, gridY ] = ndgrid(xVals,yVals);
		gridZ = gridY;
		contourf( gridX, gridY, gridZ, clev_fullScale );
		colormap( cmap_fullScale );
		set(get(gcf,"children"),"xticklabel",[]);
		grid on;
		ylabel( "value" );
		title( "colorbar" );
		return;
	end
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
	if (strcmpi(strContour,"none"))
		gridO = [];
	elseif (strcmpi(strContour,"omega"))
		if ( funchFSupportsMultiArg )
			matF = funchF(matX);
		else
			parfor n=1:numVals
				matF(:,n) = funchF(matX(:,n));
			end
		end
		rvecOmega = 0.5*sum(matF.^2,1);
		gridO = reshape( rvecOmega, [ numS1Vals, numS2Vals ] );
	elseif (strcmpi(strContour,"omegaLin"))
		matFLin = repmat(vecF0,[1,numVals]) + (matW*(matV'*matDelta));
		rvecOmegaLin = 0.5*sum(matFLin.^2,1);
		gridO = reshape( rvecOmegaLin, [ numS1Vals, numS2Vals ] );
		%
		if ( sum(sum( ((matV*(matV'*matBV)) - matBV).^2 )) > eps )
			msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
			  "Warning: Contour basis space is outside linear model space." ));
		end
	else
		error(sprintf( "Invalid value of strContour!"  ));
	end
	%
	gridZ = funchZ( gridO );
	gridZC = cap( gridZ, zCapMin, zCapMax );
	zMin = min(min(gridZC));
	zMax = max(max(gridZC));
	%
	indexC_zObsMin = 1 + round((zMin-zCapMin)*(numCol_fullScale-1.0)/(zCapMax-zCapMin));
	indexC_zObsMax = 1 + round((zMax-zCapMin)*(numCol_fullScale-1.0)/(zCapMax-zCapMin));
	cmap = cmap_fullScale(indexC_zObsMin:indexC_zObsMax,:);
	gridZV = 1 + round((gridZC-zCapMin)*(numCol_fullScale-1.0)/(zCapMax-zCapMin));
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
		if ( abs(indexB1)==n && abs(indexB2)==n )
			% Only possible if comparing a curve's end with omegaMin.
			tempLineWidth = 2;
			tempFMT = "-";
		elseif (abs(indexB1)==n)
			tempLineWidth = 2;
			tempFMT = "-";
		elseif (abs(indexB2)==n)
			tempLineWidth = 2;
			tempFMT = "-";
		else
			tempLineWidth = 2;
			tempFMT = "-";
		end
		msk = (vizDat(n).rvecR<dMax*eps^0.5);
		if (sum(msk)==max(size(vizDat(n).rvecR)))
			tempLineWidth = 5;
		end
		plot( ...
		  vizDat(n).rvecX, vizDat(n).rvecY, tempFMT, ...
		  "color", 0.7*vizDat(n).col, ...
		  "markersize", 8+4*(numCurves-n), ...
		  "linewidth", tempLineWidth );
		strLegend = [ strLegend; vizDat(n).curveName ];
	end
	legend( strLegend, "location", "northeastoutside" );
	%
	if ( ~strcmpi(strContour,"none") )
if ( DEV_VIZPT_CONTOUR_VER > 0 )
		vizPt_curve2d_contour( ...
		  gridS1, ...
		  gridS2, ...
		  gridO, ...
		  30, ...
		  0.0, ...
		  studyPtDat.contourVizDat.omegaLo, ...
		  studyPtDat.contourVizDat.omega0, ...
		  2.0 * studyPtDat.contourVizDat.omega0, ...
		  @(omega_dummy)( sqrt(omega_dummy) ), ...
		  "contourf" );
else
		if (1)
			contourf( gridS1, gridS2, gridZC, clev_fullScale );
		elseif (1)
			contourf( gridS1, gridS2, gridZV, 30 );
		else
			image( s1Vals, s2Vals, gridZV' );
			set(get(gcf,"children"),"ydir","normal");
		end
		colormap( cmap );
end
	end
	%
	plot( 0.0, 0.0, "k+", "linewidth", 2, "markersize", 30 );
	for n=1:numCurves
		m = vizDat(n).indexOfMin;
		msk = (vizDat(n).rvecR<dMax*eps^0.5);
		if (sum(msk)==max(size(vizDat(n).rvecR)))
			msk = [];
		end
		plot( ...
		  vizDat(n).rvecX(msk), vizDat(n).rvecY(msk), "o", ...
		  "color", 0.7*vizDat(n).col, ...
		  "linewidth", 3, ...
		  "markersize", 10+5*(numCurves-n) );
		plot( ...
		  vizDat(n).rvecX(m), vizDat(n).rvecY(m), "s", ...
		  "color", 0.7*vizDat(n).col, ...
		  "linewidth", 1, ...
		  "markersize", 20+5*(numCurves-n), ...
		  vizDat(n).rvecX(m), vizDat(n).rvecY(m), "x", ...
		  "color", 0.7*vizDat(n).col, ...
		  "linewidth", 1, ...
		  "markersize", 20+5*(numCurves-n) );
		  %"markersize", 10+10*(numCurves-n), ...
		  %vizDat(n).rvecX(m), vizDat(n).rvecY(m), "v", ...
		  %"color", 0.7*vizDat(n).col, ...
		  %"linewidth", tempLineWidth, ...
		  %"markersize", 10+10*(numCurves-n) );
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
