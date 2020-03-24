function vizPt_vs( studyPtDat, valFlagX, valFlagY, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commoninit;
	vizPt_init;
	thisFile = "vizPt_vs";
	%
	msg_warn( verbLev, thisFile, __LINE__, "This is TACish." );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	%
	curveDat = studyPtDat.curveDat;
	numCurves = size(curveDat,2);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	%
	if ( valFlagX == VIZPT_VALS__INDEX )
		for n=1:numCurves
			valsX(n).vals = (1:size(curveDat(n).rvecNu,2));
		end
		strNameX = "index";
	elseif ( valFlagX > VIZPT_VALS__VEC_FLIN_BASE )
		indexX = valFlagX - VIZPT_VALS__VEC_FLIN_BASE;
		for n=1:numCurves
			valsX(n).vals = curveDat(n).matY(indexX,:);
		end
		strNameX = sprintf("y(%d)",indexX);
	elseif ( valFlagX == VIZPT_VALS__NU );
		for n=1:numCurves
			valsX(n).vals = curveDat(n).rvecNu;
		end
		strNameX = "nu";
	else
		error(sprintf("Invalid value of matValFlag(%d) (%g).", n, matValFlag(n) ));
	end
	%
	if ( valFlagY == VIZPT_VALS__INDEX )
		for n=1:numCurves
			valsY(n).vals = (1:size(curveDat(n).rvecNu,2));
		end
		strNameY = "index";
	elseif ( valFlagY > VIZPT_VALS__VEC_FLIN_BASE )
		indexX = valFlagY - VIZPT_VALS__VEC_FLIN_BASE;
		for n=1:numCurves
			valsY(n).vals = curveDat(n).matY(indexX,:);
		end
		strNameY = sprintf("y(%d)",indexX);
	elseif ( valFlagY == VIZPT_VALS__NU );
		for n=1:numCurves
			valsY(n).vals = curveDat(n).rvecNu;
		end
		strNameY = "nu";
	else
		error(sprintf("Invalid value of matValFlag(%d) (%g).", n, matValFlag(n) ));
	end
	%
	for n=1:numCurves
		indexOfMin(n) = curveDat(n).indexOfMin;
		col(n,:) = curveDat(n).col;
	end
	%
	strLegend = [];
	clf;
	hold on
	%
	for n=1:numCurves
		plot( ...
		  valsX(n).vals, valsY(n).vals, ...
		  "o-", "color", 0.8*col(n,:) );
		strLegend = [ strLegend; curveDat(n).curveName ];
	end
	legend( strLegend, "location", "northeastoutside" );
	%
	plot( valsX(1).vals(1), valsY(1).vals(1), "k+", "linewidth", 2, "markersize", 20 );
	for n=1:numCurves
		m = indexOfMin(n);
		plot( ...
		  valsX(n).vals(end), valsY(n).vals(end), ...
		  "x", "color", 0.4*col(n,:), "linewidth", 2, "markersize", 10, ...
		  valsX(n).vals(m), valsY(n).vals(m), ...
		  "s", "color", 0.4*col(n,:), "linewidth", 2, "markersize", 10 );
	end
	grid on;
	hold off;
	%
	xlabel(strNameX);
	ylabel(strNameY);
	title([ strNameY " vs " strNameX ]);
	%
	%
return;
end
