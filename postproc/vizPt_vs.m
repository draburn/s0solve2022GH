function vizPt_vs( varargin )
	% 2020.03.24: May have been cleaner to use a cell:
	%  vizPts_vs( studyPtDat, "nu", { "y", 1 } ).
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	prm = [];
	datIn = [];
	commoninit;
	thisFile = "vizPt_vs";
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	argIndex = 1;
	%
	studyPtDat = varargin(argIndex){1}; argIndex++;
	numCurves = size(studyPtDat.curveDat,2);
	%
	for d=1:2
		strName = varargin(argIndex){1}; argIndex++;
		assert( ischar(strName) );
		%
		if ( strcmpi(strName,"y") )
			strField = "matY";
			indexC = varargin(argIndex){1}; argIndex++;
			assert( isposintscalar(indexC) );
		elseif ( strcmpi(strName,"delta") )
			indexC = varargin(argIndex){1}; argIndex++;
			strField = "matDelta";
			assert( isposintscalar(indexC) );
		elseif ( strcmpi(strName,"x") )
			strField = "matX";
			indexC = varargin(argIndex){1}; argIndex++;
			assert( isposintscalar(indexC) );
		elseif ( strcmpi(strName,"f") )
			strField = "matF";
			indexC = varargin(argIndex){1}; argIndex++;
			assert( isposintscalar(indexC) );
		elseif ( strcmpi(strName,"fLin") )
			strField = "matFLin";
			indexC = varargin(argIndex){1}; argIndex++;
			assert( isposintscalar(indexC) );
		elseif ( strcmpi(strName,"index") )
			strField = "rvecIndex";
			indexC = 0;
		elseif ( strcmpi(strName,"nu") )
			strField = "rvecNu";
			indexC = 0;
		elseif ( strcmpi(strName,"mu") )
			warning( "You entered \"mu\". I'll assume you meant \"nu\"." );
			strField = "rvecNu";
			indexC = 0;
		elseif ( strcmpi(strName,"deltaNorm") )
			strField = "rvecDeltaNorm";
			indexC = 0;
		elseif ( strcmpi(strName,"omega") )
			strField = "rvecOmega";
			indexC = 0;
		elseif ( strcmpi(strName,"omegaLin") )
			strField = "rvecOmegaLin";
			indexC = 0;
		elseif ( strcmpi(strName,"dac") )
			strField = "rvecDAC";
			indexC = 0;
		else
			error(sprintf( "Unsupported value of strName as argument %d.", argIndex-1 ));
		end
		%
		for n=1:numCurves
			% These values are actually indp (n).
			% Not seeing a simple alternative place to put them.
			vizT(n).strName = strName;
			vizT(n).strField = strField;
			vizT(n).indexC = indexC;
			%
			if ( 0==indexC )
				vizT(n).rvecVals = getfield(studyPtDat.curveDat(n),strField);
			else
				matTemp = getfield(studyPtDat.curveDat(n),strField);
				assert( indexC >= 1 );
				assert( indexC <= size(matTemp,1) );
				vizT(n).rvecVals = matTemp(indexC,:);
			end
			%
			if (0)
				vizT(n).rvecVals += 1.0E-3 * randn(size(vizT(n).rvecVals)) ...
				  * (max(vizT(n).rvecVals)-min(vizT(n).rvecVals));
			end
		end
		%
		%
		switch (d)
		case {1}
			vizX = vizT;
		case {2}
			vizY = vizT;
		otherwise
			error("Impossible!");
		end
	end
	%
	% Common to vizX and vizY.
	for n=1:numCurves
		vizC(n).col = studyPtDat.curveDat(n).col;
		vizC(n).indexOfMin = studyPtDat.curveDat(n).indexOfMin;
		vizC(n).curveName = studyPtDat.curveDat(n).curveName;
	end
	%
	%
	strLegend = [];
	clf;
	hold on
	%
	for n=1:numCurves
		plot( ...
		  vizX(n).rvecVals, vizY(n).rvecVals, "o-", ...
		  "color", 0.7*vizC(n).col, ...
		  "markersize", 4+2*(numCurves-n) );
		strLegend = [ strLegend; vizC(n).curveName ];
	end
	legend( strLegend, "location", "northeastoutside" );
	%
	plot( vizX(1).rvecVals(1), vizY(1).rvecVals(1), "k+", ...
	  "linewidth", 2, "markersize", 30 );
	for n=1:numCurves
		m = vizC(n).indexOfMin;
		plot( ...
		  vizX(n).rvecVals(end), vizY(n).rvecVals(end), "s", ...
		  "color", 0.5*vizC(n).col, ...
		  "linewidth", 2, ...
		  "markersize", 10+2*(numCurves-n), ...
		  vizX(n).rvecVals(m), vizY(n).rvecVals(m), "x", ...
		  "color", 0.5*vizC(n).col, ...
		  "linewidth", 2, ...
		  "markersize", 20+2*(numCurves-n) );
	end
	grid on;
	hold off;
	%
	xlabel(vizX(1).strName);
	ylabel(vizY(1).strName);
	title([ vizY(1).strName " vs " vizX(1).strName ]);
	%
	%
return;
end

%!test
%!	test_studyPt
%!	vizPt_vs( studyPtDat, "nu", "y", 1 );
