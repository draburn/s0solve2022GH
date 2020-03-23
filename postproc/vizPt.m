function [ figIndex, retCode, datOut ] = vizPt( studyPtDat, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commoninit;
	thisFile = "vizPt";
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	%
	curveDat = studyPtDat.curveDat;
	numCurves = size(curveDat,2)
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	%
	%
	figIndex++; figure(figIndex);
	hold off;
	plot( 0.0, curveDat(1).rvecOmega(1), "k+", "linewidth", 3, "markersize", 20 );
	hold on
	for n=1:numCurves
		m = curveDat(n).indexOfMin;
		plot( ...
		  curveDat(n).rvecDeltaNorm, curveDat(n).rvecOmega, "o-", ...
		  curveDat(n).rvecDeltaNorm(end), curveDat(n).rvecOmega(end), "gx", ...
		  "linewidth", 3, "markersize", 20, ...
		  curveDat(n).rvecDeltaNorm(m), curveDat(n).rvecOmega(m), "cs", ...
		  "linewidth", 3, "markersize", 20 );
	end
	grid on;
	hold off;
	%
	%
	figIndex++; figure(figIndex);
	hold off;
	plot( 0.0, 0.0, "k+", "linewidth", 3, "markersize", 20 );
	hold on;
	for n=1:numCurves
		stepType = curveDat(n).stepType;
		m = curveDat(n).indexOfMin;
		plot( ...
		  curveDat(n).matY(1,:), curveDat(n).matY(2,:), "o-", ...
		  curveDat(n).matY(1,end), curveDat(n).matY(2,end), "gx", ...
		  "linewidth", 3, "markersize", 20, ...
		  curveDat(n).matY(1,m), curveDat(n).matY(2,m), "gs", ...
		  "linewidth", 3, "markersize", 20 );
		if (STEPTYPE__SPECIFIED_VECTOR==stepType)
			plot( curveDat(n).matY(1,end), curveDat(n).matY(2,end), ...
			  "rs", "linewidth", 3, "markersize", 20 );
		end
	end
	grid on;
	hold off;
	%
	%
retCode = RETCODE__SUCCESS;
return;
end

%!test
%!	test_studyPt;
