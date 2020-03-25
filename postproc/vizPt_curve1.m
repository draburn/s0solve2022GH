function vizPt_curve1( studyPtDat, vecBU, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commoninit;
	thisFile = "vizPt_curve1";
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	numCurves = max(size(studyPtDat.curveDat));
	vecBV = myorth(vecBU);
	%
	figIndex++; figure(figIndex);
	clf();
	hold on;
	%
	for n=1:numCurves
		matDelta = studyPtDat.curveDat(n).matDelta;
		rvecVizX = vecBV' * matDelta;
		matRes = matDelta - (vecBV * rvecVizX);
		rvecVizY = sqrt(sum((matRes.^2),1));
		%
		plot( ...
		  rvecVizX, rvecVizY, 'o-' );
	end
	grid on;
	hold off;
	%
	%
return;
end
