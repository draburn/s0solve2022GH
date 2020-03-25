function vizPt_curve( studyPtDat, vecBasisIndex, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commoninit;
	thisFile = "vizPt_curve";
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	numB = max(size(vecBasisIndex));
	numCurves = max(size(studyPtDat.curveDat));
	%
	if (0==numB)
		error("vecBasisIndex appears to be empty.");
	elseif (1==numB)
		indexB = vecBasisIndex(1);
		vizPt_curve1( studyPtDat, indexB, prm, datIn );
	elseif (2==numB)
		error( "Not implemented!" );
	elseif (3<=numB)
		error( "Not implemented!" );
	end
	
	%
	%
return;
end
