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
		vizPt_curve1d( studyPtDat, indexB, prm, datIn );
	elseif (2==numB)
		indexB1 = vecBasisIndex(1);
		indexB2 = vecBasisIndex(2);
		%subplot(2,2,1);
		figIndex++; figure(figIndex);
		vizPt_curve2d( studyPtDat, indexB1, indexB2, "omegaLin", prm, datIn );
		%subplot(2,2,2);
		figIndex++; figure(figIndex);
		vizPt_curve2d( studyPtDat, indexB1, indexB2, "omega", prm, datIn );
		figIndex++; figure(figIndex);
		vizPt_curve2d( studyPtDat, indexB1, indexB2, "colorbar", prm, datIn );
	elseif (3<=numB)
		error( "Not implemented!" );
	end
	%
return;
end
