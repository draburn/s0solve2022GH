numIter = max(size(grootDatOut.iterDat));
numFigs = 0;
matI = eye(sizeX,sizeX);
vizNumPts = 20;
for n=[0,round(1.0*numIter/3.0),round(2.0*numIter/3.0),numIter]
	if (0==n)
		vecX = vecX0;
	else
		vecX = grootDatOut.iterDat(n).vecX;
	end
	matJ = demoFunc0101_evalJaco( vecX, funcPrm );
	vizPrm.numFigsOffset = numFigs;
	numNewFigs = vizLLMCurves( funchF, vecX, matI, matJ, vizNumPts, funcPrm.x0, vizPrm );
	numFigs += numNewFigs;
end
