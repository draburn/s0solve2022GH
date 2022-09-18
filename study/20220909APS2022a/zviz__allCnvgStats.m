	perAlgo_numCnvg = zeros( numAlgos, 1 );
	perAlgo_numFail = zeros( numAlgos, 1 );
	for probIndex = 1 : numProbs
	for algoIndex = 1 : numAlgos
		this_s = zcdo.prob(probIndex).grootXDatOut.s(algoIndex);
		switch ( this_s.grootFlag )
		case { GROOT_FLAG__CNVG }
			perAlgo_numCnvg(algoIndex)++;
		case { GROOT_FLAG__FAIL, GROOT_FLAG__STOP }
			perAlgo_numFail(algoIndex)++;
		otherwise
			error([ "Unsupported value of grootFlag (\"" s.grootFlag "\"." ]);
		endswitch
	endfor
	endfor
	%
	%
	for algoIndex = 1 : numAlgos
	for probIndex = 1 : numProbs
		s(algoIndex).vecGrootFlag(probIndex) = zcdo.prob(probIndex).grootXDatOut.s(algoIndex).grootFlag;
		s(algoIndex).vecFevalCount(probIndex) = zcdo.prob(probIndex).grootXDatOut.s(algoIndex).fevalCount;
	endfor
	endfor
	clear algoIndex;
	clear probIndex;
	%
	perAlgo_sumCeval = zeros( numAlgos, 1 );
	perAlgo_sumCevalSq = zeros( numAlgos, 1 );
	for algoIndex = 1 : numAlgos
		perAlgo_sumCeval(algoIndex) = sum( s(algoIndex).vecFevalCount(allCnvgMask) );
		perAlgo_sumCevalSq(algoIndex) = sum( s(algoIndex).vecFevalCount(allCnvgMask).^2 );
	endfor
	%
	% Consider reporting over just cases where ALL converge!
	perAlgo_cnvgFrac = perAlgo_numCnvg./(eps+perAlgo_numCnvg+perAlgo_numFail);
	perAlgo_cevalAvg = perAlgo_sumCeval./(eps+perAlgo_numCnvg);
	perAlgo_cevalSqAvg = perAlgo_sumCevalSq./(eps+perAlgo_numCnvg);
	perAlgo_cevalVarSq = perAlgo_cevalSqAvg - (perAlgo_cevalAvg.^2);
	for algoIndex = 1 : numAlgos
		this_cevalAvg = perAlgo_cevalAvg(algoIndex);
		this_cevalVarSq = perAlgo_cevalVarSq(algoIndex);
		if ( this_cevalVarSq > 0.0 )
			this_cevalVar = sqrt( this_cevalVarSq );
		else
			this_cevalVar = 0.0;
		endif
		this_strSolverName = zcdo.prob(1).grootXDatOut.s(algoIndex).strSolverName;
		%
		if ( prm.verbLev >= VERBLEV__MAIN )
			msg( __FILE__, __LINE__, sprintf( ...
			  "  Solver  %2d  %15s:  %10.3f +/- %10.3f @ %6.2f%%.", ...
			  algoIndex, this_strSolverName, ...
			  perAlgo_cevalAvg(algoIndex), this_cevalVar, ...
			  100.0*perAlgo_cnvgFrac(algoIndex) ) );
		endif
	endfor
