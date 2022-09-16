function scratch_zviz( csdo, prm=[] )
	mydefs;
	startTime = time();
	prm.verbLev = mygetfield( prm, "verbLev", VERBLEV__DETAILS );
	prm.valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	%
	assert( ~isempty( csdo ) );
	assert( ~isempty( mygetfield( csdo, "probSetPrm", [] ) ) );
	numProbs = mygetfield( csdo.probSetPrm, "numProbs", [] );
	assert( isposintscalar(numProbs) );
	%
	numAlgos = size( csdo.prob(1).grootXDatOut.s, 2 );
	assert( issize(csdo.prob(1).grootXDatOut.algoSetPrm.s,[1,numAlgos]) );
	%
	%for probIndex = 1 : numProbs
	%for algoIndex = 1 : numAlgos
	%	s(algoIndex).prob(probIndex) = csdo.prob(probIndex).grootXDatOut.s(algoIndex);
	%endfor
	%endfor
	%
	%
	for algoIndex = 1 : numAlgos
	for probIndex = 1 : numProbs
		s(algoIndex).vecGrootFlag(probIndex) = csdo.prob(probIndex).grootXDatOut.s(algoIndex).grootFlag;
		s(algoIndex).vecFevalCount(probIndex) = csdo.prob(probIndex).grootXDatOut.s(algoIndex).fevalCount;
	endfor
	endfor
	%
	vecPctCnvg = 100.0*(1:numProbs)/double(numProbs);
	algoIndex = 1;
	msk = (s(algoIndex).vecGrootFlag(:)==GROOT_FLAG__CNVG)';
	mymin = min(s(algoIndex).vecFevalCount(msk));
	semilogy( [ 0.0, vecPctCnvg((1:sum(msk))), vecPctCnvg(sum(msk)) ], [ mymin, sort(s(algoIndex).vecFevalCount(msk)), mymin ], 'o-' );
	cellAry_empty{algoIndex} = " ";
	cellAry_legend{algoIndex} = csdo.prob(1).grootXDatOut.s(algoIndex).strSolverName;
	hold on;
	for algoIndex = 2 : numAlgos
		msk = (s(algoIndex).vecGrootFlag(:)==GROOT_FLAG__CNVG)';
		mymin = min(s(algoIndex).vecFevalCount(msk));
		semilogy( [ 0.0, vecPctCnvg((1:sum(msk))), vecPctCnvg(sum(msk)) ], [ mymin, sort(s(algoIndex).vecFevalCount(msk)), mymin ], 'x-' );
		cellAry_empty{algoIndex} = " ";
		cellAry_legend{algoIndex} = csdo.prob(1).grootXDatOut.s(algoIndex).strSolverName;
	endfor
	ax = axis();
	axis([ 0.0, 100.0, ax(3), ax(4) ]);
	hold off;
	grid on;
	set( xlabel(""), "Interpreter", "none" );
	set( ylabel(""), "Interpreter", "none" );
	set( title(""), "Interpreter", "none" );
	set( legend( cellAry_empty, "location", "eastoutside"), "Interpreter", "none" );
	xlabel( "percentile" );
	ylabel( "feval count" );
	legend( cellAry_legend, "location", "eastoutside" );
	%
	%
	%
	perAlgo_numCnvg = zeros( numAlgos, 1 );
	perAlgo_numFail = zeros( numAlgos, 1 );
	perAlgo_sumCeval = zeros( numAlgos, 1 );
	perAlgo_sumCevalSq = zeros( numAlgos, 1 );
	for probIndex = 1 : numProbs
	for algoIndex = 1 : numAlgos
		this_s = csdo.prob(probIndex).grootXDatOut.s(algoIndex);
		switch ( this_s.grootFlag )
		case { GROOT_FLAG__CNVG }
			perAlgo_numCnvg(algoIndex)++;
			perAlgo_sumCeval(algoIndex) += this_s.fevalCount;
			perAlgo_sumCevalSq(algoIndex) += ((this_s.fevalCount)^2);
		case { GROOT_FLAG__FAIL, GROOT_FLAG__STOP }
			perAlgo_numFail(algoIndex)++;
		otherwise
			error([ "Unsupported value of grootFlag (\"" s.grootFlag "\"." ]);
		endswitch
	endfor
	endfor
	%
	perAlgo_cnvgFrac = perAlgo_numCnvg./(eps+perAlgo_numCnvg+perAlgo_numFail);
	perAlgo_cevalAvg = perAlgo_sumCeval./(eps+perAlgo_numCnvg);
	perAlgo_cevalSqAvg = perAlgo_sumCevalSq./(eps+perAlgo_numCnvg(algoIndex));
	perAlgo_cevalVarSq = perAlgo_cevalSqAvg - (perAlgo_cevalAvg.^2);
	for algoIndex = 1 : numAlgos
		this_cevalAvg = perAlgo_cevalAvg(algoIndex);
		this_cevalVarSq = perAlgo_cevalVarSq(algoIndex);
		if ( this_cevalVarSq > 0.0 )
			this_cevalVar = sqrt( this_cevalVarSq );
		else
			this_cevalVar = 0.0;
		endif
		this_strSolverName = csdo.prob(1).grootXDatOut.s(algoIndex).strSolverName;
		%
		if ( prm.verbLev >= VERBLEV__MAIN )
			msg( __FILE__, __LINE__, sprintf( ...
			  "  Solver  %2d  %15s:  %10.3f +/- %10.3f @ %6.2f%%.", ...
			  algoIndex, this_strSolverName, ...
			  perAlgo_cevalAvg(algoIndex), this_cevalVar, ...
			  100.0*perAlgo_cnvgFrac(algoIndex) ) );
		endif
	endfor
return;
endfunction
