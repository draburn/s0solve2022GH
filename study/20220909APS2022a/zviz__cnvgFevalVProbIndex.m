	for algoIndex = 1 : numAlgos
	for probIndex = 1 : numProbs
		s(algoIndex).vecGrootFlag(probIndex) = zcdo.prob(probIndex).grootXDatOut.s(algoIndex).grootFlag;
		s(algoIndex).vecFevalCount(probIndex) = zcdo.prob(probIndex).grootXDatOut.s(algoIndex).fevalCount;
	endfor
	endfor
	clear algoIndex;
	clear probIndex;
	%
	numFigs++; figure(numFigs);
	for n = 1 : numAlgos
		cnvgMask = (s(n).vecGrootFlag(:)==GROOT_FLAG__CNVG)';
		fevalCountMin = min(s(n).vecFevalCount(cnvgMask));
		semilogy( ...
		  [ (1:numProbs)(cnvgMask) ], ...
		  [ s(n).vecFevalCount(cnvgMask) ], ...
		  "linewidth", 2, "markersize", mksz{n}, mktp{n} );
		hold on;
	endfor
	hold off;
	grid on;
	set( xlabel(""), "Interpreter", "none" );
	set( ylabel(""), "Interpreter", "none" );
	set( title(""), "Interpreter", "none" );
	set( legend( cellAry_empty, "location", "eastoutside"), "Interpreter", "none" );
	xlabel( "problem index" );
	ylabel( "cnvg feval count" );
	title([ zcdo.runName ": cnvg feval V prob index" ]);
	legend( cellAry_legend, "location", "eastoutside" );
	%
	clear cnvgMask;
	clear fevalCountMin;
