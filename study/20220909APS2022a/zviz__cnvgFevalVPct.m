	for algoIndex = 1 : numAlgos
	for probIndex = 1 : numProbs
		s(algoIndex).vecGrootFlag(probIndex) = zcdo.prob(probIndex).grootXDatOut.s(algoIndex).grootFlag;
		s(algoIndex).vecFevalCount(probIndex) = zcdo.prob(probIndex).grootXDatOut.s(algoIndex).fevalCount;
	endfor
	endfor
	clear algoIndex;
	clear probIndex;
	%
	n = 1;
	cnvgMask = (s(n).vecGrootFlag(:)==GROOT_FLAG__CNVG)';
	overallCnvgFevalCountMin = min(s(n).vecFevalCount(cnvgMask));
	for n = 2 : numAlgos
		cnvgMask = (s(n).vecGrootFlag(:)==GROOT_FLAG__CNVG)';
		cnvgFevalCountMin = min(s(n).vecFevalCount(cnvgMask));
		overallCnvgFevalCountMin = min([ overallCnvgFevalCountMin, cnvgFevalCountMin ] );
	endfor
	clear n;
	clear cnvgMask;
	clear cnvgFevalCountMin;
	%
	numFigs++; figure(numFigs);
	vecPct = 100.0*(1:numProbs)/double(numProbs);
	for n = 1 : numAlgos
		cnvgMask = (s(n).vecGrootFlag(:)==GROOT_FLAG__CNVG)';
		cnvgFevalCountMin = min(s(n).vecFevalCount(cnvgMask));
		semilogy( ...
		  [ 0.0, vecPct((1:sum(cnvgMask))), vecPct(sum(cnvgMask)) ], ...
		  [ cnvgFevalCountMin, sort(s(n).vecFevalCount(cnvgMask)), overallCnvgFevalCountMin ], ...
		  "linewidth", 2, "markersize", mksz{n}, mktp{n} );
		hold on;
	endfor
	clear vecPct;
	clear n;
	clear cnvgMask;
	clear cnvgFevalCountMin;
	clear overallCnvgFevalCountMin;
	%
	ax = axis();
	axis([ 0.0, 100.0, ax(3), ax(4) ]);
	hold off;
	grid on;
	set( xlabel(""), "Interpreter", "none" );
	set( ylabel(""), "Interpreter", "none" );
	set( title(""), "Interpreter", "none" );
	set( legend( cellAry_empty, "location", "eastoutside"), "Interpreter", "none" );
	xlabel( "percentile" );
	ylabel( "cnvg feval count" );
	title([ zcdo.runName ": cnvg feval V pct" ]);
	legend( cellAry_legend, "location", "eastoutside" );
	clear fevalCountMin;
	clear ax;
