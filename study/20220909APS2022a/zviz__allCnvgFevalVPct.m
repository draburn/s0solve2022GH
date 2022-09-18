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
	allCnvgMask = (s(n).vecGrootFlag(:)==GROOT_FLAG__CNVG)';
	for n = 2 : numAlgos
		cnvgMask = (s(n).vecGrootFlag(:)==GROOT_FLAG__CNVG)';
		allCnvgMask &= cnvgMask;
	endfor
	clear cnvgMask;
	clear cnvgFevalCountMin;
	%
	n = 1;
	overallCnvgFevalCountMin = min(s(n).vecFevalCount(allCnvgMask));
	for n = 2 : numAlgos
		cnvgFevalCountMin = min(s(n).vecFevalCount(allCnvgMask));
		overallCnvgFevalCountMin = min([ overallCnvgFevalCountMin, cnvgFevalCountMin ] );
	endfor
	%
	numCnvg = sum(double(allCnvgMask));
	msg( __FILE__, __LINE__, sprintf( "All-solver convergence percent: %0.3f.", numCnvg*100.0/numProbs ) );
	%
	numFigs++; figure(numFigs);
	vecPct = 100.0*(1:numCnvg)/double(numCnvg);
	for n = 1 : numAlgos
		cnvgFevalCountMin = min(s(n).vecFevalCount(allCnvgMask));
		semilogy( ...
		  [ 0.0, vecPct, vecPct(end) ], ...
		  [ cnvgFevalCountMin, sort(s(n).vecFevalCount(allCnvgMask)), overallCnvgFevalCountMin ], ...
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
	xlabel( "all cnvg percentile" );
	ylabel( "cnvg feval count" );
	title([ zcdo.runName ": allcnvg fevalVpct" ]);
	legend( cellAry_legend, "location", "eastoutside" );
	clear fevalCountMin;
	clear ax;
