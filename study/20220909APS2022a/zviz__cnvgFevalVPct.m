	for algoIndex = 1 : numAlgos
	for probIndex = 1 : numProbs
		s(algoIndex).vecGrootFlag(probIndex) = csdo.prob(probIndex).grootXDatOut.s(algoIndex).grootFlag;
		s(algoIndex).vecFevalCount(probIndex) = csdo.prob(probIndex).grootXDatOut.s(algoIndex).fevalCount;
	endfor
	endfor
	clear algoIndex;
	clear probIndex;
	%
	vecPct = 100.0*(1:numProbs)/double(numProbs);
	for n = 1 : numAlgos
		cnvgMask = (s(n).vecGrootFlag(:)==GROOT_FLAG__CNVG)';
		fevalCountMin = min(s(n).vecFevalCount(cnvgMask));
		semilogy( ...
		  [ 0.0, vecPct((1:sum(cnvgMask))), vecPct(sum(cnvgMask)) ], ...
		  [ fevalCountMin, sort(s(n).vecFevalCount(cnvgMask)), fevalCountMin/2.0 ], ...
		  "linewidth", 2, "markersize", mksz{n}, mktp{n} );
		hold on;
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
