function groot_xviz( grootXDatOut, prm=[] );
	mydefs;
	startTime = time();
	numFigs0 = mygetfield( prm, "numFigs0", 20 );
	numFigs = numFigs0;
	%
	numFigs++; figure(numFigs);
	for n=1:grootXDatOut.algoSetPrm.n;
		semilogy( ...
		  grootXDatOut.s(n).matInfoA(:,2), ...
		  grootXDatOut.s(n).matInfoA(:,4), ...
		  "linewidth", 2, "markersize", 3+3*n, ...
		  [ STR_MARKER_TYPES(1+mod(n,length(STR_MARKER_TYPES))) "-" ] );
		hold on;
		cellAry_empty{n} = " ";
		cellAry_legend{n} = grootXDatOut.s(n).strSolverName;
	endfor
	hold off;
	grid on;
	set( xlabel(""), "Interpreter", "none" );
	set( ylabel(""), "Interpreter", "none" );
	set( title(""), "Interpreter", "none" );
	set( legend( cellAry_empty, "location", "eastoutside"), "Interpreter", "none" );
	xlabel( "feval count" );
	ylabel( "||f||" );
	legend( cellAry_legend, "location", "eastoutside" );
return;
endfunction
