function groot_xviz( grootXDatOut );
	mydefs;
	startTime = time();
	%
	n = 1;
	semilogy( grootXDatOut.s(n).matInfoA(:,2), grootXDatOut.s(n).matInfoA(:,4), 'o-' );
	cellAry_empty{n} = " ";
	cellAry_legend{n} = grootXDatOut.s(n).strSolverName;
	hold on;
	for n=2:grootXDatOut.algoSetPrm.n;
		semilogy( grootXDatOut.s(n).matInfoA(:,2), grootXDatOut.s(n).matInfoA(:,4), 'x-' );
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
