function stop = groot_fsolve__outputFcn( x, optimValues, state )
	global global_outputFcnDat;
	global_outputFcnDat.matCnvg = [ global_outputFcnDat.matCnvg; optimValues.funccount, optimValues.fval ];
	global_outputFcnDat.matX = [ global_outputFcnDat.matX, x ];
	if ( stopsignalpresent() )
		msg( __FILE__, __LINE__, "Found stop signal." );
		stop = true;
	else
		stop = false;
	endif
endfunction
