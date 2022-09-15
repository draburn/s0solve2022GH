function stop = groot_fsolve__outputFcn( x, optimValues, state )
	global global_outputFcnDat;
	global_outputFcnDat.matInfoA = [ global_outputFcnDat.matInfoA; optimValues.iter, optimValues.funccount, norm(x), optimValues.fval ];
	global_outputFcnDat.matRecordX = [ global_outputFcnDat.matRecordX, x ];
	if ( stopsignalpresent() )
		msg( __FILE__, __LINE__, "\n" );
		msg( __FILE__, __LINE__, "Found stop signal." );
		stop = true;
	else
		stop = false;
	endif
endfunction
