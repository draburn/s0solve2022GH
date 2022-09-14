function stop = outfun__global_matCnvg( x, optimValues, state )
	global global_matCnvg;
	global_matCnvg = [ global_matCnvg; optimValues.funccount, optimValues.fval ];
	if ( stopsignalpresent() )
		msg( __FILE__, __LINE__, "Found stop signal." );
		stop = true;
	else
		stop = false;
	endif
endfunction
