function val = myternary( tf, val1, val2 )
	if ( tf )
		val = val1;
	else
		val = val2;
	endif
endfunction
