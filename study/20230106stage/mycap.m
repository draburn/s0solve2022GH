function xOut = mycap( xIn, xLo, xHi )
	xOut = xIn;
	if ( isempty(xIn) )
		return;
	endif
	if ( ~isempty(xLo) )
		xOut( xOut < xLo ) = xLo;
	endif
	if ( ~isempty(xHi) )
		xOut( xOut > xHi ) = xHi;
	endif
return;
endfunction
