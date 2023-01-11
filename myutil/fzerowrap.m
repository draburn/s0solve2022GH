function [ x, fzero_x, fzero_fval, fzero_info, fzero_output ] = fzerowrap( fzero_fun, fzero_x0, fzero_options )
	if ( nargin <= 1 )
		[ fzero_x, fzero_fval, fzero_info, fzero_output ] = fzero( fzero_fun );
	elseif ( nargin <= 2 )
		[ fzero_x, fzero_fval, fzero_info, fzero_output ] = fzero( fzero_fun, fzero_x0 );
	elseif ( nargin <= 3 )
		[ fzero_x, fzero_fval, fzero_info, fzero_output ] = fzero( fzero_fun, fzero_x0, fzero_options );
	endif
	%
	xOut = fzero_output.bracketx;
	fOut = fzero_output.brackety;
	if ( abs(fOut(2)-fOut(1)) > sqrt(eps)*(abs(fOut(2))+abs(fOut(1))) )
		x = ( xOut(1)*fOut(2) - xOut(2)*fOut(1) ) / ( fOut(2) - fOut(1) );
	else
		x = fzero_x;
	endif
return;
endfunction
