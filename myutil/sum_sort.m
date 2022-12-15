function s = sum_sort( y )
	assert( 1 == nargin );
	[ foo, si ] = sort( abs(y), "descend" );
	x = y(si);
	s = sum(x);
return;
endfunction
