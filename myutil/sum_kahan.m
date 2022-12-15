function s = sum_kahan( x )
	assert( 1 == nargin );
	s = 0.0; % sum
	c = 0.0; % compensation term
	for n = 1 : numel(x)
		y = x(n) - c;
		t = s + y;
		c = ( t - s ) - y;
		s = t;
	endfor
return;
endfunction
