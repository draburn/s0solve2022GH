% Based on Kahan-Babushka-Neumaier, plus sorting.
% Is this best we can do?

function s = sum_kbn( x )
	assert( 1 == nargin );
	s = 0.0; % sum
	c = 0.0; % compensation term
	for n = 1 : numel(x)
		t = s + x(n);
		if ( abs(s) >= abs(x(n)) )
			c += (s-t) + x(n);
		else
			c += (x(n)-t) + s;
		endif
		s = t;
	endfor
	s += c;
return;
endfunction
