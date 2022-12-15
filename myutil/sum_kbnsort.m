% Based on Kahan-Babushka-Neumaier, plus sorting.
% Is this best we can do?
%  ... No:
%   we could consider all possible summation routes,
%   and each possible summation route thereafter, iteratively,
%   to minimize the compensation term.
%  But, a single sort is *reasonably* optimal... I think.

function s = sum_kbnsort( y )
	assert( 1 == nargin );
	[ foo, si ] = sort( abs(y), "descend" );
	x = y(si);
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
