function boo = isrealvector( x, sz=[] );
	if (~isreal(x))
		boo = false;
		return;
		% Compatibility with MATLAB allows logical and char to pass isreal().
	elseif (~isnumeric(x))
		boo = false;
		return;
	elseif (0==prod(isfinite(x)))
		boo = false;
		return;
	elseif (isempty(sz))
		sz = max(size(x));
	end
	boo = issize( x, sz ) || issize( x, [1,sz] ) || issize( x, [sz,1] );
	return;
end
