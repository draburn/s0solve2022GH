function boo = isrealorinfscalar( x );
	assert( 1 == nargin );
	if ( isrealscalar(x) )
		boo = true;
		return;
	elseif ( ~issize(x,[1,1]) )
		boo = false;
		return;
	elseif ( isinf(x) )
		boo = true;
		return;
	end
	boo = false;
return;
end
