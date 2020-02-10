%  Function...
%    boo = fleq( x, y, tol=10.0*eps )
%  Overview...
%    Part of myutil module.
%    Return true if the two inputs are equivalent to within tol.
%  Note: If one (but not both) of the values is exactly zero,
%    this will return FALSE.
function boo = fleq( x, y, tol=10.0*eps )
	boo = ( abs(x-y) <= abs(tol)*(abs(x)+abs(y)) );
return;
end

%!test
%!	assert( fleq(5,5) );
%!	assert( ~fleq(5,4) );
%!	assert( fleq(0,0) );
%!	assert( fleq(1+5*eps,1) );
%!	assert( ~fleq(1+100*eps,1) );
%!	assert( ~fleq(eps,0) );
%!	assert( fleq((1.0/eps)+1.0,(1.0/eps)) );
