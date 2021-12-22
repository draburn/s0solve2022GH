%  Function...
%    boo = isposintscalar( x )
%  Overview...
%    Part of myutil module.
%    Checks whether or not x is a positive integer scalar.
%    See size(), issize(), and isposintarray() for more information.
function boo = isposintscalar( x );
	boo = isposintarray( x, [1,1] );
return;
end

%!test
%!	assert( ~isposintscalar( ones(4,5,6) ) );
%!	assert( ~isposintscalar( 1+i ) );
%!	assert( ~isposintscalar( ones(4,1) ) );
%!	assert( ~isposintscalar( ones(1,4) ) );
%!	assert( isposintscalar( ones(1,1) ) );
%!	assert( ~isposintscalar( "Hello world!" ) );
%!	assert( ~isposintscalar( 1==1 ) );
%!	assert( ~isposintscalar( 1==0 ) );
%!	assert( ~isposintscalar( [1,1]==[1,1] ) );
%!	assert( ~isposintscalar( [1,1]==[1,0] ) );
%!	assert( ~isposintscalar( [1,1]==[0,0] ) );
%!	assert( ~isposintscalar( [char(1),"a"] ) );
%!	assert( isposintscalar( 1 ) );
%!	assert( ~isposintscalar( Inf ) );
%!	assert( ~isposintscalar( NaN ) );
%!	assert( ~isposintscalar(0.5) );