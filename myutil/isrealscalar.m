%  Function...
%    boo = isrealscalar( x )
%  Overview...
%    Part of myutil module.
%    Checks whether or not x is a real scalar.
%    See isrealarray(), size(), and issize() for more information.
function boo = isrealscalar( x );
	boo = isrealarray( x, [1,1] );
return;
end

%!test
%!	assert( ~isrealscalar( ones(4,5,6) ) );
%!	assert( ~isrealscalar( 1+i ) );
%!	assert( ~isrealscalar( ones(4,1) ) );
%!	assert( ~isrealscalar( ones(1,4) ) );
%!	assert( isrealscalar( ones(1,1) ) );
%!	assert( ~isrealscalar( "Hello world!" ) );
%!	assert( ~isrealscalar( 1==1 ) );
%!	assert( ~isrealscalar( 1==0 ) );
%!	assert( ~isrealscalar( [1,1]==[1,1] ) );
%!	assert( ~isrealscalar( [1,1]==[1,0] ) );
%!	assert( ~isrealscalar( [1,1]==[0,0] ) );
%!	assert( ~isrealscalar( [char(1),"a"] ) );
%!	assert( isrealscalar( 1 ) );
%!	assert( ~isrealscalar( Inf ) );
%!	assert( ~isrealscalar( NaN ) );
