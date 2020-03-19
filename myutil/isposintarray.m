%  Function...
%    boo = isposintarray( x, sizeVec=[] )
%  Overview...
%    Part of myutil module.
%    Checks whether or not x is a positive integer scalar/matrix/array.
%    Also, optionally, checks that size(x) == sizeVec.
%    See size(), issize(), and isrealarray() for more information.
function boo = isposintarray( x, sizeVec=[] );
	if ( false == isrealarray(x,sizeVec) )
		boo = false;
		return;
	elseif ( 1 ~= prod(abs(x-round(x)) < 10.0*eps*x) )
		boo = false;
		return;
	elseif ( 1 ~= prod( 0 < x ) )
		boo = false;
		return;
	else
	boo = true;
return;
end

%!test
%!	assert( isposintarray( ones(4,5,6), [4,5,6] ) );
%!	assert( ~isposintarray( ones(4,5,6)+i, [4,5,6] ) );
%!	assert( isposintarray( ones(4,5,6) ) );
%!	assert( isposintarray( ones(4,1) ) );
%!	assert( isposintarray( ones(1,4) ) );
%!	assert( ~isposintarray( "Hello world!" ) );
%!	assert( ~isposintarray( 1==1 ) );
%!	assert( ~isposintarray( 1==0 ) );
%!	assert( ~isposintarray( [1,1]==[1,1] ) );
%!	assert( ~isposintarray( [1,1]==[1,0] ) );
%!	assert( ~isposintarray( [1,1]==[0,0] ) );
%!	assert( ~isposintarray( [char(1),"a"] ) );
%!	assert( isposintarray( 1 ) );
%!	assert( ~isposintarray( [1,Inf] ) );
%!	assert( ~isposintarray( [1,NaN] ) );
%!	assert( ~isposintarray( [1,0.5], [1,2] ) );
