%  Function...
%    boo = isrealarray( x, sizeVec=[] )
%  Overview...
%    Part of myutil module.
%    Checks whether or not x is a real scalar/matrix/array.
%    Also, optionally, checks that size(x) == sizeVec.
%    See size() and issize() for more information.
function boo = isrealarray( x, sizeVec=[] );
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
	elseif (isempty(sizeVec))
		boo = true;
		return;
	else
		boo = issize( x, sizeVec );
		return;
	end
end

%!test
%!	assert( isrealarray( ones(4,5,6), [4,5,6] ) );
%!	assert( ~isrealarray( ones(4,5,6)+i, [4,5,6] ) );
%!	assert( isrealarray( ones(4,5,6) ) );
%!	assert( isrealarray( ones(4,1) ) );
%!	assert( isrealarray( ones(1,4) ) );
%!	assert( ~isrealarray( "Hello world!" ) );
%!	assert( ~isrealarray( 1==1 ) );
%!	assert( ~isrealarray( 1==0 ) );
%!	assert( ~isrealarray( [1,1]==[1,1] ) );
%!	assert( ~isrealarray( [1,1]==[1,0] ) );
%!	assert( ~isrealarray( [1,1]==[0,0] ) );
%!	assert( ~isrealarray( [char(1),"a"] ) );
%!	assert( isrealarray( 1 ) );
%!	assert( ~isrealarray( [1,Inf] ) );
%!	assert( ~isrealarray( [1,NaN] ) );
