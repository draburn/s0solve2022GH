%  Function...
%    boo = issize( x, sizeVec );
%  Overview...
%    Part of myutil module.
%    Checks whether or not size(x) == sizeVec.
%    Note that size() returns a column vector with at least two columns.
%    See size() for more information.
function boo = issize( x, sizeVec );
	szX = size(x);
	for n=1:2
		if (size(szX,n)~=size(sizeVec,n))
			boo = false;
			return;
		end
	end
	for n=1:size(sizeVec,2)
		if (szX(n)~=sizeVec(n))
			boo = false;
			return;
		end
	end
	boo = true;
	return;
end

%!test
%!	assert( issize( ones(4,5,6), [4,5,6] ) );
%!	assert( ~issize( ones(4,5,6), [4,5,7] ) );
%!	assert( ~issize( ones(4,5), [4,5,6] ) );
%!	assert( ~issize( ones(4,1), [4] ) );
%!	assert( issize( ones(4,1), [4,1] ) );
%!	assert( issize( [], [0,0] ) );
%!	assert( ~issize( [], 0 ) );
