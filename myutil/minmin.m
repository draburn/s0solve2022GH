function [ minminOfX, index1OfMinMin, index2OfMinMin ] = minmin( x )
	size1 = size(x,1);
	size2 = size(x,2);
	assert( isrealarray(x,[size1,size2]) )
	[ minOfX, index1OfMin ] = min( x );
	[ minminOfX, index2OfMinMin ] = min(minOfX);
	index1OfMinMin = index1OfMin(index2OfMinMin);
return;
end
