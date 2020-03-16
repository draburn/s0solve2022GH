function [ aryX, aryS ] = spangrid( matV, matSLoHi, rvecN, prm=[] )
	sizeX = size(matV,1);
	sizeK = size(matV,2);
	assert( isrealarray(matV,[sizeX,sizeK]) );
	assert( matV'*matV, eye(sizeK), (eps^0.75)*(sizeK^2) );
	assert( isrealarray(matSLoHi,[2,sizeK]) );
	assert( isrealarray(rvecN,[1,sizeK]) );
	assert( sizeK >= 1 );
	assert( sizeK <= sizeX );
	assert( rvecN(:) >= 1 );
	assert( round(rvecN(:)), rvecN(:), -(eps^0.75) );
	%
	numPts = prod(rvecN);
	matX = zeros([ sizeX, rvecN ]);
	%
	switch (sizeK)
	case {1}
		matS = linspace(matSLoHi(1),matSLoHi(2),rvecN(1));
	case {2}
		vecS1 = linspace(matSLoHi(1,1),matSLoHi(2,1),rvecN(1));
		vecS2 = linspace(matSLoHi(1,2),matSLoHi(2,2),rvecN(2));
		[ aryS1, aryS2 ] = ndgrid( vecS1, vecS2 );
		aryS(1,:,:) = aryS1;
		aryS(2,:,:) = aryS2;
		matS = reshape( aryS, [sizeK, numPts ] );
	case {3}
		vecS1 = linspace(matSLoHi(1,1),matSLoHi(2,1),rvecN(1));
		vecS2 = linspace(matSLoHi(1,2),matSLoHi(2,2),rvecN(2));
		vecS3 = linspace(matSLoHi(1,3),matSLoHi(2,3),rvecN(3));
		[ aryS1, aryS2, aryS3 ] = ndgrid( vecS1, vecS2, vecS3 );
		aryS(1,:,:,:) = aryS1;
		aryS(2,:,:,:) = aryS2;
		aryS(3,:,:,:) = aryS3;
		matS = reshape( aryS, [sizeK, numPts ] );
	case {4}
		vecS1 = linspace(matSLoHi(1,1),matSLoHi(2,1),rvecN(1));
		vecS2 = linspace(matSLoHi(1,2),matSLoHi(2,2),rvecN(2));
		vecS3 = linspace(matSLoHi(1,3),matSLoHi(2,3),rvecN(3));
		vecS4 = linspace(matSLoHi(1,4),matSLoHi(2,4),rvecN(4));
		[ aryS1, aryS2, aryS3, aryS4 ] = ndgrid( vecS1, vecS2, vecS3, vecS4 );
		aryS(1,:,:,:,:) = aryS1;
		aryS(2,:,:,:,:) = aryS2;
		aryS(3,:,:,:,:) = aryS3;
		aryS(4,:,:,:,:) = aryS3;
		matS = reshape( aryS, [sizeK, numPts ] );
	otherwise
		error(sprintf( "5+ dimensional case not supported." ));
		% Crib from ndgrid to get general case?
	end
	matX = matV * matS;
	aryX = reshape( matX, [sizeX, rvecN] );
	%
	return;
return;
end

%!test
%!	test_spangrid

