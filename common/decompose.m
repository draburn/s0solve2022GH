function [ rvecC, normR ] = decompose( vecX, matV )
	sizeX = size(vecX,1);
	assert( isrealarray(vecX,[sizeX,1]) );
	if (isempty(matV))
		rvecC = [];
		normR = sqrt( sum(vecX.^2) );
		return;
	end
	sizeK = size(matV,2);
	assert( sizeK <= sizeX );
	assert( isrealarray(matV,[sizeX,sizeK]) );
	assert( abs( matV'*matV - eye(sizeK,sizeK) ) < 10.0 * (sizeK^3) * eps );
	rvecC = vecX'*matV;
	vecR = vecX - (matV*(rvecC'));
	normR = sqrt( sum(vecR.^2) );
return;
end


%!test
%!	test_decompose
