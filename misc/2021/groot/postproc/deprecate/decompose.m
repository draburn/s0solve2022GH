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
	%
	% If require matV to be orthonormal...
	%assert( abs( matV'*matV - eye(sizeK,sizeK) ) < 10.0 * (sizeK^3) * eps );
	%rvecC = vecX'*matV;
	%vecR = vecX - (matV*(rvecC'));
	%
	% But, instead...
	vecR = vecX;
	rvecC = zeros(1,sizeK);
	for k=1:sizeK
		denom = sum(matV(:,k).^2);
		assert( 0.0 < denom );
		rvecC(k) = vecX' * matV(:,k) / denom;
		vecR -= rvecC(k) * matV(:,k );
	end
	%
	normR = sqrt( sum(vecR.^2) );
return;
end


%!test
%!	test_decompose
