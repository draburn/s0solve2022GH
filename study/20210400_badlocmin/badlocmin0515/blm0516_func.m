function vecF = blm0429_func( matX, matA, matB, vecX0, indexNNonLin, c1, c0 )
	numVals = size(matX,2);
	matZ = matB*(matX-repmat(vecX0,[1,numVals]));
	rvecYNL = matZ(indexNNonLin,:).^3;
	rvecZNL = rvecYNL.*( (rvecYNL.^2) + (c1*rvecYNL) + c0 );
	matZ(indexNNonLin,:) = rvecZNL;
	vecF = matA*matZ;
end
