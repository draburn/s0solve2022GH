function vecF = blm0429_func_alt( matX, matA, matB, vecX0, indexNNonLin, c1, c0 )
	numVals = size(matX,2);
	matAB = matA*matB;
	vecANL = matA(:,indexNNonLin);
	vecBNL = matB(indexNNonLin,:)';
	matD = (matX-repmat(vecX0,[1,numVals]));
	rvecYNL = vecBNL' * matD;
	rvecZNL = rvecYNL.*( (rvecYNL.^2) + (c1*rvecYNL) + c0 );
	vecF = matAB*matD + (vecANL*(rvecZNL-rvecYNL));
end
