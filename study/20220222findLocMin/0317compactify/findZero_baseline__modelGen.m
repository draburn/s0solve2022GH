function [ matJ, ary3Kappa, funchFModel, modelGen_datOut ] = findZero_baseline__modelGen( vecX, vecF, vecX_prev, vecF_prev, matJ_prev, funchF, modelGen_prm )
	modelGen_datOut = [];
	modelGen_datOut.fevalCount = 0;
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	matJ = zeros(sizeF,sizeX);
	epsFD = mygetfield( modelGen_prm, "epsFD", eps^0.25 );
	for n = 1 : sizeX
		vecXP = vecX; vecXP(n) += epsFD;
		vecXM = vecX; vecXM(n) -= epsFD;
		vecFP = funchF(vecXP); modelGen_datOut.fevalCount++;
		vecFM = funchF(vecXM); modelGen_datOut.fevalCount++;
		matJ(:,n) = (vecFP-vecFM)/(2.0*epsFD);
	endfor
	ary3Kappa = [];
	funchFModel = @(x)( vecF + matJ*(x-vecX) );
return;
endfunction
