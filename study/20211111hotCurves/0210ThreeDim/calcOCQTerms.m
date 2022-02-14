function [ vecPhiHat, vecEta, lambdaAbsMin, matPsi, matLambda ] = calcOCQTerms( vecX0, funchF, matH, epsX=eps^0.3 )
	[ matPsi, matLambda ] = eig(matH);
	[ lambdaAbsMin, nOfAbsMin ] = min(abs(diag(matLambda)));
	vecPhiHat = matPsi(:,nOfAbsMin);
	if ( nargout >= 2 )
		vecEta = (funchF( vecX0 + epsX*vecPhiHat ) - funchF( vecX0 - epsX*vecPhiHat ))/(2.0*epsX);
	end
return;
end
