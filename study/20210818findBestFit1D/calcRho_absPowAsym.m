function [ errFlag, vecRho ] = calcRho_absPowAsym( rhoArgs, vecZ )
	yLVals = (rhoArgs.xVals<vecZ(1)).*(abs(rhoArgs.xVals-vecZ(1)).^vecZ(2));
	yRVals = (rhoArgs.xVals>vecZ(1)).*(abs(rhoArgs.xVals-vecZ(1)).^vecZ(3));
	if (~isrealarray(yLVals))
		errFlag = true;
		vecRho = [];
		return;
		% For example, 47^209 = "Inf".
	end
	if (~isrealarray(yRVals))
		errFlag = true;
		vecRho = [];
		return;
		% For example, 47^209 = "Inf".
	end
	%
	matY = [ ones(size(yLVals)), yLVals, yRVals ];
	matD = diag(rhoArgs.dVals);
	matA = matD*matY;
	matB = matA'*matA;
	r = rcond(matB); % Instead, check if all yVals are nearly the same...? Slower?
	if ( ~isrealscalar(r) || eps >= r )
		errFlag = true;
		vecRho = [];
		return;
	end
	vecC = matB\(matA'*matD*(rhoArgs.fVals));
	vecRho = (rhoArgs.dVals.^2) .*( vecC(1) + (vecC(2)*yLVals) + (vecC(3)*yRVals) - rhoArgs.fVals );
	errFlag = false;
return;
end
%fExt = vecC(1);
%fL = vecC(2);
%fR = vecC(3);
