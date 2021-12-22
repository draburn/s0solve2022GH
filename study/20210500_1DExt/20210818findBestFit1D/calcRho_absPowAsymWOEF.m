function vecRho = calcRho_absPowAsymWOEF( rhoArgs, vecZ )
	vecRho = NaN;
	yLVals = (rhoArgs.xVals<vecZ(1)).*(abs(rhoArgs.xVals-vecZ(1)).^vecZ(2));
	yRVals = (rhoArgs.xVals>vecZ(1)).*(abs(rhoArgs.xVals-vecZ(1)).^vecZ(3));
	if (~isrealarray(yLVals))
		return;
		% For example, 47^209 = "Inf".
	end
	if (~isrealarray(yRVals))
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
		return;
	end
	vecC = matB\(matA'*matD*(rhoArgs.fVals));
	vecRho = (rhoArgs.dVals.^2) .*( vecC(1) + (vecC(2)*yLVals) + (vecC(3)*yRVals) - rhoArgs.fVals );
return;
end
%fExt = vecC(1);
%fL = vecC(2);
%fR = vecC(3);
