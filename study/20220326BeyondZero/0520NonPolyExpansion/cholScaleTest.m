clear;
if (0)
	setprngstates(48245168); sz = 500;
	matSF = diag(exp(0.0*randn(sz,1)));
	matSX = diag(exp(10.0*randn(sz,1)));
	matJ = matSF*( randn(sz,sz).*exp(2.0*randn(sz,sz)) )/matSX;
	vecX = matSX*( randn(sz,1).*exp(2.0*randn(sz,1)) );
else
	setprngstates(); sz = 1000;
	matSF = diag(exp(2.0*randn(sz,1)));
	matSX = diag(exp(10.0*randn(sz,1)));
	matJ = matSF*( randn(sz,sz).*exp(2.0*randn(sz,sz)) )/matSX;
	vecX = matSX*( randn(sz,1).*exp(2.0*randn(sz,1)) );
endif
vecF = matJ*vecX;
%
matH = matJ'*matJ;
vecG = matJ'*vecF;
%
matS = diag(sqrt( diag(matH) + eps*norm(diag(matH)) ));
matSInv = diag(1.0./sqrt( diag(matH) + eps*norm(diag(matH)) ));
matHHat = matSInv*matH*matSInv;
matRHat = chol( matHHat );
matR = chol( matH );
%
rc = [ rcond(matH), rcond(matHHat) ]
cr = [ max(diag(matRHat))/min(diag(matRHat)), max(diag(matR))/min(diag(matR)) ]
%
%resOrig = [ norm(matH\vecG-vecX), norm( matR\(matR'\vecG) - vecX ) ]/norm(vecX)
%resHat = [ norm(matSInv*(matHHat\(matSInv*vecG))-vecX), norm( matSInv*(matRHat\(matRHat'\(matSInv*vecG))) - vecX ) ]/norm(vecX)
res = [ norm( matR\(matR'\vecG) - vecX ), norm( matSInv*(matRHat\(matRHat'\(matSInv*vecG))) - vecX ) ]/norm(vecX)
