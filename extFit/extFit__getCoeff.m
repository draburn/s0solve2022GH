function vecC = extFit__getCoeff( bigX, bigP, xVals, fVals, wVals = [] )
	if (isempty(wVals))
		wVals = ones(size(xVals));
	end
	dVals = abs( xVals - bigX ).^bigP;
	sigma1 = sum(wVals);
	sigmaD = sum(wVals.*dVals);
	sigmaDD = sum(wVals.*dVals.*dVals);
	sigmaF = sum(wVals.*fVals);
	sigmaFD = sum(wVals.*dVals.*fVals);
	vecC = [ sigma1, sigmaD; sigmaD, sigmaDD ] \ [ sigmaF; sigmaFD ];
return;
end

%!test
%!	setprngstates(91738368);
%!	bigA = randn();
%!	bigB = randn();
%!	bigX = randn()*exp(3*randn());
%!	bigP = 0.3 + abs(randn());
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = randn(1,5);
%!	fVals = funchF(xVals);
%!	vecC = extFit__getCoeff( bigX, bigP, xVals, fVals );
%!	%echo__vecC = vecC
%!	assert(fleq(vecC(1),bigA));
%!	assert(fleq(vecC(2),bigB));

%!test
%!	setprngstates(89163680);
%!	bigA = randn();
%!	bigB = randn();
%!	bigX = randn()*exp(3*randn());
%!	bigP = 0.3 + abs(randn());
%!	numPts = 2 + round(abs(randn()*exp(3*randn())));
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = randn(1,numPts).*exp(randn()) + bigX;
%!	fVals = funchF(xVals);
%!	vecC = extFit__getCoeff( bigX, bigP, xVals, fVals );
%!	if (0)
%!		echo__numPts = numPts
%!		echo__xVals = xVals
%!		trueVals = [ bigA, bigB, bigX, bigP ]
%!		calcdVals = vecC'
%!	end
%!	assert(fleq(vecC(1),bigA));
%!	assert(fleq(vecC(2),bigB));

%!test
%!	setprngstates();
%!	bigA = randn();
%!	bigB = randn();
%!	bigX = randn()*exp(3*randn());
%!	bigP = 0.3 + abs(randn());
%!	numPts = 2 + round(abs(randn()*exp(3*randn())));
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = randn(1,numPts).*exp(randn()) + bigX;
%!	fVals = funchF(xVals);
%!	vecC = extFit__getCoeff( bigX, bigP, xVals, fVals );
%!	if (0)
%!		echo__numPts = numPts
%!		echo__xVals = xVals
%!		trueVals = [ bigA, bigB, bigX, bigP ]
%!		calcdVals = vecC'
%!	end
%!	assert(fleq(vecC(1),bigA));
%!	assert(fleq(vecC(2),bigB));
