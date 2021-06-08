function [ omega, rho, bigA, bigB ] = extFit_calcOmega( xVals, fVals, bigX, bigP, prm = [] );
	thisFile = "extFit_calcOmega";
	%
	% The logic:
	%  d_n = abs( x_n - X )^P
	%  sigma_1  = sum_n w_n
	%  sigma_d  = sum_n w_n * d_n
	%  sigma_dd = sum_n w_n * d_n^2
	%  sigma_f  = sum_n w_n * f_n
	%  sigma_fd = sum_n w_n * f_n * d_n
	%  [ A; B ] = [ sigma_1, sigma_d; sigma_d, sigma_dd ] \ [ sigma_f; sigma_df ]
	%  rho_n = A + B*d_n - f_n
	%  omega = 0.5 * ( sum_n w_n * rho_n^2 ) / sum_1
	%
	% But, do this in a way that is quick for a large number of X and P values.
	%
	numPts = size(xVals,2);
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	size1 = size(bigX,1);
	size2 = size(bigX,2);
	assert( isrealarray(bigX,[size1,size2]) );
	assert( isrealarray(bigP,[size1,size2]) );
	wVals = mygetfield( prm, "wVals", [] );
	if (isempty(wVals))
		wVals = ones(1,numPts);
	end
	assert( isrealarray(wVals,[1,numPts]) );
	%
	ary3D   = zeros(size1,size2,numPts);
	ary3WD  = zeros(size1,size2,numPts);
	ary3WDD = zeros(size1,size2,numPts);
	ary3WDF = zeros(size1,size2,numPts);
	parfor n=1:numPts
		ary3D(:,:,n) = abs( xVals(n) - bigX ).^bigP;
		ary3WD(:,:,n) = wVals(n)*ary3D(:,:,n);
		ary3WDD(:,:,n) = ary3WD(:,:,n).*ary3D(:,:,n);
		ary3WDF(:,:,n) = ary3WD(:,:,n)*fVals(n);
	end
	wfVals = wVals.*fVals;
	sigma1 = sum(wVals);
	sigmaF = sum(fVals);
	matSigmaD  = sum( ary3WD,  3 );
	matSigmaDD = sum( ary3WDD, 3 );
	matSigmaDF = sum( ary3WDF, 3 );
	%
	matDenom = sigma1 * matSigmaDD - (matSigmaD).^2;
	bigA = ( matSigmaDD * sigmaF - matSigmaD .* matSigmaDF ) ./ matDenom;
	bigB = ( sigma1 * matSigmaDF - matSigmaD * sigmaF ) ./ matDenom;
	%
	rho = zeros(size1,size2,numPts);
	parfor n=1:numPts
		rho(:,:,n) = bigA + bigB .* ary3D(:,:,n) - fVals(n);
		wRho(:,:,n) = wVals(n) * rho(:,:,n);
		wRhoSq(:,:,n) = wRho(:,:,n) .* rho(:,:,n);
	end
	omega = 0.5 * sum( wRhoSq, 3 );
	%
return;
end
