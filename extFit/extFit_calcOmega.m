function [ omega, rho, bigA, bigB, bigC ] = extFit_calcOmega( xVals, fVals, bigX, bigP, wVals = [] );
	thisFile = "extFit_calcOmega";
	
	if (1)
	% New model.
	assert( bigP < 50.0 );
	%
	[ foo, nOfPtWiseMin ] = min( fVals );
	numPts = size(fVals,2);
	nOfPtWiseMin = median([ nOfPtWiseMin, 2, numPts-1 ]);
	fOfPtWiseMin = fVals(nOfPtWiseMin);
	for i1=1:size(bigX,1)
	for i2=1:size(bigX,2)
		this_bigX = bigX(i1,i2);
		this_bigP = bigP(i1,i2);
		this_bigD = abs( xVals - this_bigX ).^this_bigP;
		%
		vecX = xVals(nOfPtWiseMin-1:nOfPtWiseMin+1)';
		vecD = this_bigD(nOfPtWiseMin-1:nOfPtWiseMin+1)';
		vecF = fVals(nOfPtWiseMin-1:nOfPtWiseMin+1)';
		matM = [ ones(3,1), vecX, vecD ];
		if (0)
		msg( thisFile, __LINE__, "..." );
		msg( thisFile, __LINE__, sprintf( "  this_bigX = %g.", this_bigX ) );
		msg( thisFile, __LINE__, sprintf( "  this_bigP = %g.", this_bigP ) );
		abs( xVals - this_bigX )
		echo__vecD = vecD'
		echo__matM = matM
		msg( thisFile, __LINE__, "..." );
		assert( isrealarray(matM,[3,3]) );
		msg( thisFile, __LINE__, "..." );
		end
		vecC = matM\vecF;
		%msg( thisFile, __LINE__, "..." );
		if (0)
		rc = rcond(matM);
		if (rc<eps)
			msg( thisFile, __LINE__, sprintf( "rcond(matM) = %g...", rc ) );
			msg( thisFile, __LINE__, sprintf( "  this_bigX = %g.", this_bigX ) );
			msg( thisFile, __LINE__, sprintf( "  this_bigP = %g.", this_bigP ) );
			%echo__matM = matM
			%echo__vecX = vecX'
			%echo__vecD = vecD'
			%echo__vecF = vecF'
			%echo__vecC = vecC'
		end
		end
		this_bigA = vecC(1);
		this_bigB = vecC(2);
		this_bigC = vecC(3);
		this_rho = this_bigA + this_bigB*xVals + this_bigC*this_bigD - fVals;	
		%
		this_omega = 0.5*sum(this_rho.^2);
		omega(i1,i2) = this_omega;
		rho(i1,i2,:) = this_rho;
		bigA(i1,i2) = this_bigA;
		bigB(i1,i2) = this_bigB;
		bigC(i1,i2) = this_bigC;
	end
	end
	if ( isrealscalar(bigX) )
		rho = this_rho;
	end
	return
	end
	
	bigC = [];
	
	if (0)
	% Least squares on 3 closes pts.
	[ foo, nOfPtWiseMin ] = min( fVals );
	numPts = size(fVals,2);
	nOfPtWiseMin = median([ nOfPtWiseMin, 2, numPts-1 ]);
	fOfPtWiseMin = fVals(nOfPtWiseMin);
	for i1=1:size(bigX,1)
	for i2=1:size(bigX,2)
		this_bigX = bigX(i1,i2);
		this_bigP = bigP(i1,i2);
		this_bigD = abs( xVals - this_bigX ).^this_bigP;
		%
		vecD = this_bigD(nOfPtWiseMin-1:nOfPtWiseMin+1)';
		vecF = fVals(nOfPtWiseMin-1:nOfPtWiseMin+1)';
		vecC = [ ones(3,1), vecD ] \ vecF;
		this_bigA = vecC(1);
		this_bigB = vecC(2);
		%
		assert( isrealscalar(this_bigB) );
		assert( isrealscalar(this_bigA) );
		this_rho = this_bigA + this_bigB*this_bigD - fVals;
		this_omega = 0.5*sum(this_rho.^2);
		omega(i1,i2) = this_omega;
		rho(i1,i2,:) = this_rho;
		bigA(i1,i2) = this_bigA;
		bigB(i1,i2) = this_bigB;
	end
	end
	if ( isrealscalar(bigX) )
		rho = this_rho;
	end
	return;
	end
	
	
	if (1)
	% Hit C exactly, least squares for L&R.
	% May be wrong.
	% Does not seem to work very well.
	[ foo, nOfPtWiseMin ] = min( fVals );
	numPts = size(fVals,2);
	nOfPtWiseMin = median([ nOfPtWiseMin, 2, numPts-1 ]);
	fOfPtWiseMin = fVals(nOfPtWiseMin);
	for i1=1:size(bigX,1)
	for i2=1:size(bigX,2)
		this_bigX = bigX(i1,i2);
		this_bigP = bigP(i1,i2);
		this_bigD = abs( xVals - this_bigX ).^this_bigP;
		fooD = this_bigD([nOfPtWiseMin-1,nOfPtWiseMin+1]) - this_bigD(nOfPtWiseMin);
		fooF = fVals([nOfPtWiseMin-1,nOfPtWiseMin+1]) - fVals(nOfPtWiseMin);
		this_bigB = sum(fooF.*fooD) / sum(fooD.^2);
		this_bigA = fOfPtWiseMin - this_bigB*abs( xVals(nOfPtWiseMin) - this_bigX ).^this_bigP;
		assert( isrealscalar(this_bigB) );
		assert( isrealscalar(this_bigA) );
		this_rho = this_bigA + this_bigB*this_bigD - fVals;
		this_omega = 0.5*sum(this_rho.^2);
		omega(i1,i2) = this_omega;
		rho(i1,i2,:) = this_rho;
		bigA(i1,i2) = this_bigA;
		bigB(i1,i2) = this_bigB;
	end
	end
	if ( isrealscalar(bigX) )
		rho = this_rho;
	end
	return;
	end
	
	
	
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
	if (isempty(wVals))
		wVals = ones(1,numPts);
	end
	assert( isrealarray(wVals,[1,numPts]) );
	%
	if ( 1==size1 && 1==size2 )
		dVals = abs(xVals-bigX).^bigP;
		wdVals = wVals.*dVals;
		sigma1 = sum(wVals);
		sigmaF = sum(wVals.*fVals);
		sigmaD  = sum(wdVals);
		sigmaDD = sum(wdVals.*dVals);
		sigmaDF = sum(wdVals.*fVals);
		%
		denom = sigma1 * sigmaDD - (sigmaD.*sigmaD);
		bigA = ( sigmaDD * sigmaF  - sigmaD * sigmaDF ) ./ denom;
		bigB = ( sigma1  * sigmaDF - sigmaD * sigmaF  ) ./ denom;
		rho = bigA + bigB * dVals - fVals;
		omega = 0.5 * sum( wVals .* rho .* rho );
	return;
	end
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
