function [ matX, datOut ] = calcOCQLevCurve( vecF0, matJ, vecEta, vecPhi, vecX0, prm=[] )
	thisFile = "calcOCQLevCurve";
	msg( thisFile, __LINE__, "WORK-IN-PROGRESS." );
	%
	% vecFModel( vecX0 + vecDelta ) = ...
	%  vecF0 + matJ * vecDelta + vecEta * ( vecPhi' * vecDelta )^2,
	% with vecPhi'*vecPhi = 1.
	%
	sizeF = size(vecF0,1);
	sizeX = size(vecX0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	assert( isrealarray(matJ,[sizeF,sizeX]) );
	assert( isrealarray(vecEta,[sizeF,1]) );
	assert( isrealarray(vecPhi,[sizeX,1]) );
	assert( abs(vecPhi'*vecPhi - 1.0) < eps*sizeX );
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	takePlus = true;
	% For production, instead of takePlus or minus,
	% we would probably take which produces a smaller xi.
	% But, this is good enough for now.
	%
	%
	% We want to expand: matJ = vecLambda * (vecPhi') + matU * matR * (matPsi'),
	% with matR invertible.
	% Note that vecPhi*(vecPhi') + matPsi*(matPsi') = matIX is NOT guaranteed,
	% since matJ could have a null space (beyond) vecPhi.
	% However, for now, we'll assume that's NOT the case,
	% and merely confirm that matU*matR is invertible.
	% Using matW == matU*matR, we'll just calc matW'*matW and ensure it's positive-definite.
	%
	matIF = eye(sizeF,sizeF);
	vecLambda = matJ*vecPhi;
	% Unfortunately, it seems neither the built-in orth nor myorth is designed for this.
	% So, here goes...
	orthoTol = eps*sizeX;
	sizeJ = 0;
	matPsi = zeros(sizeX,sizeJ);
	n = 0;
	while (1)
		n++;
		if ( n > sizeX )
			break;
		end
		%
		vecU = zeros(sizeX,1);
		vecU(n) = 1.0;
		%
		% Orthogonalize vecU against accepted vectors, and, if it hasn't vanished, add to matPsi.
		% Use two passes of whatchimacallit...
		vecU -= vecPhi * (vecPhi'*vecU);
		for k=1:sizeJ
			vecU -= matPsi(:,k) * (matPsi(:,k)'*vecU);
		end
		uNorm = norm(vecU);
		if ( uNorm < orthoTol )
			continue;
		end
		vecU /= uNorm;
		% Second pass..
		vecU -= vecPhi * (vecPhi'*vecU);
		for k=1:sizeJ
			vecU -= matPsi(:,k) * (matPsi(:,k)'*vecU);
		end
		uNorm = norm(vecU);
		if ( uNorm < orthoTol )
			continue;
		end
		vecU /= uNorm;
		%
		% Still good. Add to psi.
		matPsi(:,sizeJ+1) = vecU;
		sizeJ++;
	end
	if (0)
		msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
		echo__matPsi = matPsi
		echo__sizeJ = sizeJ
		msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
	end
	assert( sum(sum((matPsi'*matPsi-eye(sizeJ,sizeJ)).^2)) < eps*sizeJ^2 );
	%
	matIJ = eye(sizeJ,sizeJ);
	matW = matJ * matPsi;
	matH = matW' * matW;
	assert( size(matH,1) >= 1 );
	[ matR, cholFlag ] = chol( matH );
	if ( 0 ~= cholFlag )
		msg( thisFile, __LINE__, "ERROR: Part of the nullspace of matJ is in matPsi." );
		msg( thisFile, __LINE__, "This case is not (yet) supported." );
		error( "Part of the nullspace of matJ is in matPsi." );
	end
	if (0)
		msg( thisFile, __LINE__, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
		echo__matW = matW
		echo__matH = matH
		echo__matR = matR
		msg( thisFile, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
	end
	clear cholFlag;
	clear matR;
	%
	%
	%
	% Now that we have our bases...
	%   vecDelta = vecPhi * z + matPhi * vecY.
	%   xi = 0.5 * ( vecFModel' * vecFModel ) + 0.5 * mu * ( z^2 + vecY'*vecY ).
	%
	%
	matX(:,1) = vecX0;
	vecMu(1) = +Inf;
	vecZ(1) = 0.0;
	numPts = 1;
	mu = 1e5;
	breakIter = 0;
	breakLimit = 10000;
	while (1)
		breakIter++;
		if ( breakIter > breakLimit )
			msg( thisFile, __LINE__, "Hit break limit!" );
			return;
		end
		%
		matM = matH + mu*matIJ;
		%
		[ matR, cholFlag ] = chol( matM );
		if ( 0~=cholFlag )
			msg( thisFile, __LINE__, "matM is not positive definite. We'll see what happens." );
		end
		clear cholFlag;
		clear matR;
		%
		matMInv = inv(matM);
		%
		matA = matIF - ( matW * matMInv * (matW') );
		a = 2.0 * vecEta' * matA * vecEta;
		b = 3.0 * vecEta' * matA * vecLambda;
		c = mu + (2.0 * vecEta' * matA * vecF0) + (vecLambda' * matA * vecLambda);
		d = vecLambda' * matA * vecF0;
		funchRes = @(dummy)( (a*(dummy.^3)) + (b*(dummy.^2)) + (c*dummy) + d );
		funchDRes = @(dummy)( (3.0*a*(dummy.^2)) + (2.0*b*dummy) + c );
		funchD2Res = @(dummy)( (6.0*a*dummy) + (2.0*b) );
		%
		z0 = vecZ(:,numPts); % Our initial guess.
		z = z0;
		bigD = (4.0*(b^2)) - (12.0*a*c);
		if ( bigD <= 0.0 )
			% Only one solution.
			zLo = -Inf;
			zHi = +Inf;
		elseif ( takePlus )
			% Two solutions, take plus.
			msg( thisFile, __LINE__, "Model has two mins; taking plus." );
			zLo = ( -2.0*b + sqrt(bigD) ) / ( 6.0*a );
			zHi = +Inf;
			if ( z<zLo )
				z = zLo + 1.0;
			end
		else
			% Two solutions, take minus.
			msg( thisFile, __LINE__, "Model has two mins; taking minus." );
			zLo = -Inf;
			zHi = ( -2.0*b - sqrt(bigD) ) / ( 6.0*a );
			if ( z>zHi )
				z = zHi + 1.0;
			end
		end
		%
		zRes = (a*(z^3)) + (b*(z^2)) + (c*z) + d;
		zResTol = (abs(a)+abs(b)+abs(c)+abs(d))*(eps^0.75);
		while ( abs(zRes) > zResTol )
			breakIter++;
			if ( breakIter > breakLimit )
				msg( thisFile, __LINE__, "Hit break limit!" );
				return;
			end
			z -= zRes / funchDRes(z);
			zRes = (a*(z^3)) + (b*(z^2)) + (c*z) + d;
		end
		echo__z = z;
		vecZeta = vecF0 + (vecLambda*z) + (vecEta*(z^2));
		vecY = -matMInv*(matW'*vecZeta);
		vecDelta = (vecPhi*z) + (matPsi*vecY);
		vecX = vecX0 + vecDelta;
		%
		numPts++;
		matX(:,numPts) = vecX;
		vecMu(numPts) = mu;
		%
		vecA(numPts) = a;
		vecB(numPts) = b;
		vecC(numPts) = c;
		vecD(numPts) = d;
		vecBigD(numPts) = bigD;
		vecZ0(numPts) = z0;
		vecZ(numPts) = z;
		matZeta(:,numPts) = vecZeta;
		matY(:,numPts) = vecY;
		%
		if ( 0.0 == mu )
			break;
			msg( thisFile, __LINE__, "Reached mu = 0.0." );
		end
		mu /= 1.1;
		if ( mu < 5e-6 )
			mu = 0.0;
		end
	end
	%
	datOut.vecF0 = vecF0;
	datOut.matJ = matJ;
	datOut.vecEta = vecEta;
	datOut.vecPhi = vecPhi;
	datOut.vecX0 = vecX0;
	%
	datOut.vecLambda = vecLambda;
	datOut.matPsi = matPsi;
	datOut.matIF = matIF;
	datOut.matIJ = matIJ;
	datOut.matW = matW;
	datOut.matH = matH;
	%
	datOut.matX = matX;
	datOut.vecMu = vecMu;
	%
	datOut.vecA = vecA;
	datOut.vecB = vecB;
	datOut.vecC = vecC;
	datOut.vecD = vecD;
	datOut.vecBigD = vecBigD;
	datOut.vecZ0 = vecZ0;
	%
	datOut.vecZ = vecZ;
	datOut.matZeta = matZeta;
	datOut.matY = matY;
	%
return;
end
