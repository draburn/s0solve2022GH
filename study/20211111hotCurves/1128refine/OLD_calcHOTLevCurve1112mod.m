function matY = OLD_calcHOTLevCurve1112( modelFuncPrm, vecX0, prm=[] )
	thisFile = "OLD_calcHOTLevCurve1112";
	msg( thisFile, __LINE__, "This is a hacked version." );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	%
	numSteps = mygetfield( prm, "numSteps", 100 );
	assert( isposintscalar(numSteps) );
	finalStepEps = mygetfield( prm, "finalStepEps", 0.5 );
	assert( isrealscalar(finalStepEps) );
	assert( 0.0 <= finalStepEps );
	%
	matI = eye(sizeX,sizeX);
	n = 1;
	vecX = vecX0;
	while (1)
		matY(:,n) = vecX;
		n++;
		if ( n > numSteps )
			break;
		end
		s = n/(numSteps+finalStepEps);
		[ omega, vecF, matJ, vecG, matH ] = testFunc_evalDeriv( vecX, modelFuncPrm );
		vecG = s*vecG + (1.0-s)*(vecX-vecX0);
		matH = s*matH + (1.0-s)*matI;
		[ matR, cholFlag ] = chol( matH );
		if ( 0 ~= cholFlag )
			msg( thisFile, __LINE__, "Hessian is non-positive-definite." );
			%error( "Hessian is non-positive-definite." );
			%return;
			%
			% This is a test/hack. Does not work.
			[ matPsi, matLam ] = eig( matH );
			[ lamMin, indexOfLamMin ] = min(diag(matLam));
			lamAbsMax = max(abs(diag(matLam)));
			vecPsi = matPsi(:,indexOfLamMin);
			vecDelta = -(vecPsi'*vecG)*vecPsi;
			deltaNorm = norm(vecDelta);
			if (0.0 == deltaNorm)
				msg( thisFile, __LINE__, "Eigenvector is non-gradiental." );
				return
			end
			vecDelta *= 0.01/deltaNorm;
			%
			% This is another test/hack.
			%msg( thisFile, __LINE__, "Using a test/hack." );
			%thisMu = 0.1*eigs(matH,1,'lm');
			%thisMu = 0.1*lamAbsMax;%+abs(lamMin);
			%matH2 = matH + thisMu* matI;
			%vecDelta = -(matH2\vecG);
		else
			vecDelta = -( matR \ (matR'\vecG) );
		end
		maxStepSize = 0.5;
		if (norm(vecDelta)>maxStepSize)
			vecDelta *= maxStepSize/norm(vecDelta);
		end
		vecX += vecDelta;
	end
	%
return;
end
