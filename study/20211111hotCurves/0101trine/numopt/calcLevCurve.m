function [ matX, datOut ] = calcLevCurve( funchOmega, funchG, vecX0, prm=[] )
	thisFile = "calcLevCurve";
	msg( thisFile, __LINE__, "THIS IS A WORK-IN-PROGRESS!" );
	%
	deltaNormTol = 1e-6;
	vNormTol = 1e-6;
	n = 1;
	vecX = vecX0;
	matX(:,n) = vecX0;
	numSVals = 201;
	vecS = linspace( 0.0, 1.0, numSVals );
	for n=2:numSVals
		s = vecS(n);
		%
		vecX_prev = vecX;
		switch (5)
		case (1)
		funchXDotLSODE = @(x,t)( -s*funchG(x) -(1.0-s)*(x-vecX0) );
		xLSODE = lsode( funchXDotLSODE, vecX', [ 0.0, 10.0 ] );
		vecX = xLSODE(end,:)';
		case (2)
		%funchGHatIsh = @(g)( g/sqrt( deltaNormTol^2 + g'*g ) );
		funchGHatIsh = @(g)( g/sqrt( deltaNormTol + norm(g) ) );
		funchXDotLSODE = @(x,t)( -funchGHatIsh(s*funchG(x)+(1.0-s)*(x-vecX0)) );
		xLSODE = lsode( funchXDotLSODE, vecX, [ 0.0, 2.0 ] );
		vecX = xLSODE(end,:);
		case (3)
			funchXi = @(x)( s*funchOmega(x) + (1.0-s)*0.5*(x-vecX0)'*(x-vecX0) );
			for m=1:10
				xi0 = funchXi(vecX);
				assert( isrealscalar(xi0) );
				assert( 0.0 <= xi0 );
				%
				vecV = -s*funchG(vecX) -(1.0-s)*(vecX-vecX0);
				if ( norm(vecV) <= vNormTol )
					break;
				end
				epsFD = (0.1*deltaNormTol)/norm(vecV);
				%
				vecXP = vecX + (vecV*epsFD);
				xiP = funchXi(vecXP);
				assert( isrealscalar(xiP) );
				assert( 0.0 <= xiP );
				dXi = -vecV'*vecV;
				% xiModel = xi0 + dXi*t + c*(t^2);
				c = ( xiP - xi0 - (dXi*epsFD) ) / (epsFD^2);
				%
				deltaT = calcLinishRootOfQuad( c, dXi, xi0 );
				assert( 0.0 < deltaT );
				%
				vecX += vecV*deltaT;
				if ( norm(vecV*deltaT) <= deltaNormTol )
					break;
				end
				%
				% Check that omega has decreased?
			end
			%if ( norm(vecV) > vNormTol )
			%	msg( thisFile, __LINE__, "Failed to reach vNormTol." );
			%	%echo__s = s
			%	%return
			%end
		case 4
			funchBigLMod = @(x)( s*funchOmega(x) + 0.5*(1.0-s)*(x-vecX0)'*(x-vecX0) );
			funchVecGMod = @(x)( s*funchG(x) + (1.0-s)*(x-vecX0) );
			vecX = findLocMinLossG_bfgs( funchBigLMod, funchVecGMod, vecX );
		case 5
			fcn = @(x)( calcLevCurve__fnc_fminunc( x, funchOmega, funchG, s, vecX0 ) );
			%[ bigL, vecG ] = fcn(vecX0)
			opts = optimset( 'GradObj', 'on' );
			%opts = optimset( 'GradObj', 'on', 'TolX', 1e-4 );
			vecX = fminunc( fcn, vecX, opts );
		otherwise
		error( "Invalid case." );
		end
		matX(:,n) = vecX;
		%
		%vecX = (2.0*vecX) - (vecX_prev);
	end
	%
return;
end
