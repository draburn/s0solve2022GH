function [ matX, datOut ] = calcMinfordCurve_fminsurf( funchOmega, funchG, vecX0, matS=[], prm=[] )
	thisFile = "calcMinfordCurve_fminsurf";
	%
	sizeX = size(vecX0,1);
	assert(isrealarray(vecX0,[sizeX,1]));
	if (~isempty(matS))
		assert(isrealarray(matS,[sizeX,sizeX]));
		matS_nonEmpty = matS;
	else
		matS_nonEmpty = eye(sizeX,sizeX);
	end
	datOut = [];
	%
	function [ vecS, matDST ] = funchVecS( vecX, vecXC, bigR )
		vecD = vecX-vecXC;
		magD = norm(vecD);
		if (0.0==magD)
			vecD(1) = bigR;
			magD = bigR;
		end
		sizeX = size(vecX,1);
		vecS = vecXC + bigR*vecD/magD;
		matDST = ((bigR/magD)*eye(sizeX,sizeX)) - (bigR/(magD^3))*(vecD*vecD');
	end
	%
	function [ bigL, vecDL ] = funchBigL( vecX, funchOmega, funchG )
		bigL = funchOmega(vecX);
		vecDL = funchG(vecX);
	end
	%
	function [ bigP, vecDP ] = funchBigP( vecD )
		bigP = 0.5*(vecD'*vecD);
		vecDP = vecD;
	end
	%
	function [ bigF, vecDF ] = funchBigF( vecX, funchVecS, funchBigL, funchBigP )
		[ vecS, matDST ] = funchVecS(vecX);
		[ bigL, vecDL ] = funchBigL(vecS);
		[ bigP, vecDP ] = funchBigP(vecX-vecS);
		bigF = bigL + bigP;
		vecDF = (matDST*(vecDL-vecDP)) + vecDP;
	end
	%
	%
	iterLimit = 1000;
	deltaR = 0.1;
	vecXC = vecX0;
	%
	bigR = 0.0;
	iterCount = 0;
	%
	%vecX = vecX0; % Don't start at center of surface!
	vecG0 = funchG(vecX0);
	magG0 = norm(vecG0);
	if ( 0.0==magG0 )
		msg( thisFile, __LINE__, "Initial gradient is zero." );
		vecX = vecX0;
		return;
	end
	vecX = vecX0 + (0.1*deltaR)*vecG0/magG0;
	%
	onSurf = false;
	matX(:,iterCount+1) = vecX;
	while (1)
		iterCount++;
		if ( iterCount > iterLimit )
			msg( thisFile, __LINE__, "Reached iterLimit." );
			return;
		end
		%
		bigR = bigR + deltaR;
		vecX = matX(:,iterCount);
		%
		usePredictiveGuess = false;
		if (usePredictiveGuess)
		s = norm(matS_nonEmpty*(vecX-vecXC));
		if ( s >= bigR )
			msg( thisFile, __LINE__, "WARNING: Starting point lies outside surface; pulling to surface." );
			vecX = vecXC + bigR*(vecX-vecXC)/s;
		end
		switch (iterCount)
		case 1
			% vecX is fine as-is.
			onSurf = false;
		otherwise
			% Allow variable deltaR ~ find actual intersection with new surface.
			%vecX = 2*matX(:,iterCount) - matX(:,iterCount-1);
			%onSurf = [];
			% I had considered trying to extrapolate per a circle, but, I now think that wouldn't be so helpful.
			step_vecSXMXP = matS_nonEmpty * ( matX(:,iterCount) - matX(:,iterCount-1) );
			step_vecSXMXC = matS_nonEmpty * ( matX(:,iterCount) - vecXC );
			step_c2 = step_vecSXMXP'*step_vecSXMXP;
			step_c1 = 2.0*step_vecSXMXP'*step_vecSXMXC;
			step_c0 = step_vecSXMXC'*step_vecSXMXC - bigR^2;
			%
			% getQuadGoodPt() here is wrong;
			%  we always want the "forward" solution,
			%  and we have no guarantee about the sign of c1.
			%step_u = getQuadGoodPt( step_c2, step_c1, step_c0 );
			% However, the problem is simple
			%  since there should be exactly one positive solution...
			step_d = (step_c1^2) - (4.0*step_c0*step_c2);
			assert( step_d > -eps*( (step_c1^2) + abs(step_c0*step_c2) ) );
			if ( step_d > 0.0 )
				step_sqrtD = sqrt(step_d);
			else
				step_sqrtD = 0.0;
			end
			step_u = ( -step_c1 + step_sqrtD ) / ( 2.0*step_c2 );
			assert( 0.0 < step_u || 0.0 < step_c0 ); % step_u should be positive, unless we start outside the surface.
			vecX = matX(:,iterCount) + step_u * ( matX(:,iterCount) - matX(:,iterCount-1) );
			vecX = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, matS );
			onSurf = true;
			%
			%
			% This simple form seems to work better, at least some of the time.
			% But, I'll stick with the sophisticated version above until there's an explanation.
			%vecX = 2*matX(:,iterCount) - matX(:,iterCount-1);
			%onSurf = [];
		end
		end
		%
		funchVecS_at = @(x)( funchVecS( x, vecX0, bigR ) );
		funchBigL_at = @(x)( funchBigL( x, funchOmega, funchG ) );
		funchBigP_at = @(x)( funchBigP( x ) );
		funchBigF_at = @(x)( funchBigF( x, funchVecS_at, funchBigL_at, funchBigP_at ) );
		%
		if (0)
			echo__vecX0 = vecX0
			[ vecS, matDST ] = funchVecS_at( vecX0 )
			[ bigL, vecDL ] = funchBigL_at( vecS )
			[ bigP, vecDP ] = funchBigP_at( vecX0 - vecX )
			[ bigF, vecDF ] = funchBigF_at( vecX0 )
		end
		%
		switch 1
		case 1
			opts = optimset( 'GradObj', 'on' );
			vecX_next = fminunc( funchBigF_at, vecX0, opts );
		case 2
			vecX_next = fminunc( funchBigF_at, vecX0 );
		case 3
			opts = optimset( 'GradObj', 'on', 'TolX', 1e-3, 'TolFun', 1e-5 );
			vecX_next = fminunc( funchBigF_at, vecX0, opts );
		end
		%
		if (0)
			echo__vecX_next = vecX_next
			[ vecS, matDST ] = funchVecS_at( vecX_next )
			[ bigL, vecDL ] = funchBigL_at( vecS )
			[ bigP, vecDP ] = funchBigP_at( vecX_next - vecX )
			[ bigF, vecDF ] = funchBigF_at( vecX_next )
		end
		%
		assert( isrealarray(vecX_next,[sizeX,1]) );
		if ( funchG(vecX_next)'*(vecX_next-vecX0) > 0 )
			msg( thisFile, __LINE__, "Fell off surface." );
			opts = optimset( 'GradObj', 'on' );
			vecX_next = fminunc( funchBigL_at, vecX_next, opts );
			assert( isrealarray(vecX_next,[sizeX,1]) );
			s_next = norm(matS_nonEmpty*(vecX_next-vecXC));
		end
		matX(:,iterCount+1) = vecX_next;
		%s_next = norm(matS_nonEmpty*(vecX_next-vecXC));
		s_next = norm(vecX_next-vecXC);
		if ( s_next > bigR*(1.0-sqrt(eps)) )
			if ( s_next > bigR + (0.5*deltaR) )
				msg( thisFile, __LINE__, "WARNING: Generated point lies well outside surface." );
				echo__iterCount = iterCount
				echo__bigR = bigR
				return
				vecX_next = vecXC + bigR*(vecX_next-vecXC)/s_next;
			end
		end
		if ( s_next < bigR - (0.5*deltaR) )
			msg( thisFile, __LINE__, "Reached local min." );
			return;
		end
	end
	%
return;
end
