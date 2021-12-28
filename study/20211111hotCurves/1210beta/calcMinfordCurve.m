function [ matX, datOut ] = calcMinfordCurve( funchOmega, funchG, vecX0, matS=[], prm=[] )
	thisFile = "calcMinfordCurve";
	%msg( thisFile, __LINE__, "TO-DO..." );
	%msg( thisFile, __LINE__, "LATER..." );
	%msg( thisFile, __LINE__, " ~ BFGS? LSODE step? ~ for narrow valleys, esp finish." );
	%msg( thisFile, __LINE__, " ~ Obey prm." );
	%msg( thisFile, __LINE__, " ~ Review, refactor, add step limits, etc.?" );
	%msg( thisFile, __LINE__, " ~ Make funcG optional." );
	%
	sizeX = size(vecX0,1);
	assert(isrealarray(vecX0,[sizeX,1]));
	if (~isempty(matS))
		assert(isrealarray(matS,[sizeX,sizeX]));
		matS_nonEmpty = matS;
	else
		matS_nonEmpty = eye(sizeX,sizeX);
		% Use this for convenience, but, for efficiency,
		% let sub-modules know that matS is actually empty.
	end
	datOut = [];
	%
	%
	%
	iterLimit = 1000;
	deltaR = 0.1;
	vecXC = vecX0;
	%
	bigR = 0.0;
	iterCount = 0;
	%
	vecX = vecX0;
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
		switch (iterCount)
		case 1
			vecX = matX(:,iterCount);
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
			assert( 0.0 < step_u );
			vecX = matX(:,iterCount) + step_u* ( matX(:,iterCount) - matX(:,iterCount-1) );
			vecX = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, matS );
			onSurf = true;
			%
			%
			% This simple form seems to work better, at least some of the time.
			% But, I'll stick with the sophisticated version above until there's an explanation.
			%vecX = 2*matX(:,iterCount) - matX(:,iterCount-1);
			%onSurf = [];
		end
		%
		switch (3)
		case 1
		vecX_next = calcMinfordCurve__findNextPt( funchOmega, funchG, onSurf, vecX, vecXC, bigR, matS );
		case 2
		vecX_next = calcMinfordCurve__findNextPt_bfgs( funchOmega, funchG, onSurf, vecX, vecXC, bigR, matS );
		case 3
		vecX_next = calcMinfordCurve__findNextPt_mybfgs( funchOmega, funchG, onSurf, vecX, vecXC, bigR, matS );
		end
		%
		assert( isrealarray(vecX_next,[sizeX,1]) );
		if ( norm(matS_nonEmpty*(vecX_next-vecXC)) < bigR*(1.0-sqrt(eps)) )
			msg( thisFile, __LINE__, "Reached local min." );
			return;
		end
		matX(:,iterCount+1) = vecX_next;
	end
	%
return;
end
