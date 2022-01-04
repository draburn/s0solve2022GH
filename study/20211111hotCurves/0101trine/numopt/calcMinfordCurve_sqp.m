function [ matX, datOut ] = calcMinfordCurve_sqp( funchOmega, funchG, vecX0, matS=[], prm=[] )
	thisFile = "calcMinfordCurve_sqpm";
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
		%matS_nonEmpty = matS;
	else
		%matS_nonEmpty = eye(sizeX,sizeX);
	end
	datOut = [];
	hessB = -eye(sizeX,sizeX);
	%
	%
	function [ l, g ] = sqp_phi( x, dummy_funchBigL, dummy_funchVecG )
		l = dummy_funchBigL(x);
		if (nargout>=2)
			g = dummy_funchVecG(x);
		end
	end
	function [ b, gradB, hessB ] = sqp_b_sansScale( x, r, vecXC )
		bScale = 1.0;
		b = bScale*( r^2 - 0.5*(x-vecXC)'*(x-vecXC) );
		if (nargout>=2)
			gradB = -bScale*(x-vecXC);
			if ( nargout>=3 )
				sizeX = size(vecXC,1);
				hessB = -bScale*eye(sizeX,sizeX);
			end
		end
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
		vecX = matX(:,iterCount);
		%
		if (0)
		s = norm(matS_nonEmpty*(vecX-vecXC));
		if ( s >= bigR )
			msg( thisFile, __LINE__, "WARNING: Starting point lies outside surface; pulling to surface." );
			echo__s_minus_bigR = s - bigR
			vecX = vecXC + bigR*(vecX-vecXC)/s;
		end
		switch (iterCount)
		case 1
			% vecX is fine as-is.
			onSurf = false;
		otherwise
			if (0)
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
			else
			% This simple form seems to work better, at least some of the time.
			% But, I'll stick with the sophisticated version above until there's an explanation.
			vecX = 2*matX(:,iterCount) - matX(:,iterCount-1);
			%onSurf = [];
			end
		end
		end
		%
		sqp_phi_at = @(x)( sqp_phi( x, funchOmega, funchG ) );
		sqp_b_at = @(x)( sqp_b_sansScale( x, bigR, vecXC ) );
		%
		%%%vecX_next = sqp( vecX0, sqp_phi_at, [], sqp_b_at )
		lb = []
		ub = []
		maxiter = 1000
		tol = eps^0.75
		[x, obj, info, iter, nf, lambda] = sqp( vecX0, sqp_phi_at, [], sqp_b_at, lb, ub, maxiter, tol )
		return
		%
		assert( isrealarray(vecX_next,[sizeX,1]) );
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
		matX(:,iterCount+1) = vecX_next;
		if ( s_next < bigR - (0.5*deltaR) )
			msg( thisFile, __LINE__, "Reached local min." );
			return;
		end
	end
	%
return;
end
