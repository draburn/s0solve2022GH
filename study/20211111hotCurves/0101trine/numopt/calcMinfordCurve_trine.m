function [ matX, datOut ] = calcMinfordCurve_trine( funchOmega, funchG, vecX0, matS=[], prm=[] )
	thisFile = "calcMinfordCurve_trine";
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
	function [ l, g ] = octopt_fcn( x, dummy_funchBigL, dummy_funchVecG )
		l = dummy_funchBigL(x);
		if (nargout>=2)
			g = dummy_funchVecG(x);
		end
	end
	function q = octopt_inequc( x, r, vecXC )
		q = r - norm(x-vecXC);
		%q = r + (x(1)-vecXC(1));
		%q = r^2 - (x-vecXC)'*(x-vecXC)
		%q = r - (x(1)-vecXC(1));
		%q = x(1) - r
		%q = r - x(1)
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
		switch 3
		case 1
			octopt_fcn_at = @(x)( octopt_fcn( x, funchOmega, funchG ) );
			% What I want. Doesn't work.
		case 2
			octopt_fcn_at = @(x)( funchOmega( x ) );
			% Also doesn't work.
		case 3
			octopt_fcn_at = @(x)( 0.5*(x'*x) );
			% Does work.
		end
		if ( isempty(matS) )
			octopt_inequc_at = @(x)( octopt_inequc( x, bigR, vecXC ) );
		else
			error( "matS support not implemented." );
		end
		%opts = optimset( 'GradObj', 'on', 'TolX', 1e-1, 'TolFun', 1e-1, 'inequc', {octopt_inequc_at} );
		opts = optimset( 'GradObj', 'on', 'inequc', {octopt_inequc_at} );
		%opts = optimset( 'inequc', {octopt_inequc_at} );
		vecX_next = nonlin_min( octopt_fcn_at, vecX, opts );
		%
		assert( isrealarray(vecX_next,[sizeX,1]) );
		matX(:,iterCount+1) = vecX_next;
		s_next = norm(matS_nonEmpty*(vecX_next-vecXC));
		if ( s_next < bigR - (0.5*deltaR) )
			msg( thisFile, __LINE__, "Reached local min." );
			return;
		end
	end
	%
return;
end
