% Function...
%  [ vecXPts, datOut ] = calcLevCurve_alytG_minford( vecX0, funchOmegaG, prm=[] )
% Calculates points along the Minford (minimum for given distance) Levenberg curve for funcOmega using findOmegaMinWithinSurf().
% funchOmegaG must support the following interface:
%  [ omega, vecG ] = funchOmegaG( vecX )

function [ matX, datOut ] = calcLevCurve_alytG_minford( vecX0, funchOmegaG, prm=[] )
	%
	sizeX = size(vecX0,1);
	assert(isrealarray(vecX0,[sizeX,1]));
	%
	matS = mygetfield( prm, "matS", [] );
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
	iterLimit = 10000;
	deltaR = 0.05;
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
			msg( __FILE__, __LINE__, "Reached iterLimit." );
			return;
		end
		%
		bigR = bigR + deltaR;
		vecX = matX(:,iterCount);
		%
		%
		if ( isempty(matS) )
			funchSurf = @(x)( funcSurfEllip( x, vecXC, bigR ) );
		else
			funchSurf = @(x)( funcSurfEllip( x, vecXC, bigR, matS ) );
		end
		%
		vecX_next = findOmegaMinWithinSurf( vecX, funchSurf, funchOmegaG );
		%
		assert( isrealarray(vecX_next,[sizeX,1]) );
		s_next = norm(matS_nonEmpty*(vecX_next-vecXC));
		if ( s_next < bigR - (0.5*deltaR) )
			msg( __FILE__, __LINE__, "Reached local min." );
			return;
		end
		if ( s_next > bigR*(1.0-sqrt(eps)) )
			if ( s_next > bigR + (0.5*deltaR) )
				msg( __FILE__, __LINE__, "WARNING: Generated point lies well outside surface." );
			end
			vecX_next = vecXC + bigR*(vecX_next-vecXC)/s_next;
		end
		matX(:,iterCount+1) = vecX_next;
	end
	%
return;
end
