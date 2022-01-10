function [ matX, datOut ] = calcMinfordCurve( funchOmega_standalone, funchG_standalone, vecX0, matS=[], prm=[] )
	thisFile = "calcMinfordCurve";
	commondefs;
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
	funchOmega = @(x)( funchOmega_reCombo( x, funchOmega_standalone, funchG_standalone ) );
	%
	%
	%
	iterLimit = 10000;
	deltaR = 0.1;
	vecXC = vecX0;
	%[ omega0, vecNablaOmega0 ] = funchOmega( vecX0 );
	%tolMagNablaOmega = mygetfield( prm, "tolMagNablaOmega", eps025*norm(vecNablaOmega0) + eps050*sizeX );
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
		%
		if ( isempty(matS) )
			funchSurf = @(x)( funcSurf_ellip( x, bigR, vecXC ) );
		else
			funchSurf = @(x)( funcSurf_ellip( x, bigR, vecXC, matS ) );
		end
		%
		vecX_next = findObjSurfMin_simple( vecX, funchSurf, funchOmega );
		%vecX_next = findObjSurfMin_onOff( vecX, funchSurf, funchOmega );
		%
		assert( isrealarray(vecX_next,[sizeX,1]) );
		s_next = norm(matS_nonEmpty*(vecX_next-vecXC));
		if ( s_next < bigR - (0.5*deltaR) )
			msg( thisFile, __LINE__, "Reached local min." );
			return;
		end
		if ( s_next > bigR*(1.0-sqrt(eps)) )
			if ( s_next > bigR + (0.5*deltaR) )
				msg( thisFile, __LINE__, "WARNING: Generated point lies well outside surface." );
			end
			vecX_next = vecXC + bigR*(vecX_next-vecXC)/s_next;
		end
		matX(:,iterCount+1) = vecX_next;
	end
	%
return;
end
