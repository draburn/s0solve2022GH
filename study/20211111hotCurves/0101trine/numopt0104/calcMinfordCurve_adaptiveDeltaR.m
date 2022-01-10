function [ matX, datOut ] = calcMinfordCurve_adaptiveDeltaR( funchOmega_standalone, funchG_standalone, vecX0, matS=[], prm=[] )
	commondefs;
	thisFile = "calcMinfordCurve_adaptiveDeltaR";
	%msg( thisFile, __LINE__, "TO-DO..." );
	%msg( thisFile, __LINE__, "LATER..." );
	%msg( thisFile, __LINE__, " ~ BFGS? LSODE step? ~ for narrow valleys, esp finish." );
	%msg( thisFile, __LINE__, " ~ Obey prm." );
	%msg( thisFile, __LINE__, " ~ Review, refactor, add step limits, etc.?" );
	%msg( thisFile, __LINE__, " ~ Make funcG optional." );
	msg( thisFile, __LINE__, "This function may not be well-behaved. Please consider using calcMinfordCurve() instead." );
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
	desiredStepSize = 0.1;
	vecXC = vecX0;
	[ omega0, vecNablaOmega0 ] = funchOmega( vecX0 );
	tolMagNablaOmega = mygetfield( prm, "tolMagNablaOmega", eps025*norm(vecNablaOmega0) + eps050*sizeX );
	%
	bigR = 0.0;
	iterCount = 0;
	%
	vecX = vecX0;
	onSurf = false;
	matX(:,iterCount+1) = vecX;
	doFig = false;
	if (doFig)
		figure(100);
		hold off;
		plot( vecX0(1), vecX0(2), 'o-' );
		grid on;
		hold on;
	end
	while (1)
		%
		[ omega, vecNablaOmega ] = funchOmega( vecX );
		if ( norm(vecNablaOmega) < tolMagNablaOmega )
			msg( thisFile, __LINE__, "Reached tolMagNablaOmega." );
			return;
		end
		%
		iterCount++;
		if ( iterCount > iterLimit )
			msg( thisFile, __LINE__, "Reached iterLimit." );
			return;
		end
		%
		vecD_guess = -vecNablaOmega;
		vecX_guess = vecX + vecD_guess*(desiredStepSize/norm(vecD_guess));
		deltaR = norm(matS_nonEmpty*(vecX_guess-vecX0)) - norm(matS_nonEmpty*(vecX-vecX0));
		if ( 2 <= iterCount )
			%clear vecX_guess
			vecD_guess2 = vecX- matX(:,iterCount-1);
			vecX_guess2 = vecX + vecD_guess2*(desiredStepSize/norm(vecD_guess2));
			deltaR_2 = norm(matS_nonEmpty*(vecX_guess2-vecX0)) - norm(matS_nonEmpty*(vecX-vecX0));
			vecX_guess = vecX_guess2;
			deltaR = min([ deltaR, deltaR_2 ]);
		end
		assert( 0.0 < deltaR  );
		bigR += deltaR;
		%
		%
		if ( isempty(matS) )
			funchSurf = @(x)( funcSurf_ellip( x, bigR, vecXC ) );
		else
			funchSurf = @(x)( funcSurf_ellip( x, bigR, vecXC, matS ) );
		end
		%
		vecX = vecX_guess;
		%vecX_next = findObjSurfMin_simple( vecX, funchSurf, funchOmega );
		vecX_next = findObjSurfMin_onOff( vecX, funchSurf, funchOmega );
		if (doFig)
		if ( 2<= iterCount )
			plot( ...
			  [ matX(1,iterCount-1), matX(1,iterCount), vecX_guess(1), vecX_next(1) ], ...
			  [ matX(2,iterCount-1), matX(2,iterCount), vecX_guess(2), vecX_next(2) ], ...
			  '-' );
			plot( matX(1,iterCount-1), matX(2,iterCount-1), 's', 'markersize', 15 );
			plot( matX(1,iterCount), matX(2,iterCount), 'p', 'markersize', 15 );
			plot( vecX(1), vecX(2), '^', 'markersize', 15 );
			plot( vecX_next(1), vecX_next(2), 'x', 'markersize', 15 );
		end
		end
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
		vecX = vecX_next;
		matX(:,iterCount+1) = vecX_next;
	end
	%
return;
end
