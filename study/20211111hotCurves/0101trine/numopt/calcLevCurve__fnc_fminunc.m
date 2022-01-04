function [ bigL, vecG ] = calcLevCurve__fnc_fminunc( vecX, funchBigL, funchVecG, s, vecX0, matS=[] )
	if ( isempty(matS) )
		bigL = s*funchBigL(vecX) + 0.5*(1.0-s)*(vecX-vecX0)'*(vecX-vecX0);
		vecG = s*funchVecG(vecX) + (1.0-s)*(vecX-vecX0);
	else
		error( "Support for matS is not implemented." );
	end
end
