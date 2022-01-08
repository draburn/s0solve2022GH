function [ omega, vecNablaOmega ] = funchOmega_reCombo( vecX, funchOmega, funchNablaOmega )
	thisFile = "funchOmega_gradReCombo";
	omega = funchOmega( vecX );
	if ( 2 <= nargout )
		vecNablaOmega = funchNablaOmega( vecX );
	end
return;
end
