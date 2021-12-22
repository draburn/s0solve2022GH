function vecXSurf = calcMinfordCurve__evalXSurf( vecXC, bigR, vecX, matS=[], prm=[] )
	%thisFile = "calcMinfordCurve__evalXSurf";
	%
	if ( mygetfield(prm,"doChecks",true) )
		sizeX = size(vecXC,1);
		assert( isrealarray(vecXC,[sizeX,1]) );
		assert( isrealscalar(bigR) );
		assert( bigR > 0.0 );
		assert( isrealarray(vecX,[sizeX,1]) );
		if (~isempty(matS))
			assert( isrealarray(matS,[sizeX,sizeX]) );
		end
	end
	%
	vecD = vecX - vecXC;
	normD = norm(vecD);
	if ( 0.0 == normD )
		vecXSurf = vecXC;
		return;
	end
	%
	if (isempty(matS))
		vecXSurf = vecXC + (bigR*vecD/normD);
	elseif (0)
		% Old/Bad version.
		vecXSurf = vecXC + (bigR*matS*vecD/normD);
	else
		normSD = norm(matS*vecD);
		vecXSurf = vecXC + (bigR*vecD/normSD);
	end
return;
end
