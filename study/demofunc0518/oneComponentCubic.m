function matF = oneComponentCubic( matX, funcPrm )
	rvecVTX = funcPrm.vecV' * matX;
	matF = ( funcPrm.matM * matX ) ...
	  + ( funcPrm.cQuad * funcPrm.vecA * (rvecVTX.^2) ) + ( funcPrm.vecA * (rvecVTX.^3) );
return;
end
