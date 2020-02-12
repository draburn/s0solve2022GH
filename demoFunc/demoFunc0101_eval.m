function vecF = demoFunc0101_eval( vecX, funcPrm )
	vecY = vecX - funcPrm.x0;
	vecF = (funcPrm.m1 * vecY) ...
	  + ( (funcPrm.m2a*vecY) .* (funcPrm.m2b*vecY) ) ...
	  + ( (funcPrm.m3a*vecY) .* (funcPrm.m3b*vecY) .* (funcPrm.m3c*vecY) );
return;
end

%!test
%!	sizeX = 100;
%!	vecX0 = zeros(sizeX,1);
%!	seedPrm = demoFunc0101_genSeedPrm("easy");
%!	seedPrm.sizeX = sizeX;
%!	seedPrm.sizeF = sizeX;
%!	funcPrm = demoFunc0101_genFuncPrm(seedPrm);
%!	vecF0 = demoFunc0101_eval( vecX0, funcPrm );
%!	vecX4 = funcPrm.x0+(1e-4);
%!	vecX8 = funcPrm.x0+(1e-8);
%!	vecF4 = demoFunc0101_eval( vecX4, funcPrm );
%!	vecF8 = demoFunc0101_eval( vecX8, funcPrm );
%!	vecFSolu = demoFunc0101_eval( funcPrm.x0, funcPrm );
%!	res0 = sqrt(sum(vecF0.^2))
%!	res4 = sqrt(sum(vecF4.^2))
%!	res8 = sqrt(sum(vecF8.^2))
%!	resSolu = sqrt(sum(vecFSolu.^2))
