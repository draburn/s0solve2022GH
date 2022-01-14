function matJ = demoFunc0101_evalJaco( vecX, funcPrm )
	sizeX = funcPrm.sizeX;
	vecY = vecX - funcPrm.x0;
	matT2a = repmat( funcPrm.m2a * vecY, 1, sizeX );
	matT2b = repmat( funcPrm.m2b * vecY, 1, sizeX );
	matT3a = repmat( funcPrm.m3a * vecY, 1, sizeX );
	matT3b = repmat( funcPrm.m3b * vecY, 1, sizeX );
	matT3c = repmat( funcPrm.m3c * vecY, 1, sizeX );
	matJ = funcPrm.m1 ...
	  + ( funcPrm.m2a .* matT2b ) ...
	  + ( matT2a .* funcPrm.m2b ) ...
	  + ( funcPrm.m3a .* matT3b .* matT3c ) ...
	  + ( matT3a .* funcPrm.m3b .* matT3c ) ...
	  + ( matT3a .* matT3b .* funcPrm.m3c );
return;
end

%!test
%!	commondefs;
%!	sizeX = 1600;
%!	sizeF = 1000;
%!	seedPrm = demoFunc0101_genSeedPrm("easy");
%!	seedPrm.sizeX = sizeX;
%!	seedPrm.sizeF = sizeF;
%!	funcPrm = demoFunc0101_genFuncPrm(seedPrm);
%!	vecX0 = funcPrm.x0;
%!	matJ = demoFunc0101_evalJaco( vecX0, funcPrm );
%!	assert( isrealarray(matJ,[sizeF,sizeX]) );
%!	d = sqrt(sum(sum((matJ-funcPrm.m1).^2)));
%!	s = sqrt(sum(sum(matJ.^2))) + sqrt(sum(sum(funcPrm.m1)));
%!	assert( d < (eps^0.5)*sizeX*sizeF*s );


%!test
%!	commondefs;
%!	sizeX = 1234;
%!	sizeF = 555;
%!	seedPrm = demoFunc0101_genSeedPrm("moderate");
%!	seedPrm.sizeX = sizeX;
%!	seedPrm.sizeF = sizeF;
%!	funcPrm = demoFunc0101_genFuncPrm(seedPrm);
%!	vecX0 = zeros(sizeX,1);
%!	funchF = @(vecXDummy)( demoFunc0101_eval( vecXDummy, funcPrm ) );
%!	matJ1 = fdjaco( funchF, vecX0 );
%!	matJ2 = demoFunc0101_evalJaco( vecX0, funcPrm );
%!	assert( isrealarray(matJ2,[sizeF,sizeX]) );
%!	d = sqrt(sum(sum((matJ2-matJ1).^2)));
%!	s = sqrt(sum(sum(matJ2.^2))) + sqrt(sum(sum(matJ1)));
%!	assert( d < (eps^0.5)*sizeX*sizeF*s );
