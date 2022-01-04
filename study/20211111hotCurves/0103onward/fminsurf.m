function [ vecXF, retCode, datOut ] = fminsurf( funchBigL, funchVecS, funchBigP, vecX0, prm=[] )
	%
	function [ bigF, vecDF ] = fminsurf_internal( vecX, funchBigL, funchVecS, funchBigP )
		[ vecS, matDST ] = funchVecS(vecX);
		[ bigL, vecDL ] = funchBigL(vecS);
		[ bigP, vecDP ] = funchBigP(vecX-vecS);
		bigF = bigL + bigP;
		vecDF = (matDST*(vecDL-vecDP)) + vecDP;
	return;
	end
	%
	if (0)
		echo__vecX0 = vecX0
		[ vecS, matDST ] = funchVecS(vecX0)
		[ bigL, vecDL ] = funchBigL(vecS)
		[ bigP, vecDP ] = funchBigP(vecX0-vecS)
		[ bigF, vecDF ] = fminsurf_internal( vecX0, funchBigL, funchVecS, funchBigP )
	end
	%
	fcn = @(x)fminsurf_internal( x, funchBigL, funchVecS, funchBigP );
	switch 1
	case 1
		opts = optimset( 'GradObj', 'on' );
		vecXF = fminunc( fcn, vecX0, opts );
	case 2
		vecXF = fminunc( fcn, vecX0 );
	end
	%
	if (0)
		echo__vecXF = vecXF
		[ vecS, matDST ] = funchVecS(vecXF)
		[ bigL, vecDL ] = funchBigL(vecS)
		[ bigP, vecDP ] = funchBigP(vecXF-vecS)
		[ bigF, vecDF ] = fminsurf_internal( vecXF, funchBigL, funchVecS, funchBigP )
	end
	%
	retCode = [];
	datOut = [];
	return;
return
end;
