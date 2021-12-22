function matX = calcGradDescentCurve( funchGradDescent, vecX0, prm=[] )
	thisFile = "calcGradDescentCurve";
	msg( thisFile, __LINE__, "~~~ FORSAKEN! ~~~." );
	%
	stepSize = 0.01;
	vecX = vecX0;
	matX(:,1) = vecX0;
	n = 0;
	while (1)
		n++;
		if ( n > 10000 )
			break;
		end
		%
		vecGD1 = funchGradDescent( vecX );
		if ( 0.0 == norm(vecGD1) )
			break;
		end
		vecX2 = vecX + stepSize*vecGD1/norm(vecGD1);
		vecGD2 = funchGradDescent(vecX2);
		%
		vecDelta = (vecGD1/norm(vecGD1)) + (vecGD2/norm(vecGD2));
		if ( norm(vecDelta) < 0.5 )
			% I know this isn't the proper stopping criteria.
			% But, it's decent for a hack.
			break;
		end
		vecX += stepSize*( vecDelta / norm(vecDelta) );
		matX(:,n) = vecX;
		%
	end
	%
return;
end
