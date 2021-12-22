function matX = NEO_calcGradCurve( funchAntiG, vecX0, prm=[] )
	thisFile = "NEO_calcGradCurve";
	%
	matX = lsode( funchAntiG, vecX0, 1000.0*linspace(0.0,1.0,1001) );
	%
return;
end
