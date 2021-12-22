function matX = NEO_calcGradCurveX( funcPrm, vecX0, prm=[] )
	thisFile = "NEO_calcGradCurveX";
	funchAntiG = @(x)( -testFunc_gForLSODE(x,funcPrm) );
	matX = lsode( funchAntiG, vecX0, 1000.0*linspace(0.0,1.0,1001) );
return;
end
