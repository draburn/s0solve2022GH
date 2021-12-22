function matX = NEO_calcGradCurve1128( funchG, vecX0, prm=[] )
	thisFile = "NEO_calcGradCurve1128";
	%
	funchAntiG = @(x)( -funchG(x) );
	matX = lsode( funchAntiG, vecX0, 100.0*linspace(0.0,1.0,101) )';
	%
return;
end
