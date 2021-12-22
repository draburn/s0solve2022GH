prc0522_3db_init;
vecFM = funchF(vecXM);
thisFile = "prc0522_3db";
if ( norm(vecFM) < 1e-8 )
	msg(thisFile,__LINE__,"Converged!");
	return;
end
msg(thisFile,__LINE__,"Hit bad local min!");
prc0522_3db_calc;
prc0522_3db_viz;
