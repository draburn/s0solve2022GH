function [ vecX_next, vecF_next, stepGen_datOut ] = findZero_baseline__stepGen( vecX, vecF, funchXTrialOfP, funchFModel, funchF, stepGen_prm )
	stepGen_datOut = [];
	stepGen_datOut.fevalCount = 0;
	vecX_next = funchXTrialOfP(1.0);
	vecF_next = funchF(vecX_next); stepGen_datOut.fevalCount++;
return;
endfunction
