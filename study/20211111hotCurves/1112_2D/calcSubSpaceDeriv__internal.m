	if ( 0.0 == norm(vecA) || 0.0 == norm(vecB) )
		vecATKB = zeros(sizeF,1);
		return;
	end
	%
	sFD1 = norm(vecEps.*vecA)/(norm(vecA)^2);
	sFD2 = norm(vecEps.*vecB)/(norm(vecB)^2);
	vecXPP = vecX0 + (sFD1*vecA) + (sFD2*vecB);
	vecXPM = vecX0 + (sFD1*vecA) - (sFD2*vecB);
	vecXMP = vecX0 - (sFD1*vecA) + (sFD2*vecB);
	vecXMM = vecX0 - (sFD1*vecA) - (sFD2*vecB);
	clear sFD1;
	clear sFD2;
	%
	vecFPP = funchF(vecXPP);
	vecFPM = funchF(vecXPM);
	vecFMP = funchF(vecXMP);
	vecFMM = funchF(vecXMM);
	assert( isrealarray(vecFPP,[sizeF,1]) );
	assert( isrealarray(vecFPM,[sizeF,1]) );
	assert( isrealarray(vecFMP,[sizeF,1]) );
	assert( isrealarray(vecFMM,[sizeF,1]) );
	%
	vecTempPD = (vecFPP-vecFPM)/norm(vecXPP-vecXPM);
	vecTempMD = (vecFMP-vecFMM)/norm(vecXMP-vecXMM);
	vecTempDP = (vecFPP-vecFMP)/norm(vecXPP-vecXMP);
	vecTempDM = (vecFPM-vecFMM)/norm(vecXPM-vecXMM);
	%
	vecXPC = (vecXPP+vecXPM)/2.0;
	vecXMC = (vecXMP+vecXMM)/2.0;
	vecXCP = (vecXPP+vecXMP)/2.0;
	vecXCM = (vecXPM+vecXMM)/2.0;
	%
	vecTemp1 = (vecTempPD-vecTempMD)/norm(vecXPC-vecXMC);
	vecTemp2 = (vecTempDP-vecTempDM)/norm(vecXCP-vecXCM);
	vecATKB = (vecTemp1+vecTemp2)*0.5*norm(vecA)*norm(vecB);
	%
	clear vecXPP;
	clear vecXPM;
	clear vecXMP;
	clear vecXMM;
	%
	clear vecFPP;
	clear vecFPM;
	clear vecFMP;
	clear vecFMM;
	%
	clear vecTempPD;
	clear vecTempMD;
	clear vecTempDP;
	clear vecTempMD;
	%
	clear vecXPC;
	clear vecXMC;
	clear vecXCP;
	clear vecXCM;
	%
	clear vecTemp1;
	clear vecTemp2;
return;
