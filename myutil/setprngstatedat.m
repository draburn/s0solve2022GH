function setprngstatedat( statedat )
	rand("state",statedat.ofrand);
	randn("state",statedat.ofrandn);
	rande("state",statedat.ofrande);
	randp("state",statedat.ofrandp);
	randg("state",statedat.ofrandg);
return;
end
