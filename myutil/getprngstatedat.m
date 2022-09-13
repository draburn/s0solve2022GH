function statedat = getprngstatedat()
	statedat.ofrand = rand("state");
	statedat.ofrandn = randn("state");
	statedat.ofrande = rande("state");
	statedat.ofrandp = randp("state");
	statedat.ofrandg = randg("state");
return;
end
