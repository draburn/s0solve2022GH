function newState = setprngstates( newState=[], printState=true )
	thisFile = "setprngstates";
	if (isempty(newState))
		newState = mod( round(now*1E11), 1E8 );
	end
	if (printState)
		msg(thisFile,__LINE__,["Setting pRNG states to " num2str(newState) ]);
	end
	rand("state",newState);
	randn("state",newState);
	rande("state",newState);
	randp("state",newState);
	randg("state",newState);
end
