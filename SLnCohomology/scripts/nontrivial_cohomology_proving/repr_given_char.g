G := SL(2,5);
chi := Irr(G)[4];
Display(chi);
rep := IrreducibleRepresentationsDixon(G,chi);
# rep2 := IrreducibleAffordingRepresentation(chi); is this function deprecated??
Display(rep);
# Display(rep2);
reps := IrreducibleRepresentationsDixon(G:unitary);
Display(reps);