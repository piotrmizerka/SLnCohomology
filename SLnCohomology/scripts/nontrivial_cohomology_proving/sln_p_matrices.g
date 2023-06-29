OutputLogTo("./sln_p_matrices.txt") # temporary file, erased at the end
G := SL(3,3);; # these parameters have to be manually adjusted and agree with N and p from "sln_p_nontrivial_cohom.jl"
G_list := AsList(G);;
iso := IsomorphismPermGroup(G);;
for g in G_list do
    Display(g);
od;
