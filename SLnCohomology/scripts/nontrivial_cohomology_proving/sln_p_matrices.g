OutputLogTo("./sl5_2_matrices.txt") # hard-coded
G := SL(5,2);; # hard-coded: these parameters have to be manually adjusted and agree with N and p from "sln_nontrivial_cohomology.jl"
G_list := AsList(G);;
for g in G_list do
    Display(g);
od;
