# This GAP script saves all elements of SL(N,p) (N and p hard-coded, to change accordingly).

# Before running this script, run GAP with admin priviliges (sudo on mac).
OutputLogTo("./sl5_2_matrices.txt"); # hard-coded, better siubstitute with the absolute path if you don't want to look for the logs too long :)
G := SL(5,2);; # hard-coded: these parameters have to be manually adjusted and agree with N and p from "sln_nontrivial_cohomology.jl"
G_list := AsList(G);;
for g in G_list do
    Display(g);
od;
