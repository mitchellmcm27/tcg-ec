include("perplex_api.jl")
include("oxide_compositions.jl")

# Absolute path to Perple_X installation
resourcepath = "/root/resources"
perplexdir = joinpath(resourcepath,"perplex-stable")
scratchdir = "./output/" # Location of directory to store output files

phases_all = "O\nPl\nSp\nCpx\nWad\nRing\nPv\nWus\nC2/c\nOpx\nAki\nPpv\nCF\nGt\nNaAl\n"
phases_simple = "Pl\nCpx\nOpx\nGt\n"
phases_exclude = "maj\nnamaj\nnamj\n"

T_range_2d = (500+273.15, 1000+273.15) # Kelvin
P_range_2d = (5000, 25000) # bar

T_range_1d = (500+273.15,1000+273.15)
P_range_1d = (25000, 5000)

T_surf = 273.15
melt_model = "melt(G)"


force_pseudosection = false
if "-f" in ARGS
    force_pseudosection = true
end

if length(ARGS)>0
    comp_names = []
    for arg in ARGS
        if arg in keys(compositions)
            push!(comp_names, arg)
        end
    end
else 
    comp_names = keys(compositions)
end

for name in comp_names
    comp = compositions[name]
    blkpath = joinpath(scratchdir,name,name*".blk")
    print(blkpath)
    print("\n")
    pseudosection_exists = isfile(blkpath)

    if force_pseudosection || !pseudosection_exists
        print("Solving pseudosection...\n")
        perplex_build_vertex(perplexdir, scratchdir, comp["composition"], comp["elements"], P_range_2d, T_range_2d, 
            dataset="stx21ver.dat", 
            excludes=phases_exclude,
            solution_phases=phases_all, 
            composition_basis=comp["basis"], 
            mode_basis="vol", 
            name=name
        )
    end

    print("Extracting profile...\n")
    modes = perplex_werami_profile(perplexdir, scratchdir, P_range_1d, T_range_1d, name=name)

    print("Extracting density grid...\n")
    modes = perplex_werami_rho(perplexdir, scratchdir, name=name)

    print("Extracting point...\n")
    point = perplex_werami_point(perplexdir,scratchdir,1e4,800+273.15,name=name)
    print(point)

end