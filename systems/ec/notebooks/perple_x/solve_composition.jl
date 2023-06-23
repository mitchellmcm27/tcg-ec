using JSON
include("perplex_api.jl")

# Absolute path to Perple_X installation
resourcepath = "/root/resources"
perplexdir = joinpath(resourcepath,"perplex-stable")
scratchdir = "./output/" # Location of directory to store output files

dataset = "stx21ver.dat"
mode_basis = "vol"
phases_all = ["O","Pl","Sp","Cpx","Wad","Ring","Pv","Wus","C2/c","Opx","Aki","Ppv","CF","Gt","NaAl"]
phases_simple = ["Pl","Cpx","Opx","Gt"]
#phases_exclude = ["maj","namaj","namj","sp","cor","herc","neph","fo","fa",]
phases_exclude = ["maj","namaj","namj"]

T_range_2d = (500+273.15, 1000+273.15) # Kelvin
P_range_2d = (5000, 25000) # bar

T_range_1d = (500+273.15,1000+273.15) # Kelvin
P_range_1d = (25000, 5000) # bar

T_point = 800+273.15 # K
P_point = 1.0e4 # bar

T_surf = 273.15
melt_model = "melt(G)"

force_pseudosection = false
if "-f" in ARGS
    force_pseudosection = true
end

compositions = JSON.parsefile("compositions.json")

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
    oxide_comp = convert(Array{Number},comp["composition"])
    oxides = convert(Array{String},comp["elements"])
    composition_basis = comp["basis"]

    blkpath = joinpath(scratchdir,name,name*".blk")
    pseudosection_exists = isfile(blkpath)

    if force_pseudosection || !pseudosection_exists
        print("Solving pseudosection...\n")
        perplex_build_vertex(perplexdir, scratchdir, oxide_comp,oxides, P_range_2d, T_range_2d, 
            dataset=dataset, 
            excludes=join(phases_exclude,"\n")*"\n",
            solution_phases=join(phases_all,"\n")*"\n", 
            composition_basis=composition_basis, 
            mode_basis=mode_basis, 
            name=name
        )
    end

    print("Extracting profile...\n")
    modes = perplex_werami_profile(perplexdir, scratchdir, P_range_1d, T_range_1d, name=name)

    print("Extracting density grid...\n")
    modes = perplex_werami_rho(perplexdir, scratchdir, name=name)

    print("Extracting point...\n")
    point = perplex_werami_point(perplexdir,scratchdir,P_point,T_point,name=name)
    print(point)

end