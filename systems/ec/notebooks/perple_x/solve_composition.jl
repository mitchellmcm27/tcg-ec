using JSON
include("perplex_api.jl")

# Absolute path to Perple_X installation
resourcepath = "/root/resources"
perplexdir = joinpath(resourcepath,"perplex-stable")

mode_basis = "vol"

T_range_2d = (300+273.15, 1300+273.15) # Kelvin
P_range_2d = (5000, 25000) # bar

T_range_1d = (650+273.15, 850+273.15) # Kelvin
P_range_1d = (5000, 25000) # bar

T_point = 900+273.15 # K
P_point = 2.0e4 # bar

T_surf = 273.15
melt_model = "melt(G)"

force_pseudosection = false
args = ARGS
if "-f" in ARGS
    force_pseudosection = true
    args = filter(x->x!="-f", args)
end

compositions = JSON.parsefile("compositions.json")
all_comp_names = collect(keys(compositions))
comp_names = []
if length(args)>0
    for arg in args
        comp_matches = filter(contains(Regex(arg)), all_comp_names)
        if length(comp_matches) > 0
            global comp_names = vcat(comp_names, comp_matches)
        end
    end
else 
    global comp_names = all_comp_names
end

for name in comp_names
    comp = compositions[name]

    dataset = comp["dataset"] isa String ? comp["dataset"] : "stx21ver"

    if dataset == "stx21ver"
        scratchdir = "./output/"
        phases_exclude = ["maj","namaj","namj"]
        phases = ["O","Pl","Sp","Cpx","Wad","Ring","Pv","Wus","C2/c","Opx","Aki","Ppv","CF","Gt","NaAl"]
        #phases = ["Pl","Cpx","Opx","Gt"]
    else
        scratchdir = "./output"*dataset*"/"
        phases_exclude = []
        phases = [
            "Fsp(C1)",
            "melt(G)",
            "Omph(GHP)",
            "Gt(W)",
            "Opx(W)" ,
            "Mt(W)" ,
            "Ilm(WPH)",
            "cAmph(G)",
            "Bi(W)",
            "Mica(W)",
            "T"
        ]
    end
    oxide_comp = convert(Array{Number},comp["composition"])
    oxides = convert(Array{String},comp["elements"])
    composition_basis = comp["basis"]

    blkpath = joinpath(scratchdir,name,name*".blk")
    pseudosection_exists = isfile(blkpath)

    if force_pseudosection || !pseudosection_exists
        print("Solving pseudosection...\n")
        perplex_build_vertex(perplexdir, scratchdir, oxide_comp,oxides, P_range_2d, T_range_2d, 
            dataset=dataset*".dat", 
            excludes=join(phases_exclude,"\n")*"\n",
            solution_phases=join(phases,"\n")*"\n", 
            composition_basis=composition_basis, 
            mode_basis=mode_basis, 
            name=name
        )
        perplex_pssect(perplexdir, scratchdir, name=name)
    end

    print("Extracting profile...\n")
    modes = perplex_werami_profile(perplexdir, scratchdir, P_range_1d, T_range_1d, name=name)

    print("Extracting density grid...\n")
    modes = perplex_werami_rho(perplexdir, scratchdir, name=name)

    print("Extracting point...\n")
    point = perplex_werami_point(perplexdir,scratchdir,P_point,T_point,name=name)
    print(point)

end