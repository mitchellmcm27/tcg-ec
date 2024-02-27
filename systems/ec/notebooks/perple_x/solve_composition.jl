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

force_pseudosection = false
normalize_composition = false

args = ARGS
if "-f" in ARGS
    force_pseudosection = true
    args = filter(x->x!="-f", args)
end

if "-n" in ARGS
    normalize_composition = true
    args = filter(x->x!="-n", args)
end

compositions = JSON.parsefile("compositions.json")
all_comp_names = collect(keys(compositions))
comp_names = []
if length(args)>0
    for arg in args
        global comp_names = vcat(comp_names, arg)
    end
else 
    global comp_names = all_comp_names
end

for name in comp_names
    composition_name = name
    comp = compositions[name]
    dataset = comp["dataset"] isa String ? comp["dataset"] : "stx21ver"
    scratchdir = "./output/"
    if dataset == "stx21ver"
        if occursin("pyrolite",name) || occursin("harzburgite",name)
            println("excluding no phases")
            phases_exclude = []
        else
            println("excluding majoritic garnet endmembers")
            phases_exclude = ["maj","namaj","namj"]
        end
        # all solution phases
        phases = ["O","Pl","Sp","Cpx","Wad","Ring","Pv","Wus","C2/c","Opx","Aki","Ppv","CF","Gt","NaAl"]
        include_fluid="" # perplex will not ask, so leave it blank
        saturated_fluid=""
    else
        phases_exclude = [
            #melt(G)
            "sil8L",
            "ctjL",
            "fo8L",
            "wi8L",
            "fa8L",
            "q8L",
            "qjL",
            "dijL",
            "jdjL",
            "ctjL",
            "fojL",
            "fajL",
            "hmjL",
            "ekjL",
            #cAmph(G)
            #Omph(HP)
            #Gt(W)
            "andr",
            #Opx(W)
            #Bi(W)
            #Mica(W)
            "ma",
            #Ilm(WPH)
            "ilm",
            "ilm_nol",
            #T
            #Sp(WPC)
            #Ep(HP11)
            #Pl(I1,HP)
            #Fsp(C1)
            ]
        # Same as Li et al. 2022
        phases = [
            "melt(G)",
            "cAmph(G)",
            "Omph(HP)",
            "Gt(W)",
            "Opx(W)" ,
            "Bi(W)",
            "Mica(W)" ,
            "Ilm(WPH)", # ilmenite-hematite
            "T",
            "Sp(WPC)", # magnetite-spinel
            "Ep(HP11)",
            "Pl(I1,HP)",
            "Fsp(C1)"
        ]
        include_fluid="n" # perplex asks
        saturated_fluid="n"
    end
    oxide_comp = convert(Array{Number},comp["composition"])
    oxides = convert(Array{String},comp["elements"])
    if dataset == "stx21ver"
        map!((s) -> uppercase(s), oxides)
    end
    composition_basis = comp["basis"]

    println(sum(oxide_comp))
    println(oxide_comp)
    println(phases)
    println(oxides)
    if composition_basis=="wt" && dataset == "stx21ver"
        if sum(oxide_comp) < 100.0

            if normalize_composition
                composition_name = composition_name * "_norm"
            end

            iFe = findfirst(x->x=="FEO",oxides)
            iSi = findfirst(x->x=="SIO2",oxides)
            iMg = findfirst(x->x=="MGO",oxides)
            iNa = findfirst(x->x=="NA2O",oxides)
            iCa = findfirst(x->x=="CAO",oxides)
            iAl = findfirst(x->x=="AL2O3",oxides)

            iTi = findfirst(x->x=="TIO2",oxides)
            iK  = findfirst(x->x=="K2O",oxides)
            iMn = findfirst(x->x=="MNO",oxides)
            iP  = findfirst(x->x=="P2O5",oxides)

            wMg = oxide_comp[iMg]
            wFe = oxide_comp[iFe]
            Xmg = wMg/(wMg+wFe)

            if iTi != nothing
                # distribute TiO2 content to Fe and Mg
                # to make high density phases
                diff = oxide_comp[iTi]
                if normalize_composition
                    oxide_comp[iFe] += diff * (1-Xmg)
                    oxide_comp[iMg] += diff * Xmg
                end
                oxide_comp[iTi] = 0
            end
            if iK != nothing
                # distribute K2O to SiO2
                # to make more quartz, which has a similar density to kspar
                diff = oxide_comp[iK]
                if normalize_composition
                    oxide_comp[iSi] += diff
                end
                oxide_comp[iK] = 0
            end

            if iK != nothing
                # distribute K2O to SiO2
                # to make more quartz, which has a similar density to kspar
                diff = oxide_comp[iK]
                if normalize_composition
                    oxide_comp[iSi] += diff
                end
                oxide_comp[iK] = 0
            end

            if iP != nothing
                # distribute P2O5 to Al2O3
                diff = oxide_comp[iP]
                if normalize_composition
                    oxide_comp[iAl] += diff
                end
                oxide_comp[iP] = 0
            end

            
            if iMn != nothing
                # distribute MnO to MgO
                diff = oxide_comp[iMn]
                if normalize_composition
                    oxide_comp[iMg] += diff
                end
                oxide_comp[iMn] = 0
            end

      
            # any remaining deficit...
            diff = 100.0 - sum(oxide_comp)
            if diff > 0 && normalize_composition
                oxide_comp[iFe] += diff * (1-Xmg)
                oxide_comp[iMg] += diff * Xmg
            end
        end
    end

    println(oxide_comp)

    blkpath = joinpath(scratchdir,composition_name,composition_name*".blk")
    pseudosection_exists = isfile(blkpath)

    if force_pseudosection || !pseudosection_exists
        print("Solving pseudosection...\n")
        perplex_build_vertex(perplexdir, scratchdir, oxide_comp,oxides, P_range_2d, T_range_2d, 
            dataset=dataset*".dat", 
            excludes=join(phases_exclude,"\n")*"\n",
            solution_phases=join(phases,"\n")*"\n", 
            composition_basis=composition_basis, 
            mode_basis=mode_basis, 
            name=composition_name,
            saturated_fluid=saturated_fluid
        )
        perplex_pssect(perplexdir, scratchdir, name=composition_name)
    end

    print("Extracting profile...\n")
    modes = perplex_werami_profile(perplexdir, scratchdir, P_range_1d, T_range_1d,name=composition_name)

    print("Extracting density grid...\n")
    modes = perplex_werami_rho(perplexdir, scratchdir, include_fluid=include_fluid, name=composition_name)

    print("Extracting point...\n")
    point = perplex_werami_point(perplexdir,scratchdir,P_point,T_point,name=composition_name)
    print(point)
end