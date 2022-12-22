using GeneticsMakie, CairoMakie, CSV, DataFrames, SnpArrays, Arrow, StatsBase, Colors, ColorSchemes

CairoMakie.activate!(type = "png")
set_theme!(font = "Arial")

locuszoomcolors = ColorScheme([
    RGB{Float64}(213 / 255, 63 / 255, 57 / 255),
    RGB{Float64}(250 / 255, 167 / 255, 1 / 255),
    RGB{Float64}(89 / 255, 175 / 255, 88 / 255),
    RGB{Float64}(70 / 255, 185 / 255, 217 / 255),
    RGB{Float64}(54 / 255, 126 / 255, 189 / 255)
    ])

function subsetref(ref::SnpData, chr::AbstractString, range1::Real, range2::Real, path::AbstractString)
    SnpArrays.filter(ref, trues(size(ref)[1]), GeneticsMakie.findlocus(ref, chr, range1, range2); des = path)
    SnpData(path)
end

function subsetgwas(gwas, chr, range1, range2)
    dfs = Vector{DataFrame}(undef, length(gwas))
    for i in 1:length(gwas)
        dfs[i] = gwas[i][GeneticsMakie.findlocus(gwas[i], chr, range1, range2), :]
    end
    dfs
end

issig(P::AbstractVector; p = 5e-8) = any(P .< p)
issig(df::DataFrame; p = 5e-8) = issig(df.P; p = p)

@info "Load GWAS results"
delete!(GeneticsMakie.gwas, "mdd")
gwas = []
phenos = ["asd", "adhd", "anorexia", "scz", "bd", "sczvsbd", "cross", 
    "alz", "parkinson", "als", "stroke",
    "left-hand", "cp", "intelligence", "ea", "birth", "children", "risk", "neuroticism", "smoking", "alcohol", "insomnia", "icv",
    "amd", "cataract", "glaucoma", "menarche", "menopause", "pcos",
    "lupus", "ra", "t1d", "uc", "crohns", "ibd", "pbc", "asthma", "oa",
    "gerd", "barrett", "height", "weight", "tsh", 
    "cad", "afib", "dbp", "sbp", "ckd", "t2d", "covid",
    "acne", "bald", "ad", "vitiligo", "psoriasis",
    "bcc", "melanoma", "breast", "prostate", "ec",
    "rbc", "hgb", "hct", "mch", "mcv", "mchc", "rdw", "wbc", "neu", "lymph", "mono", "baso", "eos", "plt", "mpv",
    "alt", "ast", "astalt", "ap", "apoa", "apob", "crp", "phosphate", "calcium", "urate", "bun", "creatinine", "cystatinc", "ggt",
    "triglyceride", "cholesterol", "hdl", "ldl", "lpa",
    "glucose", "hba1c", "igf1", "vitamind",
    "egfr", "ucreatinine", "usodium", "upotassium", "ualbumin",
    "shbg", "testosterone", "totbilirubin", "dirbilirubin", "totprotein", "albumin", "nap"]
for key in phenos
    push!(gwas, DataFrame(Arrow.Table(joinpath.("data/gwas/processed/", key * ".tsv.arrow"))))
end
titles = [GeneticsMakie.gwas[key].title for key in phenos]
newtitles = [getindex(i, 1) for i in split.(titles, " (")]
newtitles[findfirst(isequal("Male-pattern baldness"), newtitles)] = "Baldness"
newtitles[findfirst(isequal("Eosinophil count"), newtitles)] = "Eosinophil"
newtitles[findfirst(isequal("Monocyte count"), newtitles)] = "Monocyte"
newtitles[findfirst(isequal("Mean corpuscular hemoglobin"), newtitles)] = "MCH"
newtitles[findfirst(isequal("Mean corpuscular volume"), newtitles)] = "MCV"
newtitles[findfirst(isequal("AST / ALT"), newtitles)] = "AST/ALT"
newtitles[findfirst(isequal("Type II diabetes"), newtitles)] = "T2D"
newtitles[findfirst(isequal("Red blood cell count"), newtitles)] = "RBC"
newtitles[findfirst(isequal("Hemoglobin concentration"), newtitles)] = "Hgb"
newtitles[findfirst(isequal("Platelet count"), newtitles)] = "Platelet"
newtitles[findfirst(isequal("Neutrophil count"), newtitles)] = "Neutrophil"
newtitles[findfirst(isequal("White blood cell count"), newtitles)] = "WBC"
newtitles[findfirst(isequal("RBC distribution width"), newtitles)] = "RDW"
newtitles[findfirst(isequal("Mean platelet volume"), newtitles)] = "MPV"

@info "Loading GENCODE"
@time gencode = Arrow.Table("data/gencode/gencode.v39lift37.annotation.gtf.arrow")|> DataFrame
gencode.gene_id = gencode.gene_name

@info "Loading 1000 Genomes reference panel"
@time kgp_raw = SnpData("data/1kg/kgp.eur.maf0.05")
@assert size(kgp_raw) == (503, 7_230_618)

focus = Dict(
    "WNT10A" => ["bald"],
    "CRELD2" => ["shbg"],
    "FGF10" => ["prostate"],
    "SOAT1" => ["testosterone", "igf1"],
    "CLEC16A" => ["totprotein", "nap", "eos"],
    "BCL11A" => ["scz", "nap", "astalt", "ast"],
    "CSTA" => ["phosphate", "mcv", "mch", "calcium"],
    "FCHO2" => ["nap", "mono", "mcv", "eos", "ast", "albumin"],
    "PNPLA3" => ["urate", "wbc", "totbilirubin", "testosterone", "t2d", "shbg", "rdw", "rbc", 
        "plt", "neu", "mpv", "mcv", "mch", "hgb", "hdl", "hct", "dirbilirubin", 
        "cholesterol", "astalt", "ast", "apoa", "alt"]
    )

begin
    f = Figure(resolution = (612, 612 * 2))
    axs = [Axis(f[i, j]) for i in 1:8, j in 1:8]
    window = 1e6
    for (i, gene) in enumerate(["WNT10A", "CRELD2", "FGF10", "SOAT1", "CLEC16A", "BCL11A", "CSTA", "FCHO2", "PNPLA3"])
        @info "Working on $gene gene"
        chr, start, stop = GeneticsMakie.findgene(gene, gencode)
        range1, range2 = start - window, stop + window
        range1 < 0 ? range1 = 0 : nothing
        range2 > GeneticsMakie.GRCh37_totlength[chr] ? range2 = GeneticsMakie.GRCh37_totlength[chr] : nothing
        @info "Subsetting 1kg"
        kgp = subsetref(kgp_raw, chr, range1, range2, "data/kgp.filtered")
        @info "Subsetting GWAS"
        ind = findall(in(focus[gene]), phenos)
        push!(ind, findfirst(==("acne"), phenos))
        dfs = subsetgwas(gwas[ind], chr, range1, range2)
        LD = let
            geno = convert(Matrix{Float64}, kgp.snparray)
            LD = cor(geno, dims = 1)
            LD.^2
        end
        if gene != "PNPLA3"
            for k in 1:(length(ind) - 1)
                storage = innerjoin(dfs[k], dfs[end], on = [:CHR, :BP], makeunique = true)
                storage.ind .= 0
                for j in 1:nrow(storage)
                    ind1 = findfirst(isequal(storage.BP[j]), kgp.snp_info.position)
                    isnothing(ind1) ? continue : storage.ind[j] = ind1
                end
                filter!(x -> x.ind != 0, storage)
                storage.LD1 .= 1.0
                storage.LD2 .= 1.0
                ind1 = argmin(storage.P_1)
                for j in 1:nrow(storage)
                    storage.LD1[j] = LD[storage.ind[ind1], storage.ind[j]]
                end
                ind1 = argmin(storage.P)
                for j in 1:nrow(storage)
                    storage.LD2[j] = LD[storage.ind[ind1], storage.ind[j]]
                end
                x = -log10.(storage[!, "P_1"])
                y = -log10.(storage[!, "P"])
                scatter!(axs[k, i], x, y, 
                    markersize = 1.5, color = storage[!, "LD1"], colorrange = (0, 1), 
                    colormap = cgrad(locuszoomcolors, 5, categorical = true, rev = true))
                scatter!(axs[k, i], [maximum(x)], [y[argmax(x)]], color = :purple1, markersize = 4.0, marker = '◆')
                text!(axs[k, i], "$(storage.SNP_1[argmax(x)])", position = (maximum(x), y[argmax(x)]), textsize = 5, align = (:right, :bottom))
                vlines!(axs[k, i], -log10(5e-8), color = (:gray, 0.4), linewidth = 0.5)
                hlines!(axs[k, i], -log10(5e-8), color = (:gray, 0.4), linewidth = 0.5)
                ylims!(axs[k, i], 0, maximum(y) + 2.5)
                xlims!(axs[k, i], 0, maximum(x) + 2.5) 
                # scatter!(axs[k, 2i - 1], x, y, 
                #     markersize = 1.5, color = storage[!, "LD1"], colorrange = (0, 1), 
                #     colormap = cgrad(locuszoomcolors, 5, categorical = true, rev = true))
                # scatter!(axs[k, 2i], x, y, 
                #     markersize = 1.5, color = storage[!, "LD2"], colorrange = (0, 1), 
                #     colormap = cgrad(locuszoomcolors, 5, categorical = true, rev = true))
                # scatter!(axs[k, 2i - 1], [maximum(x)], [y[argmax(x)]], color = :purple1, markersize = 4.0, marker = '◆')
                # scatter!(axs[k, 2i], [x[argmax(y)]], [maximum(y)], color = :purple1, markersize = 4.0, marker = '◆')
                # text!(axs[k, 2i - 1], "$(storage.SNP_1[argmax(x)])", position = (maximum(x), y[argmax(x)]), textsize = 5, align = (:right, :bottom))
                # text!(axs[k, 2i], "$(storage.SNP_1[argmax(y)])", position = (x[argmax(y)], maximum(y)), textsize = 5, align = (:left, :bottom))
                # [vlines!(axs[k, l], -log10(5e-8), color = (:gold, 0.5), linewidth = 0.5) for l in 2i:(2i - 1)]
                # [hlines!(axs[k, l], -log10(5e-8), color = (:gold, 0.5), linewidth = 0.5) for l in 2i:(2i - 1)]
                # [ylims!(axs[k, l], 0, maximum(y) + 2.5) for l in 2i:(2i - 1)]
                # [xlims!(axs[k, l], 0, maximum(x) + 2.5) for l in 2i:(2i - 1)]
                Label(f[k, i, Top()], "$(newtitles[ind][k]) ($(gene))", textsize = 6, valign = :bottom)
            end
        else
            for k in 1:(length(ind) - 1)
                storage = innerjoin(dfs[k], dfs[end], on = [:CHR, :BP], makeunique = true)
                storage.ind .= 0
                for j in 1:nrow(storage)
                    ind1 = findfirst(isequal(storage.BP[j]), kgp.snp_info.position)
                    isnothing(ind1) ? continue : storage.ind[j] = ind1
                end
                filter!(x -> x.ind != 0, storage)
                storage.LD1 .= 1.0
                storage.LD2 .= 1.0
                ind1 = argmin(storage.P_1)
                for j in 1:nrow(storage)
                    storage.LD1[j] = LD[storage.ind[ind1], storage.ind[j]]
                end
                ind1 = argmin(storage.P)
                for j in 1:nrow(storage)
                    storage.LD2[j] = LD[storage.ind[ind1], storage.ind[j]]
                end
                x = -log10.(storage[!, "P_1"])
                y = -log10.(storage[!, "P"])
                if k < 21
                    indx = 3 + mod1(k, 5)
                    indy = div(k - 1, 5) + 1
                    scatter!(axs[indx, indy], x, y, 
                        markersize = 1.5, color = storage[!, "LD1"], colorrange = (0, 1), 
                        colormap = cgrad(locuszoomcolors, 5, categorical = true, rev = true))
                    scatter!(axs[indx, indy], [maximum(x)], [y[argmax(x)]], color = :purple1, markersize = 4.0, marker = '◆')
                    text!(axs[indx, indy], "$(storage.SNP_1[argmax(x)])", position = (maximum(x), y[argmax(x)]), textsize = 5, align = (:right, :bottom))
                    vlines!(axs[indx, indy], -log10(5e-8), color = (:gray, 0.4), linewidth = 0.5)
                    hlines!(axs[indx, indy], -log10(5e-8), color = (:gray, 0.4), linewidth = 0.5)
                    ylims!(axs[indx, indy], 0, maximum(y) + 2.5)
                    xlims!(axs[indx, indy], 0, maximum(x) + 2.5)
                    Label(f[indx, indy, Top()], "$(newtitles[ind][k]) ($(gene))", textsize = 6, valign = :bottom)
                else
                    scatter!(axs[k - 14, 5], x, y, 
                        markersize = 1.5, color = storage[!, "LD1"], colorrange = (0, 1), 
                        colormap = cgrad(locuszoomcolors, 5, categorical = true, rev = true))
                    scatter!(axs[k - 14, 5], [maximum(x)], [y[argmax(x)]], color = :purple1, markersize = 4.0, marker = '◆')
                    text!(axs[k - 14, 5], "$(storage.SNP_1[argmax(x)])", position = (maximum(x), y[argmax(x)]), textsize = 5, align = (:right, :bottom))
                    vlines!(axs[k - 14, 5], -log10(5e-8), color = (:gray, 0.4), linewidth = 0.5)
                    hlines!(axs[k - 14, 5], -log10(5e-8), color = (:gray, 0.4), linewidth = 0.5)
                    ylims!(axs[k - 14, 5], 0, maximum(y) + 2.5)
                    xlims!(axs[k - 14, 5], 0, maximum(x) + 2.5)
                    Label(f[k - 14, 5, Top()], "$(newtitles[ind][k]) ($(gene))", textsize = 6, valign = :bottom)
                end
            end
        end
    end
    [hidedecorations!(axs[i, j], ticks = false, ticklabels = false, label = false) for i in 1:8, j in 1:8]
    [axs[i, j].spinewidth = 0.75 for i in 1:8, j in 1:8]
    [axs[i, j].ytickwidth = 0.75 for i in 1:8, j in 1:8]
    [axs[i, j].ylabelsize = 6 for i in 1:8, j in 1:8]
    [axs[i, j].yticklabelsize = 6 for i in 1:8, j in 1:8]
    [axs[i, j].yticksize = 3 for i in 1:8, j in 1:8]
    [axs[i, j].xtickwidth = 0.75 for i in 1:8, j in 1:8]
    [axs[i, j].xlabelsize = 6 for i in 1:8, j in 1:8]
    [axs[i, j].xticklabelsize = 6 for i in 1:8, j in 1:8]
    [axs[i,j ].xticksize = 3 for i in 1:8, j in 1:8]
    hidespines!(axs[3, 4])
    hidedecorations!(axs[3, 4])
    hidespines!(axs[4, 5])
    hidedecorations!(axs[4, 5])
    [hidespines!(axs[i, j]) for i in 2:3, j in 1:3]
    [hidedecorations!(axs[i, j]) for i in 2:3, j in 1:3]
    [hidespines!(axs[i, j]) for i in 5:6, j in 5:7]
    [hidedecorations!(axs[i, j]) for i in 5:6, j in 5:7]
    [hidespines!(axs[i, j]) for i in 7:8, j in 6:8]
    [hidedecorations!(axs[i, j]) for i in 7:8, j in 6:8]
    Label(f[1:end, 0], text = "Other -log[p]", textsize = 8, rotation = pi / 2, tellheight = false)
    Label(f[9, 1:end], text = "Acne -log[p]", textsize = 8)
    Colorbar(f[1:(end - 1), 10], limits = (0, 1), ticks = 0:1:1,
        colormap = cgrad(locuszoomcolors, 5, categorical = true, rev = true), 
        label = "LD", ticksize = 0, height = 25,
        tickalign = 0, ticklabelsize = 6, flip_vertical_label = true, 
        labelsize = 6, width = 5, spinewidth = 0.5, tickwidth = 0)
    [rowsize!(f.layout, i, Aspect(2, 1)) for i in 1:8]
    colgap!(f.layout, 5)
    rowgap!(f.layout, 5)
    resize_to_layout!(f)
    save("figs/coloc.png", f, px_per_unit = 4)
end