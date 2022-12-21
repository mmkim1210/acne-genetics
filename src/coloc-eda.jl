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
    "lupus", "ra", "t1d", "celiac", "uc", "crohns", "ibd", "as", "pbc", "asthma", "oa",
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

@info "Loading GENCODE"
@time gencode = Arrow.Table("data/gencode/gencode.v39lift37.annotation.gtf.arrow")|> DataFrame
gencode.gene_id = gencode.gene_name

@info "Loading 1000 Genomes reference panel"
@time kgp_raw = SnpData("data/1kg/kgp.eur.maf0.05")
@assert size(kgp_raw) == (503, 7_230_618)

ind = findfirst(isequal("acne"), phenos)
loci = GeneticsMakie.findgwasloci(gwas[ind])

filter!(x -> !(x.CHR == "6" && x.BP > 24e6 && x.BP < 34e6), loci)
filter!(x -> !(x.CHR == "5" && x.BP > 129e6 && x.BP < 132e6), loci)
filter!(x -> !(x.CHR == "8" && x.BP > 7e6 && x.BP < 13e6), loci)

loci1 = GeneticsMakie.findclosestgene(loci, gencode; proteincoding = true)

for (i, gene) in enumerate(loci1.gene)
    @info "Working on $gene gene"
    window = 1e6
    chr, start, stop = GeneticsMakie.findgene(gene, gencode)
    range1, range2 = start - window, stop + window
    range1 < 0 ? range1 = 0 : nothing
    range2 > GeneticsMakie.GRCh37_totlength[chr] ? range2 = GeneticsMakie.GRCh37_totlength[chr] : nothing
    @info "Subsetting 1kg"
    kgp = subsetref(kgp_raw, chr, range1, range2, "data/kgp.filtered")
    @info "Subsetting GWAS"
    dfs = subsetgwas(gwas, chr, range1, range2)

    LD = let
        geno = convert(Matrix{Float64}, kgp.snparray)
        LD = cor(geno, dims = 1)
        LD.^2
    end

    for i in 1:length(dfs)
        if issig(dfs[i]) && i != ind
            storage = innerjoin(dfs[i], dfs[ind], on = [:CHR, :BP], makeunique = true)
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

            f = Figure(resolution = (306, 175))
            axs = [Axis(f[1, i]) for i in 1:2]
            x = -log10.(storage[!, "P_1"])
            y = -log10.(storage[!, "P"])
            scatter!(axs[1], x, y, 
                markersize = 1.5, color = storage[!, "LD1"], colorrange = (0, 1), 
                colormap = cgrad(locuszoomcolors, 5, categorical = true, rev = true))
            scatter!(axs[2], x, y, 
                markersize = 1.5, color = storage[!, "LD2"], colorrange = (0, 1), 
                colormap = cgrad(locuszoomcolors, 5, categorical = true, rev = true))
            scatter!(axs[1], [maximum(x)], [y[argmax(x)]], color = :purple1, markersize = 4.0, marker = '◆')
            scatter!(axs[2], [x[argmax(y)]], [maximum(y)], color = :purple1, markersize = 4.0, marker = '◆')
            text!(axs[1], "$(storage.SNP_1[argmax(x)])", position = (maximum(x), y[argmax(x)]), textsize = 5, align = (:right, :bottom))
            text!(axs[2], "$(storage.SNP_1[argmax(y)])", position = (x[argmax(y)], maximum(y)), textsize = 5, align = (:left, :center))
            [vlines!(axs[i], -log10(5e-8), color = (:gold, 0.5), linewidth = 0.5) for i in 1:2]
            [hlines!(axs[i], -log10(5e-8), color = (:gold, 0.5), linewidth = 0.5) for i in 1:2]
            axs[1].spinewidth = 0.75
            axs[1].ytickwidth = 0.75
            axs[1].ylabelsize = 6
            axs[1].yticklabelsize = 6
            axs[1].yticksize = 3
            axs[1].xtickwidth = 0.75
            axs[1].xlabelsize = 6
            axs[1].xticklabelsize = 6
            axs[1].xticksize = 3
            axs[1].xlabel = "$(newtitles[ind]) -log[p]"
            axs[1].ylabel = "$(newtitles[i]) -log[p]"
            axs[2].spinewidth = 0.75
            axs[2].ytickwidth = 0.75
            axs[2].ylabelsize = 6
            axs[2].yticklabelsize = 6
            axs[2].yticksize = 3
            axs[2].xtickwidth = 0.75
            axs[2].xlabelsize = 6
            axs[2].xticklabelsize = 6
            axs[2].xticksize = 3
            axs[2].xlabel = "$(newtitles[ind]) -log[p]"
            axs[2].ylabel = "$(newtitles[i]) -log[p]"
            [hidedecorations!(axs[i], ticks = false, ticklabels = false, label = false) for i in 1:2]
            [ylims!(axs[i], 0, maximum(y) + 2.5) for i in 1:2]
            [xlims!(axs[i], 0, maximum(x) + 2.5) for i in 1:2]
            Colorbar(f[1, 3], limits = (0, 1), ticks = 0:1:1,
                colormap = cgrad(locuszoomcolors, 5, categorical = true, rev = true), 
                label = "LD", ticksize = 0, height = 25,
                tickalign = 0, ticklabelsize = 6, flip_vertical_label = true, 
                labelsize = 6, width = 5, spinewidth = 0.5, tickwidth = 0)
            rowsize!(f.layout, 1, Aspect(1, 1))
            colgap!(f.layout, 5)
            save("figs/coloc/coloc-$(gene)-$(phenos[i]).png", f, px_per_unit = 4)
        end
    end
end