using GeneticsMakie, CairoMakie, CSV, DataFrames, SnpArrays, Arrow

CairoMakie.activate!(type = "png")
set_theme!(font = "Arial")

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

@info "Loading GENCODE"
@time gencode = Arrow.Table("data/gencode/gencode.v39lift37.annotation.gtf.arrow")|> DataFrame
gencode.gene_id = gencode.gene_name

@info "Loading 1000 Genomes reference panel"
@time kgp_raw = SnpData("data/1kg/kgp.eur.maf0.05")
@assert size(kgp_raw) == (503, 7_230_618)

@info "Load GWAS results"
[delete!(GeneticsMakie.gwas, key) for key in ["mdd", "as", "celiac"]]

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

function locuszoom(genes)
    for gene in genes
        @info "Working on $gene gene"
        window = 1e6
        chr, start, stop = GeneticsMakie.findgene(gene, gencode)
        range1, range2 = start - window, stop + window
        range1 < 0 ? range1 = 0 : nothing
        range2 > GeneticsMakie.GRCh37_totlength[chr] ? range2 = GeneticsMakie.GRCh37_totlength[chr] : nothing
        
        @info "Subsetting 1000 Genomes"
        @time kgp = subsetref(kgp_raw, chr, range1, range2, "data/kgp.filtered")

        @info "Subsetting GWAS results"
        @time dfs = subsetgwas(gwas, chr, range1, range2)

        @info "Plotting phenome-wide LocusZoom"
        n = length(titles)
        m = div(length(titles), 2)
        f = Figure(resolution = (612, 1500))
        ga = f[1, 1] = GridLayout()
        gb = f[1, 2] = GridLayout()
        axs1 = [Axis(ga[i, 1]) for i in 1:(m + 1)]
        axs2 = [Axis(gb[i, 1]) for i in 1:(n - m + 1)]
        for i in 1:m
            if issig(dfs[i])
                GeneticsMakie.plotlocus!(axs1[i], chr, range1, range2, dfs[i]; ld = kgp)
                if dfs[i].BP[findmin(dfs[i].P)[2]] < (range1 + range2) / 2
                    Label(ga[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
                else
                    Label(ga[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
                end    
            else
                GeneticsMakie.plotlocus!(axs1[i], chr, range1, range2, dfs[i])
                Label(ga[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end
            rowsize!(ga, i, 30)
        end
        for i in 1:(n - m)
            if issig(dfs[i + m])
                GeneticsMakie.plotlocus!(axs2[i], chr, range1, range2, dfs[i + m]; ld = kgp)
                if dfs[i + m].BP[findmin(dfs[i + m].P)[2]] < (range1 + range2) / 2
                    Label(gb[i, 1, Top()], "$(titles[i + m])", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
                else
                    Label(gb[i, 1, Top()], "$(titles[i + m])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
                end    
            else
                GeneticsMakie.plotlocus!(axs2[i], chr, range1, range2, dfs[i + m])
                Label(gb[i, 1, Top()], "$(titles[i + m])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end
            rowsize!(gb, i, 30)
        end
        GeneticsMakie.plotgenes!(axs1[m + 1], chr, range1, range2, gencode; height = 0.1)
        rs = GeneticsMakie.plotgenes!(axs2[n - m + 1], chr, range1, range2, gencode; height = 0.1)
        rowsize!(ga, m + 1, rs)
        rowsize!(gb, n - m + 1, rs)
        GeneticsMakie.labelgenome(ga[m + 1, 1, Bottom()], chr, range1, range2)
        GeneticsMakie.labelgenome(gb[n - m + 1, 1, Bottom()], chr, range1, range2)
        Colorbar(ga[1:m, 2], limits = (0, 1), ticks = 0:1:1, height = 20,
            colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
            tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
            labelsize = 6, width = 5, spinewidth = 0.5)
        Colorbar(gb[1:m, 2], limits = (0, 1), ticks = 0:1:1, height = 20,
            colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
            tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
            labelsize = 6, width = 5, spinewidth = 0.5)
        Label(ga[1:m, 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
        Label(gb[1:m, 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
        for i in 1:(m + 1)
            vlines!(axs1[i], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs1[i], stop, color = (:gold, 0.5), linewidth = 0.5)
        end
        for i in 1:(n - m + 1)
            vlines!(axs2[i], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs2[i], stop, color = (:gold, 0.5), linewidth = 0.5)
        end
        for i in 1:m
            lines!(axs1[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        end
        for i in 1:(n - m)
            lines!(axs2[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        end
        if isodd(n)
            ax = Axis(ga[m + 2, 2])
            hidedecorations!(ax)
            hidespines!(ax)
        end
        Label(f[1, 1:2, Top()], "$(gene) locus", font = "Arial bold", textsize = 8, padding = (0, 0, 5, 0), valign = :bottom)
        rowgap!(ga, 5)
        colgap!(ga, 5)
        rowgap!(gb, 5)
        colgap!(gb, 5)
        colgap!(f.layout, 5)
        resize_to_layout!(f)
        save("figs/$(gene)-locuszoom.png", f, px_per_unit = 4)
    end
    run(`say "the job is finished"`)
end

loci = GeneticsMakie.findgwasloci(gwas[51])
filter!(x -> !(x.CHR == "6" && x.BP > 24e6 && x.BP < 34e6), loci)
filter!(x -> !(x.CHR == "5" && x.BP > 129e6 && x.BP < 132e6), loci)
filter!(x -> !(x.CHR == "8" && x.BP > 7e6 && x.BP < 13e6), loci)
loci1 = GeneticsMakie.findclosestgene(loci, gencode; proteincoding = true)
locuszoom(loci1.gene)