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
@time gencode = Arrow.Table("data/gencode/gencode.v39lift37.annotation.parsed.gtf.arrow")|> DataFrame
gencode.gene_id = gencode.gene_name

@info "Loading 1000 Genomes reference panel"
@time kgp_raw = SnpData("data/1kg/kgp.eur.maf0.05")
@assert size(kgp_raw) == (503, 7_230_618)

@info "Loading acne GWAS"
gwas = DataFrame(Arrow.Table(joinpath.("data/gwas/processed/acne.tsv.arrow")))

loci = GeneticsMakie.findgwasloci(gwas)
filter!(x -> !(x.CHR == "6" && x.BP > 26e6 && x.BP < 34e6), loci)
loci1 = GeneticsMakie.findclosestgene(loci, gencode; proteincoding = true)
loci2 = GeneticsMakie.findclosestgene(loci, gencode; proteincoding = true, start = true)
[loci1.gene loci2.gene]

begin
    @info "Plotting acne GWAS loci"
    n = nrow(loci)
    m = div(n, 2)
    f = Figure(resolution = (612, 792 * 2))
    ga = f[1, 1] = GridLayout()
    gb = f[1, 2] = GridLayout()
    axs1 = [Axis(ga[i, 1]) for i in 1:2m]
    axs2 = [Axis(gb[i, 1]) for i in 1:2(n - m)]
    window = 0.5e6
    for (i, gene) in enumerate(loci1.gene)
        @info "Working on $gene gene."
        chr, start, stop = GeneticsMakie.findgene(gene, gencode)
        range1, range2 = start - window, stop + window
        range1 < 0 ? range1 = 0 : nothing
        range2 > GeneticsMakie.GRCh37_totlength[chr] ? range2 = GeneticsMakie.GRCh37_totlength[chr] : nothing
        kgp = subsetref(kgp_raw, chr, range1, range2, "data/kgp.filtered")
        dfs = subsetgwas([gwas], chr, range1, range2)
        if 1 <= i <= m
            GeneticsMakie.plotlocus!(axs1[2i - 1], chr, range1, range2, dfs[1]; ld = kgp)
            if dfs[1].BP[findmin(dfs[1].P)[2]] < (range1 + range2) / 2
                Label(ga[2i - 1, 1, Top()], "$(gene) locus", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
            else
                Label(ga[2i - 1, 1, Top()], "$(gene) locus", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end    
            rowsize!(ga, 2i - 1, 30)
            rs = GeneticsMakie.plotgenes!(axs1[2i], chr, range1, range2, gencode; height = 0.1)
            rowsize!(ga, 2i, rs)
            GeneticsMakie.labelgenome(ga[2i, 1, Bottom()], chr, range1, range2)
            vlines!(axs1[2i - 1], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs1[2i - 1], stop, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs1[2i], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs1[2i], stop, color = (:gold, 0.5), linewidth = 0.5)
            lines!(axs1[2i - 1], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        else
            i -= m
            GeneticsMakie.plotlocus!(axs2[2i - 1], chr, range1, range2, dfs[1]; ld = kgp)
            if dfs[1].BP[findmin(dfs[1].P)[2]] < (range1 + range2) / 2
                Label(gb[2i - 1, 1, Top()], "$(gene) locus", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
            else
                Label(gb[2i - 1, 1, Top()], "$(gene) locus", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end    
            rowsize!(gb, 2i - 1, 30)
            rs = GeneticsMakie.plotgenes!(axs2[2i], chr, range1, range2, gencode; height = 0.1)
            rowsize!(gb, 2i, rs)
            GeneticsMakie.labelgenome(gb[2i, 1, Bottom()], chr, range1, range2)
            vlines!(axs2[2i - 1], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs2[2i - 1], stop, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs2[2i], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs2[2i], stop, color = (:gold, 0.5), linewidth = 0.5)    
            lines!(axs2[2i - 1], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        end
    end
    Colorbar(ga[1:2m, 2], limits = (0, 1), ticks = 0:1:1, height = 20,
        colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
        tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
        labelsize = 6, width = 5, spinewidth = 0.5)
    Colorbar(gb[1:2m, 2], limits = (0, 1), ticks = 0:1:1, height = 20,
        colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
        tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
        labelsize = 6, width = 5, spinewidth = 0.5)
    Label(ga[1:2m, 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
    Label(gb[1:2m, 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
    ga_fill = Axis(ga[2m + 1, 2])
    hidespines!(ga_fill)
    hidedecorations!(ga_fill)
    mhc = [Axis(gb[2(n - m) + i, 2]) for i in 1:2]
    @info "Plotting MHC"
    chr, range1, range2 = "6", 24e6, 36e6
    kgp = subsetref(kgp_raw, chr, range1, range2, "data/kgp.filtered")
    dfs = subsetgwas([gwas], chr, range1, range2)
    GeneticsMakie.plotlocus!(mhc[1], chr, range1, range2, dfs[1]; ld = kgp)
    rowsize!(gb, 2(n - m) + 1, 30)
    lines!(mhc[1], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
    Label(gb[2(n - m) + 1, 2, Top()], "MHC locus", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
    h = Float64[]
    _, c, _ = GeneticsMakie.findgene("H2AC1", gencode)
    push!(h, c)
    _, _, c = GeneticsMakie.findgene("MOG", gencode) # Extended class I
    _, s, _ = GeneticsMakie.findgene("ZFP57", gencode)
    push!(h, (c + s) / 2)
    _, _, c = GeneticsMakie.findgene("MICB", gencode) # Classical class I
    _, s, _ = GeneticsMakie.findgene("PPIAP9", gencode)
    push!(h, (c + s) / 2)
    _, _, c = GeneticsMakie.findgene("NOTCH4", gencode) # Classical class III
    _, s, _ = GeneticsMakie.findgene("TSBP1", gencode)
    push!(h, (c + s) / 2)
    _, _, c = GeneticsMakie.findgene("HCG24", gencode) # Classical class II
    _, s, _ = GeneticsMakie.findgene("COL11A2", gencode)
    push!(h, (c + s) / 2)
    _, c, _ = GeneticsMakie.findgene("RPL12P1", gencode) # Extended class II
    push!(h, c)
    pushfirst!(h, range1)
    push!(h, range2)
    # https://bioinformatics.stackexchange.com/questions/14649/what-are-the-coordinates-of-the-extended-major-histocompatibility-complex-xmhc
    colors = ["gray95", "#4062D8", "#CB3C33", "#389826", "#9658B2", "#4062D8", "gray95"]
    for i in 1:(length(h) - 1)
        poly!(mhc[2], Point2f[(h[i], 0), (h[i + 1], 0), (h[i + 1], 1), (h[i], 1)], color = colors[i], strokewidth = 0)
    end
    xlims!(mhc[2], range1, range2)
    ylims!(mhc[2], 0, 1)
    hidedecorations!(mhc[2])
    hidespines!(mhc[2], :t, :l, :r)
    mhc[2].spinewidth = 0.75
    rowsize!(gb, 2(n - m) + 2, 5)
    GeneticsMakie.labelgenome(gb[2(n - m) + 2, 2, Bottom()], chr, range1, range2)
    Legend(gb[2(n - m) + 3, 2], [PolyElement(color = colors[i], strokecolor = :transparent) for i in 2:5], 
        ["Extended MHC region", "Class I region", "Class III region", "Class II region"],
        rowgap = 0, labelsize = 6, tellheight = true, tellwidth = false, orientation = :horizontal, nbanks = 2,
        framevisible = false, patchsize = (3, 3), strokewidth = 0.1)
    Label(f[1, 1:2, Top()], "Acne GWAS loci", font = "Arial bold", textsize = 8, padding = (0, 0, 5, 0), valign = :bottom)
    rowgap!(ga, 5)
    colgap!(ga, 5)
    rowgap!(gb, 5)
    colgap!(gb, 5)
    colgap!(f.layout, 2.5)
    resize_to_layout!(f)
    save("figs/acne-locuszoom.png", f, px_per_unit = 4)
    run(`say "the job is finished"`)    
end