using GeneticsMakie, CairoMakie, CSV, DataFrames, SnpArrays, Arrow, StatsBase

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
@time kgp_raw = SnpData("data/1kg/kgp.eur.maf0.05.geno")
@assert size(kgp_raw) == (503, 7_230_539)

@info "Loading acne GWAS"
gwas = DataFrame(Arrow.Table(joinpath.("data/gwas/processed/acne.tsv.arrow")))

loci = GeneticsMakie.findgwasloci(gwas)
for i in eachindex(GeneticsMakie.GRCh37_highld)
    storage = deepcopy(loci)
    chr, range1, range2 = GeneticsMakie.GRCh37_highld[i].chr, GeneticsMakie.GRCh37_highld[i].range1, GeneticsMakie.GRCh37_highld[i].range2
    storage.region .= "$(chr):$(range1)-$(range2)"
    ld = filter(x -> x.CHR == GeneticsMakie.GRCh37_highld[i].chr && 
        x.BP > GeneticsMakie.GRCh37_highld[i].range1 && 
        x.BP < GeneticsMakie.GRCh37_highld[i].range2, storage)
    nrow(ld) > 0 ? println(ld) : nothing
end

filter!(x -> !(x.CHR == "6" && x.BP > 24e6 && x.BP < 34e6), loci)
filter!(x -> !(x.CHR == "5" && x.BP > 129e6 && x.BP < 132e6), loci)
filter!(x -> !(x.CHR == "8" && x.BP > 7e6 && x.BP < 13e6), loci)

loci1 = GeneticsMakie.findclosestgene(loci, gencode; proteincoding = true)
loci2 = GeneticsMakie.findclosestgene(loci, gencode; proteincoding = true, start = true)
[loci1.gene loci2.gene]

begin
    @info "Plotting acne GWAS loci"
    n = nrow(loci)
    m = div(n, 4)
    f = Figure(resolution = (612 * 2, 792 * 2))
    ga = f[1, 1] = GridLayout()
    gb = f[1, 2] = GridLayout()
    gc = f[1, 3] = GridLayout()
    gd = f[1, 4] = GridLayout()
    axs1 = [Axis(ga[i, 1]) for i in 1:(2 * 10)]
    axs2 = [Axis(gb[i, 1]) for i in 1:(2 * 12)]
    axs3 = [Axis(gc[i, 1]) for i in 1:(2 * 11)]
    axs4 = [Axis(gd[i, 1]) for i in 1:(2 * 6)]
    window = 0.5e6
    for (i, gene) in enumerate(loci1.gene)
        @info "Working on $gene gene."
        chr, start, stop = GeneticsMakie.findgene(gene, gencode)
        range1, range2 = start - window, stop + window
        range1 < 0 ? range1 = 0 : nothing
        range2 > GeneticsMakie.GRCh37_totlength[chr] ? range2 = GeneticsMakie.GRCh37_totlength[chr] : nothing
        kgp = subsetref(kgp_raw, chr, range1, range2, "data/kgp.filtered")
        dfs = subsetgwas([gwas], chr, range1, range2)
        if 1 <= i <= 10
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
        elseif 10 < i <= 22
            i -= 10
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
        elseif 22 < i <= 33
            i -= 22
            GeneticsMakie.plotlocus!(axs3[2i - 1], chr, range1, range2, dfs[1]; ld = kgp)
            if dfs[1].BP[findmin(dfs[1].P)[2]] < (range1 + range2) / 2
                Label(gc[2i - 1, 1, Top()], "$(gene) locus", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
            else
                Label(gc[2i - 1, 1, Top()], "$(gene) locus", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end    
            rowsize!(gc, 2i - 1, 30)
            rs = GeneticsMakie.plotgenes!(axs3[2i], chr, range1, range2, gencode; height = 0.1)
            rowsize!(gc, 2i, rs)
            GeneticsMakie.labelgenome(gc[2i, 1, Bottom()], chr, range1, range2)
            vlines!(axs3[2i - 1], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs3[2i - 1], stop, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs3[2i], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs3[2i], stop, color = (:gold, 0.5), linewidth = 0.5)    
            lines!(axs3[2i - 1], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        else
            i -= 33
            GeneticsMakie.plotlocus!(axs4[2i - 1], chr, range1, range2, dfs[1]; ld = kgp)
            if dfs[1].BP[findmin(dfs[1].P)[2]] < (range1 + range2) / 2
                Label(gd[2i - 1, 1, Top()], "$(gene) locus", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
            else
                Label(gd[2i - 1, 1, Top()], "$(gene) locus", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end    
            rowsize!(gd, 2i - 1, 30)
            rs = GeneticsMakie.plotgenes!(axs4[2i], chr, range1, range2, gencode; height = 0.1)
            rowsize!(gd, 2i, rs)
            GeneticsMakie.labelgenome(gd[2i, 1, Bottom()], chr, range1, range2)
            vlines!(axs4[2i - 1], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs4[2i - 1], stop, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs4[2i], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs4[2i], stop, color = (:gold, 0.5), linewidth = 0.5)    
            lines!(axs4[2i - 1], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        end        
    end
    Colorbar(gd[1:2m, 2], limits = (0, 1), ticks = 0:1:1, height = 20,
        colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
        tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
        labelsize = 6, width = 5, spinewidth = 0.5)
    Label(ga[1:2(m + 1), 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
    long = [Axis(gd[i, 1]) for i in 13:18]
    @info "Plotting 5q31"
    chr, range1, range2 = "5", 129e6, 132e6
    kgp = subsetref(kgp_raw, chr, range1, range2, "data/kgp.filtered")
    dfs = subsetgwas([gwas], chr, range1, range2)
    GeneticsMakie.plotlocus!(long[1], chr, range1, range2, dfs[1]; ld = kgp)
    gene = "SLC22A4"
    _, start, stop = GeneticsMakie.findgene(gene, gencode)
    if dfs[1].BP[findmin(dfs[1].P)[2]] < (range1 + range2) / 2
        Label(gd[13, 1, Top()], "5q31 locus", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
    else
        Label(gd[13, 1, Top()], "5q31 locus", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
    end    
    rowsize!(gd, 13, 30)
    rs = GeneticsMakie.plotgenes!(long[2], chr, range1, range2, gencode; height = 0.1)
    rowsize!(gd, 14, rs)
    GeneticsMakie.labelgenome(gd[14, 1, Bottom()], chr, range1, range2)
    vlines!(long[1], start, color = (:gold, 0.5), linewidth = 0.5)
    vlines!(long[1], stop, color = (:gold, 0.5), linewidth = 0.5)
    vlines!(long[2], start, color = (:gold, 0.5), linewidth = 0.5)
    vlines!(long[2], stop, color = (:gold, 0.5), linewidth = 0.5)    
    lines!(long[1], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
    @info "Plotting 8p23"
    chr, range1, range2 = "8", 8011830, 11727042
    kgp = subsetref(kgp_raw, chr, range1, range2, "data/kgp.filtered")
    dfs = subsetgwas([gwas], chr, range1, range2)
    GeneticsMakie.plotlocus!(long[3], chr, range1, range2, dfs[1]; ld = kgp)
    gene = "PINX1"
    _, start, stop = GeneticsMakie.findgene(gene, gencode)
    if dfs[1].BP[findmin(dfs[1].P)[2]] < (range1 + range2) / 2
        Label(gd[15, 1, Top()], "8p23 locus", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
    else
        Label(gd[15, 1, Top()], "8p23 locus", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
    end    
    rowsize!(gd, 15, 30)
    rs = GeneticsMakie.plotgenes!(long[4], chr, range1, range2, gencode; height = 0.1)
    rowsize!(gd, 16, rs)
    GeneticsMakie.labelgenome(gd[16, 1, Bottom()], chr, range1, range2)
    vlines!(long[3], start, color = (:gold, 0.5), linewidth = 0.5)
    vlines!(long[3], stop, color = (:gold, 0.5), linewidth = 0.5)
    vlines!(long[4], start, color = (:gold, 0.5), linewidth = 0.5)
    vlines!(long[4], stop, color = (:gold, 0.5), linewidth = 0.5)    
    lines!(long[3], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
    @info "Plotting MHC"
    chr, range1, range2 = "6", 24e6, 36e6
    kgp = subsetref(kgp_raw, chr, range1, range2, "data/kgp.filtered")
    dfs = subsetgwas([gwas], chr, range1, range2)
    GeneticsMakie.plotlocus!(long[5], chr, range1, range2, dfs[1]; ld = kgp)
    rowsize!(gd, 17, 30)
    lines!(long[5], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
    Label(gd[17, 1, Top()], "MHC locus", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
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
        poly!(long[6], Point2f[(h[i], 0), (h[i + 1], 0), (h[i + 1], 1), (h[i], 1)], color = colors[i], strokewidth = 0)
    end
    xlims!(long[6], range1, range2)
    ylims!(long[6], 0, 1)
    hidedecorations!(long[6])
    hidespines!(long[6], :t, :l, :r)
    long[6].spinewidth = 0.75
    rowsize!(gd, 18, 5)
    GeneticsMakie.labelgenome(gd[18, 1, Bottom()], chr, range1, range2)
    Legend(gd[19, 1], [PolyElement(color = colors[i], strokecolor = :transparent) for i in 2:5], 
        ["Extended MHC region", "Class I region", "Class III region", "Class II region"],
        rowgap = 0, labelsize = 6, tellheight = true, tellwidth = false, orientation = :horizontal, nbanks = 2,
        framevisible = false, patchsize = (3, 3), strokewidth = 0.1)
    gb_fill = Axis(gb[25, 1])
    hidespines!(gb_fill)
    hidedecorations!(gb_fill)
    gc_fill = Axis(gc[23, 1])
    hidespines!(gc_fill)
    hidedecorations!(gc_fill)
    gd_fill = Axis(gd[20, 1])
    hidespines!(gd_fill)
    hidedecorations!(gd_fill)
    Label(f[1, 1:4, Top()], "42 acne GWAS loci", font = "Arial bold", textsize = 8, padding = (0, 0, 5, 0), valign = :bottom)
    rowgap!(ga, 5)
    colgap!(ga, 5)
    rowgap!(gb, 5)
    colgap!(gb, 5)
    rowgap!(gc, 5)
    colgap!(gc, 5)
    rowgap!(gd, 5)
    colgap!(gd, 5)
    colgap!(f.layout, 10)
    resize_to_layout!(f)
    save("figs/acne-locuszoom.png", f, px_per_unit = 4)
    run(`say "the job is finished"`)    
end