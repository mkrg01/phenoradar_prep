#!/usr/bin/env Rscript
# MonoPhy supplies all monophyly and intruder/outlier decisions. This wrapper
# handles missing rank annotations, exports, and display-only tree folding.
suppressPackageStartupMessages(library(MonoPhy))
request <- jsonlite::fromJSON(commandArgs(trailingOnly=TRUE)[1])
if (as.character(packageVersion("MonoPhy")) != request$version)
    stop("This workflow requires MonoPhy ", request$version, "; run monophy.post-deploy.sh")
settings <- request$settings
outdir <- request$outdir
taxonomy <- read.delim(file.path(outdir, "taxonomy.tsv"), check.names=FALSE,
                       stringsAsFactors=FALSE, na.strings="", quote='"')
tree <- ape::read.tree(request$tree)
if (!ape::is.rooted(tree)) stop("MonoPhy requires a rooted species tree")
if (!setequal(tree$tip.label, taxonomy$tip)) stop("R tree/taxonomy tip mismatch")

write_tsv <- function(x, path) {
    write.table(x, path, sep="\t", quote=TRUE, row.names=FALSE, na="")
}
empty <- function(fields) {
    as.data.frame(setNames(rep(list(character()), length(fields)), fields),
                  stringsAsFactors=FALSE)
}
tips_for <- function(x, taxon) {
    sort(unique(as.character(unlist(x[[taxon]], use.names=FALSE))))
}
events <- empty(c("candidate_id", "rank", "focal_taxon", "species", "role"))
groups <- empty(c("rank", "taxon", "monophyly", "mrca_node", "tip_count",
                  "intruder_tips", "outlier_tips", "missing_rank_tips_omitted"))
members <- empty(c("rank", "focal_taxon", "species", "role"))
states <- empty(c("species", "rank", "registered_taxon", "status", "intruder", "outlier"))
plot_members <- empty(c("rank", "species", "display_tip", "registered_taxon", "color", "symbol"))
warnings_seen <- character()

# Publication figures: native MonoPhy roles, measured text, vector graphics.
# Tree geometry is a cladogram. No titles, subtitles or explanatory prose are
# drawn; role definitions and full taxids remain in the report and documentation.
plot_rank <- function(rank, native, rank_events, folder) {
    labels <- setNames(taxonomy[[rank]], taxonomy$tip)
    taxa <- sort(unique(labels[!is.na(labels)]))
    palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#332288", "#882255", "#117733")
    if (length(taxa) > length(palette))
        palette <- grDevices::hcl(h=seq(15, 375, length.out=length(taxa)+1)[seq_along(taxa)], c=60, l=45)
    colors <- setNames(palette[seq_along(taxa)], taxa)
    keep <- tree$tip.label
    representative <- setNames(tree$tip.label, tree$tip.label)
    counts <- setNames(rep(1L, length(keep)), keep)
    if (settings$collapse_monophyletic && !is.null(native)) {
        for (taxon in rownames(native$result)) {
            species <- sort(names(labels)[!is.na(labels) & labels == taxon])
            if (native$result[taxon, "Monophyly"] != "Yes" || length(species) < 2L ||
                any(species %in% rank_events$species)) next
            first <- species[1]
            keep <- setdiff(keep, species[-1])
            representative[species] <- first
            counts[first] <- length(species)
        }
    }
    display <- if (length(keep) == length(tree$tip.label)) tree else ape::keep.tip(tree, keep)
    color <- setNames(rep("#777777", length(labels)), names(labels))
    color[!is.na(labels)] <- colors[labels[!is.na(labels)]]
    intruder <- names(labels) %in% rank_events$species[rank_events$role == "intruder"]
    outlier <- names(labels) %in% rank_events$species[rank_events$role == "outlier"]
    symbol <- setNames(ifelse(intruder & outlier, "both", ifelse(intruder, "intruder",
                      ifelse(outlier, "outlier", "none"))), names(labels))
    plot_members <<- rbind(plot_members, data.frame(rank=rank, species=names(labels),
        display_tip=unname(representative[names(labels)]), registered_taxon=unname(labels),
        color=unname(color), symbol=unname(symbol)))
    display$node.label <- NULL
    display$edge.length <- NULL
    ape::write.tree(display, file=file.path(folder, "display_tree.nwk"))

    shown <- display$tip.label
    n <- length(shown)
    collapsed <- counts[shown] > 1L
    taxon_text <- sub(" \\[[0-9]+\\]$", "", labels[shown])
    taxon_text[is.na(taxon_text)] <- ""
    scientific_names <- unlist(request$scientific_names, use.names=TRUE)
    main_text <- ifelse(collapsed, taxon_text, scientific_names[shown])
    count_text <- ifelse(collapsed, paste0(" (n = ", counts[shown], ")"), "")
    taxon_text[collapsed] <- ""
    fonts <- ifelse(collapsed & !rank %in% c("genus", "subgenus", "species", "subspecies"), 1L, 3L)
    font_family <- "DejaVu Sans"
    font_size <- 9
    if (!capabilities("cairo")) stop("Publication figures require R Cairo graphics support")
    # Measure with the same font/device used for the final vector figures.
    measurement <- file.path(folder, ".font-metrics.pdf")
    grDevices::cairo_pdf(measurement, width=7, height=7, family=font_family, pointsize=font_size)
    main_widths <- vapply(seq_len(n), function(i)
        strwidth(main_text[i], units="inches", font=fonts[i]), numeric(1))
    count_widths <- strwidth(count_text, units="inches", font=1, cex=0.9)
    taxon_width <- max(strwidth(taxon_text, units="inches", font=1, cex=0.9))
    grDevices::dev.off()
    unlink(measurement)

    margin <- 0.10
    row_pitch <- 0.19
    depth <- ape::node.depth(display, method=2)
    tree_width <- max(1.5, min(2.5, (max(depth)-1)*0.3))
    marker_x <- margin + tree_width + 0.07
    label_x <- marker_x + 0.12
    group_x <- label_x + max(main_widths + count_widths) + 0.22
    width <- if (all(taxon_text == "")) group_x - 0.22 + margin else group_x + taxon_width + margin
    height <- 2*margin + (n-1)*row_pitch + font_size/72
    x <- margin + (max(depth) - depth)/max(1, max(depth)-1)*tree_width
    y <- rep(NA_real_, length(depth))
    y[seq_len(n)] <- margin + font_size/144 + rev(seq(0, by=row_pitch, length.out=n))
    children <- split(display$edge[, 2], display$edge[, 1])
    postorder <- unique(ape::reorder.phylo(display, "postorder")$edge[, 1])
    for (node in postorder) y[node] <- mean(range(y[children[[as.character(node)]]]))
    draw <- function() {
        par(mai=rep(0, 4), family=font_family, ps=font_size, cex=1, xpd=NA,
            lend="butt", ljoin="mitre")
        plot.new()
        plot.window(xlim=c(0, width), ylim=c(0, height), xaxs="i", yaxs="i")
        for (node in postorder) {
            child <- children[[as.character(node)]]
            segments(x[node], min(y[child]), x[node], max(y[child]), col="#333333", lwd=0.75)
            segments(x[node], y[child], x[child], y[child], col="#333333", lwd=0.75)
        }
        flags <- symbol[shown] != "none"
        points(rep(marker_x, n), y[seq_len(n)],
               pch=c(none=21L, intruder=24L, outlier=22L, both=23L)[symbol[shown]],
               bg=color[shown], col=ifelse(flags, "#222222", color[shown]),
               cex=ifelse(flags, 0.85, 0.45), lwd=0.65)
        for (i in seq_len(n)) {
            text(label_x, y[i], main_text[i], adj=c(0, 0.5), font=fonts[i],
                 col=if (collapsed[i]) color[shown[i]] else "#111111")
            text(label_x+main_widths[i], y[i], count_text[i], adj=c(0, 0.5),
                 font=1, cex=0.9, col="#555555")
            text(group_x, y[i], taxon_text[i], adj=c(0, 0.5), font=1, cex=0.9, col=color[shown[i]])
        }
    }
    for (format in c("pdf", "svg")) {
        device <- if (format == "pdf") grDevices::cairo_pdf else grDevices::svg
        device(file.path(folder, paste0("tree.", format)), width=width, height=height,
               pointsize=font_size, family=font_family, bg="white")
        tryCatch(draw(), finally=grDevices::dev.off())
    }
    jsonlite::write_json(list(width_inches=width, height_inches=height,
        font_family=font_family, label_size_pt=font_size, taxonomy_size_pt=font_size*0.9,
        row_pitch_inches=row_pitch, topology="cladogram", collapsed_tip_counts=as.list(counts[shown]),
        titles=FALSE, subtitles=FALSE, footer=FALSE), file.path(folder, "figure.json"),
        pretty=TRUE, auto_unbox=TRUE)
}

for (rank in settings$ranks) {
    message("MonoPhy species-tree assessment: ", rank)
    folder <- file.path(outdir, "ranks", rank)
    dir.create(folder, recursive=TRUE)
    known <- taxonomy$tip[!is.na(taxonomy[[rank]])]
    rank_events <- events[FALSE, ]
    native <- NULL
    status <- setNames(rep("missing_taxonomy_rank", nrow(taxonomy)), taxonomy$tip)
    if (length(known) < 2L) {
        status[known] <- "monotypic"
        if (length(known) == 1L) {
            taxon <- taxonomy[[rank]][match(known, taxonomy$tip)]
            groups <- rbind(groups, data.frame(rank=rank, taxon=taxon, monophyly="Monotypic",
                mrca_node=NA, tip_count=1L, intruder_tips=0L, outlier_tips=0L,
                missing_rank_tips_omitted=nrow(taxonomy)-1L))
            members <- rbind(members, data.frame(rank=rank, focal_taxon=taxon,
                                                species=known, role="member"))
        }
        writeLines("Fewer than two tips have this rank; MonoPhy was not run.",
                   file.path(folder, "not_assessed.txt"))
    } else {
        assessed_tree <- if (length(known) == length(tree$tip.label)) tree else ape::keep.tip(tree, known)
        if (!ape::is.rooted(assessed_tree))
            stop("Rank pruning left an unresolved root for ", rank,
                 "; inspect the source root before assessing monophyly")
        rank_taxonomy <- data.frame(tip=known, group=taxonomy[[rank]][match(known, taxonomy$tip)])
        names(rank_taxonomy)[2] <- rank
        solution <- withCallingHandlers(
            MonoPhy::AssessMonophyly(assessed_tree, taxonomy=rank_taxonomy,
                outliercheck=TRUE, outlierlevel=settings$outlierlevel),
            warning=function(w) { warnings_seen <<- c(warnings_seen, paste(rank, conditionMessage(w))) })
        native <- solution[[1]]
        saveRDS(solution, file.path(folder, "monophy.rds"))
        ape::write.tree(assessed_tree, file=file.path(folder, "assessment_tree.nwk"), digits=17)
        write_tsv(data.frame(taxon=rownames(native$result), native$result, check.names=FALSE),
                  file.path(folder, "native_results.tsv"))
        write_tsv(native$TipStates, file.path(folder, "native_tip_states.tsv"))
        for (taxon in rownames(native$result)) {
            own <- known[rank_taxonomy[[rank]] == taxon]
            outsiders <- tips_for(native$IntruderTips, taxon)
            outliers <- tips_for(native$OutlierTips, taxon)
            mono <- as.character(native$result[taxon, "Monophyly"])
            status[own] <- switch(mono, Yes="monophyletic", No="non_monophyletic", Monotypic="monotypic")
            groups <- rbind(groups, data.frame(rank=rank, taxon=taxon, monophyly=mono,
                mrca_node=native$result[taxon, "MRCA"], tip_count=length(own),
                intruder_tips=length(outsiders), outlier_tips=length(outliers),
                missing_rank_tips_omitted=nrow(taxonomy)-length(known)))
            members <- rbind(members, data.frame(rank=rank, focal_taxon=taxon, species=own,
                            role=ifelse(own %in% outliers, "outlier", "member")))
            if (length(outsiders))
                members <- rbind(members, data.frame(rank=rank, focal_taxon=taxon,
                                                     species=outsiders, role="intruder"))
            for (role in c("intruder", "outlier")) {
                flagged <- if (role == "intruder") outsiders else outliers
                if (length(flagged)) {
                    ids <- sprintf("candidate_%06d", nrow(events) + nrow(rank_events) + seq_along(flagged))
                    rank_events <- rbind(rank_events, data.frame(candidate_id=ids, rank=rank,
                        focal_taxon=taxon, species=flagged, role=role))
                }
            }
        }
    }
    states <- rbind(states, data.frame(species=taxonomy$tip, rank=rank,
        registered_taxon=taxonomy[[rank]], status=unname(status[taxonomy$tip]),
        intruder=taxonomy$tip %in% rank_events$species[rank_events$role == "intruder"],
        outlier=taxonomy$tip %in% rank_events$species[rank_events$role == "outlier"]))
    events <- rbind(events, rank_events)
    plot_rank(rank, native, rank_events, folder)
}
write_tsv(events, file.path(outdir, "events.tsv"))
write_tsv(groups, file.path(outdir, "taxon_results.tsv"))
write_tsv(members, file.path(outdir, "group_members.tsv"))
write_tsv(states, file.path(outdir, "rank_status.tsv"))
write_tsv(plot_members, file.path(outdir, "plot_members.tsv"))
capture.output(sessionInfo(), file=file.path(outdir, "R_session.txt"))
packages <- c("MonoPhy", "ape", "phytools", "phangorn", "RColorBrewer", "jsonlite")
engine <- list(name="MonoPhy", version=as.character(packageVersion("MonoPhy")),
               R=R.version.string, packages=setNames(lapply(packages, function(p)
               as.character(packageVersion(p))), packages), warnings=unique(warnings_seen),
               outliercheck=TRUE, missing_rank_policy="omit separately at each rank",
               plots="ape; registered taxonomy colors and native MonoPhy role symbols")
jsonlite::write_json(engine, file.path(outdir, "engine.json"), pretty=TRUE, auto_unbox=TRUE)
