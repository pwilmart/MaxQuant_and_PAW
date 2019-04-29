
# load libraries
library("tidyverse")
library("psych")
library("gridExtra")
library("scales")
library("limma") 
library("edgeR") 

# read the grouped protein summary file
paw_raw <- read_tsv("grouped_protein_summary_TMT_8_sorted.txt", skip = 4,
                    n_max = 5427, guess_max = 5427)

# extract protein accession column and TMT data
# need to exclude any rows with an entry in "Filter" column
# there are also 4 unused channels (of the 11-plex slots)
paw_tmt <- paw_raw %>% 
  filter(., is.na(Filter)) %>%
  select(Accession, starts_with("TotInt_")) %>%
  select(-contains("127N"), -contains("128N"), -contains("130N"), -contains("131C"))

# separate accessions from the data
accession <- paw_tmt$Accession
paw_tmt <- paw_tmt %>% select(-Accession)

# rename the columns and gather the conditions
colnames(paw_tmt) <- c("Media_2.1", "Exo_2.1", "Exo_2.2", "Media_2.2",
                      "Exo_3.1", "Exo_3.2", "Media_3.2")
paw_tmt <- paw_tmt %>% select(contains("Media"), contains("Exo"))

head(paw_tmt)
nrow(paw_tmt)

# load data into DGEList object
group <- c(rep("media", 3), rep("exosome", 4))
y <- DGEList(counts = paw_tmt, group = group, genes = accession)
y$samples

# run the TMM normalization
y <- calcNormFactors(y)

# set colors for plotting
color <- c(rep("red", 3), rep("blue", 4))

# function to compute the normalized intensities
apply_tmm_factors <- function(y, color = NULL, plot = TRUE) {
    # computes the tmm normalized data from the DGEList object
        # y - DGEList object
        # returns a dataframe with normalized intensities
    
    # compute grand total (library size) scalings
    lib_facs <- mean(y$samples$lib.size) / y$samples$lib.size

    # the TMM factors are library adjustment factors (so divide by them)
    norm_facs <- lib_facs / y$samples$norm.factors
    cat("Overall Factors (lib.size+TMM):\n", sprintf("%-5s -> %f\n", 
                                                     colnames(y$counts), norm_facs))

    # compute the normalized data as a new data frame
    tmt_tmm <- as.data.frame(sweep(y$counts, 2, norm_facs, FUN = "*"))
    colnames(tmt_tmm) <- str_c(colnames(y$counts), "_tmm")
    
    # visualize results and return data frame
    if(plot == TRUE) {
        boxplot(log10(tmt_tmm), col = color, notch = TRUE, main = "TMM Normalized data")
    }
    tmt_tmm
}

# get the normalized data values
paw_tmt_tmm <- apply_tmm_factors(y, color)

# we need to get dispersion estimates
y <- estimateDisp(y)
plotBCV(y, main = "Dispersion trends")

collect_results <- function(df, tt, x, xlab, y, ylab) {
    # Computes new columns and extracts some columns to make results frame
        # df - data in data.frame
        # tt - top tags from edgeR test
        # x - columns for first condition
        # xlab - label for x
        # y - columns for second condition
        # ylab - label for y
        # returns a new dataframe
    
    # condition average vectors
    ave_x <- rowMeans(df[x])
    ave_y <- rowMeans(df[y])
    
    # FC, direction, candidates
    fc <- ifelse(ave_y > ave_x, (ave_y / ave_x), (-1 * ave_x / ave_y))
    direction <- ifelse(ave_y > ave_x, "up", "down")
    candidate <- cut(tt$FDR, breaks = c(-Inf, 0.01, 0.05, 0.10, 1.0), 
                     labels = c("high", "med", "low", "no"))

    # make data frame
    temp <- cbind(df[c(x, y)], data.frame(logFC = tt$logFC, FC = fc, 
                                          PValue = tt$PValue, FDR = tt$FDR, 
                                          ave_x = ave_x, ave_y = ave_y, 
                                          direction = direction, candidate = candidate, 
                                          Acc = tt$genes)) 
    
    # fix column headers for averages
    names(temp)[names(temp) %in% c("ave_x", "ave_y")]  <- str_c("ave_", c(xlab, ylab))    
    
    temp # return the data frame
}

# compute the exact test models, p-values, FC, etc.
et <- exactTest(y, pair = c("media", "exosome"))

# define variables for the columns in each condition
M <- 1:3
E <- 4:7

# make the results table 
tt <- topTags(et, n = Inf, sort.by = "none")$table
med_exo <- collect_results(paw_tmt_tmm, tt, M, "med", E, "exo")

pvalue_plot <- function(results, title) {
    # Makes p-value distribution plots
        # results - results data frame
        # title - plot title
    ggplot(results, aes(PValue)) + 
        geom_histogram(bins = 100, fill = "white", color = "black") +
        geom_hline(yintercept = mean(hist(results$PValue, breaks = 100, 
                                     plot = FALSE)$counts[26:100])) +
        ggtitle(str_c(title, " p-value distribution"))
}

# check the p-value distrubution
pvalue_plot(med_exo, "Media vs Exosome-dosed")

# see how many up and down candidates (10% FDR)
summary(decideTests(et, p.value = 0.10))

# see which proteins have the smallest p-values
topTags(et)$table

# see how many candidates are in each category
med_exo %>% count(candidate)

log2FC_plots <- function(results, range, title) {
    # Makes faceted log2FC plots by candidate
        # results - results data frame
        # range - plus/minus log2 x-axis limits
        # title - plot title
    ggplot(results, aes(x = logFC, fill = candidate)) +
        geom_histogram(binwidth=0.1, color = "black") +
        facet_wrap(~candidate) +
        ggtitle(title) + 
        coord_cartesian(xlim = c(-range, range))
}

# can look at log2FC distributions as a check
log2FC_plots(med_exo, 3, "LogFC by candidate for Media vs Exosome-dosed")

transform <- function(results, x, y) {
    # Make data frame with some transformed columns
        # results - results data frame
        # x - columns for x condition
        # y - columns for y condition
        # return new data frame
    df <- data.frame(log10((results[x] + results[y])/2), 
                     log2(results[y] / results[x]), 
                     results$candidate,
                     -log10(results$FDR))
    colnames(df) <- c("A", "M", "candidate", "P")
    df # return the frame
}

MA_plots <- function(results, x, y, title, make_facet = TRUE) {
    # makes MA-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # make_facet - flag to plot facet views
        # returns a list of plots 
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # 2-fold change lines
    ma_lines <- list(geom_hline(yintercept = 0.0, color = "black"),
                     geom_hline(yintercept = 1.0, color = "black", linetype = "dotted"),
                     geom_hline(yintercept = -1.0, color = "black", linetype = "dotted"))

    # make main MA plot
    ma <- ggplot(temp, aes(x = A, y = M)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scale_y_continuous(paste0("logFC (", y, "/", x, ")")) +
        scale_x_continuous("Ave_intensity") +
        ggtitle(title) + 
        ma_lines
    
    # make separate MA plots
    if (make_facet == TRUE) {
        ma_facet <- ggplot(temp, aes(x = A, y = M)) +
            geom_point(aes(color = candidate, shape = candidate)) +
            scale_y_continuous(paste0("log2 FC (", y, "/", x, ")")) +
            scale_x_continuous("log10 Ave_intensity") +
            ma_lines +
            facet_wrap(~ candidate) +
            ggtitle(str_c(title, " (separated)"))
    }

    # make the plots visible
    print(ma)
    if (make_facet == TRUE) {
        print(ma_facet)
    }
}    

# MA plots of DE candidates
MA_plots(med_exo, "ave_med", "ave_exo", "Media versus exosome-dosed")

scatter_plots <- function(results, x, y, title, make_facet = TRUE) {
    # makes scatter-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # make_facet - flag to plot facet views
        # returns a list of plots
    
    # 2-fold change lines
    scatter_lines <- list(geom_abline(intercept = 0.0, slope = 1.0, color = "black"),
                          geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          scale_y_log10(),
                          scale_x_log10())

    # make main scatter plot
    scatter <- ggplot(results, aes_string(x, y)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        ggtitle(title) + 
        scatter_lines

    # make separate scatter plots
    if (make_facet == TRUE) {
        scatter_facet <- ggplot(results, aes_string(x, y)) +
            geom_point(aes(color = candidate, shape = candidate)) +
            scatter_lines +
            facet_wrap(~ candidate) +
            ggtitle(str_c(title, " (separated)")) 
    }

    # make the plots visible
    print(scatter)
    if (make_facet == TRUE) {
        print(scatter_facet)
    }
}

scatter_plots(med_exo, "ave_med", "ave_exo", "Media versus exosome-dosed (edgeR)")

volcano_plot <- function(results, x, y, title, ymax) {
    # makes a volcano plot
        # results - a data frame with edgeR results
        # x - string for the x-axis column
        # y - string for y-axis column
        # ymax - upper limit for y-axis
        # title - plot title string
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # build the plot
    ggplot(temp, aes(x = M, y = P)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        xlab("log2 FC") +
        ylab("-log10 FDR") +
        coord_cartesian(xlim = c(-5, 5), ylim = c(0, ymax)) + 
        ggtitle(str_c(title, " Volcano Plot"))
}

# finally, a volcano plot
volcano_plot(med_exo, "ave_med", "ave_exo", "Media versus exosome-dosed (edgeR)", 50)

# function to extract the identifier part of the accesssion
get_identifier <- function(accession) {
    identifier <- str_split(accession, "\\|", simplify = TRUE)
    identifier[,3]
}

set_plot_dimensions <- function(width_choice, height_choice) {
    options(repr.plot.width=width_choice, repr.plot.height=height_choice)
}

plot_top_tags <- function(results, nleft, nright, top_tags) {
    # results should have data first, then test results (two condition summary table)
    # nleft, nright are number of data points in each condition
    # top_tags is number of up and number of down top DE candidates to plot
    # get top up-regulated
    up <- results %>% 
        filter(logFC >= 0) %>%
        arrange(FDR)
    up <- up[1:top_tags, ]
    
    # get top down-regulated
    down <- results %>% 
        filter(logFC < 0) %>%
        arrange(FDR)
    down <- down[1:top_tags, ]
    
    # pack them
    proteins <- rbind(up, down)
        
    color = c(rep("red", nleft), rep("blue", nright))
    for (row_num in 1:nrow(proteins)) {
        row <- proteins[row_num, ]
        vec <- as.vector(unlist(row[1:(nleft + nright)]))
        names(vec) <- colnames(row[1:(nleft + nright)])
        title <- str_c(get_identifier(row$Acc), ", int: ", scientific(mean(vec), 2), 
                       ", p-val: ", scientific(row$FDR, digits = 3), 
                       ", FC: ", round(row$FC, digits = 1))
        barplot(vec, col = color, main = title)
    }    
}

# plot the top 10 up and 10 down proteins
set_plot_dimensions(7, 4)
plot_top_tags(med_exo, 3, 4, 10)
set_plot_dimensions(7, 7)

# copy the data
ttest_PAW <- paw_tmt_tmm

# add average ratio columns (non-logged ratios), fold-change column, and row names
ttest_PAW$ave_med <- rowMeans(paw_tmt_tmm[1:3])
ttest_PAW$ave_exo  <- rowMeans(paw_tmt_tmm[4:7])
ttest_PAW$logFC <- log2(ttest_PAW$ave_exo / ttest_PAW$ave_med)
ttest_PAW$FC <- ifelse(ttest_PAW$ave_exo > ttest_PAW$ave_med, 
                          (ttest_PAW$ave_exo / ttest_PAW$ave_med), 
                          (-1 * ttest_PAW$ave_med / ttest_PAW$ave_exo))

# apply the basic two-sample t-test (we will pool variance)
t.result <- apply(log2(paw_tmt_tmm), 1, function(x) t.test(x[1:3], x[4:7], 
                                                           var.equal = TRUE))
                  
# extract the p-value column from the t-test thingy 
ttest_PAW$PValue <- unlist(lapply(t.result, function(x) x$p.value))
                                  
# do a Benjamini-Hochberg multiple testing correction
ttest_PAW$FDR <- p.adjust(ttest_PAW$PValue, method = "BH")

# add a DE candidate status column
ttest_PAW$candidate <- cut(ttest_PAW$FDR, breaks = c(-Inf, 0.01, 0.05, 0.10, 1.0), 
                           labels = c("high", "med", "low", "no"))
ttest_PAW$Acc <- accession
    
# count up, down and the rest (FDR less than 0.05)
all <- dim(ttest_PAW)[1]
up <- dim(ttest_PAW[(ttest_PAW$FDR <= 0.10) & (ttest_PAW$logFC > 0.0), ])[1]
down <- dim(ttest_PAW[(ttest_PAW$FDR <= 0.10) & (ttest_PAW$logFC <= 0.0), ])[1]
print("This is like decideTest in edgeR - 10% FDR cut:")
up 
all - up - down
down
print("Candidate Counts:")
summary(ttest_PAW$candidate)
                                  
# what does the test p-value distribution look like?
pvalue_plot(ttest_PAW, "PAW data with t-test")

# can look at log2FC distributions as a check
log2FC_plots(med_exo, 3, "LogFC by candidate for edgeR")
log2FC_plots(ttest_PAW, 3, "LogFC by candidate for t-test")

MA_plots(med_exo, "ave_med", "ave_exo", "PAW edgeR results")
MA_plots(ttest_PAW, "ave_med", "ave_exo", "PAW t-test results")

# the scatter plots
#scatter_plots(med_exo, "ave_med", "ave_exo", "PAW edgeR results")
#scatter_plots(ttest_PAW, "ave_med", "ave_exo", "PAW t-test results")

# compare volcano plots
volcano_plot(med_exo, "ave_med", "ave_exo", "PAW edgeR results", 40)
volcano_plot(med_exo, "ave_med", "ave_exo", "PAW edgeR results (expanded y-axis)", 4)
volcano_plot(ttest_PAW, "ave_med", "ave_exo", "PAW t-test results", 4)

set_plot_dimensions(7, 4)
plot_top_tags(ttest_PAW, 3, 4, 10)
set_plot_dimensions(7, 7)

for (cutoff in c(0.10, 0.05, 0.01)) {
    edgeR <- get_identifier(filter(med_exo, FDR < cutoff)$Acc)
    ttest <- get_identifier(filter(ttest_PAW, FDR < cutoff)$Acc)
    cat("cutoff:", cutoff, "\n")
    cat("number of candidates (edgeR, t-test):", length(edgeR), length(ttest), "\n")

    x_union <- union(edgeR, ttest)
    x_intersect <- intersect(edgeR, ttest)
    unique_edgeR <- setdiff(edgeR, ttest)
    unique_ttest <- setdiff(ttest, edgeR)

    cat("union and intersection:", length(x_union), length(x_intersect), "\n")
    cat("intersection out of union:", round(length(x_intersect)/length(x_union), 2), "\n")
    cat("intersection out of edgeR:", round(length(x_intersect)/length(edgeR), 2), "\n")
    cat("intersection out of t-test:", round(length(x_intersect)/length(ttest), 2), "\n")
    cat("unique to each (edgeR, ttest):", length(unique_edgeR), 
        length(unique_ttest), "\n\n")
}

edgeR_candidates <- get_identifier(filter(med_exo, FDR < 0.10)$Acc)
ttest_candidates <- get_identifier(filter(ttest_PAW, FDR < 0.10)$Acc)
unique_edgeR_candidates <- setdiff(edgeR_candidates, ttest_candidates)
unique_ttest_candidates <- setdiff(ttest_candidates, edgeR_candidates)

cat("Candidates in t-test but not candidate (<0.1) in edgeR")
interesting_ttest <- setdiff(unique_ttest_candidates, edgeR_candidates)
length(interesting_ttest)
cat("Candidates in edgeR but not candidate (<0.1) in t-test")
interesting_edgeR <- setdiff(unique_edgeR_candidates, ttest_candidates)
length(interesting_edgeR)


cat("High candidates in t-test but not candidate (<0.1) in edgeR")
interesting <- setdiff(unique_ttest, edgeR_candidates) # high t-test but not candidate in edgeR
not_interesting <- setdiff(unique_ttest, interesting) # high in t-test and also candidate in edgeR
length(interesting)

plot_selected <- function(results, nleft, nright, selected) {
    # results should have data first, then test results (two condition summary table)
    # nleft, nright are number of data points in each condition
    # selected is list of identifiers to plot
    proteins <- results
    proteins$ident <- get_identifier(results$Acc)
    proteins <- filter(proteins, ident %in% selected)
        
    color = c(rep("red", nleft), rep("blue", nright))
    for (row_num in 1:nrow(proteins)) {
        row <- proteins[row_num, ]
        vec <- as.vector(unlist(row[1:(nleft + nright)]))
        names(vec) <- colnames(row[1:(nleft + nright)])
        title <- str_c(row$ident, ", int: ", scientific(mean(vec), 2), 
                       ", p-val: ", scientific(row$FDR, digits = 3), 
                       ", FC: ", round(row$FC, digits = 1))
        barplot(vec, col = color, main = title)
    }    
}

# plot the t-test candidates that were not significant in edgeR
set_plot_dimensions(7, 4)
plot_selected(ttest_PAW, 3, 4, interesting)
set_plot_dimensions(7, 7)

volcano_plot_facet <- function(results, x, y, other, title, ymax) {
    # makes a volcano plot
        # results - a data frame with edgeR results
        # x - string for the x-axis column
        # y - string for y-axis column
        # ymax - upper limit for y-axis
        # title - plot title string
        # other - candidate flags from the other test
    
    # uses transformed data
    temp <- transform(results, x, y)
    temp$other <- other
    
    # build the plot
    ggplot(temp, aes(x = M, y = P)) +
        geom_point(aes(color = candidate, shape = candidate), alpha = 0.5) +
        xlab("log2 FC") +
        ylab("-log10 FDR") +
        coord_cartesian(xlim = c(-5, 5), ylim = c(0, ymax)) + 
        facet_wrap(~ other) +
        ggtitle(str_c(title, " Volcano Plot"))
}

# cross-candidate faceted volcano plots
volcano_plot_facet(med_exo, "ave_med", "ave_exo", str_c("t-test_", ttest_PAW$candidate), 
                   "edgeR faceted by t-test candidate", 10)

# with the tests switched
volcano_plot_facet(ttest_PAW, "ave_med", "ave_exo", str_c("edgeR_", med_exo$candidate), 
                   "t-test faceted by edgeR candidate", 4)

# save the testing results
#write.table(med_exo, file = "KUR1502_results.txt", sep = "\t",
#            row.names = FALSE, na = " ")

# log the session
sessionInfo()
