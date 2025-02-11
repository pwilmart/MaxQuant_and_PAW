
# load libraries
library("tidyverse")
library("psych")
library("gridExtra")
library("scales")
library("limma") 
library("edgeR") 

# read the grouped protein summary file
MQ_raw <- read_tsv("KUR1502_data-export.txt")

# separate accessions from the data
accession <- MQ_raw$Accession

# get MS1 quantities
intensities <- MQ_raw$Intensity
ibaqs <- MQ_raw$iBAQ

# get the reporter ion intensities
MQ_tmt <- MQ_raw %>% select(-Accession, -Intensity, -iBAQ)

head(MQ_tmt)
nrow(MQ_tmt)

counts <- 1:nrow(MQ_tmt) # create vector of appropriate length
for(i in 1:nrow(MQ_tmt)){
    # TRUE will be coerced to 1 for the summing
    counts[i] <- sum(MQ_tmt[i, ] == 20)
}
table(counts) # create the summary

# let's see what the starting data look like
color = c(rep("red", 3), rep("blue", 4))
boxplot(log10(MQ_tmt), col = color, notch = TRUE, main = "Starting MaxQuant data")

SL_Norm <- function(df, color = NULL, plot = TRUE) {
    # This makes each channel sum to the average grand total
        # df - data frame of TMT intensities
        # returns a new data frame with normalized values
    
    # compute scaling factors to make colsums match the average sum
    norm_facs <- mean(c(colSums(df))) / colSums(df)
    cat("SL Factors:\n", sprintf("%-5s -> %f\n", colnames(df), norm_facs))
    df_sl  <- sweep(df, 2, norm_facs, FUN = "*")
    
    # visualize results and return data frame
    if(plot == TRUE) {
        boxplot(log10(df_sl), col = color, notch = TRUE, main = "SL Normalized data")
    }
    df_sl
}

# SL norm the tmt data
MQ_tmt_sl <- SL_Norm(MQ_tmt, color)

# compute column sums before and after normalization
print("Before:")
format(round(colSums(MQ_tmt), digits = 0), big.mark = ",")
print("After:")
format(round(colSums(MQ_tmt_sl), digits = 0), big.mark = ",")

# load data into DGEList object
group <- c(rep("media", 3), rep("exosome", 4))
y <- DGEList(counts = MQ_tmt, group = group, genes = accession)

# we can see what y looks like
#y
y$samples # this is shorter to view

# run the TMM normalization
y <- calcNormFactors(y)

# check what changed
y$samples

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
MQ_tmt_tmm <- apply_tmm_factors(y, color)

# check returned table
head(MQ_tmt_tmm)

# check column totals
print("After TMM:")
format(round(colSums(MQ_tmt_tmm), digits = 0), big.mark = ",")

# define variables for the columns in each condition
M <- 1:3
E <- 4:7

CV <- function(df) {
    # Computes CVs of data frame rows
        # df - data frame, 
        # returns vector of CVs (%)
    ave <- rowMeans(df)    # compute averages
    sd <- apply(df, 1, sd) # compute standard deviations
    cv <- 100 * sd / ave   # compute CVs in percent (last thing gets returned)
}

labeled_boxplot <- function(df, ylim, title) {
    # Makes a box plot with the median value labeled
        # df - data frame with data to compute CVs of
        # ylim - upper limit for y-axis
        # title - plot title
    cv = CV(df)
    boxplot(cv, ylim = c(0, ylim), notch = TRUE, main = title)
    text(x = 0.65, y = boxplot.stats(cv)$stats[3], 
         labels = round(boxplot.stats(cv)$stats[3], 1))
}

# compare CV distributions
par(mfrow = c(2, 2))
labeled_boxplot(MQ_tmt[M], ylim = 150, title = "Starting Media CVs")
labeled_boxplot(MQ_tmt[E], ylim = 150, title = "Starting Exosome CVs")
labeled_boxplot(MQ_tmt_tmm[M], ylim = 75, title = "Media CVs after TMM")
labeled_boxplot(MQ_tmt_tmm[E], ylim = 75, title = "Exosome CVs after TMM")
par(mfrow = c(1, 1))

# see how things cluster after we have normalized data
plotMDS(log2(MQ_tmt_tmm), col = c(rep("red", 3), rep("blue", 4)), main = "MaxQuant TMM")

# scatter plots within groups and betwen groups
pairs.panels(log10(MQ_tmt_tmm), lm = TRUE, main = "After TMM")

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

# make the results table 
tt <- topTags(et, n = Inf, sort.by = "none")$table
med_exo <- collect_results(MQ_tmt_tmm, tt, M, "med", E, "exo")

pvalue_plots <- function(results, ylim, title) {
    # Makes p-value distribution plots
        # results - results data frame
        # ylim - ymax for expanded view
        # title - plot title
    p_plot <- ggplot(results, aes(PValue)) + 
        geom_histogram(bins = 100, fill = "white", color = "black") +
        geom_hline(yintercept = mean(hist(results$PValue, breaks = 100, 
                                     plot = FALSE)$counts[26:100]))

    # we will need an expanded plot
    p1 <- p_plot + ggtitle(str_c(title, " p-value distribution"))
    p2 <- p_plot + coord_cartesian(xlim = c(0, 1.0), ylim = c(0, ylim)) + 
        ggtitle("p-values expanded")
    grid.arrange(p1, p2, nrow = 2) # from gridExtra package
}

# check the p-value distrubution
pvalue_plots(med_exo, 100, "Media vs Exosome-dosed")

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
    
    df # return the data frame
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

scatter_plots(med_exo, "ave_med", "ave_exo", "Media versus exosome-dosed")

volcano_plot <- function(results, x, y, title, ymax) {
    # makes a volcano plot
        # results - a data frame with edgeR results
        # x - string for the x-axis column
        # y - string for y-axis column
        # title - plot title string
        # ymax - upper limit for y-axis
    
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
volcano_plot(med_exo, "ave_med", "ave_exo", "Media versus exosome-dosed", 50)

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
    # get top ipregulated
    up <- results %>% 
        filter(logFC >= 0) %>%
        arrange(FDR)
    up <- up[1:top_tags, ]
    
    # get top down regulated
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
# plot the top 20 up and 20 down proteins
set_plot_dimensions(7, 4)
plot_top_tags(med_exo, 3, 4, 20)
set_plot_dimensions(7, 7)

# column totals from the dilution series (25, 20, 15, 10, 5, 2.5) data
MQ <- c(252042061, 194468567, 147076787, 101203918, 49272033, 26071287)
PAW <- c(1996324081, 1562525113, 1207642244, 835353557, 400841992, 221591661)
dilution <- data.frame(PAW = PAW, MQ = MQ)

# from https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
ggplotRegression <- function (fit, title) {
    require(ggplot2)
    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
        geom_point(size = 3) +
        stat_smooth(method = "lm", col = "red") +
        labs(title = title, 
             subtitle = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5), 
                              "Intercept =",signif(fit$coef[[1]],5 ),
                              " Slope =",signif(fit$coef[[2]], 5),
                              " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(lm(MQ ~ PAW, data = dilution), "MQ reporter ions versus peak heights")

# get summ of reporter ions for each protein
totals  <- rowSums(MQ_tmt_tmm)

# add to data frame with MS1 abundances (drop zeros and take logs)
abundances <- data.frame(Total_TMT = totals, MS1_Intensity = intensities, iBAQ_Intensity = ibaqs)
abundances <- log10(abundances[abundances$MS1_Intensity > 0, ])

# make some scatter plots
pairs.panels(abundances, lm = TRUE, 
             main = "Protein relative abundance measures")
ggplotRegression(lm(MS1_Intensity ~ Total_TMT, data = abundances), 
                 "Correlation between MS1 and reporter ion intensities")

# save the testing results
write.table(med_exo, file = "KUR1502_MQ-results.txt", sep = "\t",
            row.names = FALSE, na = " ")

# log the session
sessionInfo()

# do the t-test on log transformed intensities to be safe
ttest_MQ <- log2(MQ_tmt_tmm)
# add average ratio columns (non-logged ratios), fold-change column, and row names
ttest_MQ$ave_med <- rowMeans(MQ_tmt_tmm[1:3])
ttest_MQ$ave_exo  <- rowMeans(MQ_tmt_tmm[4:7])
ttest_MQ$logFC <- log2(ttest_MQ$ave_exo / ttest_MQ$ave_med)
row.names(ttest_MQ) <- accession

# apply the basic two-sample t-test (we will pool variance)
t.result <- apply(ttest_MQ, 1, function(x) t.test(x[1:3], x[4:7], var.equal = TRUE))
# extract the p-value column from the t-test thingy 
ttest_MQ$PValue <- unlist(lapply(t.result, function(x) x$p.value))
# do a Benjamini-Hochberg multiple testing correction
ttest_MQ$FDR <- p.adjust(ttest_MQ$PValue, method = "BH")

# add a DE candidate status column
ttest_MQ$candidate <- cut(ttest_MQ$FDR, breaks = c(-Inf, 0.01, 0.05, 0.10, 1.0), 
                          labels = c("high", "med", "low", "no"))
    
# count up, down and the rest (FDR less than 0.05)
all <- dim(ttest_MQ)[1]
up <- dim(ttest_MQ[(ttest_MQ$FDR <= 0.10) & (ttest_MQ$logFC > 0.0), ])[1]
down <- dim(ttest_MQ[(ttest_MQ$FDR <= 0.10) & (ttest_MQ$logFC <= 0.0), ])[1]
print("This is like decideTest in edgeR - 10% FDR cut:")
up 
all - up - down
down
print("Candidate Counts:")
summary(ttest_MQ$candidate)
    
# what does the test p-value distribution look like?
ggplot(ttest_MQ, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(ttest_MQ$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("MQ data with t-test p-value distribution")

MA_plots(med_exo, "ave_med", "ave_exo", "MQ edgeR results", FALSE)
MA_plots(ttest_MQ, "ave_med", "ave_exo", "MQ t-test results")

scatter_plots(med_exo, "ave_med", "ave_exo", "MQ edgeR results", FALSE)
scatter_plots(ttest_MQ, "ave_med", "ave_exo", "MQ t-test results")

# compare volcano plots
volcano_plot(med_exo, "ave_med", "ave_exo", "EdgeR results", 4)
volcano_plot(ttest_MQ, "ave_med", "ave_exo", "t-test results", 4)


