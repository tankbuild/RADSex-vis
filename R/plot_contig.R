library(ggplot2)
library(cowplot)
library(readr)


transformM = function(x) {
    x/10^6
}


plot_contig <- function(data, contig, signif.threshold = 0.05, color.sex.bias = TRUE, sex.bias.palette = c("firebrick1", "black", "dodgerblue2"),
                        contig_names = NULL) {


    if (contig %in% names(data$lengths)) {
        contig_name = contig
    } else if (!is.null(contig_names) & contig %in% contig_names) {
        contig_name = names(contig_names)[which(contig_names == contig)]
    } else {
        stop("Contig not found")
    }

    contig_length = data$lengths[contig_name]
    contig_data = subset(data$data, data$data$Contig == contig_name)

    threshold = -log(signif.threshold / dim(data$data)[1], 10)

    g = ggplot(contig_data, aes(x = Position, y = SexBias, color = SexBias)) +
        geom_point(size=1.5) +
        scale_y_continuous(name = "Sex Bias", limits=c(-1, 1)) +
        theme(axis.title.x = element_blank(), legend.position = "none") +
        scale_x_continuous(name = paste("Position on ", contig, " (Mb)", sep=""), limits = c(0, contig_length), labels = transformM, minor_breaks = waiver())

    if (color.sex.bias) {
        g = g + scale_color_gradientn(name = "Sex Bias", colours = colorRampPalette(sex.bias.palette)(20))
    }

    h = ggplot(contig_data, aes(x = Position, y = P)) +
        geom_point(size=1.5, color = "grey20") +
        scale_y_continuous(name = expression(paste('-log'['10'], '(p'['Association with sex'], ')', sep='')), limits=c(0, max(contig_data$P, 1.25*threshold))) +
        geom_hline(yintercept = threshold, color="red3", lty=1, lwd=1) +
        annotate("text", x=0.05*max(contig_data$Position), y=threshold + 0.075*max(contig_data$P), label="p = 0.05", color="red3", size=5) +
        scale_x_continuous(name = paste("Position on ", contig, " (Mb)", sep=""), limits = c(0, contig_length), labels = transformM, minor_breaks = waiver())

    combined = plot_grid(g, h, ncol=1)

    return(combined)
}
