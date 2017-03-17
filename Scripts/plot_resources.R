library(ggplot2)
d = read.csv(file="resources.tsv", sep="\t", header=T)
xbreaks = c("E.faecium.mlst", "E.faecium", "S.sonnei")
xlabels = c(expression(paste(italic("E. faecium"), " MLST")), expression(italic("E. faecium")), expression(italic("S. sonnei")))

ggplot(d, aes(factor(Dataset), Wall_clock/60, color=Tool)) +
  geom_boxplot() +
  xlab("Data set") +
  ylab("Wall clock time (m)") +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  scale_x_discrete(breaks=xbreaks, labels=xlabels)
ggsave(filename="resources.time.pdf", width=6, height=3, scale=0.85)


ggplot(d, aes(factor(Dataset), RAM, color=Tool)) +
  geom_boxplot() +
  xlab("Data set") +
  ylab("Peak RAM (GB)") +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent") +
  scale_x_discrete(breaks=xbreaks, labels=xlabels)
ggsave(filename="resources.ram.pdf", width=6, height=3, scale=0.85)

file.remove("Rplots.pdf")
