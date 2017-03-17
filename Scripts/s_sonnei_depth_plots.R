library(ggplot2)

a=read.csv(file="summary.calls.tsv", sep="\t", header=TRUE)

ariba_positions = which(a$ARIBA.read.depth > 0)
kmer_positions = which(a$kmerResistance.read.depth > 0)
srst2_positions = which(a$SRST2.read.depth > 0)


ggplot() +
  geom_freqpoly(aes(a$ARIBA.read.depth[ariba_positions], colour="ARIBA"), binwidth=5) +
  geom_freqpoly(aes(a$kmerResistance.read.depth[kmer_positions], colour="KmerResistance"), binwidth=5) +
  geom_freqpoly(aes(a$SRST2.read.depth[srst2_positions], colour="SRST2"), binwidth=5) +
  scale_color_manual("", values=c("ARIBA"="#7fc97f", "KmerResistance"="#beaed4", "SRST2"="#fdc086")) +
  xlab("Read depth") +
  ylab("Number of calls")

ggsave("summary.s_sonnei.called_depth.pdf", scale=0.4, width=12, height=5)


ggplot() +
  geom_freqpoly(aes(a$ARIBA.read.depth[ariba_positions], colour="ARIBA"), binwidth=5) +
  geom_freqpoly(aes(a$kmerResistance.read.depth[kmer_positions], colour="KmerResistance"), binwidth=5) +
  geom_freqpoly(aes(a$SRST2.read.depth[srst2_positions], colour="SRST2"), binwidth=5) +
  scale_color_manual("", values=c("ARIBA"="#7fc97f", "KmerResistance"="#beaed4", "SRST2"="#fdc086")) +
  xlim(c(0,150)) +
  xlab("Read depth") +
  ylab("Number of calls")

ggsave("summary.s_sonnei.called_depth.0-150.pdf", scale=0.4, width=12, height=5)

file.remove("Rplots.pdf")
