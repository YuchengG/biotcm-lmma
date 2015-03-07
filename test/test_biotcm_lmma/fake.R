# Load
tab <- read.table('../../edge.tab', header = TRUE, sep = '\t')
genes <- union(tab[,1], tab[,2])
mset <- data.frame(row.names = genes)

# Fake
means <- rnorm(genes, mean = 5, sd = 1)
for(i in 1:10) { # index of current sample
  mset[genes, paste('Sample', i)] <- rnorm(genes, mean = means, sd = 1)
}

# Export
write(paste(c('Symbol', colnames(mset)), collapse = '\t'), 'fake.tab')
write.table(mset, 'fake.tab', quote = FALSE, sep = '\t', col.names = FALSE, append = TRUE)
