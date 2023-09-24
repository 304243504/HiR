Args <- commandArgs()
if (length(Args) != 7)
{
	stop("the commans is like:\n\t\t\t\tRscript hypergeometric.r inputfile outputfile\n")
	
}

if (length(grep(".gz$",Args[6])) || length(grep(".GZ$",Args[6])) || length(grep(".gzip$",Args[6])) || length(grep(".GZIP$",Args[6]))) {
  gf = gzfile(Args[6],'rt')
} else {
  gf = Args[6]
}

data <- read.table(gf,header=F)
data[,16] <- phyper(data[,8]-1, data[,11], data[,15]*2-data[,11], data[,13], lower.tail = FALSE)
data[[16]] <- as.numeric(format(data[[16]], digits = 3))
data[,17] = data[,8]/sqrt(data[,11])/sqrt(data[,13])
write.table(data, file=Args[7], sep = "\t", row.names = FALSE, col.names = FALSE,quote=F)
