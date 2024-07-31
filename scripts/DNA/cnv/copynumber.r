library("optparse")

option_list = list(
    make_option(c('-s','--sample'),default=NULL,help='sample name'),
    make_option(c('-o','--output'),default=NULL,help='output directory'),
    make_option(c('-i','--input'),default=NULL,help='input directory')
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

peaks <- function(series, span=3, ties.method = "first") {
###  This peaks function from Petr Pikal https://stat.ethz.ch/pipermail/r-help/2007-February/125241.html
	if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
	z <- embed(series, span)
	s <- span%/%2
	v <- max.col(z, ties.method=ties.method) == 1 + s
	pad <- rep(FALSE, s)
	result <- c(pad, v, pad)
	result
}
data <- paste0(opt$input,opt$sample,'.hg38.50k.k50.nobad.varbin.data.txt')
df <- read.table(file=data, header=T)
short <- paste0(opt$input,opt$name,'.hg38.50k.k50.nobad.varbin.short.txt')
dfs <- read.table(short, header=T)

starts <- c()
ends <- c()
prevEnd <- 0
len <- nrow(dfs)
for (j in 1:len) {
	thisStart = prevEnd + 1
	thisEnd = thisStart + dfs$num.mark[j] - 1
	starts <- c(starts, thisStart)
	ends <- c(ends, thisEnd)
	prevEnd = thisEnd
}

amat <- matrix(data=0, nrow=1500000, ncol=1)
counter <- 1
for (j in 1:(len-1)) {
	for (k in (j+1):len) {
		N <- round((starts[j] - ends[j] + 1) * (starts[k] - ends[k] + 1)/1000)
		D <- abs(2^dfs$seg.mean[j] - 2^dfs$seg.mean[k])
		cat(N, "\t")
		if (N > 0) {
			amat[(counter:(counter+N-1)), 1] <- rep.int(D, N)
			counter <- counter+N
		}
	}
}
a3 <- amat[(1:counter),1]
a3.95 <- sort(a3)[round(.95*counter)]
a3d <- density(a3[which(a3 < a3.95)], n=1000)
cn0 <- a3d$x[which(peaks(as.vector(a3d$y), span=59))][1]
cn1 <- a3d$x[which(peaks(as.vector(a3d$y), span=59))][2]

df$cn.ratio <- df$lowratio / cn1
df$cn.seg <- df$seg.mean.LOWESS / cn1
df$copy.number <- round(df$cn.seg)

out = paste0(opt$output,opt$name,'.hg38.50k.k50.varbin.data.copynumber.txt')
write.table(df, sep="\t", file=out, quote=F, row.names=F)

postscript(paste0(opt$name,'.wg.cn.density.ps'), height=400, width=600)
par(mar=c(5.1,4.1,4.1,4.1))
plot(a3d, main=paste0(opt$name,' seg.mean difference density'))
dev.off()

for (a in 1:24) {
	postscript(paste0(opt$name,"wg.cn.chr", a, ".ps", sep=""), height=400, width=600)
	par(mar=c(5.1,4.1,4.1,4.1))
	plot(df$cn.ratio[df$chrom==a], main=paste0(opt$name," chr", a, sep=""), xlab="Bin", ylab="Ratio", col="#CCCCCC")
	lines(df$cn.ratio[df$chrom==a], col="#CCCCCC")
	points(df$cn.seg[df$chrom==a], col="#0000DD")
	lines(df$cn.seg[df$chrom==a], col="#0000DD")
	points(df$copy.number[df$chrom==a], col="#DD0000")
	lines(df$copy.number[df$chrom==a], col="#DD0000")
	dev.off()
}

