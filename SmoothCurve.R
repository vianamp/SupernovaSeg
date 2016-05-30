filename <- '/Users/viana/Desktop/SupernovaSeg/build/curve.txt'
A <- data.frame(read.table(filename,header = F))
names(A) <- c('r','d')
#A <- within(A,ds<-smooth.spline(A$r,A$d,spar=0.1)$y)
A <- within(A,ds<-runmed(A$d, 15))
write.table(x = A[,c(1,3)],filename,sep = ' ',col.names = F,row.names = F)
