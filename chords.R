library(chorddiag)

all.dat <- read.csv("data/gene_gene_sl_relationships.csv", comment="#")
pair.data <- all.dat[all.dat$SL==1,c(1,2)]

# converting everything to characters (as opposed to factors) because R insists that
# vectors of factors have all the same levels in order for == to work properly below
u.gene1 <- as.character(unique(pair.data$gene1))
u.gene2 <- as.character(unique(pair.data$gene2))
# remove u.gene1 things from u.gene2
u.gene2 <- as.character(setdiff(u.gene2, u.gene1))

# final matrix must be square for chorddiag,
# but we only care about the first 5 rows (each of u.gene1)
xy.matrix.axis = as.character(as.factor(union(u.gene1, u.gene2)))

# iterate through u.gene1, and see if there is a SL pair with each item in xy.matrix.axis
#
data.by.gene1 <- lapply(u.gene1,
    function(thisgene1){
        message("processing ", thisgene1)
        col <- unlist(
            lapply(xy.matrix.axis,
              function(thisgene2){
                length(pair.data$gene2[pair.data$gene1 == thisgene1 & pair.data$gene2==thisgene2])
                })
        )
        return(col)
    })
tdata <- as.matrix(do.call(rbind, data.by.gene1))

# add zero'd out rows for all things in u.gene2
newRows <- lapply(1:length(u.gene2), function(x){rep(0, dim(tdata)[2])})

tdata <- rbind(tdata, do.call(rbind, newRows))
dimnames(tdata) <- list(
    gene1 = xy.matrix.axis,
    gene2 = xy.matrix.axis
)


