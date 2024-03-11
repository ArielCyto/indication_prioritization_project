#signature from IO
IO_signatures <- readRDS(get_workflow_outputs("wf-d2ad9d9a50"))
sign <- IO_signatures@.Data
names(sign)<- IO_signatures@names
IO_signatures <- sign

genes_to_check <- sign$`mature natural killer cell`
genes_to_check <- genes_to_check[genes_to_check %in% colnames(expression_matrix)]
dfi <- rank_groups_by_genes(genes_to_check)

install.packages('bipartite')
library(bipartite)
library(igraph)

rankchange <- function(list.1, list.2){
  grp = c(rep(0,length(list.1)),rep(1,length(list.2)))
  m = match(list.1, list.2)
  m = m + length(list.1)
  pairs = cbind(1:length(list.1), m)
  pairs = pairs[!is.na(pairs[,1]),]
  pairs = pairs[!is.na(pairs[,2]),]
  g = graph.bipartite(grp, as.vector(t(pairs)), directed=TRUE)
  V(g)$color =  c("red","green")[grp+1]
  V(g)$label = c(list.1, list.2)
  V(g)$x = grp
  V(g)$y = c(length(list.1):1, length(list.2):1)
  V(g)$size = 4
  g
}

# manual excel file
list.1 = c('lung cancer - LUSC_inflamed',
           'breast cancer_inflamed',
           'lung cancer - LUAD_inflamed',
           'head and neck cancer_inflamed',
           'lung cancer - LUSC_inflamed_mild_stromal',
           'prostate cancer_inflamed',
           'bladder cancer_inflamed',
           'colorectal cancer_inflamed_mild_stromal',
           'head and neck cancer_inflamed_mild_stromal',
           'lung cancer - LUAD_inflamed_stromal',
           'breast cancer_inflamed_mild_stromal',
           'lung cancer - LUSC_inflamed_stromal',
           'colorectal cancer_inflamed',
           'lung cancer - LUAD_inflamed_mild_stromal',
           'bladder cancer_inflamed_mild_stromal',
           'colorectal cancer_inflamed_stromal',
           'lung cancer - LUAD_stromal_epithelial',
           'bladder cancer_inflamed_stromal',
           'head and neck cancer_inflamed_stromal',
           'breast cancer_inflamed_stromal',
           'lung cancer - LUSC_stromal_epithelial',
           'lung cancer - LUAD_epithelial',
           'prostate cancer_stromal_epithelial',
           'colorectal cancer_stromal_epithelial',
           'head and neck cancer_stromal_epithelial',
           'lung cancer - LUSC_epithelial',
           'breast cancer_stromal_epithelial',
           'colorectal cancer_epithelial',
           'bladder cancer_stromal_epithelial',
           'head and neck cancer_epithelial',
           'bladder cancer_epithelial'
)
list.2 = c('head and neck cancer_inflamed_mild_stromal',
           'head and neck cancer_inflamed',
           'bladder cancer_inflamed_mild_stromal',
           'lung cancer - LUSC_inflamed',
           'colorectal cancer_inflamed_mild_stromal',
           'lung cancer - LUAD_inflamed',
           'breast cancer_inflamed',
           'bladder cancer_inflamed',
           'colorectal cancer_inflamed',
           'lung cancer - LUSC_inflamed_mild_stromal',
           'breast cancer_inflamed_mild_stromal',
           'head and neck cancer_inflamed_stromal',
           'lung cancer - LUSC_inflamed_stromal',
           'head and neck cancer_stromal_epithelial',
           'lung cancer - LUAD_inflamed_stromal',
           'prostate cancer_inflamed',
           'colorectal cancer_inflamed_stromal',
           'lung cancer - LUAD_inflamed_mild_stromal',
           'lung cancer - LUAD_stromal_epithelial',
           'bladder cancer_inflamed_stromal',
           'lung cancer - LUSC_stromal_epithelial',
           'head and neck cancer_epithelial',
           'colorectal cancer_epithelial',
           'lung cancer - LUAD_epithelial',
           'breast cancer_inflamed_stromal',
           'colorectal cancer_stromal_epithelial',
           'lung cancer - LUSC_epithelial',
           'prostate cancer_stromal_epithelial',
           'bladder cancer_stromal_epithelial',
           'breast cancer_stromal_epithelial',
           'bladder cancer_epithelial')

g = rankchange(list.1, list.2)
plot(g)
