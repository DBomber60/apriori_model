library(arules)
library(ggplot2)
nItems = 20

############ SIMULATE TRANSACTIONS (using Agrawal method) ############
nTrans = 2000
# generate random itemsets and give them a probability distribution
patterns <- random.patterns(nItems = nItems, nPats = 20, corr = .1, 
                            lPats = 2, cmean = .5, cvar = 0.1)

data = list()
for (i in 1:nTrans) {
  j = rpois(1,1)+1 # how many itemsets does the transaction contain
  k = findInterval(runif(j), cumsum(patterns@quality$pWeights)) + 1 #which itemsets
  trans1 = c()
    for (q in 1:length(k)) {
      trans1 = c(trans1,patterns@items@itemInfo$labels[which(patterns@items@data[,k[q]]==TRUE)])
    }
  data[[i]] = unique(trans1)
}
data <- as(data, "transactions")


### Generate design matrix and response vector for full model
toBits <- function (x, nBits = nodes) {tail(rev(as.numeric(intToBits(x))),nBits)}
bitsToInt<-function(x) {packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")}
nodes = 20
X = t(sapply(seq(1:(2^nodes)), function(x) {toBits(x, nBits=nodes)}))[-2^nodes,] # not sure we need this, ever
Y = array(0, dim=c(2^nodes-1,1)) # initialize response vector
for (i in 1:data@data@Dim[2]) { # and add observations
  Y[bitsToInt(as.numeric(data@data[,i]))] = Y[bitsToInt(as.numeric(data@data[,i]))] + 1
}
str(data)

##### FULL GLM ####
m0 = glm(Y~X, family=poisson)

#### REDUCED GLM (cut out large edges) ####
gamma = apply(X,1,sum)
X = X[-which(gamma>6),]
Y = Y[-which(gamma>6)]
m1 = glm(Y~X, family=poisson)


##### Reduced GLM (Using APriori Algorithm) ####

# generate frequent subsets using apriori algorithm
rules <- apriori(data, parameter = list(supp = 0.01, conf = 0.90, minlen=2, maxlen=4))

# addtrans takes a transactions data structure and a vector, adds vector as new transaction
addtrans = function(dat, vec) {
  dat@data@p[as.integer(length(dat@data@p)+1)] = as.integer(dat@data@p[length(dat@data@p)] + length(vec))
  dat@data@i = c(dat@data@i, as.integer(vec))
  dat@data@Dim[2] = as.integer(dat@data@Dim[2] + 1)
  return(dat@data)
}

# which items to break up - let's say those with edge size > 6
bigtrans = which(size(data)>6)

# loops through rule list - adds frequently occurring subsets as observations
for (k in 1: length(bigtrans)) {
  for (i in 1:100) {
    ruleitems = c(rules@lhs@itemInfo$labels[which(rules@lhs@data[,i]==TRUE)],rules@rhs@itemInfo$labels[which(rules@rhs@data[,i]==TRUE)])
    vec = which(data@itemInfo$labels %in% ruleitems)-1 # vector of item indices
    l = ruleitems %in% data@itemInfo$labels[which(data@data[,bigtrans[k]]==TRUE)]
    if (sum(l)/length(l)==1) {
      data@data = addtrans(data, vec)
      }
  }
}

Y = array(0, dim=c(2^nodes-1,1))

data@data@Dim[2]
for (i in 1:data@data@Dim[2]) {
  Y[bitsToInt(as.numeric(data@data[,i]))] = Y[bitsToInt(as.numeric(data@data[,i]))] + 1
}

Y = Y[-which(gamma>6)]
m2 = glm(Y ~ X, family = poisson())

####### VISUALIZATION #######
model1= m0
model2 = m1
model3 = m2

model1Frame <- data.frame(Variable = rownames(summary(model1)$coef),
                          Coefficient = summary(model1)$coef[, 1],
                          SE = summary(model1)$coef[, 2],
                          modelName = "Full Model")
model2Frame <- data.frame(Variable = rownames(summary(model2)$coef),
                          Coefficient = summary(model2)$coef[, 1],
                          SE = summary(model2)$coef[, 2],
                          modelName = "Reduced Model")
model3Frame <- data.frame(Variable = rownames(summary(model3)$coef),
                          Coefficient = summary(model3)$coef[, 1],
                          SE = summary(model3)$coef[, 2],
                          modelName = "Apriori Reduced Model")
allModelFrame <- data.frame(rbind(model1Frame, model2Frame, model3Frame))
allModelFrame = allModelFrame[-which(allModelFrame$SE>100),]
interval1 = 1.96
zp1 <- ggplot(allModelFrame, aes(colour = modelName))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
zp1 <- zp1 + geom_linerange(aes(x = Variable, ymin = Coefficient - SE*interval1,
                                ymax = Coefficient + SE*interval1),
                            lwd = 1, position = position_dodge(width = 1/2))
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + ggtitle("Comparing Model Estimates")
print(zp1)  # The trick to these is position_dodge().

