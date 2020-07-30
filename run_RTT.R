library(ape)
library(hash)

args = commandArgs(TRUE)

tree = args[1] # first input
tDates = args[2] # second input
outputTree = args[3] # output

t = read.tree(tree)
d <- read.table(tDates)

Dict = hash()
# build dictionary
for (j in 1:length(d$V2)) {
  Dict[d$V1[j]] = d$V2[j]
}

# match t$tip.label to d
d_ordered = rep(0,length(t$tip.label))
for (j in 1:length(t$tip.label)) {
  d_ordered[j] = Dict[[t$tip.label[j]]]
}
  
start_time <- Sys.time()
r = rtt(t,d_ordered,objective="rms")
end_time <- Sys.time()
time = end_time - start_time

write.tree(r,outputTree)
print(time)
