library(CompQuadForm)
library(SKAT)
library(MASS)
library(Matrix)
library(mvtnorm)

source("sim.R")
source("gen.R")
source("mtd.R")
source("imp.R")
source("err.R")
source("cor.R")

for(. in dir("R", "[.]R$", ful=TRUE))
{
    source(.)
}
