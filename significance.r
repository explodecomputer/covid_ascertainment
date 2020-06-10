library(data.table)
library(dplyr)
library(parallel)
library(ieugwasr)
library(ggplot2)

# read in 
covid <- fread("covid19_result.txt", he=T)
linker <- fread("linker_app15825.csv")
table(duplicated(covid$eid))
table(covid$eid %in% linker$app)
linker$covid <- as.numeric(linker$app %in% covid$eid)
co <- subset(linker, select=c(ieu, covid))

covid2 <- fread("covid19_result_2020_06_05.txt", he=T)
linker$covid2 <- as.numeric(linker$app %in% covid2$eid)
co2 <- subset(linker, select=c(ieu, covid2))
co2$random <- co2$ieu %in% sample(co2$ieu, sum(co2$covid2) * 3, replace=FALSE) %>% as.numeric()
table(co2$random)
table(co2$covid2)

phenpath <- "/mnt/storage/private/mrcieu/research/mr-eve/UKBB_replication2/replication/results"
load("covid_ascertainment.rdata")
res$fdr <- p.adjust(res$`Pr`, "fdr")
res <- subset(res, grepl("ukb-b", exposure) & fdr < 0.05)
phens <- res$exposure
phens <- list.files(phenpath)



covidieu <- subset(co2, covid2==1)$ieu
randomieu <- subset(co2, random==1)$ieu

phenout <- mclapply(phens, function(x)
{
	message(x)
	fn <- file.path(phenpath, x, "phen.txt")
	file.exists(fn)
	a <- fread(fn)
	a$discovery[is.na(a$discovery)] <- a$replication[is.na(a$discovery)]
	names(a)[names(a) == "discovery"] <- x

	a1 <- a %>% subset(., FID %in% covidieu, select=c("FID", x))
	a2 <- a %>% subset(., FID %in% randomieu, select=c("FID", x))
	a1$what <- "random"
	a2$what <- "covid"
	rbind(a1, a2)
}, mc.cores=16)



table(phenout[[1]]$FID == phenout[[2]]$FID)

rand <- mclapply(phenout, function(x){
	subset(x, what=="random")[[2]]
}) %>% do.call("cbind", .)
covi <- mclapply(phenout, function(x){
	subset(x, what=="covid")[[2]]
}) %>% do.call("cbind", .)
all(sapply(phenout, function(x) names(x)[2]) == phens)

colnames(covi) <- colnames(rand) <- phens

covicor <- cor(covi, use="pair")
covicor[1:10,1:10]
covicor[lower.tri(covicor)] <- 10
covicorl <- reshape2::melt(covicor) %>% subset(., Var1 != Var2 & value != 10) 

randcor <- cor(rand, use="pair")
randcor[1:10,1:10]
randcor[lower.tri(randcor)] <- 10
randcorl <- reshape2::melt(randcor) %>% subset(., Var1 != Var2 & value != 10) 

dim(randcorl)
dim(covicorl)

dat <- merge(randcorl, covicorl, by=c("Var1", "Var2"))
table(sign(dat$value.x) == sign(dat$value.y))

dat$dif <- abs(dat$value.x - dat$value.y)

iscont <- colnames(covi)[apply(covi, 2, function(x) sum(is.na(x)) > 100 & length(unique(x)) > (sum(is.na(x))/2))]
table(iscont)

dat$cont <- dat$Var1 %in% iscont & dat$Var2 %in% iscont
dat <- arrange(dat, desc(dif))
head(dat)

ao <- gwasinfo(phens) %>% dplyr::select(id, trait)
index <- match(dat$Var1, ao$id)
dat$t1 <- ao$trait[index]

index <- match(dat$Var2, ao$id)
dat$t2 <- ao$trait[index]


subset(dat, cont) %>% arrange(desc(dif)) %>% write.csv(., file="dif.csv")

save(dat, file="significance.rdata")


fn <- function(id1, id2)
{
	example <- rbind(
		tibble(id1 = rand[,id1], id2 = rand[,id2], what="Random"),
		tibble(id1 = covi[,id1], id2 = covi[,id2], what="COVID test")
	)
	ggplot(example, aes(x=id1, y=id2)) +
	geom_point(aes(colour=what), size=0.1) +
	geom_smooth(method="lm", aes(colour=what))
}



p <- fn("ukb-b-13686","ukb-b-14540") + labs(x="Age", y="Adiposity", colour="Sample\nselection") + ylim(c(-2.5, 2.5)) + xlim(c(-2.5, 2.5))
ggsave(p, file="significance_plot.pdf")



