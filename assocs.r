library(data.table)
library(dplyr)
library(parallel)

# read in 
covid <- fread("~/covid19_result.txt", he=T)
linker <- fread("~/linker_app15825.csv")
table(duplicated(covid$eid))
table(covid$eid %in% linker$app)
linker$covid <- as.numeric(linker$app %in% covid$eid)
co <- subset(linker, select=c(ieu, covid))


# Analyse sex and chip
results <- list()
covs <- fread("/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/covariates/data.covariates_ieu.bolt.txt")
temp <- merge(co, covs, by.x="ieu", by.y="FID")
results$sex <- summary(glm(covid ~ sex, temp, family="binomial"))$coefficients
results$chip <- summary(glm(covid ~ chip, temp, family="binomial"))$coefficients


# Analyse PCs
pcs <- fread("/mnt/storage/private/mrcieu/research/mr-eve/UKBB_replication/replication/data/pcs.txt")
nom <- paste0("V", 3:42)

temp <- merge(co, pcs, by.x="ieu", by.y="FID")
out <- mclapply(nom, function(x)
{
	message(x)
	form <- paste0("covid ~ ", x) %>% as.formula
	return(summary(glm(form, temp, family="binomial"))$coefficients)
}, mc.cores=5)

phenpath <- "/mnt/storage/private/mrcieu/research/mr-eve/UKBB_replication2/replication/results"
phens <- list.files(phenpath)

phenout <- mclapply(phens, function(x)
{
	message(x)
	fn <- file.path(phenpath, x, "phen.txt")
	file.exists(fn)
	a <- fread(fn)
	a$discovery[is.na(a$discovery)] <- a$replication[is.na(a$discovery)]
	a <- merge(co, a, by.x="ieu", by.y="FID")
	form <- paste0("covid ~ discovery")
	return(try(summary(glm(form, a, family="binomial"))$coefficients))
}, mc.cores=6)

names(out) <- paste0("pc", 1:40)
names(phenout) <- phens

organise <- function(l)
{
	nom <- names(l)
	lapply(nom, function(x)
	{
		l[[x]] %>% as.data.frame() %>% mutate(exposure = x) %>% {.[2,]}
	}) %>% bind_rows()
}

res <- bind_rows(
	organise(results),
	organise(out),
	organise(phenout)
)

library(ieugwasr)
ao <- gwasinfo(phens) %>% dplyr::select(id, trait)
res$trait <- res$exposure
index <- match(ao$id, res$exposure)
length(index)
all(ao$id == res$trait[index])
res$trait[index] <- ao$trait

res %>% arrange(`Pr(>|z|)`) %>% head(30)


save(res, file="covid_ascertainment.rdata")

