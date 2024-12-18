# This file was downloaded from here:
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5061
dat <- read.table("../data/E-MTAB-5061.sdrf.txt.gz",sep = "\t",header = TRUE,
                  check.names = FALSE,stringsAsFactors = FALSE)
cols <- c("Source Name",
          "Characteristics [sex]",
          "Characteristics [disease]",
          "Characteristics [single cell quality]",
          "Characteristics [submitted single cell quality]",
          "Characteristics [age]",
          "Characteristics [body mass index]",
          "Characteristics [clinical information]",
          "Characteristics [inferred cell type]")
dat <- dat[cols]
names(dat) <- c("id","sex","disease_status","quality","submitted_quality",
                "age","BMI","clinical_info","cell_type")
dat <- transform(dat,
                 sex               = factor(sex),
                 disease_status    = factor(disease_status),
                 quality           = factor(quality),
                 submitted_quality = factor(submitted_quality),
                 clinical_info     = factor(clinical_info),
                 cell_type         = factor(cell_type))
rows <- match(sample_info$id,dat$id)
dat  <- dat[rows,]
sample_info <- cbind(sample_info,dat[-1])

k <- ncol(L)
colnames(L) <- paste0("k",1:k)
L   <- as.data.frame(L)
pdat <- cbind(sample_info,L)
fit <- lm(k1 ~ sex + disease_status + quality +
               submitted_quality + age + BMI,
          pdat)
ggplot(pdat,aes(x = sex,y = k1)) +
  geom_boxplot() +
  theme_cowplot()
