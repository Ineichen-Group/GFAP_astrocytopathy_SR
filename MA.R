###call libraries
library(metafor)
library(meta)
library(readxl)
library(weightr)

###read file
GFAP_MA <- read_excel("C:/Users/benja/Dropbox/Research/Experiments/Karolinska/NMOSD-review/Meta_analysis/MA_masterfile.xlsx", sheet = 1)

#Remove children studies
GFAP_MA <- subset(GFAP_MA, MA_included!=0)

############Clinical MA
### tumor
GFAP_MA_neoplasm <- subset(GFAP_MA, N_neoplasm>=0)
GFAP_MA_neoplasm$neoplasm_percentage <- (GFAP_MA_neoplasm$N_neoplasm/GFAP_MA_neoplasm$N_total_neoplasm)*100

GFAP_MA_neoplasm.ordered <- GFAP_MA_neoplasm[order(GFAP_MA_neoplasm$neoplasm_percentage), ]
pes.forest <- metaprop(N_neoplasm, N_total_neoplasm, Author_year, data = GFAP_MA_neoplasm.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_neoplasm.png", width = 2600, height = 1250, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### viral prodrome
GFAP_MA_viral <- subset(GFAP_MA, N_viral>=0)
GFAP_MA_viral$GFAP_MA_viral_percentage <- (GFAP_MA_viral$N_viral/GFAP_MA_viral$N_total_viral)*100

GFAP_MA_viral.ordered <- GFAP_MA_viral[order(GFAP_MA_viral$GFAP_MA_viral_percentage), ]
pes.forest <- metaprop(N_viral, N_total_viral, Author_year, data = GFAP_MA_viral.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_prodrome.png", width = 2600, height = 1000, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### coexisting antibodies
GFAP_MA_coexistab <- subset(GFAP_MA, N_coexistab>=0)
GFAP_MA_coexistab$GFAP_MA_coexistab_percentage <- (GFAP_MA_coexistab$N_coexistab/GFAP_MA_coexistab$N_total_coexistab)*100

GFAP_MA_coexistab.ordered <- GFAP_MA_coexistab[order(GFAP_MA_coexistab$GFAP_MA_coexistab_percentage), ]
pes.forest <- metaprop(N_coexistab, N_total_coexistab, Author_year, data = GFAP_MA_coexistab.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_coexistab.png", width = 2600, height = 1200, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()


### immunotherapy response
GFAP_MA_immunotherapyresponse <- subset(GFAP_MA, N_immunotherapyresponse>=0)
GFAP_MA_immunotherapyresponse$GFAP_MA_immunotherapyresponse_percentage <- (GFAP_MA_immunotherapyresponse$N_immunotherapyresponse/GFAP_MA_immunotherapyresponse$N_total_immunotherapyresponse)*100

GFAP_MA_immunotherapyresponse.ordered <- GFAP_MA_immunotherapyresponse[order(GFAP_MA_immunotherapyresponse$GFAP_MA_immunotherapyresponse_percentage), ]
pes.forest <- metaprop(N_immunotherapyresponse, N_total_immunotherapyresponse, Author_year, data = GFAP_MA_immunotherapyresponse.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_immunotherapyresponse.png", width = 2600, height = 1050, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

############Imaging MA
### LME
GFAP_MA_LME <- subset(GFAP_MA, N_LME>=0)
GFAP_MA_LME$LME_percentage <- (GFAP_MA_LME$N_LME/GFAP_MA_LME$N_total_LME)*100

GFAP_MA_LME.ordered <- GFAP_MA_LME[order(GFAP_MA_LME$LME_percentage), ]
pes.forest <- metaprop(N_LME, N_total_LME, Author_year, data = GFAP_MA_LME.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_LME.png", width = 2600, height = 1100, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### radial
GFAP_MA_radial <- subset(GFAP_MA, N_radial>=0)
GFAP_MA_radial$radial_percentage <- (GFAP_MA_radial$N_radial/GFAP_MA_radial$N_total_radial)*100

GFAP_MA_radial.ordered <- GFAP_MA_radial[order(GFAP_MA_radial$radial_percentage), ]
pes.forest <- metaprop(N_radial, N_total_radial, Author_year, data = GFAP_MA_radial.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_radial.png", width = 2600, height = 1250, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### gd
GFAP_MA_gd <- subset(GFAP_MA, N_gd>=0)
GFAP_MA_gd$gd_percentage <- (GFAP_MA_gd$N_gd/GFAP_MA_gd$N_total_gd)*100

GFAP_MA_gd.ordered <- GFAP_MA_gd[order(GFAP_MA_gd$gd_percentage), ]
pes.forest <- metaprop(N_gd, N_total_gd, Author_year, data = GFAP_MA_gd.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_gd.png", width = 2600, height = 1200, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### t2
GFAP_MA_t2 <- subset(GFAP_MA, N_t2>=0)
GFAP_MA_t2$t2_percentage <- (GFAP_MA_t2$N_t2/GFAP_MA_t2$N_total_t2)*100

GFAP_MA_t2.ordered <- GFAP_MA_t2[order(GFAP_MA_t2$t2_percentage), ]
pes.forest <- metaprop(N_t2, N_total_t2, Author_year, data = GFAP_MA_t2.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_t2.png", width = 2600, height = 1100, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### t2perivent
GFAP_MA_t2perivent <- subset(GFAP_MA, N_t2perivent>=0)
GFAP_MA_t2perivent$t2perivent_percentage <- (GFAP_MA_t2perivent$N_t2perivent/GFAP_MA_t2perivent$N_total_t2perivent)*100

GFAP_MA_t2perivent.ordered <- GFAP_MA_t2perivent[order(GFAP_MA_t2perivent$t2perivent_percentage), ]
pes.forest <- metaprop(N_t2perivent, N_total_t2perivent, Author_year, data = GFAP_MA_t2perivent.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_t2perivent.png", width = 2600, height = 900, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### t2wm
GFAP_MA_t2wm <- subset(GFAP_MA, N_t2wm>=0)
GFAP_MA_t2wm$t2wm_percentage <- (GFAP_MA_t2wm$N_t2wm/GFAP_MA_t2wm$N_total_t2wm)*100

GFAP_MA_t2wm.ordered <- GFAP_MA_t2wm[order(GFAP_MA_t2wm$t2wm_percentage), ]
pes.forest <- metaprop(N_t2wm, N_total_t2wm, Author_year, data = GFAP_MA_t2wm.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_t2wm.png", width = 2600, height = 800, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### sc
GFAP_MA_sc <- subset(GFAP_MA, N_sc>=0)
GFAP_MA_sc$sc_percentage <- (GFAP_MA_sc$N_sc/GFAP_MA_sc$N_total_sc)*100

GFAP_MA_sc.ordered <- GFAP_MA_sc[order(GFAP_MA_sc$sc_percentage), ]
pes.forest <- metaprop(N_sc, N_total_sc, Author_year, data = GFAP_MA_sc.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_sc.png", width = 2600, height = 950, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### brainnormal
GFAP_MA_brainnormal <- subset(GFAP_MA, N_brainnormal>=0)
GFAP_MA_brainnormal$brainnormal_percentage <- (GFAP_MA_brainnormal$N_brainnormal/GFAP_MA_brainnormal$N_total_brainnormal)*100

GFAP_MA_brainnormal.ordered <- GFAP_MA_brainnormal[order(GFAP_MA_brainnormal$brainnormal_percentage), ]
pes.forest <- metaprop(N_brainnormal, N_total_brainnormal, Author_year, data = GFAP_MA_brainnormal.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_brainnormal.png", width = 2600, height = 1100, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()

### scnormal
GFAP_MA_scnormal <- subset(GFAP_MA, N_scnormal>=0)
GFAP_MA_scnormal$scnormal_percentage <- (GFAP_MA_scnormal$N_scnormal/GFAP_MA_scnormal$N_total_scnormal)*100

GFAP_MA_scnormal.ordered <- GFAP_MA_scnormal[order(GFAP_MA_scnormal$scnormal_percentage), ]
pes.forest <- metaprop(N_scnormal, N_total_scnormal, Author_year, data = GFAP_MA_scnormal.ordered, sm = "PLO", method.ci = "NAsm", method.tau = "DL", incr = 0.5, allincr = F, addincr = F, title = "")
png("forestplot_meta_scnormal.png", width = 2600, height = 900, res = 300)
forest(pes.forest,
       xlim = c(0,1), pscale = 1,
       rightcols = c("effect", "ci", "w.random"),
       rightlabs = c("Proportion", "95%-CI", "Weights"),
       leftcols = c("studlab", "event", "n"),
       leftlabs = c("Study", "Cases", "Total"),
       clab = "Prevalence", 
       xlab = "Proportion",
       fs.xlab = 12,
       fs.study = 12,
       fs.study.lables = 12,
       fs.heading = 12,
       squaresize = 0.5, col.square = "navy", col.square.lines = "navy",
       col.diamond = "maroon", col.diamond.lines = "maroon",
       comb.fixed = FALSE,
       lty.fixed = 0,
       lty.random = 2,
       type.study = "square",
       type.random = "diamond",
       ff.fixed = "bold.italic",
       ff.random = "bold.italic",
       hetlab = "Heterogeneity",
       fs.hetstat = 10,
       smlab = "",
       printQ = TRUE,
       print.pval.Q = TRUE,
       print.I2 = TRUE,
       print.tau2 = TRUE,
       col.by = "black",
       digits = 3)
dev.off()





