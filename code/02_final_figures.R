
library(ggplot2)

# 1. Load results
setwd("~/Desktop/biomedical-kpca-benchmark")
out <- "results/benchmark_results_final"
res <- readRDS(file.path(out, "complete_benchmark_results_6datasets_FINAL.rds"))
outer <- res$outer_results
summ <- res$summary_results

stopifnot(nrow(outer) == 720, nrow(summ) == 48)

ds <- c("Colon", "Gravier", "Leukemia", "Pima", "Prostate", "Wisconsin")
reps <- c(raw="Raw", pca="PCA", rbf_kpca="RBF-KPCA", local_kpca="Local-KPCA")
cols <- c(Raw="#00BFC4", PCA="#7CAE00",
          "RBF-KPCA"="#C77CFF", "Local-KPCA"="#F8766D")
shapes <- c(Linear=16, Radial=17)

label_data <- function(d) {
  d$Dataset <- factor(d$Dataset, levels=ds)
  d$Rep <- factor(reps[d$Representation], levels=unname(reps))
  d$Class <- factor(d$Classifier, levels=c("linear","radial"),
                    labels=c("Linear","Radial"))
  d
}
outer <- label_data(outer)
summ <- label_data(summ)

save_fig <- function(p, name, h) {
  for (ext in c("pdf","png"))
    ggsave(file.path(out, paste0(name,".",ext)),
           p, width=11, height=h, units="in", dpi=600)
}

base_theme <- theme_bw(base_size=11) +
  theme(strip.background=element_rect(fill="grey82", color="grey40"),
        strip.text=element_text(face="bold", size=11),
        axis.text=element_text(size=9),
        panel.grid.major=element_line(color="grey88", linewidth=.35),
        panel.grid.minor=element_blank(),
        panel.spacing=grid::unit(.25,"cm"),
        legend.position="bottom")

# 2. Figure 1: Outer-test ROC-AUC
methods <- as.vector(rbind(paste(unname(reps),"Linear",sep=" + "),
                           paste(unname(reps),"Radial",sep=" + ")))
outer$Method <- factor(paste(outer$Rep,outer$Class,sep=" + "),
                       levels=rev(methods))
means <- aggregate(ROC_AUC ~ Dataset+Method+Class+Rep, outer, mean)

fig1 <- ggplot(outer, aes(ROC_AUC,Method,color=Rep,shape=Class)) +
  geom_point(alpha=.45,size=2,
             position=position_jitter(width=0,height=.045)) +
  geom_point(data=means,aes(ROC_AUC,Method,shape=Class),
             inherit.aes=FALSE,color="black",size=3.6) +
  facet_wrap(~Dataset,ncol=2,scales="free_x") +
  scale_color_manual(values=cols,
                     labels=c("Raw","PCA","RBF-KPCA","Locally scaled KPCA")) +
  scale_shape_manual(values=shapes) +
  labs(x="Outer-test ROC-AUC",y=NULL,
       color="Representation",shape="Classifier") +
  base_theme +
  theme(strip.background=element_rect(fill="grey88",color="grey40"))

print(fig1)
save_fig(fig1,"Figure1_Outer_Test_ROC_AUC_6Datasets",7.4)

# 3. Figure 2: Matched ROC-AUC differences
keys <- c("Dataset","Outer_Repeat","Outer_Fold","Classifier")
raw <- outer[outer$Representation=="raw",c(keys,"ROC_AUC")]
trans <- outer[outer$Representation!="raw",
               c(keys,"Representation","ROC_AUC")]
names(raw)[5] <- "Raw_AUC"
names(trans)[6] <- "Trans_AUC"

matched <- merge(trans,raw,by=keys)
matched$Delta <- matched$Trans_AUC-matched$Raw_AUC
matched$Contrast <- paste0(reps[matched$Representation],
                           " vs Raw: ",
                           ifelse(matched$Classifier=="linear","Linear","Radial"))
contrast_order <- as.vector(rbind(
  paste0(unname(reps[-1])," vs Raw: Linear"),
  paste0(unname(reps[-1])," vs Raw: Radial")))
matched$Contrast <- factor(matched$Contrast,levels=rev(contrast_order))
stopifnot(nrow(matched)==540)

mean_delta <- aggregate(Delta ~ Dataset+Contrast,matched,mean)

fig2 <- ggplot(matched,aes(Delta,Contrast)) +
  geom_vline(xintercept=0,color="grey30",linewidth=.4) +
  geom_point(color="grey50",alpha=.35,size=1.4,
             position=position_jitter(width=0,height=.045,seed=2026)) +
  geom_point(data=mean_delta,
             aes(Delta,Contrast,
                 color=sub(" vs Raw:.*$","",as.character(Contrast))),
             inherit.aes=FALSE,shape=18,size=2.8) +
  scale_color_manual(values=cols[-1],guide="none") +
  facet_wrap(~Dataset,ncol=2,scales="free_x") +
  labs(title="Change in ROC-AUC Relative to Raw Predictors",
       x=expression(Delta*"ROC-AUC"),y=NULL) +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face="bold",size=13,
                                margin=margin(b=8)),
        strip.background=element_rect(fill="grey90",
                                      color="grey70",linewidth=.35),
        strip.text=element_text(face="bold",size=9),
        axis.text=element_text(size=8),
        panel.grid.major=element_line(color="grey92",linewidth=.25),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="grey70",fill=NA,linewidth=.35),
        panel.spacing=grid::unit(.25,"cm"))

print(fig2)
save_fig(fig2,"Figure2_Matched_Delta_ROC_AUC_6Datasets",6.4)

# 4. Figure 3: Performance versus computation time
summ$Rep <- factor(summ$Rep,
                   levels=c("Local-KPCA","PCA","Raw","RBF-KPCA"))

fig3 <- ggplot(summ,
               aes(Mean_Elapsed_Seconds,ROC_AUC_Mean,color=Rep,shape=Class)) +
  geom_point(size=4.6,alpha=.95) +
  facet_wrap(~Dataset, ncol=2, scales="free_y") +
  scale_y_continuous(expand=expansion(mult=c(0.08, 0.20))) +
  scale_x_log10() +
  scale_color_manual(values=cols,
                     breaks=c("Local-KPCA","PCA","Raw","RBF-KPCA")) +
  scale_shape_manual(values=shapes) +
  labs(title="Predictive performance versus computational cost",
       x="Mean computation time per outer evaluation (seconds, log scale)",
       y="Mean ROC-AUC",color="Representation",shape="Classifier") +
  base_theme +
  theme(plot.title=element_text(face="bold",size=13,
                                margin=margin(b=8)))

print(fig3)
save_fig(fig3,"Figure3_Performance_vs_Cost_6Datasets",7.8)
