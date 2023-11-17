load("sctype_time.rds")
load("seurat_time.rds")
load("rctd_time.rds")
load("spotlight_time.rds")

#proc.time 
#The first two entries are the total user and system CPU times of the current R 
#process and any child processes on which it has waited, and the third entry is
#the ‘real’ elapsed time since the process was started.

times <- data.frame(method = c("sctype", "seurat", "rctd", "spotlight"),
                    first = c(sctype_time_1[["elapsed"]],
                              seurat_time_1[["elapsed"]],
                              rctd_time_1[["elapsed"]],
                              spotlight_time_1[["elapsed"]]),
                    second = c(sctype_time_2[["elapsed"]],
                              seurat_time_2[["elapsed"]],
                              rctd_time_2[["elapsed"]],
                              spotlight_time_2[["elapsed"]]),
                    third = c(sctype_time_3[["elapsed"]],
                              seurat_time_3[["elapsed"]],
                              rctd_time_3[["elapsed"]],
                              spotlight_time_3[["elapsed"]]),
                    mean = c(mean(c(sctype_time_1[["elapsed"]],
                                sctype_time_2[["elapsed"]],
                                sctype_time_3[["elapsed"]])),
                             mean(c(seurat_time_1[["elapsed"]],
                                  seurat_time_2[["elapsed"]],
                                  seurat_time_3[["elapsed"]])),
                             mean(c(rctd_time_1[["elapsed"]],
                                  rctd_time_2[["elapsed"]],
                                  rctd_time_3[["elapsed"]])),
                             mean(c(spotlight_time_1[["elapsed"]],
                                  spotlight_time_2[["elapsed"]],
                                  spotlight_time_3[["elapsed"]]))),
                    sd = c(sd(c(sctype_time_1[["elapsed"]],
                                sctype_time_2[["elapsed"]],
                                sctype_time_3[["elapsed"]])),
                           sd(c(seurat_time_1[["elapsed"]],
                                seurat_time_2[["elapsed"]],
                                seurat_time_3[["elapsed"]])),
                           sd(c(rctd_time_1[["elapsed"]],
                                rctd_time_2[["elapsed"]],
                                rctd_time_3[["elapsed"]])),
                           sd(c(spotlight_time_1[["elapsed"]],
                                spotlight_time_2[["elapsed"]],
                                spotlight_time_3[["elapsed"]]))))

ggplot(times, aes(x=factor(method, levels = c("sctype", "seurat", "rctd", "spotlight")), y=mean)) + 
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  geom_errorbar( aes(x=factor(method, levels = c("sctype", "seurat", "rctd", "spotlight")), ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  xlab("method")
