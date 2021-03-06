# Author: Lianqing Zhang
# Date: 2022.03.30
# Version:1
# The purpose of this R script is to analyze hippocampal/amygdala subfields/subregions volume automatically,
# including linear model with input designs, generate tables and plots ready for publish.

# -------------------- 1.read inputs ------------------------
# inputs should include:
# 1 data, generated by hippocampal-qa pipeline, with id, volumes(generated 
#   automatically), clinical data(should already be merged previously)
# 2.filepath to design table

# set analysis settings
args <- commandArgs(trailingOnly=TRUE)
#group = 'Diagnosis'
group <- args[1]
#covariates = c("ICV","Age","Sex","Education")
covariates <- unlist(strsplit(args[2],"__"))
#x_list = c("HAMD_Total_17","HAMA_Total") #used for correlation analysis
x_list <- unlist(strsplit(args[3],"__"))
#inter_terms = c('Diagnosis*Age','Diagnosis*Sex','Diagnosis*Education')  # interaction terms
inter_terms <- unlist(strsplit(args[4],"__"))

#set data input and output
#out_dir = "D:/hmrrc/Lmri/SampleData/statisitcs"
out_dir <- args[5]
#data_path = ""
data_path <- args[6]

# config pic size and legend position
#pic_size_box = c(3,7) # width,height, should be less than 10, otherwise the front would be too small
pic_size_box <- as.numeric(unlist(strsplit(args[7],"__")))
#legend.position.box = "none"
legend.position.box <- args[8]
#BoxColorSetting = c("darkolivegreen4","darkorchid3")
BoxColorSetting <- args[9]
if (BoxColorSetting != "none"){BoxColorSetting <- unlist(strsplit(args[9],"__"))}
#pic_size_cor=c(7,3)  # width,height, should be less than 10, otherwise the front would be too small
pic_size_cor <- as.numeric(unlist(strsplit(args[10],"__")))
#legend.position.cor = "top"
legend.position.cor <- args[11]
################### should not change after this ###################

data = read.csv(data_path)
setwd(out_dir)

design = paste(c(group,covariates),collapse ="+")
interaction_ds = paste(c(design,inter_terms),collapse ="+")

# set log file
sink("log.txt",split = F)
print(paste0('Reading data from: ',data_path))
print(paste0("saving outputs to ",out_dir))
print(paste0("Analyzing GLM design: ",design,",with group:",group))


# -------------- end of 1.read inputs ------------------------


## -------------------- 2.setting of functions ------------------------
library(ggplot2)
library(ppcor)
library(car)
library(sjstats)
library(DescTools)
### -------------------- 2.1 setting of vectors ------------------------
adj.HBT=c("Whole_hippocampal_head","Whole_hippocampal_body","Hippocampal_tail")
adj.sub=c("CA1", "CA3", "CA4", "molecular_layer_HP","GC.ML.DG","subiculum","presubiculum","parasubiculum","fimbria","HATA")
adj.Amyg = c("Lateral.nucleus","Basal.nucleus","Accessory.Basal.nucleus" ,
             "Central.nucleus","Medial.nucleus" ,"Cortical.nucleus" , 
             "Anterior.amygdaloid.area.AAA", "Corticoamygdaloid.transitio","Paralaminar.nucleus")

hippocampal_subfields <- c("Whole_hippocampus",adj.HBT,adj.sub)
amygdala_subregions <- c("Whole_amygdala",adj.Amyg )
all_subfield = c(hippocampal_subfields,amygdala_subregions)

L_sub_name=paste("lh",all_subfield,sep="_")
R_sub_name=paste("rh",all_subfield,sep="_")
AI_sub_name=paste("AI",all_subfield,sep="_")
AVG_sub_name=paste("AVG",all_subfield,sep="_")

flags_lst = c('lh','rh','AI','AVG')
all_scalar = lapply(flags_lst, paste,all_subfield,sep="_")

# calculate hippocampal scalar   #ToDo:move this to python
cal_hippocampal_scalar <- function(data){
  data[,AI_sub_name] <- 100*(data[,L_sub_name]-data[,R_sub_name]) / (data[,L_sub_name]/2+data[,R_sub_name]/2)
  data[,AVG_sub_name] <- (data[,L_sub_name]/2+data[,R_sub_name]/2)
  data
}
data = cal_hippocampal_scalar(data)

transf_str2factor <- function(df,columnnames){
  for (columnname in columnnames){
    if (typeof(df[[columnname]])=="character"){
      df[[columnname]] = factor(df[[columnname]])
    } else {
      df[[columnname]] = df[[columnname]]}
  }
  df
}
data = transf_str2factor(data,c(group,covariates,x_list))
### --------------------end of  2.1 setting of vectors ------------------------

### -------------------- 2.2 description function ------------------------
#setting description function
mydescrib <- function(x) {
  m <- mean(x,na.rm = T)
  n <- length(na.omit(x))
  s <- sd(x,na.rm = T)
  Mean.sd <- paste(round(m,2)," (", round(s,1),")",sep="")
  return(c(n=n, Mean.sd=Mean.sd ))}

mydescribByGroup <- function(i,group,data){
  describe <- aggregate(data[[i]], by = list(data[[group]]), FUN=mydescrib)
  describe_colnames = describe$Group.1
  describe_table<- as.data.frame(describe$x)
  describe_table <- as.data.frame(sapply(list(describe$x)[[1]],cbind.data.frame))
  colnames(describe_table) = c(paste('n',describe_colnames),paste(describe_colnames,'mean (S.D.)'))
  describe_table
}
### --------------------end of 2.2 description function ------------------------

### -------------------- 2.3 main function ------------------------

Do_MANCOVA <- function(region_list,group,design,data){
  out_df <- data.frame(matrix(NA,0,2*length(levels(data[[group]]))))
  for (i in region_list){
    describe_table = mydescribByGroup(i,group,data)
    y<- data[[i]]
    lm_model = lm(eval(parse(text = paste0("y ~ ", design))), data)
    fit=anova_stats(Anova(lm_model,type=3))
    describe_table[1,c("F","EffectSize","p.value")] = fit[group,c("statistic","etasq","p.value")]
    if (describe_table[1,"EffectSize"] < 0.001){describe_table[1,"EffectSize"]="<0.001"}
    row.names(describe_table)[1]<- i
    out_df <- rbind.data.frame(out_df,describe_table)
  }
  out_df
}
Do_FDR <- function(flag,in_df){
  for (adj in list(adj.HBT,adj.sub,adj.Amyg)){
    full_name = paste(flag,adj,sep='_')
    in_df[full_name,"P.FDR"]=round(p.adjust(in_df[full_name,"p.value"],method = "fdr"),3)
  }
  in_df
}
Do_post_hoc <- function(in_res_table,data,group,p.col,adj='fdr'){
  n_levels=length(levels(data[[group]]))
  for (i in rownames(in_res_table)){
    pairwise_res = pairwise.t.test(data[[i]],data[[group]],p.adjust.method = adj)
    post_hoc_p = round(as.data.frame(pairwise_res$p.value),3)
    for (colvar in colnames(post_hoc_p)){
      for (rowvar in rownames(post_hoc_p)){
        if (is.na(post_hoc_p[rowvar,colvar]) == FALSE){
          contrast_var = paste(colvar,rowvar,sep = ' vs. ')
          if (in_res_table[i,p.col]<0.05){
            in_res_table[i,contrast_var]=post_hoc_p[rowvar,colvar]
          }else{
            #print '-'
            in_res_table[i,contrast_var]='-'
          }
        }
      }
    }
  }
  in_res_table
}
find_significance <- function(in_res_table,alpha=0.05,p_col){
  level2 = row.names(in_res_table[in_res_table[,p_col]<0.001,])
  in_res_table[in_res_table[p_col]<alpha & in_res_table[p_col]>=0.001,p_col] = paste(in_res_table[in_res_table[p_col]<alpha & in_res_table[p_col]>=0.001,p_col],"*",sep = '')
  in_res_table[level2,p_col] = "<0.001*"
  in_res_table
}
Wrap_MANCOVA <- function(flag,group,design,data){
  print(paste0("flag:",flag))
  out_df_name = paste('MANCOVA',flag,'csv',sep = '.')
  region_list = paste(flag,all_subfield,sep="_")
  MANCOVA_out_df = Do_MANCOVA(region_list,group,design,data)
  MANCOVA_out_df = Do_FDR(flag,MANCOVA_out_df)
  raw_p_names = row.names(MANCOVA_out_df[is.na(MANCOVA_out_df$P.FDR),])
  MANCOVA_out_df[raw_p_names,'P.FDR']=MANCOVA_out_df[raw_p_names,'p.value']
  if (length(levels(data[[group]])) > 2){MANCOVA_out_df = Do_post_hoc(MANCOVA_out_df,data,group,p.col = 'P.FDR')}
  
  #MANCOVA_out_df[raw_p_names,'P.FDR'] = '-'
  MANCOVA_out_df = find_significance(MANCOVA_out_df,0.05,'P.FDR')
  print(paste0('finished, saving outputs to ',out_df_name))
  write.csv(MANCOVA_out_df,out_df_name)
  MANCOVA_out_df
}


Do_interaction <- function(region_list,design,data){
  inter_terms = strsplit(design,"+",fixed = T)[[1]]
  ResTable=data.frame(matrix(NA,0,length(inter_terms)))
  for (i in region_list ){
    y=data[[i]]
    lm_model = lm(eval(parse(text = paste0("y ~ ", design))), data)
    fit=Anova(lm_model,type=3)
    ResTable[i,]=round(fit[1:length(inter_terms)+1,"Pr(>F)"],3)}
  colnames(ResTable)=rownames(fit)[1:length(inter_terms)+1]
  ResTable
}
Wrap_interaction <- function(flag,design,data){
  print(paste0("flag:",flag))
  out_df_name = paste('interaction',flag,'csv',sep = '.')
  region_list = paste(flag,all_subfield,sep="_")
  int_out_df = Do_interaction(region_list,design,data)
  for (c in colnames(int_out_df)){
    int_out_df = find_significance(int_out_df,p_col = c)
  }
  print(paste0('finished, saving outputs to ',out_df_name))
  write.csv(int_out_df ,out_df_name)
  int_out_df 
}

# Partial Corelation
Marix2DF <- function(in_matrix, colnames.list, rownames.list){
  out_df = as.data.frame(in_matrix)
  out_df = data.frame(lapply(out_df,as.numeric),stringsAsFactors = FALSE)
  rownames(out_df) <- rownames.list
  colnames(out_df) <- colnames.list
  out_df
}

Do_PartialCorr <- function(x.list,y.list,covariates,data = data){
  
  p.matrix <- list()
  r.matrix <- list()
  n.matrix <- list()
  for (x in x.list){
    pvalues <- list()
    rvalues <- list()
    nvalues <- list()
    for (y  in y.list){
      y.cordata <- na.omit(data[,c(x,y,covariates)])
      y.cordata <- as.data.frame(lapply(y.cordata, as.numeric))
      pcor <- pcor.test(y.cordata[,1],y.cordata[,2],y.cordata[,-c(1:2)])
      pvalues <- append(pvalues,round(pcor$p.value,3))
      rvalues <- append(rvalues,round(pcor$estimate,2))
      nvalues <- append(nvalues,pcor$n)}
    p.matrix <- cbind(p.matrix, pvalues)
    r.matrix <- cbind(r.matrix, rvalues)
    n.matrix <- cbind(n.matrix, nvalues)
  }
  out = lapply(list(r = r.matrix, p = p.matrix ,n = n.matrix), Marix2DF,x.list,y.list)
  return(out)
}
Do_Compare_Correlation <- function(r1.matrix, r2.matrix, n1.matrix, n2.matrix){
  
  z1 <- FisherZ(r1.matrix)
  z2 <- FisherZ(r2.matrix)
  z.observed <- (z1-z2)/sqrt( 1/(n1.matrix-3)+1/(n2.matrix-3) )
  pvalue2sided=2*pnorm(-abs(as.matrix(z.observed)))
  pvalue2sided
}

Wrap_PartialCorr<- function(flag,x_list,covariates,data, out_df_flag){
  print(paste0("saving outputs to:"))
  print(paste('ParCorr',out_df_flag,flag,'*.csv',sep = '.'))
  region_list = paste(flag,all_subfield,sep="_")
  Corr_res = Do_PartialCorr(region_list, x_list,covariates,data)
  
  for (i in c('r','p','n')){
    write.csv(Corr_res[[i]],paste('ParCorr',out_df_flag,flag,i,'csv',sep = '.'))
  }
  Corr_res
}

### --------------------end of 2.3 main function ------------------------

### -------------------- 2.4 plot function ------------------------

# box plot
Draw_boxplot <- function(x,data,group){
  pic_name = paste("BoxPlot",x,"png",sep = '.')
  print(paste0("saving plot to ",pic_name))
  out_p = ggplot(data,aes_string(x=group,y=x,color=group))+ 
   geom_boxplot(aes_string(group=group), fill=NA)+ 
    geom_jitter(aes_string(color = group),shape=16, position=position_jitter(0.2))+
    labs(y = x)+
    theme(legend.position=legend.position.box)
  if (BoxColorSetting!= "none"){out_p = out_p+scale_color_manual(values=BoxColorSetting)}
  ggsave(pic_name,out_p,width = pic_size_box[1],height = pic_size_box[2])
  out_p
}
Wrap_Draw_boxplot <-function(region_list,data,group){
  lapply(region_list,Draw_boxplot,data=data,group=group)
}

# correlation plot with residual
Dw_cor_residual_by_group <- function(x,y,covariates,group,data){
  pic_name = paste("CorPlot",paste0(x,'~',y),"png",sep = '.')
  print(paste0("saving plot to ",pic_name))
  plot_data = na.omit(data[,c(group,covariates,x,y)])
  Y_lm = lm(eval(parse(text = paste0(y, " ~ ", paste(covariates,collapse ="+")))), plot_data)
  Y_resid<-resid(Y_lm) 
  X_lm = lm(eval(parse(text = paste0(x, " ~ ", paste(covariates,collapse ="+")))), plot_data)
  X_resid<-resid(X_lm) 
  out_p = ggplot(plot_data,aes_string(X_resid,Y_resid,group=group,color=group)) +
    geom_point()+geom_smooth(formula = y ~ x, method = "lm")+
    labs(x = paste0('residual of ',x),y = paste0('residual of ',y))+
    theme(legend.position=legend.position.cor)
  ggsave(pic_name,out_p,width = pic_size_cor[1],height = pic_size_cor[2])
  out_p
}
Wrap_Dw_cor_residual_by_group <-function(x_list, y_list,covariates,group,data ){
  lapply(x_list, function(x){
    lapply(y_list,function(y){
      Dw_cor_residual_by_group(x,y,covariates=covariates,group=group,data=data)})
  })
}

### --------------------end of 2.4 plot function ------------------------

## ------------- end of 2.setting of functions ------------------------

## -------------------- 3.interaction ------------------------
print(" ")
print("start: interaction analysis")
print('testing interaction with design:')
print(interaction_ds)
out_int_res = lapply(flags_lst, Wrap_interaction,design=interaction_ds,data=data)
print("finished: interaction analysis ")
print(" ")

## ------------- end of 3.interaction ------------------------

## -------------------- 4.group difference ------------------------
print(" ")
print("start: MANCOVA analysis")
print("testing design:")
print(design)
out_MANCOVA_res = lapply(flags_lst, Wrap_MANCOVA,group=group,design=design,data=data)

print("Draw box plots:")
BoxPlots=lapply(unlist(all_scalar),Draw_boxplot,data=data,group=group)

print("finished: MANCOVA analysis")
print(" ")

## ------------- end of 4.group difference ------------------------

## -------------------- 4.partial correlation ------------------------
print(" ")
print("start: partial correlation analysis")
print("using covariates:")
print(paste(covariates,collapse = ', '))

print("1.Test partial correlation in the whole sample ... ")
out_df_flag = 'Whole'
ParCorr_Whole = lapply(flags_lst,Wrap_PartialCorr,x_list,covariates,data, out_df_flag)

print(" ")
print("2.Test partial correlation in each group ... ")
# update x_list to see if data were collected in both groups
update_covariates_lst = function(x_list){
  x_list_updated = c()
  for (x in x_list){
    x_levels = length(levels(factor(as.character(na.omit(data[,c(group,x)])[[group]]))))
    if (x_levels >1){
      x_list_updated = c(x_list_updated,x)
    }}
  return(x_list_updated)}
x_list_updated = update_covariates_lst(x_list)

if (length(x_list_updated) >0){
  
  print("testing correlation variable with multiple levels...")
  
  # 2.Test partial correlation in each group ...
  ParCorr_by_group = list()
  for (gp in levels(data[[group]])){
    print(paste0("Analyzing group: ",gp))
    data_gp = data[data[,group]==gp,]
    out_df_flag = gp
    ParCorr = lapply(flags_lst,Wrap_PartialCorr,x_list_updated,covariates,data_gp, out_df_flag)
    names(ParCorr) = paste(gp,flags_lst,sep = '_')
    ParCorr_by_group = append(ParCorr_by_group,ParCorr)
    print("Done.")
  }
  
  # 3.Compare correlation coefficient differences between groups
  print(" ")
  print("3.Compare correlation coefficient differences between groups")
  #ToDo:how about 3 groups?
  for (flag in flags_lst) {
    gps = levels(data[[group]])
    out_p = Do_Compare_Correlation(ParCorr_by_group[[paste(gps[1],flag,sep = "_")]]$r, 
                                   ParCorr_by_group[[paste(gps[2],flag,sep = "_")]]$r, 
                                   ParCorr_by_group[[paste(gps[1],flag,sep = "_")]]$n, 
                                   ParCorr_by_group[[paste(gps[2],flag,sep = "_")]]$n)
    out_df_flag = paste(gps,collapse = '~')
    out_p_name = paste('ParCorr',out_df_flag,flag,'p','csv',sep = '.')
    write.csv(out_p,out_p_name)
  }
  
} else{
  print("correlation variables only have 1 level, nothing to do.")
}


print("Drawing correlation plots with residuals")
print("x:")
print(paste(x_list,collapse = ', '))
print("y: volumes")

CorPlots = lapply(all_scalar,function(y){
  Wrap_Dw_cor_residual_by_group(x_list,y,covariates,group,data)
})
print("finished: partial correlation analysis")
print(" ")
## ------------- end of 4.partial correlation ------------------------