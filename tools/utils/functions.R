
## -------------------- 2.setting of functions ------------------------
library(ggplot2)
library(ppcor)
library(car)
library(sjstats)
library(DescTools)
library(multcomp)
library(ggsignif)
### -------------------- 2.1 setting of vectors ------------------------
adj.HBT=c("Whole_hippocampal_head","Whole_hippocampal_body","Hippocampal_tail")
adj.sub=c("CA1", "CA3", "CA4", "molecular_layer_HP","GC.ML.DG","subiculum","presubiculum","parasubiculum","fimbria","HATA")
adj.Amyg = c("Lateral.nucleus","Basal.nucleus","Accessory.Basal.nucleus" ,
             "Central.nucleus","Medial.nucleus" ,"Cortical.nucleus" , 
             "Anterior.amygdaloid.area.AAA", "Corticoamygdaloid.transitio","Paralaminar.nucleus")
flags_lst = c('lh','rh','AI','AVG')

hippocampal_subfields <- c("Whole_hippocampus",adj.HBT,adj.sub)
amygdala_subregions <- c("Whole_amygdala",adj.Amyg )
all_subfield = c(hippocampal_subfields,amygdala_subregions)

all_scalar = lapply(flags_lst, paste,all_subfield,sep="_")
names(all_scalar) = flags_lst

all_scalar = apply(flags_lst, paste,all_subfield,sep="_")



transf_str2factor <- function(df,columnnames){
  for (columnname in columnnames){
    if (typeof(df[[columnname]])=="character"){
      df[[columnname]] = factor(df[[columnname]])
    } else {
      df[[columnname]] = df[[columnname]]}
  }
  df
}

update_covariates_lst = function(x_list){
  # update x_list to see if data were collected in both groups
  x_list_updated = c()
  for (x in x_list){
    x_levels = length(levels(factor(as.character(na.omit(data[,c(group,x)])[[group]]))))
    if (x_levels >1){
      x_list_updated = c(x_list_updated,x)
    }}
  return(x_list_updated)}

### --------------------end of  2.1 setting of vectors ------------------------

### -------------------- 2.2 description function ------------------------
#setting description function
mydescrib <- function(x) {
  m <- mean(x,na.rm = T)
  n <- length(na.omit(x))
  s <- sd(x,na.rm = T)
  Mean.sd <- paste(round(m,2)," (", round(s,1),")",sep="")
  return(c(n=n, mean=m,sd=s,Mean.sd=Mean.sd ))}

mydescribByGroup <- function(y_list,group,data){
  describe <- aggregate(data[,y_list], by = list(data[[group]]), FUN=mydescrib)
  
  describe_colnames = describe$Group.1
  
  describe_table <- t(as.data.frame(sapply(colnames(describe[2:length(describe)]), function(x){
    as.data.frame(sapply(list(describe[x][[1]])[[1]],cbind.data.frame))
  })))
  
  colnames(describe_table) = unlist(lapply(c('n','mean','sd','Mean(S.D.)'),paste,describe_colnames))
  as.data.frame(describe_table)
}

### --------------------end of 2.2 description function ------------------------




### -------------------- 2.3 main function ------------------------

Do_MANCOVA <- function(y_list,group,design,data){
  describe_table = mydescribByGroup(y_list ,group,data)
  for (y in y_list){
    lm_model = lm(eval(parse(text = paste0(y," ~ ", design))), data)
    fit=anova_stats(Anova(lm_model,type=3))
    describe_table[y,c("F","EffectSize","p.value")] = fit[group,c("statistic","etasq","p.value")]
    if (length(levels(data[[group]])) >2 ){
      # Do post-hoc if has multiple levels 
      post_hoc <- summary(multcomp::glht(lm_model, linfct = do.call(mcp, setNames(list('Tukey'), group))))
      describe_table[y,gsub('-','vs',names(post_hoc$test$coefficients))]=post_hoc$test$pvalues
    }
  }
  describe_table$EffectSize = round(describe_table$EffectSize,3)
  describe_table[describe_table$EffectSize< 0.001,"EffectSize"] = "<0.001"
  return(describe_table)
}

Do_FDR <- function(flag,in_df){
  for (adj in list(adj.HBT,adj.sub,adj.Amyg)){
    full_name = paste(flag,adj,sep='_')
    in_df[full_name,"P.FDR"]=round(p.adjust(in_df[full_name,"p.value"],method = "fdr"),3)
  }
  in_df
}

find_significance <- function(in_res_table,alpha=0.05,p_col){
  level2 = row.names(in_res_table[in_res_table[,p_col]<0.001,])
  in_res_table[[p_col]] = round( in_res_table[[p_col]],3)
  in_res_table[in_res_table[p_col]<alpha & in_res_table[p_col]>=0.001,p_col] = paste(in_res_table[in_res_table[p_col]<alpha & in_res_table[p_col]>=0.001,p_col],"*",sep = '')
  in_res_table[level2,p_col] = "<0.001*"
  in_res_table
}



# Correlation
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

Do_PartialCorr <- function(x,y,covariates,in_data = data){
  y.cordata = na.omit(in_data[,c(x,y,covariates)])
  y.cordata <- as.data.frame(lapply(y.cordata, as.numeric))
  pcor <- pcor.test(y.cordata[,1],y.cordata[,2],y.cordata[,-c(1:2)])
  out_df = data.frame(
    p = pcor$p.value,
    r = pcor$estimate,
    z = FisherZ(pcor$estimate),
    n = pcor$n, 
    p.print = format_p(pcor$p.value)
  )
  
  return(y = out_df)
}

Do_Compare_Correlation <- function(z1,z2,n1,n2){
  z.observed <- (z1-z2)/sqrt( 1/(n1-3)+1/(n2-3) )
  pvalue2sided=2*pnorm(-abs(as.matrix(z.observed)))
  pvalue2sided
}

Wrap_Compare_Correlation <- function(ParRes1, ParRes2){
  z = lapply(list(ParRes1, ParRes2),function(ParRes){
    Par_res = as.data.frame(t(ParRes))
    z = unlist(Par_res$z,use.names = F)
    z
  })
  n = lapply(list(ParRes1, ParRes2),function(ParRes){
    Par_res = as.data.frame(t(ParRes))
    unlist(Par_res$n,use.names = F)
  })
  res = do.call(Do_Compare_Correlation,c(z,n))
  res = as.list(res)
  names(res) = colnames(ParRes1)
  return(res)
}



#2022/7/27 改到这里


Do_interaction <- function(region_list,design,data){
  inter_terms = strsplit(design,"+",fixed = T)[[1]]
  ResTable=data.frame(matrix(NA,0,length(inter_terms)))
  for (i in region_list ){
    y=data[[i]]
    lm_model = lm(eval(parse(text = paste0("y ~ ", design))), data)
    fit=Anova(lm_model,type=3)
    ResTable[i,]=round(fit[gsub("\\*",':',inter_terms),'Pr(>F)'],3)}
  colnames(ResTable)=gsub("\\*",':',inter_terms)
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
#Marix2DF <- function(in_matrix, colnames.list, rownames.list){
#  out_df = as.data.frame(in_matrix)
#  out_df = data.frame(lapply(out_df,as.numeric),stringsAsFactors = FALSE)
#  rownames(out_df) <- rownames.list
#  colnames(out_df) <- colnames.list
#  out_df
#}





### --------------------end of 2.3 main function ------------------------

### -------------------- 2.4 plot function ------------------------
# box plot
Draw_boxplot <- function(y_name,data,group,
                         p_values = 'none',
                         legend.position.box='none'){
  
  out_p = ggplot(data,aes_string(x=group,y=y_name,color=group))+ 
    geom_boxplot(aes_string(group=group), fill=NA)+ 
    geom_jitter(aes_string(color = group),shape=16, position=position_jitter(0.2))+
    labs(y = y_name,x="")+
    theme(legend.position=legend.position.box)
  
  if (BoxColorSetting!= "none"){out_p = out_p+scale_color_manual(values=BoxColorSetting)}
  if (class(p_values) == "list"){
    out_p = out_p+do.call(geom_signif,p_values)
  }
  
  return(out_p)
}

format_p <- function(p){
  p=round(p,3)
  if (p<0.001){ p = 'p < 0.001*'} else {
    if (p<0.05){ p = paste0('p = ',p)} else{
      p = 'N.S.'
    }
  }
  return(p)
}
geneate_p_ls <- function(in_df,data,group){
  gps = levels(data[[group]])
  n = length(gps)
  p.adj = round(in_df[['P.FDR']],3)
  max_value = max(data[[row.names(in_df)[1]]],na.rm = T)
  if (n==2) {
    p_values = list(
      xmin = c(gps[1]),
      xmax = c(gps[n]),
      annotation = c(paste("p =",p.adj)),
      y_position = c(1.05*max_value),
      color = 'black'
    )
  }
  else if ( p.adj >= 0.05){
    p_values = list(
      xmin = c(gps[1]),
      xmax = c(gps[n]),
      annotation = c(paste("All Groups p =",p.adj)),
      y_position = c(1.05*max_value),
      color = 'black',
      tip_length = 0
    )
    
  } else{
    # adding p values for multiple gps
    post_hoc_cols = colnames(in_df)[grep('vs',colnames(in_df))]
    post_hoc_cols_sp = strsplit(post_hoc_cols,' vs ')
    
    p_values = list(
      xmin = c(sapply(post_hoc_cols_sp, function(g){g[1]}),gps[1]),
      xmax = c(sapply(post_hoc_cols_sp, function(g){g[2]}),gps[n]),
      annotation = c(sapply(in_df[1, post_hoc_cols],format_p),paste("All Groups p =",p.adj)),
      y_position = sapply((0:n)+1, function(x){(1+x*.05)*max_value}),
      tip_length = c(rep(0.03,n),0),
      color = 'black'
    )
  }
  return(p_values)
}


# correlation plot with residual

generate_plot_data <- function(x,y,covariates,group,data){
  plot_data = na.omit(data[,c(group,covariates,x,y)])
  Y_lm = lm(eval(parse(text = paste0(y, " ~ ", paste(covariates,collapse ="+")))), plot_data)
  Y_resid<-resid(Y_lm) 
  X_lm = lm(eval(parse(text = paste0(x, " ~ ", paste(covariates,collapse ="+")))), plot_data)
  X_resid<-resid(X_lm) 
  out_df = data.frame(
    X_resid=X_resid,
    Y_resid = Y_resid)
  out_df = cbind.data.frame(out_df,plot_data[,c(group,covariates)])
  return(out_df)
}

Dw_cor_residual_whole_group <- function(x,y,covariates,group,data,ParCorr_Whole){
  plot_data = generate_plot_data(x,y,covariates,group,data)
  out_p = ggplot(plot_data,aes_string('X_resid','Y_resid'))+geom_point(aes_string(color=group))+
    geom_smooth(formula = y ~ x, method = "lm")+
    labs(x = paste0('residual of ',x),y = paste0('residual of ',y))+
    theme(legend.position=legend.position.cor)
  
  anno = ParCorr_Whole[[x]][c('r','p.print'),y]
  label = paste0('r = ',round(anno$r,3),', ',anno$p.print)
  annotation <- data.frame(
    x = c(max(plot_data$X_resid)),
    y = c(min(plot_data$Y_resid)),
    label = c(label)
  )
  out_p = out_p+ geom_text(data=annotation, aes( x=x, y=y, label=label),hjust='right',vjust ='bottom',check_overlap = TRUE)
  
  return(out_p)
  
}



generate_annotation_labels <- function(x,y,ParCorr_by_group,post_hoc_corr){
  post_hoc_cols = names(post_hoc_corr)
  gps = names(ParCorr_by_group)
  gp_labels = sapply(gps,function(gp){
    anno = ParCorr_by_group[[gp]][[x]][c('r','p.print'),y]
    paste0(gp,': r = ',round(anno$r,3),', ',anno$p.print)
  })
  post_hoc_labels = sapply(post_hoc_cols, function(post_hoc){
    p = round(post_hoc_corr[[post_hoc]][y,x][[1]],3)
    p = format_p(p)
    paste0(post_hoc,": ",p)
  })

  label = c(post_hoc_labels,gp_labels)
  return(label = label)

}

Dw_cor_residual_by_group <- function(x,y,covariates,group,data,annotation='none'){
  #pic_name = paste("CorPlot",paste0(x,'~',y),"png",sep = '.')
  #print(paste0("saving plot to ",pic_name))
  plot_data = generate_plot_data(x,y,covariates,group,data)
  X_resid = plot_data$X_resid
  Y_resid = plot_data$Y_resid
  out_p = ggplot(plot_data,aes_string(X_resid,Y_resid,group=group,color=group)) +
    geom_point()+geom_smooth(formula = y ~ x, method = "lm")+
    labs(x = paste0('residual of ',x),y = paste0('residual of ',y))+
    theme(legend.position=legend.position.cor)
  if (annotation != 'none'){
    n=length(annotation)
    anno_df <- data.frame(
      x = rep(max(X_resid),n),
      y = seq(1,by=-0.12,to = (1-0.12*(n-1)))*min(Y_resid),
      label = c(annotation)
    )
    out_p = out_p+ geom_text(data=anno_df, aes( x=x, y=y, label=label),hjust='right', inherit.aes = FALSE)
  }
  
  #ggsave(pic_name,out_p,width = pic_size_cor[1],height = pic_size_cor[2])
  out_p
}



### --------------------end of 2.4 plot function ------------------------

## ------------- end of 2.setting of functions ------------------------

