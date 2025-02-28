
###Functions to use with for plotting and modeling clusterlong

##Moved here to help keep the ClusterLong_slopeModelingAndPrecision_240606.Rmd more organized
# /path/to/data/ClusterLong/Scripts



###Ploting function to make group plots of interest

organize_LME<- function(brainVar,templateName,df,sessions = c(1:7)){ ##Given a brainVar name that is a column in plot data, this fits a basic lme and saves intercept and slope in a convenient format. df should be in long format with a column called region and measurements in column called value
  plotData.filtered <- df %>% filter(region == brainVar) %>% filter(template == templateName)  %>% filter(session %in% sessions)
  plotData.filtered<- plotData.filtered %>% pivot_wider(names_from = region, values_from = value)
  model_formula <- as.formula(paste(brainVar, "~ yearsFromBaseline + (yearsFromBaseline | ID)"))
  model.fullSamp.lmer <- lmer(model_formula, data = plotData.filtered)
  coef_CS1mm<- as_tibble(coef(model.fullSamp.lmer)$ID,rownames = "ID") %>% select(-yearsFromBaseline)
  random_effects_tidy <- tidy(model.fullSamp.lmer, effects = "ran_vals", conf.int = TRUE)
  fixed_Slope<-as.numeric(tidy(model.fullSamp.lmer,effects="fixed")[2,3])
  random_slopes_with_ci <- random_effects_tidy %>% filter(term == "yearsFromBaseline") %>%
    mutate(estimate = estimate + fixed_Slope, conf.low = conf.low + fixed_Slope, conf.high = conf.high + fixed_Slope) %>%
    select(-effect,-group,-term)
  colnames(random_slopes_with_ci)<-c("ID","slope","std.error","conf.low","conf.high")
  lmer.Ind.Data<- left_join(random_slopes_with_ci,coef_CS1mm)
  colnames(lmer.Ind.Data)<-c(colnames(random_slopes_with_ci),"intercept")
  lmer.Ind.Data<-lmer.Ind.Data %>% mutate(perAnnChange = slope/intercept*100,per.conf.low = conf.low/intercept*100,per.conf.high = conf.high/intercept*100)
  lmer.Ind.Data.plot<-left_join(lmer.Ind.Data,demo)
  return(lmer.Ind.Data.plot)
}

brainVar="Left-Hippocampus"
templateName="CSx6_inclusive"
#df<-plotData.6inclusive
sessions = c(1:6)

#### Need to check to see how much the more complicated LME model is influencing your results.... So fit a simple lm model in each ind separately and extract slope CI and intercept from those models to compare

organize_LM_simple<- function(brainVar,templateName=NULL,df,sessions = NULL,scans = NULL){ ##Given a brainVar name that is a column in plot data, this fits a basic lme and saves intercept and slope in a convenient format. df should be in long format with a column called region and measurements in column called value
  
  # If templateName is not provided, include all unique values from the template column
  if (is.null(templateName)) {
    templateName <- unique(df$template)
  }
  
  # If scans is not provided, include all unique values from the scanID column
  if (is.null(scans)) {
    scans <- unique(df$scanID)
  }
  
  # If sessopms is not provided, include all unique values from the sessions column
  if (is.null(sessions)) {
    scans <- unique(df$session)
  }
  
  plotData.filtered <- df %>% filter(region == brainVar) %>% filter(template %in% templateName)  %>% filter(session %in% sessions) %>% 
    filter(scanID %in% scans)
  plotData.filtered<- plotData.filtered %>% pivot_wider(names_from = region, values_from = value) %>% na.omit()
  model_formula <- as.formula(paste0("`",brainVar,"`", "~ yearsFromBaseline"))
  tidy.lm.full<-tibble()
  for(sub in unique(plotData.filtered$ID)){
    plotData.filtered.sub <- plotData.filtered %>% filter(ID == sub)
    model.lm <- lm(model_formula, data = plotData.filtered.sub)
    tidy.lm<-tidy(model.lm,conf.int = T) %>% filter(term == "yearsFromBaseline") %>% select(estimate,std.error,conf.low,conf.high)
    tidy.lm$intercept <- model.lm$coefficients[1]
    tidy.lm <- tidy.lm %>% rename(slope = estimate)
    tidy.lm$ID<-sub
    tidy.lm$template<-paste(templateName,collapse = ".")
    tidy.lm$nSesh<-length(unique(plotData.filtered.sub$session))
    tidy.lm$nScans<-length(plotData.filtered.sub$scanID)
    tidy.lm <- tidy.lm %>% mutate(perAnnChange = slope/intercept*100,per.conf.low = conf.low/intercept*100,per.conf.high = conf.high/intercept*100)
    tidy.lm.full <-bind_rows(tidy.lm.full,tidy.lm)
  }
  dem<-plotData.filtered %>% select(ID,Age,Sex,dx) %>% unique()
  tidy.lm.full<-left_join(tidy.lm.full,dem)
  return(tidy.lm.full)
}

brainVar="Left-Hippocampus"
templateName="CSx6_inclusive"
#df=plotData.ID
sessions = c(1,2,5,6)

organize_changeScore_simple<- function(brainVar,templateName=NULL,df,sessions = NULL,scans = NULL){ ##Given a brainVar name that is a column in plot data, this fits a basic lme and saves intercept and slope in a convenient format. df should be in long format with a column called region and measurements in column called value
  
  # If templateName is not provided, include all unique values from the template column
  if (is.null(templateName)) {
    templateName <- unique(df$template)
  }
  
  # If scans is not provided, include all unique values from the scanID column
  if (is.null(scans)) {
    scans <- unique(df$scanID)
  }
  
  # If sessopms is not provided, include all unique values from the sessions column
  if (is.null(sessions)) {
    scans <- unique(df$session)
  }
  
  plotData.filtered <- df %>% filter(region == brainVar) %>% filter(template %in% templateName)  %>% filter(session %in% sessions) %>% 
    filter(scanID %in% scans)
  plotData.filtered<- plotData.filtered %>% pivot_wider(names_from = region, values_from = value) %>% na.omit()
  model_formula <- as.formula(paste0("`",brainVar,"`", "~ yearsFromBaseline"))
  tidy.lm.full<-tibble()
  tidy.changeScore <-tibble()
  tidy.changeScore.full <-tibble()
  for(sub in unique(plotData.filtered$ID)){
    plotData.filtered.sub <- plotData.filtered %>% filter(ID == sub)
    plotData.filtered.sub <- plotData.filtered.sub %>% mutate(value = .[[brainVar]])
    
    model.lm <- lm(model_formula, data = plotData.filtered.sub)
    tidy.lm<-tidy(model.lm,conf.int = T) %>% filter(term == "yearsFromBaseline") %>% select(estimate,std.error,conf.low,conf.high)
    tidy.lm$intercept <- model.lm$coefficients[1]
    tidy.lm <- tidy.lm %>% rename(slope = estimate)
    tidy.lm$ID<-sub
    tidy.lm$template<-paste(templateName,collapse = ".")
    tidy.lm$nSesh<-length(unique(plotData.filtered.sub$session))
    tidy.lm$nScans<-length(plotData.filtered.sub$scanID)
    tidy.lm <- tidy.lm %>% mutate(perAnnChange = slope/intercept*100,per.conf.low = conf.low/intercept*100,per.conf.high = conf.high/intercept*100)
    tidy.lm.full <-bind_rows(tidy.lm.full,tidy.lm)
    
    ##Calculate the change score by averaging Visit 3 and Visit 1 and taking the difference divided by the mean to get percent change
    tmp.tidy.changeScore <- plotData.filtered.sub %>% group_by(Visit) %>% summarise(mean = mean(value),meanGap = mean(yearsFromBaseline))
    len<-nrow(tmp.tidy.changeScore)
    overallMean = (tmp.tidy.changeScore$mean[len] + tmp.tidy.changeScore$mean[1])/2
    yearsFromBaseline = (tmp.tidy.changeScore$meanGap[len] - tmp.tidy.changeScore$meanGap[1])
    tidy.changeScore <- tibble(slope = (tmp.tidy.changeScore$mean[len] - tmp.tidy.changeScore$mean[1])/yearsFromBaseline) ## Annualized Rate of Change
    tidy.changeScore$intercept <- overallMean
    tidy.changeScore$ID<-sub
    tidy.changeScore$template<-paste(templateName,collapse = ".")
    tidy.changeScore$nSesh<-length(unique(plotData.filtered.sub$session))
    tidy.changeScore$nScans<-length(plotData.filtered.sub$scanID)
    tidy.changeScore <- tidy.changeScore %>% mutate(perAnnChange = slope/intercept*100)
    tidy.changeScore.full <-bind_rows(tidy.changeScore.full,tidy.changeScore)

  }
  dem<-plotData.filtered %>% select(ID,Age,Sex,dx) %>% unique()
  tidy.lm.full<-left_join(tidy.lm.full,dem)
  tidy.changeScore.full<-left_join(tidy.changeScore.full,dem)
  return(tidy.changeScore.full)
}


templateName_test="CSx6_inclusive"
templateName_retest="CSx6_inclusive"
scans_test = NULL
scans_retest = NULL
sessions_test=c(1,3,5)
sessions_retest=c(2,4,6)

slopeErrors_LM_simple <- function(brainVar,df,templateName_test=NULL,templateName_retest=NULL,
                                  sessions_test=c(1,3,5),sessions_retest=c(2,4,6), scans_test = NULL,scans_retest = NULL,label = "leftBlank"){
  ### This is a wrapper around organize_LM_simple that allows you to calculate test retest error for different combinations of templates and/or scans. You can specificy both templates and scans or just templates or just scans to filter out the data you want to estimate slopes from.
  
  test<-organize_LM_simple(brainVar,df = df,templateName = templateName_test, sessions = sessions_test,scans = scans_test) %>% select(ID,Age,Sex,dx,slope,nSesh,nScans,perAnnChange)
  retest<-organize_LM_simple(brainVar,df = df, templateName = templateName_retest,sessions = sessions_retest,scans = scans_retest) %>% select(ID,Age,Sex,dx,slope,nSesh,nScans,perAnnChange)
  test$condition<-"test"
  retest$condition<-"retest"
  slope_testRetest_IDsum<-bind_rows(test,retest) %>% group_by(ID) %>% summarise(numSesh = mean(nSesh),numScans = mean(nScans))
  slope_testRetest<-bind_rows(test,retest) %>% select(-nSesh,-nScans) %>% pivot_wider(names_from = condition,values_from = c(slope,perAnnChange))
  slope_testRetest <- slope_testRetest %>% mutate(absError = abs(slope_test - slope_retest),absError.percent=abs(perAnnChange_test - perAnnChange_retest))
  slope_testRetest <- left_join(slope_testRetest,slope_testRetest_IDsum)
  slope_testRetest$label <- label
  slope_testRetest$brainVar <-brainVar
  slopes_TR <- slope_testRetest %>% select(slope_test,slope_retest)
  sum<-slope_testRetest %>% summarise(meanSlopeError = mean(absError),meanSlopeError.percent = mean(absError.percent),medianSlopeError = median(absError), sampSize=n(), meanNumSesh = mean(numSesh), 
                                      meanNumScans = mean(numScans),testRetest_cor = cor(slope_test,slope_retest),
                                      sdAbsError = sd(absError), seAbsError = sdAbsError / sqrt(sampSize), ciAbsError = qt(0.975, df = sampSize - 1) * seAbsError,
                                      sdAbsError.percent = sd(absError.percent), seAbsError.percent = sdAbsError.percent / sqrt(sampSize), ciAbsError.percent = qt(0.975, df = sampSize - 1) * seAbsError.percent)
  sum$label <- label
  sum$brainVar <- brainVar
  return(list(df = slope_testRetest,summary = sum))
}

###Median regression using quantile regression

organize_LM_median<- function(brainVar,templateName,df,sessions = c(1:7)){ ##Given a brainVar name that is a column in plot data, this fits a basic lme and saves intercept and slope in a convenient format. df should be in long format with a column called region and measurements in column called value
  plotData.filtered <- df %>% filter(region == brainVar) %>% filter(template == templateName)  %>% filter(session %in% sessions)
  plotData.filtered<- plotData.filtered %>% pivot_wider(names_from = region, values_from = value) %>% na.omit()
  model_formula <- as.formula(paste(brainVar, "~ yearsFromBaseline"))
  tidy.lm.full<-tibble()
  for(sub in unique(plotData.filtered$ID)){
    plotData.filtered.sub <- plotData.filtered %>% filter(ID == sub)
    model.lm <- rq(model_formula, data = plotData.filtered.sub)
    tidy.lm<-tidy(model.lm,conf.int = T) %>% filter(term == "yearsFromBaseline") %>% select(estimate,std.error,conf.low,conf.high)
    tidy.lm$intercept <- model.lm$coefficients[1]
    tidy.lm <- tidy.lm %>% rename(slope = estimate)
    tidy.lm$ID<-sub
    tidy.lm$template<-templateName
    tidy.lm$nSesh<-length(unique(plotData.filtered.sub$session))
    tidy.lm$nScans<-length(plotData.filtered.sub$scanID)
    tidy.lm <- tidy.lm %>% mutate(perAnnChange = slope/intercept*100,per.conf.low = conf.low/intercept*100,per.conf.high = conf.high/intercept*100)
    tidy.lm.full <-bind_rows(tidy.lm.full,tidy.lm)
  }
  dem<-plotData.filtered %>% select(ID,Age,Sex,dx) %>% unique()
  tidy.lm.full<-left_join(tidy.lm.full,dem)
  return(tidy.lm.full)
}

organize_bm_simple<- function(brainVar,templateName,df,sessions = c(1:7)){ ##Given a brainVar name that is a column in plot data, this fits a basic lme and saves intercept and slope in a convenient format. df should be in long format with a column called region and measurements in column called value
  plotData.filtered <- df %>% filter(region == brainVar) %>% filter(template == templateName)  %>% filter(session %in% sessions)
  plotData.filtered<- plotData.filtered %>% pivot_wider(names_from = region, values_from = value) %>% na.omit()
  model_formula <- as.formula(paste(brainVar, "~ yearsFromBaseline"))
  tidy.lm.full<-tibble()
  for(sub in unique(plotData.filtered$ID)){
    plotData.filtered.sub <- plotData.filtered %>% filter(ID == sub)
    model.bm <- stan_glm(model_formula , iter = 10000, data = plotData.filtered.sub, prior_intercept = student_t(df = 7))
    #+ session/SessionWithinVisit/OrderWithinSession
    #+ # Extract posterior samples for fixed effects
    
    
    posterior_samples <- as.data.frame(posterior_samples(model.bm))
    
    # Calculate 95% credible intervals for the slope (annual rate of change)
    annual_rate_change <- posterior_samples %>%
      select(contains("yearsFromBaseline")) %>%
      summarise(
        Estimate = mean(yearsFromBaseline),
        CI_Lower = quantile(yearsFromBaseline, 0.025),
        CI_Upper = quantile(yearsFromBaseline, 0.975)
      )
    
    tidy.lm<-tidy(model.lm,conf.int = T) %>% filter(term == "yearsFromBaseline") %>% select(estimate,std.error,conf.low,conf.high)
    tidy.lm$intercept <- model.lm$coefficients[1]
    tidy.lm <- tidy.lm %>% rename(slope = estimate)
    tidy.lm$ID<-sub
    tidy.lm$template<-templateName
    tidy.lm$nSesh<-length(unique(plotData.filtered.sub$session))
    tidy.lm$nScans<-length(plotData.filtered.sub$scanID)
    tidy.lm <- tidy.lm %>% mutate(perAnnChange = slope/intercept*100,per.conf.low = conf.low/intercept*100,per.conf.high = conf.high/intercept*100)
    tidy.lm.full <-bind_rows(tidy.lm.full,tidy.lm)
  }
  return(tidy.lm.full)
}

organize_LM_robust<- function(brainVar,templateName,df,sessions = c(1:7)){ ##Given a brainVar name that is a column in plot data, this fits a basic lme and saves intercept and slope in a convenient format. df should be in long format with a column called region and measurements in column called value
  plotData.filtered <- df %>% filter(region == brainVar) %>% filter(template == templateName)  %>% filter(session %in% sessions)
  plotData.filtered<- plotData.filtered %>% pivot_wider(names_from = region, values_from = value) %>% na.omit()
  model_formula <- as.formula(paste(brainVar, "~ yearsFromBaseline"))
  tidy.lm.full<-tibble()
  for(sub in unique(plotData.filtered$ID)){
    plotData.filtered.sub <- plotData.filtered %>% filter(ID == sub)
    model.lm <- rlm(model_formula, data = plotData.filtered.sub)
    tidy.lm<-tidy(model.lm,conf.int = T) %>% filter(term == "yearsFromBaseline") %>% select(estimate,std.error,conf.low,conf.high)
    tidy.lm$intercept <- model.lm$coefficients[1]
    tidy.lm <- tidy.lm %>% rename(slope = estimate)
    tidy.lm$ID<-sub
    tidy.lm$template<-templateName
    tidy.lm$nSesh<-length(unique(plotData.filtered.sub$session))
    tidy.lm$nScans<-length(plotData.filtered.sub$scanID)
    tidy.lm <- tidy.lm %>% mutate(perAnnChange = slope/intercept*100,per.conf.low = conf.low/intercept*100,per.conf.high = conf.high/intercept*100)
    tidy.lm.full <-bind_rows(tidy.lm.full,tidy.lm)
  }
  return(tidy.lm.full)
}

sessions<-c(1:6)
###Get order from HippChange
#hippOrderModel<- organize_LME(brainVar = brainVar,templateName ="CSx6_inclusive",plotData )
#lmer.Ind.Data.plot <-hippOrderModel %>% mutate(ID = fct_reorder(ID, perAnnChange,.desc = T))
#hippOrder<-levels(lmer.Ind.Data.plot$ID)
brainVar<-"WMH"
templateName="CSx6_inclusive"
#data=plotData.6sesh
title="Test"

generate_plot_ChangeSummary <- function(brainVar,templateName,data,sessions = c(1:7),title=""){
  
  ###
  plotData.filtered <- data %>% filter(region == brainVar) %>% filter(template == templateName) %>% filter(session %in% sessions)
  
  lm.Ind.Data.plot <-organize_LM_simple(brainVar,templateName,plotData.filtered,sessions = sessions,scans = NULL)
  ###Set yaxis limits based on all templates so that you can compare
  #lm.Ind.Data.plot.1 <-organize_LM_simple(brainVar,"CSx6_16",data,sessions = sessions,scans = NULL)
  #lm.Ind.Data.plot.2 <-organize_LM_simple(brainVar,"ADNI_2",data,sessions = sessions,scans = NULL)
  #lm.Ind.Data.plot.3 <-organize_LM_simple(brainVar,"CSx6_crossRes",data,sessions = sessions,scans = NULL)
  
  ###For now set limits solely on CSx6_inclusive because this is most important and comparisons will simply be partially cut off
  yMin=min(c(lm.Ind.Data.plot$conf.low))*1.03
  yMax=max(c(lm.Ind.Data.plot$conf.high))*1.05
  
  ##Order by Age for plotting
  #order<-intersect(hippOrder, lmer.Ind.Data.plot$ID)
  lm.Ind.Data.plot <-lm.Ind.Data.plot %>% mutate(ID = fct_reorder(ID, Age)) 
  #lmer.Ind.Data.plot <-lmer.Ind.Data.plot %>% mutate(ID = factor(ID, levels = order))
  
  lm.Ind.Data.plot <-lm.Ind.Data.plot %>% mutate(ID = fct_reorder(ID, Age),demID = paste0(Age,"-",Sex)) 
  lm.Ind.Data.plot <-lm.Ind.Data.plot %>% ungroup() %>% arrange(demID) %>% group_by(demID) %>% mutate(demIDNum=dense_rank(hemi))
  
  
  ####Make average decline plot by group
  
  lm.Ind.Data.plot.YA <- lm.Ind.Data.plot %>% filter(dx == "YA")
  lm.Ind.Data.plot.YA.label<-lm.Ind.Data.plot.YA %>% group_by(ID) %>% summarize(demID = first(demID), .groups = 'drop')
  p1<-ggplot(data = lm.Ind.Data.plot.YA,aes(y = slope, x = ID)) +
    geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
    geom_segment(
      x = -Inf, y = 0, xend = Inf, yend = 0,  # Extend line beyond x-axis limits
      color = "black", size = 1, linetype = "dashed"
    ) +
    geom_point(position = position_dodge(width = 0.25), size = 3, color = "#8cb369") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.25), width = 0.2, color = "#8cb369") +
    ylab("Annualized Change") +
    xlab("Young Adults") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    scale_x_discrete(labels = lm.Ind.Data.plot.YA.label$demID) + # Set x-axis labels
    theme_classic() +
    theme(text = element_text(family = "Helvetica",hjust = 0.5,size = 16),
                              axis.line = element_line(size = 1),
                              axis.ticks = element_line(colour = "black", size = 1),
                              axis.text =  element_text(color="black",size=14),
                              axis.title= element_text(size=18), # Remove the x-axis title
                              axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
                              axis.text.y = element_text(size = 14,color = "black"),
                              legend.text = element_text(size=14),
                              # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
                              legend.position = "none",
                              legend.justification=c(0,1),
                              legend.background = element_rect(fill = "transparent"), # Make legend background transparent
                              legend.key = element_rect(fill = "transparent", color = "transparent"),
                              legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  lm.Ind.Data.plot.CU <- lm.Ind.Data.plot %>% filter(dx == "CU")
  lm.Ind.Data.plot.CU.label<-lm.Ind.Data.plot.CU %>% group_by(ID) %>% summarize(demID = first(demID), .groups = 'drop')
  p2<-ggplot(data = lm.Ind.Data.plot.CU,aes(y = slope, x = ID)) +
    geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
    geom_segment(
      x = -Inf, y = 0, xend = Inf, yend = 0,  # Extend line beyond x-axis limits
      color = "black", size = 1, linetype = "dashed"
    ) +
    geom_point(position = position_dodge(width = 0.25), size = 3, color = "#8cb369") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.25), width = 0.2, color = "#8cb369") +
    ylab("Annualized Change") +
    xlab("Older Adults (CDR 0)") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    scale_x_discrete(labels = lm.Ind.Data.plot.CU.label$demID) + # Set x-axis labels
    theme(text = element_text(family = "Helvetica",hjust = 0.5,size = 16),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.text =  element_text(color="black",size=14),
          axis.title.x = element_text(size=18), # Remove the x-axis title
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
          axis.text.y = element_blank(),
          legend.text = element_text(size=14),
          # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
          legend.position = "none",
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  lm.Ind.Data.plot.AD <- lm.Ind.Data.plot %>% filter(dx == "AD")
  lm.Ind.Data.plot.AD.label<-lm.Ind.Data.plot.AD %>% group_by(ID) %>% summarize(demID = first(demID), .groups = 'drop')
  p3<-ggplot(data = lm.Ind.Data.plot.AD,aes(y = slope, x = ID)) +
    geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
    geom_segment(
      x = -Inf, y = 0, xend = Inf, yend = 0,  # Extend line beyond x-axis limits
      color = "black", size = 1, linetype = "dashed"
    ) +
    geom_point(position = position_dodge(width = 0.25), size = 3, color = "#8cb369") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.25), width = 0.2, color = "#8cb369") +
    ylab("Annualized Change") +
    xlab("MCI and AD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    scale_x_discrete(labels = lm.Ind.Data.plot.AD.label$demID) + # Set x-axis labels
    theme(text = element_text(family = "Helvetica",hjust = 0.5,size = 16),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.text =  element_text(color="black",size=14),
          axis.title.x = element_text(size=18), # Remove the x-axis title
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
          axis.text.y = element_blank(),
          legend.text = element_text(size=14),
          # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
          legend.position = "none",
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  lm.Ind.Data.plot.FTLD <- lm.Ind.Data.plot %>% filter(dx == "FTLD")
  lm.Ind.Data.plot.FTLD.label<-lm.Ind.Data.plot.FTLD %>% group_by(ID) %>% summarize(demID = first(demID), .groups = 'drop')
  p4<-ggplot(data = lm.Ind.Data.plot.FTLD,aes(y = slope, x = ID)) +
    geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
    geom_segment(
      x = -Inf, y = 0, xend = Inf, yend = 0,  # Extend line beyond x-axis limits
      color = "black", size = 1, linetype = "dashed"
    ) +
    geom_point(position = position_dodge(width = 0.25), size = 3, color = "#8cb369") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.25), width = 0.2, color = "#8cb369") +
    ylab("Annualized Change") +
    xlab("FTLD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    scale_x_discrete(labels = lm.Ind.Data.plot.FTLD$demID) + # Set x-axis labels
    theme(text = element_text(family = "Helvetica",hjust = 0.5,size = 16),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.text =  element_text(color="black",size=14),
          axis.title.x = element_text(size=18), # Remove the x-axis title
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
          axis.text.y = element_blank(),
          legend.text = element_text(size=14),
          # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
          legend.position = "none",
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  ### Save out combined figure
  relWidths=c(nrow(lm.Ind.Data.plot.YA)+4.5,nrow(lm.Ind.Data.plot.CU),nrow(lm.Ind.Data.plot.AD),nrow(lm.Ind.Data.plot.FTLD))
  figure <- ggarrange(p1, p2,p3,p4, ncol = 4, nrow = 1,widths = relWidths)
  
  title <- ggdraw() + draw_label(title, size = 20, fontface = 'bold', hjust = 0.5, x = 0.5, y = .8, vjust = 1)
  # Combine the title and the combined figure
  outFigure <- plot_grid(title, figure, ncol = 1, rel_heights = c(0.08, 1))
  
  return(outFigure)
  
  #ggsave("/path/to/data//Visualizations/Plots/groupSummary_perAnnChange/group.perAnnChange.Hippocampal.240529.png",figure,dpi = 500,width = 12,height = 6,units = "in")
  
}



generate_plot_percentChangeSummary <- function(brainVar,templateName,data,sessions = c(1:7),title=""){
  
  ###
  plotData.filtered <- data %>% filter(region == brainVar) %>% filter(template == templateName) %>% filter(session %in% sessions)
  
  lm.Ind.Data.plot <-organize_LM_simple(brainVar,templateName,plotData.filtered,sessions = sessions,scans = NULL)
  ###Set yaxis limits based on all templates so that you can compare
  #lm.Ind.Data.plot.1 <-organize_LM_simple(brainVar,"CSx6_16",data,sessions = sessions,scans = NULL)
  #lm.Ind.Data.plot.2 <-organize_LM_simple(brainVar,"ADNI_2",data,sessions = sessions,scans = NULL)
  #lm.Ind.Data.plot.3 <-organize_LM_simple(brainVar,"CSx6_crossRes",data,sessions = sessions,scans = NULL)
  
  ###For now set limits solely on CSx6_inclusive because this is most important and comparisons will simply be partially cut off
  yMin=min(c(lm.Ind.Data.plot$per.conf.low))*1.03
  yMax=max(c(lm.Ind.Data.plot$per.conf.high))*1.05
  
  ##Order by Age for plotting
  #order<-intersect(hippOrder, lmer.Ind.Data.plot$ID)
  lm.Ind.Data.plot <-lm.Ind.Data.plot %>% mutate(ID = fct_reorder(ID, Age)) 
  #lmer.Ind.Data.plot <-lmer.Ind.Data.plot %>% mutate(ID = factor(ID, levels = order))
  
  ####Make average decline plot by group
  
  lm.Ind.Data.plot.YA <- lm.Ind.Data.plot %>% filter(dx == "YA")
  p1<-ggplot(data = lm.Ind.Data.plot.YA,aes(y = perAnnChange, x = ID)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.25), size = 3, color = "#8cb369") +
    geom_errorbar(aes(ymin = per.conf.low, ymax = per.conf.high), position = position_dodge(width = 0.25), width = 0.2, color = "#8cb369") +
    ylab("% Annualized Change") +
    xlab("Young Adults") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
                            axis.title.x = element_text(size = 18),axis.title.y = element_text(size = 18),legend.position = "none") 
  
  lm.Ind.Data.plot.CU <- lm.Ind.Data.plot %>% filter(dx == "CU")
  p2<-ggplot(data = lm.Ind.Data.plot.CU,aes(y = perAnnChange, x = ID)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.25), size = 3, color = "#8cb369") +
    geom_errorbar(aes(ymin = per.conf.low, ymax = per.conf.high), position = position_dodge(width = 0.25), width = 0.2, color = "#8cb369") +
    ylab("% Annualized Change") +
    xlab("Older Adults (CDR 0)") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.title.y = element_blank(),axis.title.x = element_text(size = 18),legend.position = "none") 
  
  lm.Ind.Data.plot.AD <- lm.Ind.Data.plot %>% filter(dx == "AD")
  p3<-ggplot(data = lm.Ind.Data.plot.AD,aes(y = perAnnChange, x = ID)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.25), size = 3, color = "#8cb369") +
    geom_errorbar(aes(ymin = per.conf.low, ymax = per.conf.high), position = position_dodge(width = 0.25), width = 0.2, color = "#8cb369") +
    ylab("% Annualized Change") +
    xlab("MCI and AD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.title.y = element_blank(), axis.title.x = element_text(size = 18),legend.position = "none") 
  
  lm.Ind.Data.plot.FTLD <- lm.Ind.Data.plot %>% filter(dx == "FTLD")
  p4<-ggplot(data = lm.Ind.Data.plot.FTLD,aes(y = perAnnChange, x = ID)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.25), size = 3, color = "#8cb369") +
    geom_errorbar(aes(ymin = per.conf.low, ymax = per.conf.high), position = position_dodge(width = 0.25), width = 0.2, color = "#8cb369") +
    ylab("% Annualized Change") +
    xlab("FTLD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.title.y = element_blank(), axis.title.x = element_text(size = 18),legend.position = "none") 
  
  ### Save out combined figure
  relWidths=c(nrow(lm.Ind.Data.plot.YA),nrow(lm.Ind.Data.plot.CU),nrow(lm.Ind.Data.plot.AD),nrow(lm.Ind.Data.plot.FTLD))
  figure <- ggarrange(p1, p2,p3,p4, ncol = 4, nrow = 1,widths = relWidths)
  
  title <- ggdraw() + draw_label(title, size = 20, fontface = 'bold', hjust = 0.5, x = 0.5, y = .8, vjust = 1)
  # Combine the title and the combined figure
  outFigure <- plot_grid(title, figure, ncol = 1, rel_heights = c(0.08, 1))
  
  return(outFigure)
  
  #ggsave("/path/to/data//Visualizations/Plots/groupSummary_perAnnChange/group.perAnnChange.Hippocampal.240529.png",figure,dpi = 500,width = 12,height = 6,units = "in")
  
}


##################### Figure with Left and Right on the same plot
##Just use "Hippocampus" for brain Var

brainVar<-"Hippocampus"

generate_plot_Bilateral_ChangeSummary <- function(brainVar,templateName,data,sessions = c(1:6),title=""){
  
  
  # Define Unicode characters for half circles
  left_half_circle <- "\u25D6"
  right_half_circle <- "\u25D7"
  
  ###
  brainVar.left<- paste0("Left-",brainVar)
  brainVar.right <- paste0("Right-",brainVar)
  plotData.filtered.left <- data %>% filter(region == brainVar.left) %>% filter(template == templateName) %>% filter(session %in% sessions)
  plotData.filtered.right <- data %>% filter(region == brainVar.right) %>% filter(template == templateName) %>% filter(session %in% sessions)
  
  lm.Ind.Data.plot.left <-organize_LM_simple(brainVar.left,templateName,plotData.filtered.left,sessions = sessions)
  lm.Ind.Data.plot.left$hemi<-"Left"
  lm.Ind.Data.plot.right <-organize_LM_simple(brainVar.right,templateName,plotData.filtered.right,sessions = sessions)
  lm.Ind.Data.plot.right$hemi<-"Right"
  lm.Ind.Data.plot.comb<-bind_rows(lm.Ind.Data.plot.left,lm.Ind.Data.plot.right)
  yMin=min(c(lm.Ind.Data.plot.comb$conf.low))*1.03
  yMax=max(c(lm.Ind.Data.plot.comb$conf.high))*1.05
  
  
  ##Order by perAnnChange for plotting
  #order<-intersect(hippOrder, lm.Ind.Data.plot.comb$ID)
  lm.Ind.Data.plot.comb <-lm.Ind.Data.plot.comb %>% mutate(ID = fct_reorder(ID, Age),demID = paste0(Age,"-",Sex)) 
  lm.Ind.Data.plot.comb <-lm.Ind.Data.plot.comb %>% ungroup() %>% arrange(demID) %>% group_by(demID) %>% mutate(demIDNum=dense_rank(hemi))
  #lm.Ind.Data.plot <-lm.Ind.Data.plot %>% mutate(ID = factor(ID, levels = order))
  
  ####Make average decline plot by group
  
  #colorPal<-c("firebrick","black")
  custom_colors <- c("Left" = "#8cb369", "Right" = "#bc4b51")
  
  
  # Create a column for the Unicode symbols based on the Hemi column
  lm.Ind.Data.plot.comb <- lm.Ind.Data.plot.comb %>% mutate(shape = ifelse(hemi == "Left", left_half_circle, right_half_circle))
  
  lm.Ind.Data.plot.YA <- lm.Ind.Data.plot.comb %>% filter(dx == "YA")
  lm.Ind.Data.plot.YA.label<-lm.Ind.Data.plot.YA %>% group_by(ID) %>% summarize(demID = first(demID), .groups = 'drop')
  p1<-ggplot(data = lm.Ind.Data.plot.YA,aes(y = slope, x = ID,colour = hemi,label = shape)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_text(position = position_dodge(width = 0.15), size = 10, family = "Arial Unicode MS") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = custom_colors) + 
    ylab("Annualized Change") +
    xlab("Young Adults") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    scale_x_discrete(labels = lm.Ind.Data.plot.YA.label$demID) + # Set x-axis labels
    theme_classic() + theme(plot.title = element_text(hjust = 0.5,size = 16), 
                            axis.title.y = element_text(size = 18,color="black"), 
                            axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
                            axis.text.y = element_text(size = 14,color = "black"),
                            axis.title.x = element_text(size = 18),legend.position = "none") 
  
  lm.Ind.Data.plot.CU <- lm.Ind.Data.plot.comb %>% filter(dx == "CU")
  lm.Ind.Data.plot.CU.label<-lm.Ind.Data.plot.CU %>% group_by(ID) %>% summarize(demID = first(demID), .groups = 'drop')
  p2<-ggplot(data = lm.Ind.Data.plot.CU,aes(y = slope, x = ID,colour = hemi,label = shape)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_text(position = position_dodge(width = 0.15), size = 10, family = "Arial Unicode MS") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = custom_colors) + 
    ylab("Annualized Change (%)") +
    xlab("Older Adults (CDR 0)") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    scale_x_discrete(labels = lm.Ind.Data.plot.CU.label$demID) + # Set x-axis labels
    theme(plot.title = element_text(hjust = 0.5,size = 16), 
          axis.title.y = element_blank(),
                            axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 18),legend.position = "none") 
  
  lm.Ind.Data.plot.AD <- lm.Ind.Data.plot.comb %>% filter(dx == "AD")
  lm.Ind.Data.plot.AD.label<-lm.Ind.Data.plot.AD %>% group_by(ID) %>% summarize(demID = first(demID), .groups = 'drop')
  p3<-ggplot(data = lm.Ind.Data.plot.AD,aes(y = slope, x = demID,colour = hemi,label = shape)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_text(position = position_dodge(width = 0.15), size = 10, family = "Arial Unicode MS") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = custom_colors) + 
    ylab("Annualized Change (%)") +
    xlab("MCI and AD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    scale_x_discrete(labels = lm.Ind.Data.plot.AD.label$demID) + # Set x-axis labels
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5,size = 16), axis.title.y = element_blank(),
                            axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 18),legend.position = "none") 
  
  
  lm.Ind.Data.plot.FTLD <- lm.Ind.Data.plot.comb %>% filter(dx == "FTLD")
  lm.Ind.Data.plot.FTLD.label<-lm.Ind.Data.plot.FTLD %>% group_by(ID) %>% summarize(demID = first(demID), .groups = 'drop')
  p4<-ggplot(data = lm.Ind.Data.plot.FTLD,aes(y = slope, x = demID,colour = hemi,label = shape)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_text(position = position_dodge(width = 0.15), size = 10, family = "Arial Unicode MS") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = custom_colors) + 
    ylab("Annualized Change") +
    xlab("FTLD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    scale_x_discrete(labels = lm.Ind.Data.plot.FTLD.label$demID) + # Set x-axis labels
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5,size = 16),
          axis.title.y = element_blank(), 
          axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 18),legend.position = "none") 
  
  ### Save out combined figure
  relWidths=c(nrow(lm.Ind.Data.plot.YA)+4.5,nrow(lm.Ind.Data.plot.CU),nrow(lm.Ind.Data.plot.AD),nrow(lm.Ind.Data.plot.FTLD))
  figure <- ggarrange(p1, p2,p3,p4, ncol = 4, nrow = 1,widths = relWidths)
  
  title <- ggdraw() + draw_label(title, size = 20, fontface = 'bold', hjust = 0.5, x = 0.5, y = .8, vjust = 1)
  # Combine the title and the combined figure
  outFigure <- plot_grid(title, figure, ncol = 1, rel_heights = c(0.08, 1))
  
  return(outFigure)
  #ggsave("/path/to/data//Visualizations/Plots/groupSummary_perAnnChange/group.perAnnChange.Hippocampal.240529.png",figure,dpi = 500,width = 12,height = 6,units = "in")
  
}


generate_plot_Bilateral_percentChangeSummary <- function(brainVar,templateName,data,sessions = c(1:6),title=""){
  
  # Define Unicode characters for half circles
  left_half_circle <- "\u25D6"
  right_half_circle <- "\u25D7"
  
  ###
  brainVar.left<- paste0("Left-",brainVar)
  brainVar.right <- paste0("Right-",brainVar)
  plotData.filtered.left <- data %>% filter(region == brainVar.left) %>% filter(template == templateName) %>% filter(session %in% sessions)
  plotData.filtered.right <- data %>% filter(region == brainVar.right) %>% filter(template == templateName) %>% filter(session %in% sessions)
  
  lm.Ind.Data.plot.left <-organize_LM_simple(brainVar.left,templateName,plotData.filtered.left,sessions = sessions)
  lm.Ind.Data.plot.left$hemi<-"Left"
  lm.Ind.Data.plot.right <-organize_LM_simple(brainVar.right,templateName,plotData.filtered.right,sessions = sessions)
  lm.Ind.Data.plot.right$hemi<-"Right"
  lm.Ind.Data.plot.comb<-bind_rows(lm.Ind.Data.plot.left,lm.Ind.Data.plot.right)
  yMin=min(c(lm.Ind.Data.plot.comb$per.conf.low))*1.03
  yMax=max(c(lm.Ind.Data.plot.comb$per.conf.high))*1.05
  
  
  ##Order by perAnnChange for plotting
  #order<-intersect(hippOrder, lm.Ind.Data.plot.comb$ID)
  lm.Ind.Data.plot.comb <-lm.Ind.Data.plot.comb %>% mutate(ID = fct_reorder(ID, Age),demID = paste0(Age,"-",Sex))
  lm.Ind.Data.plot.comb <- lm.Ind.Data.plot.comb %>% mutate(Age_bin = cut(Age,
                         breaks = seq(0, 100, by = 5),  # Adjust the sequence as needed
                         right = F,
                         labels = paste(seq(0, 95, by = 5), seq(4, 99, by = 5), sep = "-"))) %>%
                        mutate(Age_bin = ifelse(Age_bin %in% c("15-19", "20-24"), "18-24", as.character(Age_bin)))
  lm.Ind.Data.plot.comb <-lm.Ind.Data.plot.comb %>% mutate(demID_AgeBin = paste0(Age_bin," ",Sex))
  
  lm.Ind.Data.plot.comb <-lm.Ind.Data.plot.comb %>% ungroup() %>% arrange(demID_AgeBin) %>% group_by(demID_AgeBin) %>% mutate(demID_AgeBinNum=dense_rank(hemi))
  #lm.Ind.Data.plot <-lm.Ind.Data.plot %>% mutate(ID = factor(ID, levels = order))
  
  ####Make average decline plot by group
  
  #colorPal<-c("firebrick","black")
  custom_colors <- c("Left" = "#8cb369", "Right" = "#bc4b51")
  
    
  # Create a column for the Unicode symbols based on the Hemi column
  lm.Ind.Data.plot.comb <- lm.Ind.Data.plot.comb %>% mutate(shape = ifelse(hemi == "Left", left_half_circle, right_half_circle))
  
  lm.Ind.Data.plot.YA <- lm.Ind.Data.plot.comb %>% filter(dx == "YA")
  lm.Ind.Data.plot.YA.label<-lm.Ind.Data.plot.YA %>% group_by(ID) %>% summarize(demID_AgeBin = first(demID_AgeBin), .groups = 'drop')
  p1<-ggplot(data = lm.Ind.Data.plot.YA,aes(y = perAnnChange, x = ID,colour = hemi,label = shape)) +
    geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
    geom_segment(
      x = -Inf, y = 0, xend = Inf, yend = 0,  # Extend line beyond x-axis limits
      color = "black", size = 1, linetype = "dashed"
    ) +
    geom_text(position = position_dodge(width = 0.15), size = 10, family = "Arial Unicode MS") +
    geom_errorbar(aes(ymin = per.conf.low, ymax = per.conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = custom_colors) + 
    ylab("Annualized Change (%)") +
    xlab("Adults") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    scale_x_discrete(labels = lm.Ind.Data.plot.YA.label$demID_AgeBin) + # Set x-axis labels
    theme_classic() + 
  theme(text = element_text(family = "Helvetica",hjust = 0.5,size = 16),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text =  element_text(color="black",size=14),
        axis.title= element_text(size=18), # Remove the x-axis title
        axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
        axis.text.y = element_text(size = 14,color = "black"),
        legend.text = element_text(size=14),
        # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
        legend.position = "none",
        legend.justification=c(0,1),
        legend.background = element_rect(fill = "transparent"), # Make legend background transparent
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  lm.Ind.Data.plot.CU <- lm.Ind.Data.plot.comb %>% filter(dx == "CU")
  lm.Ind.Data.plot.CU.label<-lm.Ind.Data.plot.CU %>% group_by(ID) %>% summarize(demID_AgeBin = first(demID_AgeBin), .groups = 'drop')
  p2<-ggplot(data = lm.Ind.Data.plot.CU,aes(y = perAnnChange, x = ID,colour = hemi,label = shape)) +
    geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
    geom_segment(
      x = -Inf, y = 0, xend = Inf, yend = 0,  # Extend line beyond x-axis limits
      color = "black", size = 1, linetype = "dashed"
    ) +
    geom_text(position = position_dodge(width = 0.15), size = 10, family = "Arial Unicode MS") +
    geom_errorbar(aes(ymin = per.conf.low, ymax = per.conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = custom_colors) + 
    ylab("Annualized Change (%)") +
    xlab("Older Adults (CDR 0)") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    scale_x_discrete(labels = lm.Ind.Data.plot.CU.label$demID_AgeBin) + # Set x-axis labels
    theme(text = element_text(family = "Helvetica",hjust = 0.5,size = 16),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.text =  element_text(color="black",size=14),
          axis.title.x = element_text(size=18), # Remove the x-axis title
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
          axis.text.y = element_blank(),
          legend.text = element_text(size=14),
          # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
          legend.position = "none",
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  lm.Ind.Data.plot.AD <- lm.Ind.Data.plot.comb %>% filter(dx == "AD")
  lm.Ind.Data.plot.AD.label<-lm.Ind.Data.plot.AD %>% group_by(ID) %>% summarize(demID_AgeBin = first(demID_AgeBin), .groups = 'drop')
  p3<-ggplot(data = lm.Ind.Data.plot.AD,aes(y = perAnnChange, x = ID,colour = hemi,label = shape)) +
    geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
    geom_segment(
      x = -Inf, y = 0, xend = Inf, yend = 0,  # Extend line beyond x-axis limits
      color = "black", size = 1, linetype = "dashed"
    ) +
    geom_text(position = position_dodge(width = 0.15), size = 10, family = "Arial Unicode MS") +
    geom_errorbar(aes(ymin = per.conf.low, ymax = per.conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = custom_colors) + 
    ylab("Annualized Change (%)") +
    xlab("MCI and AD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    scale_x_discrete(labels = lm.Ind.Data.plot.CU.label$demID_AgeBin) + # Set x-axis labels
    theme_classic() + 
    theme(text = element_text(family = "Helvetica",hjust = 0.5,size = 16),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.text =  element_text(color="black",size=14),
          axis.title.x = element_text(size=18), # Remove the x-axis title
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
          axis.text.y = element_blank(),
          legend.text = element_text(size=14),
          # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
          legend.position = "none",
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  lm.Ind.Data.plot.FTLD <- lm.Ind.Data.plot.comb %>% filter(dx == "FTLD")
  lm.Ind.Data.plot.FTLD.label<-lm.Ind.Data.plot.FTLD %>% group_by(ID) %>% summarize(demID_AgeBin = first(demID_AgeBin), .groups = 'drop')
  p4<-ggplot(data = lm.Ind.Data.plot.FTLD,aes(y = perAnnChange, x = demID_AgeBin,colour = hemi,label = shape)) +
    geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
    geom_segment(
      x = -Inf, y = 0, xend = Inf, yend = 0,  # Extend line beyond x-axis limits
      color = "black", size = 1, linetype = "dashed"
    ) +
    geom_text(position = position_dodge(width = 0.15), size = 10, family = "Arial Unicode MS") +
    geom_errorbar(aes(ymin = per.conf.low, ymax = per.conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = custom_colors) + 
    ylab("Annualized Change (%)") +
    xlab("FTLD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    scale_x_discrete(labels = lm.Ind.Data.plot.FTLD.label$demID_AgeBin) + # Set x-axis labels
    theme_classic() + 
    theme(text = element_text(family = "Helvetica",hjust = 0.5,size = 16),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text =  element_text(color="black",size=14),
        axis.title.x = element_text(size=18), # Remove the x-axis title
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14, angle = 40, hjust = 1, vjust = 1,color = "black"),
        axis.text.y = element_blank(),
        legend.text = element_text(size=14),
        # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
        legend.position = "none",
        legend.justification=c(0,1),
        legend.background = element_rect(fill = "transparent"), # Make legend background transparent
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  ### Save out combined figure
  relWidths=c(nrow(lm.Ind.Data.plot.YA)+3.5,nrow(lm.Ind.Data.plot.CU),nrow(lm.Ind.Data.plot.AD),nrow(lm.Ind.Data.plot.FTLD))
  figure <- ggarrange(p1, p2,p3,p4, ncol = 4, nrow = 1,widths = relWidths)
  
  title <- ggdraw() + draw_label(title, size = 20, fontface = 'bold', hjust = 0.5, x = 0.5, y = .8, vjust = 1)
  # Combine the title and the combined figure
  outFigure <- plot_grid(title, figure, ncol = 1, rel_heights = c(0.08, 1))
  
  return(outFigure)
  
  #ggsave("/path/to/data//Visualizations/Plots/groupSummary_perAnnChange/group.perAnnChange.Hippocampal.240529.png",figure,dpi = 500,width = 12,height = 6,units = "in")
  
}


##################### Figure with ADNI and CS on same plot

brainVar="Left-Hippocampus"

generate_plot_ChangeSummaryWithADNI <- function(brainVar,templateName,data,sessions = c(1:6),title=""){
  
  ###
  plotData.filtered <- data %>% filter(region == brainVar) %>% filter(template == templateName) %>% filter(session %in% sessions)
  
  lm.Ind.Data.plot <-organize_LM_simple(brainVar,templateName,plotData.filtered,sessions = sessions)
  ###Set yaxis limits based on all templates so that you can compare
  lm.Ind.Data.plot.1 <-organize_LM_simple(brainVar,"CSx6_inclusive",data,sessions = sessions)
  lm.Ind.Data.plot.2 <-organize_LM_simple(brainVar,"ADNI_inclusive",data,sessions = sessions)
  #lm.Ind.Data.plot.3 <-organize_LM_simple(brainVar,"CSx6_crossRes",data,sessions = sessions)
  
  ###For now set limits solely on CSx6_inclusive because this is most important and comparisons will simply be partially cut off
  yMin=min(c(lm.Ind.Data.plot.1$conf.low*1.2))
  yMax=max(c(lm.Ind.Data.plot.1$conf.high*1.2))
  
  
  
  ##Combine ADNI And CS for plotting
  lm.Ind.Data.plot$scanType = "CS 1mm"
  lm.Ind.Data.plot.2$scanType = "ADNI"
  lm.Ind.Data.plot.comb <- full_join(lm.Ind.Data.plot,lm.Ind.Data.plot.2)
  
  ##Order by perAnnChange for plotting
  #order<-intersect(hippOrder, lm.Ind.Data.plot.comb$ID)
  lm.Ind.Data.plot.comb <-lm.Ind.Data.plot.comb %>% mutate(ID = fct_reorder(ID, Age)) 
  
  #lm.Ind.Data.plot <-lm.Ind.Data.plot %>% mutate(ID = factor(ID, levels = order))
  
  ####Make average decline plot by group
  
  colorPal<-c("firebrick","black")
  
  lm.Ind.Data.plot.YA <- lm.Ind.Data.plot.comb %>% filter(dx == "YA")
  p1<-ggplot(data = lm.Ind.Data.plot.YA,aes(y = slope, x = ID,colour = scanType)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.15), size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = colorPal) + 
    ylab("Annualized Change") +
    xlab("Young Adults") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
                            axis.title.x = element_blank(),legend.position = "none")
  
  lm.Ind.Data.plot.CU <- lm.Ind.Data.plot.comb %>% filter(dx == "CU")
  p2<-ggplot(data = lm.Ind.Data.plot.CU,aes(y = slope, x = ID,colour = scanType)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.15), size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = colorPal) + 
    ylab("Annualized Change") +
    xlab("Older Adults (CDR 0)") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.title.y = element_blank(),axis.title.x = element_blank(),legend.position = "none")
  
  lm.Ind.Data.plot.AD <- lm.Ind.Data.plot.comb %>% filter(dx == "AD")
  p3<-ggplot(data = lm.Ind.Data.plot.AD,aes(y = slope, x = ID,colour = scanType)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.15), size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = colorPal) + 
    ylab("Annualized Change") +
    xlab("MCI and AD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.title.y = element_blank(),axis.title.x = element_blank(),legend.position = "none")
  
  lm.Ind.Data.plot.FTLD <- lm.Ind.Data.plot.comb %>% filter(dx == "FTLD")
  p4<-ggplot(data = lm.Ind.Data.plot.FTLD,aes(y = slope, x = ID,colour = scanType)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.15), size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.15), width = 0.2) +
    scale_color_manual(values = colorPal) + 
    ylab("Annualized Change") +
    xlab("FTLD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.title.y = element_blank(),axis.title.x = element_blank(),legend.position = "none")
  
  ### Save out combined figure
  relWidths=c(nrow(lm.Ind.Data.plot.YA),nrow(lm.Ind.Data.plot.CU),nrow(lm.Ind.Data.plot.AD),nrow(lm.Ind.Data.plot.FTLD))
  figure <- ggarrange(p1, p2,p3,p4, ncol = 4, nrow = 1,widths = relWidths)
  
  title <- ggdraw() + draw_label(title, size = 14, fontface = 'bold', hjust = 0.5, x = 0.5, y = .8, vjust = 1)
  # Combine the title and the combined figure
  outFigure <- plot_grid(title, figure, ncol = 1, rel_heights = c(0.08, 1))
  
  return(outFigure)
  
  #ggsave("/path/to/data//Visualizations/Plots/groupSummary_perAnnChange/group.perAnnChange.Hippocampal.240529.png",figure,dpi = 500,width = 12,height = 6,units = "in")
  
}


brainVar="Left-Hippocampus"
templateName1<-"CSx6_8_test"
templateName2<-"CSx6_8_retest"
#data<-plotData.all6sesh

###Probably use this for test retest, but made so you can use any two tamples

generate_plot_ChangeSummaryComparison <- function(brainVar,templateName1,templateName2,data,sessions = c(1:6),title=""){
  
  ###
  plotData.filtered <- data %>% filter(region == brainVar) %>% filter(template == templateName) %>% filter(session %in% sessions)
  
  lm.Ind.Data.plot.1 <-organize_LM_simple(brainVar,templateName1,data,sessions = sessions)
  lm.Ind.Data.plot.2 <-organize_LM_simple(brainVar,templateName2,data,sessions = sessions)
  
  ###set axis limits
  yMin=min(c(lm.Ind.Data.plot.1$conf.low,lm.Ind.Data.plot.2$conf.low))
  yMax=max(c(lm.Ind.Data.plot.1$conf.high,lm.Ind.Data.plot.2$conf.high))
  
  
  
  ##Combine ADNI And CS for plotting
  lm.Ind.Data.plot.1$scanType = templateName1
  lm.Ind.Data.plot.2$scanType = templateName2
  lm.Ind.Data.plot.comb <- full_join(lm.Ind.Data.plot.1,lm.Ind.Data.plot.2)
  
  ##Order by perAnnChange for plotting
  #order<-intersect(hippOrder, lm.Ind.Data.plot.comb$ID)
  lm.Ind.Data.plot.comb <-lm.Ind.Data.plot.comb %>% mutate(ID = fct_reorder(ID, Age)) 
  
  #lm.Ind.Data.plot <-lm.Ind.Data.plot %>% mutate(ID = factor(ID, levels = order))
  
  ####Make average decline plot by group
  
  colorPal<-c("steelblue4","olivedrab4")
  
  lm.Ind.Data.plot.YA <- lm.Ind.Data.plot.comb %>% filter(dx == "YA")
  p1<-ggplot(data = lm.Ind.Data.plot.YA,aes(y = slope, x = ID,colour = scanType)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.35), size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.35), width = 0.2) +
    scale_color_manual(values = colorPal) + 
    ylab("Annualized Change") +
    xlab("Young Adults") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
                            axis.title.x = element_blank(),legend.position = "none")
  
  lm.Ind.Data.plot.CU <- lm.Ind.Data.plot.comb %>% filter(dx == "CU")
  p2<-ggplot(data = lm.Ind.Data.plot.CU,aes(y = slope, x = ID,colour = scanType)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.35), size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.35), width = 0.2) +
    scale_color_manual(values = colorPal) + 
    ylab("Annualized Change") +
    xlab("Older Adults (CDR 0)") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.title.y = element_blank(),axis.title.x = element_blank(),legend.position = "none")
  
  lm.Ind.Data.plot.AD <- lm.Ind.Data.plot.comb %>% filter(dx == "AD")
  p3<-ggplot(data = lm.Ind.Data.plot.AD,aes(y = slope, x = ID,colour = scanType)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.35), size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.35), width = 0.2) +
    scale_color_manual(values = colorPal) + 
    ylab("Annualized Change") +
    xlab("MCI and AD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.title.y = element_blank(),axis.title.x = element_blank(),legend.position = "none")
  
  lm.Ind.Data.plot.FTLD <- lm.Ind.Data.plot.comb %>% filter(dx == "FTLD")
  p4<-ggplot(data = lm.Ind.Data.plot.FTLD,aes(y = slope, x = ID,colour = scanType)) +
    geom_hline(yintercept = 0, color = "slategray", size = 3, alpha = 0.3) +
    geom_point(position = position_dodge(width = 0.35), size = 3) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.35), width = 0.2) +
    scale_color_manual(values = colorPal) + 
    ylab("Annualized Change") +
    xlab("FTLD") +
    #scale_y_continuous(limits = c(yMin,yMax), expand = c(0, 0)) +
    coord_cartesian(ylim = c(yMin, yMax)) + 
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
          axis.title.y = element_blank(),axis.title.x = element_blank(),legend.position = "none")
  
  ### Save out combined figure
  relWidths=c(nrow(lm.Ind.Data.plot.YA),nrow(lm.Ind.Data.plot.CU),nrow(lm.Ind.Data.plot.AD),nrow(lm.Ind.Data.plot.FTLD))
  figure <- ggarrange(p1, p2,p3,p4, ncol = 4, nrow = 1,widths = relWidths)
  
  title <- ggdraw() + draw_label(title, size = 14, fontface = 'bold', hjust = 0.5, x = 0.5, y = .8, vjust = 1)
  # Combine the title and the combined figure
  outFigure <- plot_grid(title, figure, ncol = 1, rel_heights = c(0.08, 1))
  
  return(outFigure)
  
  #ggsave("/path/to/data//Visualizations/Plots/groupSummary_perAnnChange/group.perAnnChange.Hippocampal.240529.png",figure,dpi = 500,width = 12,height = 6,units = "in")
  
}

#######################################
##### Dashboard Style plotting functions
######################################



### Adapted from "/path/to/data//ClusterLong/Scripts/clusterIndDashboards_Functions_240613.Rmd"

id=""
brainVar="ADsignature_CT"
type="Thick"
templateName="CSx6_inclusive"
title=""
yaxisMargin=0
#df=plotData.6sesh.only6sessions


generate_plot_testRetest <- function(df,id, brainVar,type,templateName,title) { #type has to be either vol, Thick or GWcon or Area # templateName is either ADNI, crossRes or inclusive
  
  plotData.ID.allTemplate.left <- df %>% filter(ID == id, grepl(brainVar,region), origData == type) %>% filter(grepl("Left-",region))
  plotData.ID.allTemplate.right <- df %>% filter(ID == id, grepl(brainVar,region), origData == type) %>% filter(grepl("Right-",region))
  
  plotData <- df %>% filter(ID == id, grepl(brainVar,region), origData == type) 
  if(templateName == "inclusive"){
    plotData <- plotData %>% filter(grepl("CSx6",template)) %>% filter(!grepl("crossRes",template))
  }else{
    plotData <- plotData %>% filter(grepl(templateName,template))
  }
  plotData <- plotData %>% mutate(visit = as_factor(ceiling(session/2)))
  plotData <- plotData %>% mutate(visit = as_factor(if_else(visit == 1, "Baseline", if_else(visit == 2, "Visit 2",if_else(visit == 3, "Visit 3","Visit 4")))))
  
  ###Estimate time from baseline for each scan
  plotData$Date <- ymd(plotData$Date)
  baselineDate<-min(plotData$Date)
  plotData<-plotData %>% mutate(monthsFromBaseline = as.numeric(Date - baselineDate)/(365.25/12))
  
  
  #Calculate average number of months from baseline for each visit
  plotData <- plotData %>% group_by(visit) %>% mutate(visitMonths = round(mean(monthsFromBaseline),2) )
  allPlotData <- plotData %>% group_by(ID) %>% mutate(visit = as_factor(ceiling(session/2))) %>% ungroup()
  
  allPlotData.left <- allPlotData %>% filter(grepl("Left-",region))
  allPlotData.right <- allPlotData %>% filter(grepl("Right-",region))
  
  ###Separate ID data
  plotData.ID <- allPlotData %>% filter(ID == id)
  
  ###Grab demo for plot
  dx=plotData.ID$dx[1]
  id=plotData.ID$ID[1]
  age=plotData.ID$Age[1]
  sex=plotData.ID$Sex[1]
  
  ###Split into hemisperes and handle % change
  plotData.ID.left <- plotData.ID %>% filter(grepl("Left-",region))
  plotData.ID.right <- plotData.ID %>% filter(grepl("Right-",region))
  
  #Adjust x-axis to be bigger only if 7 sessions
  plotData.allIDs <- df %>% filter(grepl(brainVar,region), origData == type) 
  if(templateName == "inclusive"){
    plotData.allIDs <- plotData.allIDs %>% filter(grepl("CSx6",template)) %>% filter(!grepl("crossRes",template))
  }else{
    plotData.allIDs <- plotData.allIDs %>% filter(grepl(templateName,template))
  }
  
  #xmax<-(max(allPlotData.left$monthsFromBaseline) + 1) ## Use this if you want a different x-axis for each ID
  xmax<-(max(plotData.allIDs$yearsFromBaseline*12) +1) ## Use this if you want the same x-axis for each ID
  xmin <- 15-xmax # so that the spacing between y-axis and 0 is the same as the end of the x-axis and last tick spacing
  ###Simplify to a linear regression in each person - used to have lme model but the complexity caused some issues...
  
  model_lm.left <- lm(value ~ monthsFromBaseline, data = plotData.ID.left)
  lineMax<-xmax
  new_data.ID.left <- data.frame(ID = id, monthsFromBaseline = seq(0, lineMax, length.out = 100))
  pred_intervals.left <- predict(model_lm.left, newdata = new_data.ID.left)# , level = 0.80, n.sims = 500, seed = 123, include.resid.var = F,type = "linear.prediction")
  new_data.ID.left<-bind_cols(new_data.ID.left,pred_intervals.left)
  new_data.ID.left <- new_data.ID.left %>% rename(value = ...3) ##To make it plot right
  ##For the plot estimate the annual %change from the lme regression intercept and slope
  model.ID.left<-coef(model_lm.left)
  annChangePer.left <- round((model.ID.left[2]*12)/model.ID.left[1]*100,2)
  CI.model.ID.left<-confint(model_lm.left)
  CI.annChangePer.left <- round((CI.model.ID.left[2,]*12)/model.ID.left[1]*100,2)
  
  # Assuming new_data_specific_id contains 'monthsFromBaseline' from 0 to 15 and 'predicted_value'
  # Calculate slope (m) using the first and last fit values
  x1 <- 0  # Starting point of x
  x2 <- xmax  # Ending point of x
  
  # Ensure new_data_specific_id is ordered by monthsFromBaseline if it's not already
  new_data.ID.left <- new_data.ID.left[order(new_data.ID.left$monthsFromBaseline), ]
  y1 <- predict(model_lm.left, newdata = data.frame(monthsFromBaseline = x1, ID = id), re.form = NULL)
  y2 <- predict(model_lm.left, newdata = data.frame(monthsFromBaseline = x2), re.form = NULL)
  # Define the line endpoints based on the calculated slope
  line_data.left <- data.frame(monthsFromBaseline = c(x1, x2), value = c(y1, y2))
  
  
  model_lm.right <- lm(value ~ monthsFromBaseline, data = plotData.ID.right)
  new_data.ID.right <- data.frame(ID = id, monthsFromBaseline = seq(0, lineMax, length.out = 100))
  pred_intervals.right <- predict(model_lm.right, newdata = new_data.ID.right)# , level = 0.80, n.sims = 500, seed = 123, include.resid.var = F,type = "linear.prediction")
  new_data.ID.right<-bind_cols(new_data.ID.right,pred_intervals.right)
  new_data.ID.right <- new_data.ID.right %>% rename(value = ...3) ##To make it plot right
  ##For the plot estimate the annual %change from the lme regression intercept and slope
  model.ID.right<-coef(model_lm.right)
  annChangePer.right <- round((model.ID.right[2]*12)/model.ID.right[1]*100,2)
  CI.model.ID.right<-confint(model_lm.right)
  CI.annChangePer.right <- round((CI.model.ID.right[2,]*12)/model.ID.right[1]*100,2)
  
  # Ensure new_data_specific_id is ordered by monthsFromBaseline if it's not already
  new_data.ID.right <- new_data.ID.right[order(new_data.ID.right$monthsFromBaseline), ]
  y1 <- predict(model_lm.right, newdata = data.frame(monthsFromBaseline = x1, ID = id), re.form = NULL)
  y2 <- predict(model_lm.right, newdata = data.frame(monthsFromBaseline = x2), re.form = NULL)
  # Define the line endpoints based on the calculated slope
  line_data.right <- data.frame(monthsFromBaseline = c(x1, x2), value = c(y1, y2))
  
  ##Left
  mean.base.plotData.ID.left <- plotData.ID.left %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(mean=mean(value))
  median.base.plotData.ID.left <- plotData.ID.left %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(median=median(value))
  median.test.plotData.ID.left <- plotData.ID.left %>% ungroup() %>% filter(session == 1) %>% summarise(median=median(value))
  median.retest.plotData.ID.left <- plotData.ID.left %>% ungroup() %>% filter(session == 2) %>% summarise(median=median(value))
  
  plotData.ID.left <- plotData.ID.left %>% mutate(perDiff=(value-median.base.plotData.ID.left$median)/median.base.plotData.ID.left$median*100)
  plotData.ID.left <- plotData.ID.left %>% mutate(perDiffperYear = perDiff/(monthsFromBaseline/12))
  medians.left <- plotData.ID.left %>% group_by(monthsFromBaseline,session) %>% summarize(median = median(value),monthsFromBaseline = mean(monthsFromBaseline),se = sd(value, na.rm = TRUE)/sqrt(n()),n=n(), vis = median(as.numeric(visit)))
  medians.left <- medians.left %>% mutate(CI95 = 1.96*se)
  ###Estimate of the benefits of reducing error with many samples
  sum.sdMed.left<- medians.left %>% group_by(vis) %>% summarise(sdMed = sd(median))
  meanMedSD.left <- mean(sum.sdMed.left$sdMed)
  
  visit_Avgs_long.left <- plotData.ID.left %>% group_by(visit) %>% summarise(perDiff2 = mean(perDiff),se=sd(perDiff)/sqrt(n()))
  visit_Avgs_long.left <- visit_Avgs_long.left %>% mutate(CI95 = 1.96*se)
  
  df.summary.left <- plotData.ID.left %>% group_by(visit) %>% summarise(median = round(median(value),0),sd = sd(value),se = sd(value, na.rm = TRUE)/sqrt(n()), perDiff = round(median(perDiff),2), perDiffAnn = round(median(perDiffperYear),2))
  df.summary.left <- df.summary.left %>% mutate(CI95 = 1.96*se)
  meanSD.left<-mean(df.summary.left$sd)
  
  ##Right
  mean.base.plotData.ID.right <- plotData.ID.right %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(mean=mean(value))
  median.base.plotData.ID.right <- plotData.ID.right %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(median=median(value))
  median.test.plotData.ID.right <- plotData.ID.right %>% ungroup() %>% filter(session == 1) %>% summarise(median=median(value))
  median.retest.plotData.ID.right <- plotData.ID.right %>% ungroup() %>% filter(session == 2) %>% summarise(median=median(value))
  
  plotData.ID.right <- plotData.ID.right %>% mutate(perDiff=(value-median.base.plotData.ID.right$median)/median.base.plotData.ID.right$median*100)
  plotData.ID.right <- plotData.ID.right %>% mutate(perDiffperYear = perDiff/(monthsFromBaseline/12))
  medians.right <- plotData.ID.right %>% group_by(monthsFromBaseline,session) %>% summarize(median = median(value),monthsFromBaseline = mean(monthsFromBaseline),se = sd(value, na.rm = TRUE)/sqrt(n()), vis = median(as.numeric(visit)))
  medians.right <- medians.right %>% mutate(CI95 = 1.96*se)
  sum.sdMed.right<- medians.right %>% group_by(vis) %>% summarise(sdMed = sd(median))
  meanMedSD.right <- mean(sum.sdMed.right$sdMed)
  
  plotData.ID.right.test <- plotData.ID.right %>% filter(session %in% c(1,3,5)) %>% mutate(perDiff=(value - median.test.plotData.ID.right$median) / median.test.plotData.ID.right$median * 100)
  plotData.ID.right.test <- plotData.ID.right.test %>% mutate(perDiffperYear = perDiff/(monthsFromBaseline/12))
  medians.right.test <- plotData.ID.right.test %>% group_by(monthsFromBaseline,session) %>% summarize(value = median(value))
  plotData.ID.right.retest <- plotData.ID.right %>% filter(session %in% c(2,4,6)) %>% mutate(perDiff=(value - median.retest.plotData.ID.right$median) / median.retest.plotData.ID.right$median * 100)
  plotData.ID.right.retest <- plotData.ID.right.retest %>% mutate(perDiffperYear = perDiff/(monthsFromBaseline/12))
  medians.right.retest <- plotData.ID.right.retest %>% group_by(monthsFromBaseline,session) %>% summarize(value = median(value))
  
  plotData.ID.right.testRetest<-bind_rows(plotData.ID.right.test,plotData.ID.right.retest)
  
  #session_averages_long.right <- plotData.ID.right.testRetest %>% group_by(session) %>%
  # summarise(avg_perDiff = mean(perDiff), sd_perDiff = sd(perDiff), se = sd(perDiff) / sqrt(n()), med_perDiff = median(perDiff)) %>% 
  #   mutate(visit = rep(1:(n() %/% 2), each = 2))
  
  # visit_Avgs.right <- session_averages_long.right %>% group_by(visit) %>% summarise(visitAvg = mean(avg_perDiff),se=mean(se))
  # visit_Avgs_long.right <- plotData.ID.right %>% group_by(visit) %>% summarise(perDiff2 = mean(perDiff),se=sd(perDiff)/sqrt(n()))
  # visit_Avgs_long.right <- visit_Avgs_long.right %>% mutate(CI95 = 1.96*se)
  
  df.summary.right <- plotData.ID.right %>% group_by(visit) %>% summarise(median = round(median(value),0),se = sd(value, na.rm = TRUE)/sqrt(n()), perDiff = round(median(perDiff),2), perDiffAnn = round(median(perDiffperYear),2),sd=sd(value))
  df.summary.right <- df.summary.right %>% mutate(CI95 = 1.96*se)
  meanSD.right<-mean(df.summary.right$sd)
  
  medians.right <- medians.right %>% rename(value = median) %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  medians.left <- medians.left %>% rename(value = median) %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  ###Make Plots
  
  plotData.ID.left<- plotData.ID.left %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  plotData.ID.right<- plotData.ID.right %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  
  # Define the colors for odd and even visits
  colors <- c("Even" = "steelblue4", "Odd" = "red4")
  
  #colors <- c("red4", "steelblue4","red4", "steelblue4","red4","steelblue4","red4")
  
  ### Make range so that it's consistent except in extreme cases
  ### Make range so that it is the same size for left and right 
  
  ## use the larger size of the range as the range size for both left and right plots
  range.left <- range(plotData.ID.left$value)[2] - range(plotData.ID.left$value)[1] 
  range.right <- range(plotData.ID.right$value)[2] - range(plotData.ID.right$value)[1]
  rangeSize<-max(340,range.left,range.right)*1.2 ## Range should be at least 300 so that there aren't severe differences between participants range
  
  ##Need to make sure there is enough space for a label in the bottom left corner of the plot.
  plotData.ID.left.12 <- plotData.ID.left %>% filter(session %in% c(1,2))
  
  minRange.left<-min(c(min(plotData.ID.left$value)*.99,mean(plotData.ID.left$value)*.97,min(plotData.ID.left.12$value)-100))
  maxRange.left<-minRange.left + rangeSize
  
  plot_title.left<-paste0("Left ",brainVar,"\n", annChangePer.left,"% [",CI.annChangePer.left[1],",  ",CI.annChangePer.left[2],"]")
  
  
  # Define a function to compute tick marks that are multiples of 50
  compute_fixed_breaks <- function(min_val, max_val, step = 50) {
    start_val <- ceiling(min_val / step) * step
    end_val <- floor(max_val / step) * step
    return(seq(start_val, end_val, by = step))
  }
  
  # Calculate y-axis breaks
  breaks_y_left <- compute_fixed_breaks(minRange.left-50, maxRange.left)
  
  print(plot_title.left)
  
  p1<-ggplot(plotData.ID.left, aes(monthsFromBaseline, value, color = session_type)) + theme_classic() + scale_color_manual(values = colors) + 
    guides(color = FALSE) + 
    geom_line(data = line_data.left, aes(x = monthsFromBaseline, y = value), color = "black",linetype = "dashed",size = 1) +
    geom_pointrange(aes(ymin = value-CI95, ymax = value+CI95),alpha = 1,data = medians.left,shape = 5,size = .7,linewidth = .65,stroke = 1.4) +
    geom_jitter(position = position_jitter(0.005),alpha = .45,size = 2) + 
    #ylab(expression(Hippocampal~Volume~(mm^3))) + 
    xlab("Months") +
    ylab(NULL) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_y_continuous(limits = c(min(breaks_y_left+25), max(breaks_y_left)+1), breaks = breaks_y_left, expand = c(0, 0)) +
    scale_x_continuous(limits = c(xmin,xmax), breaks = seq(0,21,3), expand = c(0, 0)) +
    theme(text = element_text(family = "Helvetica"),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.text =  element_text(color="black",size=14, margin = margin(b = 0)),
          axis.title.x = element_text(size=16), # Remove the x-axis title
          axis.title.y = element_text(size=16),
          legend.text = element_text(size=14),
          # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
          legend.position=c(0.01,1), 
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  ##Need to make sure there is enough space for a label in the bottom left corner of the plot.
  plotData.ID.right.12 <- plotData.ID.right %>% filter(session %in% c(1,2))
  minRange.right<-min(c(min(plotData.ID.right$value)*.99,mean(plotData.ID.right$value)*.97,min(plotData.ID.right.12$value)-100))
  maxRange.right<-minRange.right + rangeSize
  
  breaks_y_right <- compute_fixed_breaks(minRange.right-50, maxRange.right)
  plot_title.right<-paste0("Right ",brainVar,"\n", annChangePer.right,"% [",CI.annChangePer.right[1],",  ",CI.annChangePer.right[2],"]")
  print(plot_title.right)
  p2<-ggplot(plotData.ID.right, aes(monthsFromBaseline, value, color = session_type)) + theme_classic() + scale_color_manual(values = colors) + 
    guides(color = FALSE) + 
    geom_line(data = line_data.right, aes(x = monthsFromBaseline, y = value), color = "black",linetype = "dashed",size = 1) +
    geom_pointrange(aes(ymin = value-CI95, ymax = value+CI95),alpha = 1,data = medians.right,shape = 5,size = .7,linewidth = .65,stroke = 1.4) +
    geom_jitter(position = position_jitter(0.005),alpha = .45,size = 2) + 
    #ylab(expression(Volume~(mm^3))) +
    #xlab("Months After Baseline") +
    xlab("Months") +
    ylab(NULL) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_y_continuous(limits = c(min(breaks_y_right+25), max(breaks_y_right)+1), breaks = breaks_y_right, expand = c(0, 0)) +
    scale_x_continuous(limits = c(xmin,xmax), breaks = seq(0,21,3), expand = c(0, 0)) +
    theme(text = element_text(family = "Helvetica"),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.text =  element_text(color="black",size=14, margin = margin(b = 0)),
          axis.title.x = element_text(size=16), # Remove the x-axis title
          axis.title.y = element_text(size=16),
          legend.text = element_text(size=14),
          # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
          legend.position=c(0.01,1), 
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  figure <- ggarrange(p1, p2,ncol = 2, nrow = 1)
  #figure_with_xlabel <- annotate_figure(figure,vbottom = text_grob("Months After Baseline", family = "Helvetica", size = 16,vjust = -.10))
  plotTemplateName<-as_tibble(templateName) %>% mutate(newVar = case_when(
    value == "ADNI" ~ "ADNI",
    value == "inclusive" ~ "All 1.0 mm CS",
    value == "crossRes" ~ "Multi-resolution CS",
    TRUE ~ "other"  # default case
  ))
  
  title <- ggdraw() + 
    draw_label(title, size = 14, fontface = 'bold', hjust = 0.5, x = 0.5, y = .8, vjust = 1)
  # Combine the title and the combined figure
  outFigure <- plot_grid(title, figure, ncol = 1, rel_heights = c(0.08, 1))
  
  ## Return outFigure if you want a title to be added
  return(figure)
}

##Similar to above but adopted for non-hippocampus ranges
## Yaxis Margin is used to help line up plots in dashboard
## Yaxis breaks is used to fix weird y-axes that don't automate well. This should be a vector of the ticks you want on the yaxis.

generate_plot_testRetest.dash <- function(df,id, brainVar,type,templateName,title,yaxisMargin,yAxisBreaks.left,yAxisBreaks.right) { #type has to be either vol, Thick or GWcon or Area # templateName is either ADNI, crossRes or inclusive
  
  plotData.ID.allTemplate.left <- df %>% filter(ID == id, grepl(brainVar,region), origData == type) %>% filter(grepl("Left-",region))
  plotData.ID.allTemplate.right <- df %>% filter(ID == id, grepl(brainVar,region), origData == type) %>% filter(grepl("Right-",region))
  
  plotData <- df %>% filter(ID == id, grepl(brainVar,region), origData == type) 
  if(templateName == "inclusive"){
    plotData <- plotData %>% filter(grepl("CSx6",template)) %>% filter(!grepl("crossRes",template))
  }else{
    plotData <- plotData %>% filter(grepl(templateName,template))
  }
  plotData <- plotData %>% mutate(visit = as_factor(ceiling(session/2)))
  plotData <- plotData %>% mutate(visit = as_factor(if_else(visit == 1, "Baseline", if_else(visit == 2, "Visit 2",if_else(visit == 3, "Visit 3","Visit 4")))))
  
  ###Estimate time from baseline for each scan
  plotData$Date <- ymd(plotData$Date)
  baselineDate<-min(plotData$Date)
  plotData<-plotData %>% mutate(monthsFromBaseline = as.numeric(Date - baselineDate)/(365.25/12))
  
  
  #Calculate average number of months from baseline for each visit
  plotData <- plotData %>% group_by(visit) %>% mutate(visitMonths = round(mean(monthsFromBaseline),2) )
  allPlotData <- plotData %>% group_by(ID) %>% mutate(visit = as_factor(ceiling(session/2))) %>% ungroup()
  
  allPlotData.left <- allPlotData %>% filter(grepl("Left-",region))
  allPlotData.right <- allPlotData %>% filter(grepl("Right-",region))
  
  ###Separate ID data
  plotData.ID <- allPlotData %>% filter(ID == id)
  
  ###Grab demo for plot
  dx=plotData.ID$dx[1]
  id=plotData.ID$ID[1]
  age=plotData.ID$Age[1]
  sex=plotData.ID$Sex[1]
  
  ###Split into hemisperes and handle % change
  plotData.ID.left <- plotData.ID %>% filter(grepl("Left-",region))
  plotData.ID.right <- plotData.ID %>% filter(grepl("Right-",region))
  
  #Adjust x-axis to be bigger only if 7 sessions
  plotData.allIDs <- df %>% filter(grepl(brainVar,region), origData == type) 
  if(templateName == "inclusive"){
    plotData.allIDs <- plotData.allIDs %>% filter(grepl("CSx6",template)) %>% filter(!grepl("crossRes",template))
  }else{
    plotData.allIDs <- plotData.allIDs %>% filter(grepl(templateName,template))
  }
  
  #xmax<-(max(allPlotData.left$monthsFromBaseline) + 1) ## Use this if you want a different x-axis for each ID
  xmax<-(max(plotData.allIDs$yearsFromBaseline*12) +1) ## Use this if you want the same x-axis for each ID
  xmin <- 15-xmax # so that the spacing between y-axis and 0 is the same as the end of the x-axis and last tick spacing
  ###Simplify to a linear regression in each person - used to have lme model but the complexity caused some issues...
  
  model_lm.left <- lm(value ~ monthsFromBaseline, data = plotData.ID.left)
  lineMax<-xmax
  new_data.ID.left <- data.frame(ID = id, monthsFromBaseline = seq(0, lineMax, length.out = 100))
  pred_intervals.left <- predict(model_lm.left, newdata = new_data.ID.left)# , level = 0.80, n.sims = 500, seed = 123, include.resid.var = F,type = "linear.prediction")
  new_data.ID.left<-bind_cols(new_data.ID.left,pred_intervals.left)
  new_data.ID.left <- new_data.ID.left %>% rename(value = ...3) ##To make it plot right
  ##For the plot estimate the annual %change from the lme regression intercept and slope
  model.ID.left<-coef(model_lm.left)
  annChangePer.left <- round((model.ID.left[2]*12)/model.ID.left[1]*100,2)
  CI.model.ID.left<-confint(model_lm.left)
  CI.annChangePer.left <- round((CI.model.ID.left[2,]*12)/model.ID.left[1]*100,2)
  
  # Assuming new_data_specific_id contains 'monthsFromBaseline' from 0 to 15 and 'predicted_value'
  # Calculate slope (m) using the first and last fit values
  x1 <- 0  # Starting point of x
  x2 <- xmax  # Ending point of x
  
  # Ensure new_data_specific_id is ordered by monthsFromBaseline if it's not already
  new_data.ID.left <- new_data.ID.left[order(new_data.ID.left$monthsFromBaseline), ]
  y1 <- predict(model_lm.left, newdata = data.frame(monthsFromBaseline = x1, ID = id), re.form = NULL)
  y2 <- predict(model_lm.left, newdata = data.frame(monthsFromBaseline = x2), re.form = NULL)
  # Define the line endpoints based on the calculated slope
  line_data.left <- data.frame(monthsFromBaseline = c(x1, x2), value = c(y1, y2))
  
  
  model_lm.right <- lm(value ~ monthsFromBaseline, data = plotData.ID.right)
  new_data.ID.right <- data.frame(ID = id, monthsFromBaseline = seq(0, lineMax, length.out = 100))
  pred_intervals.right <- predict(model_lm.right, newdata = new_data.ID.right)# , level = 0.80, n.sims = 500, seed = 123, include.resid.var = F,type = "linear.prediction")
  new_data.ID.right<-bind_cols(new_data.ID.right,pred_intervals.right)
  new_data.ID.right <- new_data.ID.right %>% rename(value = ...3) ##To make it plot right
  ##For the plot estimate the annual %change from the lme regression intercept and slope
  model.ID.right<-coef(model_lm.right)
  annChangePer.right <- round((model.ID.right[2]*12)/model.ID.right[1]*100,2)
  CI.model.ID.right<-confint(model_lm.right)
  CI.annChangePer.right <- round((CI.model.ID.right[2,]*12)/model.ID.right[1]*100,2)
  
  # Ensure new_data_specific_id is ordered by monthsFromBaseline if it's not already
  new_data.ID.right <- new_data.ID.right[order(new_data.ID.right$monthsFromBaseline), ]
  y1 <- predict(model_lm.right, newdata = data.frame(monthsFromBaseline = x1, ID = id), re.form = NULL)
  y2 <- predict(model_lm.right, newdata = data.frame(monthsFromBaseline = x2), re.form = NULL)
  # Define the line endpoints based on the calculated slope
  line_data.right <- data.frame(monthsFromBaseline = c(x1, x2), value = c(y1, y2))
  
  ##Left
  mean.base.plotData.ID.left <- plotData.ID.left %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(mean=mean(value))
  median.base.plotData.ID.left <- plotData.ID.left %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(median=median(value))
  median.test.plotData.ID.left <- plotData.ID.left %>% ungroup() %>% filter(session == 1) %>% summarise(median=median(value))
  median.retest.plotData.ID.left <- plotData.ID.left %>% ungroup() %>% filter(session == 2) %>% summarise(median=median(value))
  
  plotData.ID.left <- plotData.ID.left %>% mutate(perDiff=(value-median.base.plotData.ID.left$median)/median.base.plotData.ID.left$median*100)
  plotData.ID.left <- plotData.ID.left %>% mutate(perDiffperYear = perDiff/(monthsFromBaseline/12))
  medians.left <- plotData.ID.left %>% group_by(monthsFromBaseline,session) %>% summarize(median = median(value),monthsFromBaseline = mean(monthsFromBaseline),se = sd(value, na.rm = TRUE)/sqrt(n()),n=n(), vis = median(as.numeric(visit)))
  medians.left <- medians.left %>% mutate(CI95 = 1.96*se)
  ###Estimate of the benefits of reducing error with many samples
  sum.sdMed.left<- medians.left %>% group_by(vis) %>% summarise(sdMed = sd(median))
  meanMedSD.left <- mean(sum.sdMed.left$sdMed)
  
  visit_Avgs_long.left <- plotData.ID.left %>% group_by(visit) %>% summarise(perDiff2 = mean(perDiff),se=sd(perDiff)/sqrt(n()))
  visit_Avgs_long.left <- visit_Avgs_long.left %>% mutate(CI95 = 1.96*se)
  
  df.summary.left <- plotData.ID.left %>% group_by(visit) %>% summarise(median = round(median(value),0),sd = sd(value),se = sd(value, na.rm = TRUE)/sqrt(n()), perDiff = round(median(perDiff),2), perDiffAnn = round(median(perDiffperYear),2))
  df.summary.left <- df.summary.left %>% mutate(CI95 = 1.96*se)
  meanSD.left<-mean(df.summary.left$sd)
  
  ##Right
  mean.base.plotData.ID.right <- plotData.ID.right %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(mean=mean(value))
  median.base.plotData.ID.right <- plotData.ID.right %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(median=median(value))
  median.test.plotData.ID.right <- plotData.ID.right %>% ungroup() %>% filter(session == 1) %>% summarise(median=median(value))
  median.retest.plotData.ID.right <- plotData.ID.right %>% ungroup() %>% filter(session == 2) %>% summarise(median=median(value))
  
  plotData.ID.right <- plotData.ID.right %>% mutate(perDiff=(value-median.base.plotData.ID.right$median)/median.base.plotData.ID.right$median*100)
  plotData.ID.right <- plotData.ID.right %>% mutate(perDiffperYear = perDiff/(monthsFromBaseline/12))
  medians.right <- plotData.ID.right %>% group_by(monthsFromBaseline,session) %>% summarize(median = median(value),monthsFromBaseline = mean(monthsFromBaseline),se = sd(value, na.rm = TRUE)/sqrt(n()), vis = median(as.numeric(visit)))
  medians.right <- medians.right %>% mutate(CI95 = 1.96*se)
  sum.sdMed.right<- medians.right %>% group_by(vis) %>% summarise(sdMed = sd(median))
  meanMedSD.right <- mean(sum.sdMed.right$sdMed)
  
  plotData.ID.right.test <- plotData.ID.right %>% filter(session %in% c(1,3,5)) %>% mutate(perDiff=(value - median.test.plotData.ID.right$median) / median.test.plotData.ID.right$median * 100)
  plotData.ID.right.test <- plotData.ID.right.test %>% mutate(perDiffperYear = perDiff/(monthsFromBaseline/12))
  medians.right.test <- plotData.ID.right.test %>% group_by(monthsFromBaseline,session) %>% summarize(value = median(value))
  plotData.ID.right.retest <- plotData.ID.right %>% filter(session %in% c(2,4,6)) %>% mutate(perDiff=(value - median.retest.plotData.ID.right$median) / median.retest.plotData.ID.right$median * 100)
  plotData.ID.right.retest <- plotData.ID.right.retest %>% mutate(perDiffperYear = perDiff/(monthsFromBaseline/12))
  medians.right.retest <- plotData.ID.right.retest %>% group_by(monthsFromBaseline,session) %>% summarize(value = median(value))
  
  plotData.ID.right.testRetest<-bind_rows(plotData.ID.right.test,plotData.ID.right.retest)
  
  #session_averages_long.right <- plotData.ID.right.testRetest %>% group_by(session) %>%
  # summarise(avg_perDiff = mean(perDiff), sd_perDiff = sd(perDiff), se = sd(perDiff) / sqrt(n()), med_perDiff = median(perDiff)) %>% 
  #   mutate(visit = rep(1:(n() %/% 2), each = 2))
  
  # visit_Avgs.right <- session_averages_long.right %>% group_by(visit) %>% summarise(visitAvg = mean(avg_perDiff),se=mean(se))
  # visit_Avgs_long.right <- plotData.ID.right %>% group_by(visit) %>% summarise(perDiff2 = mean(perDiff),se=sd(perDiff)/sqrt(n()))
  # visit_Avgs_long.right <- visit_Avgs_long.right %>% mutate(CI95 = 1.96*se)
  
  df.summary.right <- plotData.ID.right %>% group_by(visit) %>% summarise(median = round(median(value),0),se = sd(value, na.rm = TRUE)/sqrt(n()), perDiff = round(median(perDiff),2), perDiffAnn = round(median(perDiffperYear),2),sd=sd(value))
  df.summary.right <- df.summary.right %>% mutate(CI95 = 1.96*se)
  meanSD.right<-mean(df.summary.right$sd)
  
  medians.right <- medians.right %>% rename(value = median) %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  medians.left <- medians.left %>% rename(value = median) %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  ###Make Plots
  
  plotData.ID.left<- plotData.ID.left %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  plotData.ID.right<- plotData.ID.right %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  
  # Define the colors for odd and even visits
  colors <- c("Even" = "steelblue4", "Odd" = "red4")
  
  #colors <- c("red4", "steelblue4","red4", "steelblue4","red4","steelblue4","red4")
  
  ### Make range so that it's consistent except in extreme cases
  ### Make range so that it is the same size for left and right 
  
  ## use the larger size of the range as the range size for both left and right plots
  range.left <- range(plotData.ID.left$value)[2] - range(plotData.ID.left$value)[1] 
  range.right <- range(plotData.ID.right$value)[2] - range(plotData.ID.right$value)[1]
  #max(c(plotData.ID.left$value,plotData.ID.right$value))
  #rangeSize<-max(340,range.left,range.right)*1.2 ## Range should be at least 300 so that there aren't severe differences between participants range
  rangeSize<-max(range.left,range.right)*1.4 
  
  ##Need to make sure there is enough space for a label in the bottom left corner of the plot.
  plotData.ID.left.12 <- plotData.ID.left %>% filter(session %in% c(1,2))
  
  minRange.left<-min(c(min(plotData.ID.left$value)*.99,mean(plotData.ID.left$value)*.97,min(plotData.ID.left.12$value)))
  maxRange.left<-minRange.left + rangeSize
  
  plot_title.left<-paste0("Left ",brainVar,"\n", annChangePer.left,"% [",CI.annChangePer.left[1],",  ",CI.annChangePer.left[2],"]")
  
  
  # Define a function to compute 5 tick marks
  compute_fixed_breaks <- function(min_val, max_val) {
    # Calculate the total range
    total_range <- max_val - min_val
    
    # Calculate the step size to ensure we get 5 ticks
    # 5 ticks means 4 steps between them
    step <- total_range / 4
    
    # Generate the 5 tick marks
    return(seq(from = min_val, to = max_val, by = step))
  }
  
  # Calculate y-axis breaks
  breaks_y_left <- compute_fixed_breaks(minRange.left, maxRange.left)
  
  print(plot_title.left)
  
  ### Check if user wants to manually set y-axis ticks and range and use those if supplied
  if (missing(yAxisBreaks.left)) {
    plot.yAxisBreaks.left <- extended(dmin = min(breaks_y_left), dmax = max(breaks_y_left), m = 6)
  } else {
    plot.yAxisBreaks.left <- yAxisBreaks.left
  }
  
  p1<-ggplot(plotData.ID.left, aes(monthsFromBaseline, value, color = session_type)) + theme_classic() + scale_color_manual(values = colors) + 
    guides(color = FALSE) + 
    geom_line(data = line_data.left, aes(x = monthsFromBaseline, y = value), color = "black",linetype = "dashed",size = 1) +
    geom_pointrange(aes(ymin = value-CI95, ymax = value+CI95),alpha = 1,data = medians.left,shape = 5,size = .7,linewidth = .65,stroke = 1.4) +
    geom_jitter(position = position_jitter(0.005),alpha = .45,size = 2) + 
    xlab("Months") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_y_continuous(
      breaks = plot.yAxisBreaks.left) +
    scale_x_continuous(
      breaks = seq(0, 21, 3), 
      expand = c(0, 0)) +
    coord_cartesian(
      ylim = c(min(plot.yAxisBreaks.left), max(plot.yAxisBreaks.left)),
      xlim = c(xmin, xmax)) + 
    theme(text = element_text(family = "Helvetica"),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          axis.text.y = element_text(color = "black", size = 14, margin = margin(t = 0, r = 1, b = 0, l = yaxisMargin)),
          axis.text.x = element_text(color = "black", size = 14, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(size=16), # Remove the x-axis title
          axis.title.y = element_text(size=16),
          legend.text = element_text(size=14),
          plot.margin = margin(0.1,0,0.1,0, "in"),
          legend.position=c(0.01,1), 
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  ##Need to make sure there is enough space for a label in the bottom left corner of the plot.
  plotData.ID.right.12 <- plotData.ID.right %>% filter(session %in% c(1,2))
  minRange.right<-min(c(min(plotData.ID.right$value)*.99,mean(plotData.ID.right$value)*.97,min(plotData.ID.right.12$value)))
  maxRange.right<-minRange.right + rangeSize
  
  breaks_y_right <- compute_fixed_breaks(minRange.right, maxRange.right)
  plot_title.right<-paste0("Right ",brainVar,"\n", annChangePer.right,"% [",CI.annChangePer.right[1],",  ",CI.annChangePer.right[2],"]")
  print(plot_title.right)
  
  ### Check if user wants to manually set y-axis ticks and range and use those if supplied
  if (missing(yAxisBreaks.right)) {
    plot.yAxisBreaks.right <- extended(dmin = min(breaks_y_right), dmax = max(breaks_y_right), m = 6)
  } else {
    plot.yAxisBreaks.right<- yAxisBreaks.right
  }
  
  p2<-ggplot(plotData.ID.right, aes(monthsFromBaseline, value, color = session_type)) + theme_classic() + scale_color_manual(values = colors) + 
    guides(color = FALSE) + 
    geom_line(data = line_data.right, aes(x = monthsFromBaseline, y = value), color = "black",linetype = "dashed",size = 1) +
    geom_pointrange(aes(ymin = value-CI95, ymax = value+CI95),alpha = 1,data = medians.right,shape = 5,size = .7,linewidth = .65,stroke = 1.4) +
    geom_jitter(position = position_jitter(0.005),alpha = .45,size = 2) + 
    xlab("Months") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_y_continuous(
      breaks = plot.yAxisBreaks.right) +
    scale_x_continuous(
      breaks = seq(0, 21, 3), 
      expand = c(0, 0)) +
    coord_cartesian(
      ylim = c(min(plot.yAxisBreaks.right), max(plot.yAxisBreaks.right)),
      xlim = c(xmin, xmax)) + 
    theme(text = element_text(family = "Helvetica"),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          #axis.text =  element_text(color="black",size=14, margin = margin(b = 0)),
          axis.text.y = element_text(color = "black", size = 14, margin = margin(t = 0, r = 1, b = 0, l = 0)),
          axis.text.x = element_text(color = "black", size = 14, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(size=16), # Remove the x-axis title
          axis.title.y = element_text(size=16),
          legend.text = element_text(size=14),
          plot.margin = margin(0.1,0,0.1,0, "in"),
          legend.position=c(0.01,1), 
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  figure <- ggarrange(p1, p2,ncol = 2, nrow = 1)

  title <- ggdraw() + 
    draw_label(title, size = 14, fontface = 'bold', hjust = 0.5, x = 0.5, y = .8, vjust = 1)
  # Combine the title and the combined figure
  outFigure <- plot_grid(title, figure, ncol = 1, rel_heights = c(0.08, 1))
  
  ## Return outFigure if you want a title to be added
  return(outFigure)
}

### Unilateral companion to function above

id=""
brainVar="WMH"
type="vol"
templateName="CSx6_inclusive"
title=""


generate_plot_global <- function(df, id, brainVar,type,templateName,title,yaxisMargin,yAxisBreaks) { #type has to be either vol, Thick or GWcon or Area # templateName is either ADNI, crossRes or inclusive
  
  plotData.ID.allTemplate<- df %>% filter(ID == id, region == brainVar, origData == type)
  
  plotData <- df %>% filter(ID == id, region == brainVar, origData == type) 
  if(templateName == "inclusive"){
    plotData <- plotData %>% filter(grepl("CSx6",template)) %>% filter(!grepl("crossRes",template))
  }else{
    plotData <- plotData %>% filter(grepl(templateName,template))
  }
  plotData <- plotData %>% mutate(visit = as_factor(ceiling(session/2)))
  plotData <- plotData %>% mutate(visit = as_factor(if_else(visit == 1, "Baseline", if_else(visit == 2, "Visit 2",if_else(visit == 3, "Visit 3","Visit 4")))))
  
  ###Estimate time from baseline for each scan
  plotData$Date <- ymd(plotData$Date)
  baselineDate<-min(plotData$Date)
  plotData<-plotData %>% mutate(monthsFromBaseline = as.numeric(Date - baselineDate)/(365.25/12))
  
  
  #Calculate average number of months from baseline for each visit
  plotData <- plotData %>% group_by(visit) %>% mutate(visitMonths = round(mean(monthsFromBaseline),2) )
  allPlotData <- plotData %>% group_by(ID) %>% mutate(visit = as_factor(ceiling(session/2))) %>% ungroup()
  
  ###Separate ID data
  plotData.ID <- allPlotData %>% filter(ID == id)
  
  ###Grab demo for plot
  dx=plotData.ID$dx[1]
  id=plotData.ID$ID[1]
  age=plotData.ID$Age[1]
  sex=plotData.ID$Sex[1]
  
  #Adjust x-axis to be bigger only if 7 sessions
  plotData.allIDs <- df %>% filter(grepl(brainVar,region), origData == type) 
  if(templateName == "inclusive"){
    plotData.allIDs <- plotData.allIDs %>% filter(grepl("CSx6",template)) %>% filter(!grepl("crossRes",template))
  }else{
    plotData.allIDs <- plotData.allIDs %>% filter(grepl(templateName,template))
  }
  #Adjust x-axis to be bigger only if 7 sessions
  xmax<-(max(allPlotData$monthsFromBaseline) + 1)
  xmax<-(max(plotData.allIDs$yearsFromBaseline*12) +1)
  ###Simplify to a linear regression in each person - used to have lme model but the complexity caused some issues...
  
  model_lm <- lm(value ~ monthsFromBaseline, data = plotData.ID)
  lineMax<-xmax #max(plotData.ID$monthsFromBaseline) +0.3
  new_data.ID <- data.frame(ID = id, monthsFromBaseline = seq(0, lineMax, length.out = 100))
  pred_intervals <- predict(model_lm, newdata = new_data.ID)# , level = 0.80, n.sims = 500, seed = 123, include.resid.var = F,type = "linear.prediction")
  
  new_data.ID<-bind_cols(new_data.ID,pred_intervals)
  new_data.ID <- new_data.ID %>% rename(value = ...3) ##To make it plot right
  ##For the plot estimate the annual %change from the lme regression intercept and slope
  model.ID<-coef(model_lm)
  annChangePer <- round((model.ID[2]*12)/model.ID[1]*100,2)
  CI.model.ID<-confint(model_lm)
  CI.annChangePer <- round((CI.model.ID[2,]*12)/model.ID[1]*100,2)
  
  # Assuming new_data_specific_id contains 'monthsFromBaseline' from 0 to 15 and 'predicted_value'
  # Calculate slope (m) using the first and last fit values
  x1 <- 0  # Starting point of x
  x2 <- lineMax  # Ending point of x
  
  # Ensure new_data_specific_id is ordered by monthsFromBaseline if it's not already
  new_data.ID <- new_data.ID[order(new_data.ID$monthsFromBaseline), ]
  y1 <- predict(model_lm, newdata = data.frame(monthsFromBaseline = x1, ID = id), re.form = NULL)
  y2 <- predict(model_lm, newdata = data.frame(monthsFromBaseline = x2), re.form = NULL)
  # Define the line endpoints based on the calculated slope
  line_data <- data.frame(monthsFromBaseline = c(x1, x2), value = c(y1, y2))
  
  
  ##Left
  mean.base.plotData.ID <- plotData.ID %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(mean=mean(value))
  median.base.plotData.ID <- plotData.ID %>% ungroup() %>% filter(session %in% c(1,2)) %>% summarise(median=median(value))
  median.test.plotData.ID <- plotData.ID %>% ungroup() %>% filter(session == 1) %>% summarise(median=median(value))
  median.retest.plotData.ID <- plotData.ID %>% ungroup() %>% filter(session == 2) %>% summarise(median=median(value))
  
  plotData.ID <- plotData.ID %>% mutate(perDiff=(value-median.base.plotData.ID$median)/median.base.plotData.ID$median*100)
  plotData.ID <- plotData.ID %>% mutate(perDiffperYear = perDiff/(monthsFromBaseline/12))
  medians <- plotData.ID %>% group_by(monthsFromBaseline,session) %>% summarize(median = median(value),monthsFromBaseline = mean(monthsFromBaseline),se = sd(value, na.rm = TRUE)/sqrt(n()),n=n(), vis = median(as.numeric(visit)))
  medians <- medians %>% mutate(CI95 = 1.96*se)
  ###Estimate of the benefits of reducing error with many samples
  sum.sdMed<- medians %>% group_by(vis) %>% summarise(sdMed = sd(median))
  meanMedSD <- mean(sum.sdMed$sdMed)
  
  visit_Avgs_long <- plotData.ID %>% group_by(visit) %>% summarise(perDiff2 = mean(perDiff),se=sd(perDiff)/sqrt(n()))
  visit_Avgs_long <- visit_Avgs_long %>% mutate(CI95 = 1.96*se)
  
  df.summary <- plotData.ID %>% group_by(visit) %>% summarise(median = round(median(value),0),sd = sd(value),se = sd(value, na.rm = TRUE)/sqrt(n()), perDiff = round(median(perDiff),2), perDiffAnn = round(median(perDiffperYear),2))
  df.summary <- df.summary %>% mutate(CI95 = 1.96*se)
  meanSD<-mean(df.summary$sd)
  
  medians <- medians %>% rename(value = median) %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  ###Make Plots
  
  plotData.ID <- plotData.ID %>% mutate(session_type = ifelse(as.numeric(session) %% 2 == 0, "Even", "Odd"))
  
  # Define the colors for odd and even visits
  colors <- c("Even" = "steelblue4", "Odd" = "red4")
  
  #Adjust x-axis to be bigger only if 7 sessions
  plotData.allIDs <- df %>% filter(grepl(brainVar,region), origData == type) 
  if(templateName == "inclusive"){
    plotData.allIDs <- plotData.allIDs %>% filter(grepl("CSx6",template)) %>% filter(!grepl("crossRes",template))
  }else{
    plotData.allIDs <- plotData.allIDs %>% filter(grepl(templateName,template))
  }
  
  #xmax<-(max(allPlotData.left$monthsFromBaseline) + 1) ## Use this if you want a different x-axis for each ID
  xmax<-(max(plotData.allIDs$yearsFromBaseline*12) +1) ## Use this if you want the same x-axis for each ID
  xmin <- 15-xmax # so that the spacing between y-axis and 0 is the same as the end of the x-axis and last tick spacing
  
  #colors <- c("red4", "steelblue4","red4", "steelblue4","red4","steelblue4","red4")
  
  ###Make range so that it's consistent except in extreme cases
  minRange<-min(c(min(plotData.ID.allTemplate$value)*.99,mean(plotData.ID.allTemplate$value)*.95))
  maxRange<-max(c(max(plotData.ID.allTemplate$value)*1.01,mean(plotData.ID.allTemplate$value)*1.05))
  
  plot_title<-paste0(brainVar,"\n", annChangePer,"% [",CI.annChangePer[1],",  ",CI.annChangePer[2],"]")
  
  print(plot_title)
  
  #breaks_y_right <- compute_fixed_breaks(minRange.right, maxRange.right)
  #plot_title.right<-paste0("Right ",brainVar,"\n", annChangePer.right,"% [",CI.annChangePer.right[1],",  ",CI.annChangePer.right[2],"]")
  #print(plot_title.right)
  
  ### Check if user wants to manually set y-axis ticks and range and use those if supplied
  if (missing(yAxisBreaks)) {
    plot.yAxisBreaks <- extended(dmin =  minRange,dmax = maxRange,m = 6)
  } else {
    plot.yAxisBreaks<- yAxisBreaks
  }
  
  p1<-ggplot(plotData.ID, aes(monthsFromBaseline, value, color = session_type)) + theme_classic() + scale_color_manual(values = colors) + 
    guides(color = FALSE) + 
    geom_pointrange(aes(ymin = value-CI95, ymax = value+CI95),alpha = 1,data = medians,shape = 5,size = .7,linewidth = .65,stroke = 1.4) +
    geom_jitter(position = position_jitter(0.005),alpha = .45,size = 2) + 
    geom_line(data = line_data, aes(x = monthsFromBaseline, y = value), color = "black",linetype = "dashed",size = 1) +
    xlab("Months") +
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_y_continuous(limits = c(min(plot.yAxisBreaks), max(plot.yAxisBreaks)),
                       breaks = plot.yAxisBreaks) +
    scale_x_continuous(limits = c(xmin,xmax), breaks = seq(0,21,3), expand = c(0, 0)) +
    theme(text = element_text(family = "Helvetica"),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          #axis.text =  element_text(color="black",size=14, margin = margin(b = 0)),
          axis.text.y = element_text(color = "black", size = 14, margin = margin(t = 0, r = 1, b = 0, l = yaxisMargin)),
          axis.text.x = element_text(color = "black", size = 14, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(size=16), # Remove the x-axis title
          axis.title.y = element_text(size=16),
          legend.text = element_text(size=14),
          # plot.margin = margin(0.15,0.1,0.1,0.1, "in"),
          plot.margin = margin(0.1,0,0.1,0, "in"),
          legend.position=c(0.01,1), 
          legend.justification=c(0,1),
          legend.background = element_rect(fill = "transparent"), # Make legend background transparent
          legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.text.align = 0, plot.title=element_text(size=16,hjust=0.5,margin=margin(b=10))) 
  
  #figure <- ggarrange(p1, p2,ncol = 2, nrow = 1)
  plotTemplateName<-as_tibble(templateName) %>% mutate(newVar = case_when(
    value == "ADNI" ~ "ADNI",
    value == "inclusive" ~ "All 1.0 mm CS",
    value == "crossRes" ~ "Multi-resolution CS",
    TRUE ~ "other"  # default case
  ))
  
  title <- ggdraw() + 
    draw_label(title, size = 14, fontface = 'bold', hjust = 0.5, x = 0.5, y = .8, vjust = 1)
  # Combine the title and the combined figure
  outFigure <- plot_grid(title, p1, ncol = 1, rel_heights = c(0.08, 1))
  
  return(outFigure)
}




