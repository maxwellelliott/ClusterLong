library(conflicted)
library(RColorBrewer)
library(boot)
library(corrplot)
library(psych)
library(svglite)
library(magick)
library(png)
library(ggpubr)
library(lme4)
library(merTools)
library(ggdist)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(rstanarm)
library(ggthemes)
library(broom.mixed)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")


set.seed(221219)

demo<-read_csv("/path/to/data/ClusterLong/Data/ClusterDemographics_240129.csv")


volData<-read_table("/path/to/data/ClusterLong/Data/aseg.vol.long.table.240710.allTemplates")
lhThickData<-read_table("/path/to/data/ClusterLong/Data/aparc.lh.thickness.long.table.240710.allTemplates")
rhThickData<-read_table("/path/to/data/ClusterLong/Data/aparc.rh.thickness.long.table.240710.allTemplates")
lhGWconData<-read_table("/path/to/data/ClusterLong/Data/aparc.lh.GWratio.long.table.240710.allTemplates")
rhGWconData<-read_table("/path/to/data/ClusterLong/Data/aparc.rh.GWratio.long.table.240710.allTemplates")

ADsigData<-read_csv("/path/to/data/ClusterLong/Data/ADsignatureMask_meanCT_inclusive_240722.csv",col_names = c("ID","template","longScanID","hemisphere","value"))
ADsigData_GWR<-read_csv("/path/to/data/ClusterLong/Data/ADsignatureMask_meanGWR_inclusive_240722.csv",col_names = c("ID","template","longScanID","hemisphere","value"))
FTDsigData<-read_csv("/path/to/data/ClusterLong/Data/FTDsignatureMask_meanCT_inclusive_240722.csv",col_names = c("ID","template","longScanID","hemisphere","value"))
FTDsigData_GWR<-read_csv("/path/to/data/ClusterLong/Data/FTDsignatureMask_meanGWR_inclusive_240722.csv",col_names = c("ID","template","longScanID","hemisphere","value"))

##Grab data from scan names
volData<-volData %>% separate(c("Date", "ID","scanID"),col = "Measure:volume",sep = "_",extra = "merge",convert = T)  
lhThickData <- lhThickData %>% separate(c("Date", "ID","scanID"),col = "lh.aparc.thickness",sep = "_",extra = "merge",convert = T)
rhThickData <- rhThickData %>% separate(c("Date", "ID","scanID"),col = "rh.aparc.thickness",sep = "_",extra = "merge",convert = T)
lhGWconData <- lhGWconData %>% separate(c("Date", "ID","scanID"),col = "Measure:mean",sep = "_",extra = "merge",convert = T)
rhGWconData <- rhGWconData %>% separate(c("Date", "ID","scanID"),col = "Measure:mean",sep = "_",extra = "merge",convert = T)
ADsigData <- ADsigData %>% separate(c("Date", "ID","scanID"),col = "longScanID",sep = "_",extra = "merge",convert = T)
ADsigData_GWR <- ADsigData_GWR %>% separate(c("Date", "ID","scanID"),col = "longScanID",sep = "_",extra = "merge",convert = T)
FTDsigData <- FTDsigData %>% separate(c("Date", "ID","scanID"),col = "longScanID",sep = "_",extra = "merge",convert = T)
FTDsigData_GWR <- FTDsigData_GWR %>% separate(c("Date", "ID","scanID"),col = "longScanID",sep = "_",extra = "merge",convert = T)

##Clean up ADsig
ADsigDataWide<-ADsigData %>% pivot_wider(names_from = hemisphere,values_from = value,names_glue = "{hemisphere}-ADsignature") %>% rename(`Left-ADsignature_CT` = `lh-ADsignature`, `Right-ADsignature_CT` = `rh-ADsignature`) %>% select(-template) %>% 
  mutate(meanADsignature_CT = (`Left-ADsignature_CT` + `Right-ADsignature_CT`)/2)
ADsigData_GWRWide<-ADsigData_GWR %>% pivot_wider(names_from = hemisphere,values_from = value,names_glue = "{hemisphere}-ADsignature") %>% rename(`Left-ADsignature_GWR` = `lh-ADsignature`, `Right-ADsignature_GWR` = `rh-ADsignature`) %>% select(-template)  %>% 
  mutate(meanADsignature_GWR = (`Left-ADsignature_GWR` + `Right-ADsignature_GWR`)/2)
FTDsigDataWide<-FTDsigData %>% pivot_wider(names_from = hemisphere,values_from = value,names_glue = "{hemisphere}-FTDsignature") %>% rename(`Left-FTDsignature_CT` = `lh-FTDsignature`, `Right-FTDsignature_CT` = `rh-FTDsignature`) %>% select(-template) %>% 
  mutate(meanFTDsignature_CT = (`Left-FTDsignature_CT` + `Right-FTDsignature_CT`)/2)
FTDsigData_GWRWide<-FTDsigData_GWR %>% pivot_wider(names_from = hemisphere,values_from = value,names_glue = "{hemisphere}-FTDsignature") %>% rename(`Left-FTDsignature_GWR` = `lh-FTDsignature`, `Right-FTDsignature_GWR` = `rh-FTDsignature`) %>% select(-template)  %>% 
  mutate(meanFTDsignature_GWR = (`Left-FTDsignature_GWR` + `Right-FTDsignature_GWR`)/2)

#ADD mean GWcon for each Hemispere
lhGWconData <- lhGWconData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(`MeanGWCon` = mean(c_across(5:36), na.rm = TRUE)) %>% ungroup()
rhGWconData <- rhGWconData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(`MeanGWCon` = mean(c_across(5:36), na.rm = TRUE)) %>% ungroup()

###ADD in combined ventricle variables for each Hemisphere
volData <- volData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(allVentricles = sum(`Left-Inf-Lat-Vent`,`Left-Lateral-Ventricle`,`Right-Inf-Lat-Vent`,`Right-Lateral-Ventricle`,`3rd-Ventricle`,`4th-Ventricle`,`5th-Ventricle`, na.rm = TRUE)) %>% ungroup()
volData <- volData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(`Left-totalLateralVentricle` = sum(`Left-Inf-Lat-Vent`,`Left-Lateral-Ventricle`, na.rm = TRUE), `Right-totalLateralVentricle` = sum(`Right-Inf-Lat-Vent`,`Right-Lateral-Ventricle`, na.rm = TRUE)) %>% ungroup()

##Add in left and right total brain volume measures for later plotting
volData <- volData %>% rowwise() %>%
  mutate(`Left-totalBrainVolume` = sum(`lhCortexVol`, `lhCerebralWhiteMatterVol`, `Left-Accumbens-area`, `Left-Amygdala`, `Left-Caudate`, 
                                       `Left-Hippocampus`, `Left-Pallidum`, `Left-Putamen`, `Left-Thalamus-Proper`, `Left-VentralDC`, 
                                       `Left-Cerebellum-Cortex`, `Left-Cerebellum-White-Matter`, na.rm = TRUE),
         `Right-totalBrainVolume` = sum(`rhCortexVol`, `rhCerebralWhiteMatterVol`, `Right-Accumbens-area`, `Right-Amygdala`, `Right-Caudate`, 
                                        `Right-Hippocampus`, `Right-Pallidum`, `Right-Putamen`, `Right-Thalamus-Proper`, `Right-VentralDC`, 
                                        `Right-Cerebellum-Cortex`, `Right-Cerebellum-White-Matter`, na.rm = TRUE)) %>% ungroup()

#total hippocampus
volData <- volData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(totalHippocampus = sum(`Left-Hippocampus`,`Right-Hippocampus`)) %>% ungroup()

### Convert big variables to CM for plotting purposes
volData <- volData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(BrainSegVolNotVentSurf_cm = BrainSegVolNotVentSurf/1000, `Left-totalLateralVentricle_cm` = `Left-totalLateralVentricle`/1000, 
         `Right-totalLateralVentricle_cm` = `Right-totalLateralVentricle`/1000, `WM-hypointensities_cm` = `WM-hypointensities`/1000,
         `Left-Hippocampus_cm` = `Left-Hippocampus`/1000,`Right-Hippocampus_cm` = `Right-Hippocampus`/1000) %>% ungroup()

###Make dataframes longer
volDataLong<-volData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
lhThickLong<-lhThickData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
rhThickLong<-rhThickData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
lhGWconLong<-lhGWconData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
rhGWconLong<-rhGWconData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
ADsigDataLong <- ADsigDataWide %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
ADsigDataLong$Date <- as.numeric(ADsigDataLong$Date)
ADsigData_GWRLong <- ADsigData_GWRWide %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
ADsigData_GWRLong$Date <- as.numeric(ADsigData_GWRLong$Date)
FTDsigDataLong <- FTDsigDataWide %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
FTDsigDataLong$Date <- as.numeric(FTDsigDataLong$Date)
FTDsigData_GWRLong <- FTDsigData_GWRWide %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
FTDsigData_GWRLong$Date <- as.numeric(FTDsigData_GWRLong$Date)

##Add in session variable
volDataLong<-volDataLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
lhThickLong<-lhThickLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
rhThickLong<-rhThickLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
lhGWconLong<-lhGWconLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
rhGWconLong<-rhGWconLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
ADsigDataLong <- ADsigDataLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
ADsigData_GWRLong <- ADsigData_GWRLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
FTDsigDataLong <- FTDsigDataLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
FTDsigData_GWRLong <- FTDsigData_GWRLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))

##Can run this to check if dense_rank work properly there should be 6 sessions, 
table(volDataLong$session)
#table(ADsigDataLong$session)
table(lhGWconLong$session)

##Fix naming so combining data and plotting is easier
volDataLong<-volDataLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
lhThickLong<-lhThickLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
rhThickLong<-rhThickLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
lhGWconLong<-lhGWconLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
rhGWconLong<-rhGWconLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
ADsigDataLong <- ADsigDataLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
ADsigData_GWRLong <- ADsigData_GWRLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
FTDsigDataLong <- FTDsigDataLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
FTDsigData_GWRLong <- FTDsigData_GWRLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID))

lhThickLong<-lhThickLong %>% mutate(region = gsub("lh_", "Left-", region))
rhThickLong<-rhThickLong %>% mutate(region = gsub("rh_", "Right-", region))
lhGWconLong<-lhGWconLong %>% mutate(region = paste0("Left-", region))
rhGWconLong<-rhGWconLong %>% mutate(region = paste0("Right-", region))

###Combine all data into big DF
allFSdataLong<-bind_rows(list(vol=volDataLong,Thick=lhThickLong,Thick=rhThickLong,GWcon=lhGWconLong,GWcon=rhGWconLong, 
                              Thick = ADsigDataLong, Thick = FTDsigDataLong,GWcon = ADsigData_GWRLong, GWcon = FTDsigData_GWRLong), .id = "origData")
allFSdataLong<- allFSdataLong %>% mutate(ID = gsub(".ANAT", "", ID))
allFSdataLong<- allFSdataLong %>% mutate(region = gsub("_thickness", "", region))
#rename CortexVol
allFSdataLong<-allFSdataLong %>% mutate(region = gsub("rhCortexVol", "Right-CortexVol", region))
allFSdataLong<-allFSdataLong %>% mutate(region = gsub("lhCortexVol", "Left-CortexVol", region))

allFSdataLong <- allFSdataLong %>% filter(str_detect(ID,"",negate=T))

##Add in column for Dx
###Edited to remove IDs
allFSdataLong <- allFSdataLong %>% add_column(dx)

###Add in column for exclusion
demo2<-inner_join(demo,allFSdataLong %>% select(ID,dx) %>% unique())

###Edited to remove IDs
allFSdataLong$excluded <- ifelse(allFSdataLong$ID %in% c(""),1,0)

#add demographics
allFSdataLong<-left_join(allFSdataLong,demo) 

##Fix issue with group tibble instead of regular tibble
allFSdataLong<-as_tibble(allFSdataLong)

##Filter out repeated rows
allFSdataLong <- allFSdataLong %>% unique()

####Make variables to prepare for plotting

summaryTable<-allFSdataLong %>% group_by(ID, session,origData) %>% summarise(n=n())
numADNIsummaryTable<-allFSdataLong %>% filter(grepl("ADNI",scanID)) %>% group_by(ID) %>% filter(region == "Left-Hippocampus") %>% summarise(n=n())
nSesh<-allFSdataLong %>% group_by(ID,session) %>% summarise(n=n()) %>% summarise(nSessions = max(session)) ###How many sessions does each person have
#nSesh<-allFSdataLong %>% group_by(ID, template,session) %>% summarise(n=n()) %>% summarise(nSessions = n(session)) ###How many sessions does each person have
IDs_allData <- numADNIsummaryTable$ID


##Finish up reformatting here so you don't have to have a messy plotting script

plotData<-drop_na(allFSdataLong)
###Set up templates so that they make sense
plotData <- plotData %>% mutate(template = str_remove(template, "^[^_]*_"))

##Count number of scans present for each template and participant
numScansPerTemplatePre<-plotData %>% filter(region == "Left-Hippocampus") %>% group_by(ID,template) %>% summarise(numScans = n())

##Change template names if there is not already a template called inclusive
plotData<-plotData %>% group_by(ID) %>%
  mutate(hasCSx6_inclusive = if_else(any(template == "CSx6_inclusive"), 1, 0),hasADNI_inclusive = if_else(any(template == "ADNI_inclusive"), 1, 0)) %>%
  mutate(templateInclusive = case_when(
    template == "CSx6_16" & hasCSx6_inclusive == 0 ~ "CSx6_inclusive",
    template == "ADNI_2" & hasADNI_inclusive == 0 ~ "ADNI_inclusive",
    TRUE ~ template))

templateCountTable<-table(plotData$ID,plotData$templateInclusive)
numScansPerTemplatePost<-plotData %>% filter(region == "Left-Hippocampus") %>% group_by(ID,templateInclusive) %>% summarise(numScans = n())
numCSx6_inclusive<-numScansPerTemplatePost %>% filter(templateInclusive == "CSx6_inclusive")
##Replace the old template with the new template column
plotData<-plotData %>% select(-template) %>% rename(template = templateInclusive)


###Estimate time from baseline for each scan
plotData$Date <- ymd(plotData$Date)
plotData <- plotData %>% group_by(ID) %>% mutate(baselineDate = min(Date))
plotData<-plotData %>% mutate(yearsFromBaseline = as.numeric(Date - baselineDate)/(365.25))

plotData<-plotData %>% mutate(region = gsub("WM-hypointensities", "WMH", region)) 

###Remove bigVen templateInclusive
plotData<- plotData %>% filter(!grepl("bigVen",template)) %>% filter(!grepl("adStats",template)) 

#######Remove people who don't have as many sessions
##How many sessions per person
numSesh<-plotData %>% group_by(ID) %>% summarize(numSesh = max(session))
SixSeshIDs<-numSesh %>% filter(numSesh > 5) %>% select(ID) %>% as_vector()

plotData <- plotData %>% group_by(ID) %>% mutate(maxSesh = max(session)) %>% ungroup()

###Add in info about sessions and timing
plotData <- plotData %>% mutate(Visit = (session + 1) %/% 2)
plotData <- plotData %>% mutate(SessionWithinVisit = if_else(session %% 2 == 1, "test", "retest"))
#plotData <- plotData %>% mutate(newScanType = str_split(scanID, "_") %>% map_chr(~ tail(.x, 1))) %>% 
#  mutate(OrderWithinSession = match(str_to_lower(newScanType), letters)) %>% mutate(OrderWithinBreak = ifelse(OrderWithinSession > 4,yes = OrderWithinSession - 4,no = OrderWithinSession))
#Add in break info
beforeBreakScans<- c("0.8_CSx6_a","1.2_CSx6_a","1.0_CSx6_a","1.0_CSx6_b","1.0_CSx6_c","1.0_CSx6_d","0.9_CSx6_a","1.1_CSx6_a","1.0_ADNI_a")
plotData <- plotData %>% mutate(breakInfo = if_else(scanID %in% beforeBreakScans, "before", "after"))

###Clean up - remove atypical scans and NAs
regScans = c("0.8_CSx6_a","1.2_CSx6_a","1.0_CSx6_a","1.0_CSx6_b","1.0_CSx6_c","1.0_CSx6_d","0.9_CSx6_a","1.1_CSx6_a","1.0_ADNI_a", "0.8_CSx6_b","1.2_CSx6_b","1.0_CSx6_e","1.0_CSx6_f","1.0_CSx6_g","1.0_CSx6_h","0.9_CSx6_b","1.1_CSx6_b")
CS1mmScans <- regScans[c(3:6,12:15)]
CScrossResScans <-regScans[c(1:2,7:8,10:11,16:17)]
plotData <- plotData %>% filter(scanID %in% regScans)


#save out data for use here and other places
write_rds(plotData,"/path/to/data/ClusterLong/Data/ClusterLong_240803_allFSdataLong.rds")

################################
##### Generate Errors
################################

##Use if you just edit error generation below and don't want to rerun above
#plotData <- read_rds("/path/to/data/ClusterLong/Data/ClusterLong_240803_allFSdataLong.rds")

altPath="/path/"

##Source the functions that you use below
source(paste0(altPath,"Manuscripts/ClusterLong/Scripts/ClusterLong_Functions_240712.R"))

###Save out baseline test retest errors for each paper region 
regScans = c("0.8_CSx6_a","1.2_CSx6_a","1.0_CSx6_a","1.0_CSx6_b","1.0_CSx6_c","1.0_CSx6_d","0.9_CSx6_a","1.1_CSx6_a","1.0_ADNI_a", "0.8_CSx6_b","1.2_CSx6_b","1.0_CSx6_e","1.0_CSx6_f","1.0_CSx6_g","1.0_CSx6_h","0.9_CSx6_b","1.1_CSx6_b")
CS1mmScans <- regScans[c(3:6,12:15)]
CScrossResScans <-regScans[c(1:2,7:8,10:11,16:17)]

allTestRetest.IDs <- as_vector(read_csv(paste0("/path/to/data/ClusterLong/Data/allTestRetestID.csv")))
paperRegions<-read_csv("/path/to/data/ClusterLong/Data/paperRegions.csv")
plotRegions<-c("Left-ADsignature_CT","Right-ADsignature_CT","Left-FTDsignature_CT","Right-FTDsignature_CT", "Left-totalLateralVentricle","Right-totalLateralVentricle", "Left-totalBrainVolume","Right-totalBrainVolume", "WMH", "Left-Hippocampus","Right-Hippocampus", "allVentricles","totalHippocampus","BrainSegVolNotVentSurf","Right-MeanGWCon","Left-MeanThickness","Right-MeanThickness","Left-Amygdala","Right-Amygdala","Left-parahippocampal","Right-parahippocampal")
testRetestRegions <-unique(c(as_vector(paperRegions),plotRegions))

plotData.allTestRetest <- plotData %>% filter(ID %in% allTestRetest.IDs & dx != "MCInos") %>% select(-breakInfo,-SessionWithinVisit,-Visit,-maxSesh,-hasADNI_inclusive,-hasCSx6_inclusive)


#### Baseline Errors for all paper regions for test retest subs
# ADNI 1 and CS 1 2 4 8

###Generate long Errors and summary for use in plotting scripts
allMorph.SlopeError.allData<-tibble()
allMorph.SlopeError.SummaryData<-tibble()

for(brainVar in testRetestRegions){
  print(brainVar)
  for(temp in c("ADNI_1","CSx6_1","CSx6_2","CSx6_4","CSx6_8")){
    ##Get data types to separate GWcon and thickness if needed
    dataTypes<-plotData.allTestRetest %>% filter(region == brainVar) %>% ungroup() %>% select(origData) %>% unique() %>% as_vector()
    for(type in dataTypes){
      errorData<-plotData.allTestRetest %>% filter(origData == type & region == brainVar & str_detect(template,temp))
      template.test <-paste0(temp,"_test")
      template.retest <-paste0(temp,"_retest")
      errorList<-slopeErrors_LM_simple(brainVar = brainVar,df = errorData,templateName_test = template.test,templateName_retest = template.retest,sessions_test = c(1,3,5),sessions_retest = c(2,4,6),label = temp)
      errorSum<-bind_cols(errorList$summary,tibble(origData = type))
      errorAll<-bind_cols(errorList$df,tibble(origData = type))
      allMorph.SlopeError.SummaryData<-bind_rows(allMorph.SlopeError.SummaryData,errorSum)
      allMorph.SlopeError.allData<-bind_rows(allMorph.SlopeError.allData,errorAll)
    }
  }
}

#Save summary for future use if you don't want to rerun above
write_csv(allMorph.SlopeError.SummaryData,paste0(altPath,"/Manuscripts/ClusterLong/Data/allMorph.SlopeError.SummaryData.240803.csv"))
write_csv(allMorph.SlopeError.allData,paste0(altPath,"/Manuscripts/ClusterLong/Data/allMorph.SlopeError.AllData.240803.csv"))




##################################################
####### Cross-sectional errors for test retest 
####### Goal: get true estimates of baseline error consistent with past papers, not biased by long pipeline
##################################################

##Read in cross sectional data



volData<-read_csv("/path/to/data/ClusterLong/Data/ClusterOrig_cs_FS_aseg_vol_allCrossSectional_240710.csv")
lhThickData<-read_csv("/path/to/data/ClusterLong/Data/ClusterOrig_cs_FS_aparc_thk_lh_allCrossSectional_240710.csv")
rhThickData<-read_csv("/path/to/data/ClusterLong/Data/ClusterOrig_cs_FS_aparc_thk_rh_allCrossSectional_240710.csv")
lhGWconData<-read_csv("/path/to/data/ClusterLong/Data/ClusterOrig_cs_FS_aseg_w-g.pct.stats_lh_allCrossSectional_240710.csv")
rhGWconData<-read_csv("/path/to/data/ClusterLong/Data/ClusterOrig_cs_FS_aseg_w-g.pct.stats_rh_allCrossSectional_240710.csv")

ADsigData<-read_csv("/path/to/data/ClusterLong/Data/ADsignatureMask_meanCT_crossSectional_240803.csv",col_names = c("ID","longScanID","hemisphere","value"))
ADsigData_GWR<-read_csv("/path/to/data/ClusterLong/Data/ADsignatureMask_meanGWR_crossSectional_240803.csv",col_names = c("ID","longScanID","hemisphere","value"))
FTDsigData<-read_csv("/path/to/data/ClusterLong/Data/FTDsignatureMask_meanCT_crossSectional_240803.csv",col_names = c("ID","longScanID","hemisphere","value"))
FTDsigData_GWR<-read_csv("/path/to/data/ClusterLong/Data/FTDsignatureMask_meanGWR_crossSectional_240803.csv",col_names = c("ID","longScanID","hemisphere","value"))

##Grab data from scan names
volData<-volData %>% separate(c("Date", "ID","scanID"),col = "Measure:volume",sep = "_",extra = "merge",convert = T)  
lhThickData <- lhThickData %>% separate(c("Date", "ID","scanID"),col = "lh.aparc.thickness",sep = "_",extra = "merge",convert = T)
rhThickData <- rhThickData %>% separate(c("Date", "ID","scanID"),col = "rh.aparc.thickness",sep = "_",extra = "merge",convert = T)
lhGWconData <- lhGWconData %>% separate(c("Date", "ID","scanID"),col = "Measure:mean",sep = "_",extra = "merge",convert = T)
rhGWconData <- rhGWconData %>% separate(c("Date", "ID","scanID"),col = "Measure:mean",sep = "_",extra = "merge",convert = T)
ADsigData <- ADsigData %>% separate(c("Date", "ID","scanID"),col = "longScanID",sep = "_",extra = "merge",convert = T)
ADsigData_GWR <- ADsigData_GWR %>% separate(c("Date", "ID","scanID"),col = "longScanID",sep = "_",extra = "merge",convert = T)
FTDsigData <- FTDsigData %>% separate(c("Date", "ID","scanID"),col = "longScanID",sep = "_",extra = "merge",convert = T)
FTDsigData_GWR <- FTDsigData_GWR %>% separate(c("Date", "ID","scanID"),col = "longScanID",sep = "_",extra = "merge",convert = T)

##Clean up ADsig
ADsigDataWide<-ADsigData %>% pivot_wider(names_from = hemisphere,values_from = value,names_glue = "{hemisphere}-ADsignature") %>% rename(`Left-ADsignature_CT` = `lh-ADsignature`, `Right-ADsignature_CT` = `rh-ADsignature`) %>% 
  mutate(meanADsignature_CT = (`Left-ADsignature_CT` + `Right-ADsignature_CT`)/2)
ADsigData_GWRWide<-ADsigData_GWR %>% pivot_wider(names_from = hemisphere,values_from = value,names_glue = "{hemisphere}-ADsignature") %>% rename(`Left-ADsignature_GWR` = `lh-ADsignature`, `Right-ADsignature_GWR` = `rh-ADsignature`) %>% 
  mutate(meanADsignature_GWR = (`Left-ADsignature_GWR` + `Right-ADsignature_GWR`)/2)
FTDsigDataWide<-FTDsigData %>% pivot_wider(names_from = hemisphere,values_from = value,names_glue = "{hemisphere}-FTDsignature") %>% rename(`Left-FTDsignature_CT` = `lh-FTDsignature`, `Right-FTDsignature_CT` = `rh-FTDsignature`) %>% 
  mutate(meanFTDsignature_CT = (`Left-FTDsignature_CT` + `Right-FTDsignature_CT`)/2)
FTDsigData_GWRWide<-FTDsigData_GWR %>% pivot_wider(names_from = hemisphere,values_from = value,names_glue = "{hemisphere}-FTDsignature") %>% rename(`Left-FTDsignature_GWR` = `lh-FTDsignature`, `Right-FTDsignature_GWR` = `rh-FTDsignature`) %>% 
  mutate(meanFTDsignature_GWR = (`Left-FTDsignature_GWR` + `Right-FTDsignature_GWR`)/2)

#ADD mean GWcon for each Hemispere
lhGWconData <- lhGWconData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(`MeanGWCon` = mean(c_across(5:36), na.rm = TRUE)) %>% ungroup()
rhGWconData <- rhGWconData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(`MeanGWCon` = mean(c_across(5:36), na.rm = TRUE)) %>% ungroup()

###ADD in combined ventricle variables for each Hemisphere
volData <- volData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(allVentricles = sum(`Left-Inf-Lat-Vent`,`Left-Lateral-Ventricle`,`Right-Inf-Lat-Vent`,`Right-Lateral-Ventricle`,`3rd-Ventricle`,`4th-Ventricle`,`5th-Ventricle`, na.rm = TRUE)) %>% ungroup()
volData <- volData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(`Left-totalLateralVentricle` = sum(`Left-Inf-Lat-Vent`,`Left-Lateral-Ventricle`, na.rm = TRUE), `Right-totalLateralVentricle` = sum(`Right-Inf-Lat-Vent`,`Right-Lateral-Ventricle`, na.rm = TRUE)) %>% ungroup()

##Add in left and right total brain volume measures for later plotting
volData <- volData %>% rowwise() %>%
  mutate(`Left-totalBrainVolume` = sum(`lhCortexVol`, `lhCerebralWhiteMatterVol`, `Left-Accumbens-area`, `Left-Amygdala`, `Left-Caudate`, 
                                       `Left-Hippocampus`, `Left-Pallidum`, `Left-Putamen`, `Left-Thalamus-Proper`, `Left-VentralDC`, 
                                       `Left-Cerebellum-Cortex`, `Left-Cerebellum-White-Matter`, na.rm = TRUE),
         `Right-totalBrainVolume` = sum(`rhCortexVol`, `rhCerebralWhiteMatterVol`, `Right-Accumbens-area`, `Right-Amygdala`, `Right-Caudate`, 
                                        `Right-Hippocampus`, `Right-Pallidum`, `Right-Putamen`, `Right-Thalamus-Proper`, `Right-VentralDC`, 
                                        `Right-Cerebellum-Cortex`, `Right-Cerebellum-White-Matter`, na.rm = TRUE)) %>% ungroup()

#total hippocampus
volData <- volData %>% rowwise() %>%  # Apply the following operations row by row
  mutate(totalHippocampus = sum(`Left-Hippocampus`,`Right-Hippocampus`)) %>% ungroup()

###Make dataframes longer
volDataLong<-volData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
lhThickLong<-lhThickData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
rhThickLong<-rhThickData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
lhGWconLong<-lhGWconData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
rhGWconLong<-rhGWconData %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
ADsigDataLong <- ADsigDataWide %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
ADsigDataLong$Date <- as.numeric(ADsigDataLong$Date)
ADsigData_GWRLong <- ADsigData_GWRWide %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
ADsigData_GWRLong$Date <- as.numeric(ADsigData_GWRLong$Date)
FTDsigDataLong <- FTDsigDataWide %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
FTDsigDataLong$Date <- as.numeric(FTDsigDataLong$Date)
FTDsigData_GWRLong <- FTDsigData_GWRWide %>% pivot_longer(!(Date:scanID),names_to = "region", values_to = "value")
FTDsigData_GWRLong$Date <- as.numeric(FTDsigData_GWRLong$Date)

##Add in session variable
volDataLong<-volDataLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
lhThickLong<-lhThickLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
rhThickLong<-rhThickLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
lhGWconLong<-lhGWconLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
rhGWconLong<-rhGWconLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
ADsigDataLong <- ADsigDataLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
ADsigData_GWRLong <- ADsigData_GWRLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
FTDsigDataLong <- FTDsigDataLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))
FTDsigData_GWRLong <- FTDsigData_GWRLong %>% arrange(Date,ID) %>% group_by(ID) %>% mutate(session=dense_rank(Date))

##Can run this to check if dense_rank work properly there should be 6 sessions, 
table(volDataLong$session)
#table(ADsigDataLong$session)
table(lhGWconLong$session)

##Fix naming so combining data and plotting is easier
volDataLong<-volDataLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
lhThickLong<-lhThickLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
rhThickLong<-rhThickLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
lhGWconLong<-lhGWconLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
rhGWconLong<-rhGWconLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
ADsigDataLong <- ADsigDataLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
ADsigData_GWRLong <- ADsigData_GWRLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
FTDsigDataLong <- FTDsigDataLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID)) 
FTDsigData_GWRLong <- FTDsigData_GWRLong %>% separate(c("scanID","template"),col = "scanID",sep = ".long.",extra = "merge") %>% mutate(scanID = gsub("_ADNI", "_ADNI_a", scanID))

lhThickLong<-lhThickLong %>% mutate(region = gsub("lh_", "Left-", region))
rhThickLong<-rhThickLong %>% mutate(region = gsub("rh_", "Right-", region))
lhGWconLong<-lhGWconLong %>% mutate(region = paste0("Left-", region))
rhGWconLong<-rhGWconLong %>% mutate(region = paste0("Right-", region))

###Combine all data into big DF
allFSdataLong<-bind_rows(list(vol=volDataLong,Thick=lhThickLong,Thick=rhThickLong,GWcon=lhGWconLong,GWcon=rhGWconLong,
                              Thick = ADsigDataLong, Thick = FTDsigDataLong,GWcon = ADsigData_GWRLong, GWcon = FTDsigData_GWRLong), .id = "origData")
allFSdataLong<- allFSdataLong %>% mutate(ID = gsub(".ANAT", "", ID))
allFSdataLong<- allFSdataLong %>% mutate(region = gsub("_thickness", "", region))
#rename CortexVol
allFSdataLong<-allFSdataLong %>% mutate(region = gsub("rhCortexVol", "Right-CortexVol", region))
allFSdataLong<-allFSdataLong %>% mutate(region = gsub("lhCortexVol", "Left-CortexVol", region))

allFSdataLong <- allFSdataLong %>% filter(str_detect(ID,"",negate=T))


##Add in column for Dx
###Edited to remove IDs
allFSdataLong <- allFSdataLong %>% add_column(dx)

###Add in column for exclusion
demo2<-inner_join(demo,allFSdataLong %>% select(ID,dx) %>% unique())

###Edited to remove IDs
allFSdataLong$excluded <- ifelse(allFSdataLong$ID %in% c(""),1,0)

#add demographics
allFSdataLong<-left_join(allFSdataLong,demo) 

##Fix issue with group tibble instead of regular tibble
allFSdataLong<-as_tibble(allFSdataLong)

##Filter out repeated rows
allFSdataLong <- allFSdataLong %>% select(-template) %>% unique()

plotData<-drop_na(allFSdataLong)


###Estimate time from baseline for each scan
plotData$Date <- ymd(plotData$Date)
plotData <- plotData %>% group_by(ID) %>% mutate(baselineDate = min(Date))
plotData<-plotData %>% mutate(yearsFromBaseline = as.numeric(Date - baselineDate)/(365.25))

plotData<-plotData %>% mutate(region = gsub("WM-hypointensities", "WMH", region)) %>% mutate(region = gsub("WM-hypointensities_cm", "WMH_cm", region)) 

###Remove bigVen templateInclusive 
plotData <- plotData %>% group_by(ID) %>% mutate(maxSesh = max(session)) %>% ungroup()

###Add in info about sessions and timing
plotData <- plotData %>% mutate(Visit = (session + 1) %/% 2)
plotData <- plotData %>% mutate(SessionWithinVisit = if_else(session %% 2 == 1, "test", "retest"))
#plotData <- plotData %>% mutate(newScanType = str_split(scanID, "_") %>% map_chr(~ tail(.x, 1))) %>% 
#  mutate(OrderWithinSession = match(str_to_lower(newScanType), letters)) %>% mutate(OrderWithinBreak = ifelse(OrderWithinSession > 4,yes = OrderWithinSession - 4,no = OrderWithinSession))
#Add in break info
beforeBreakScans<- c("0.8_CSx6_a","1.2_CSx6_a","1.0_CSx6_a","1.0_CSx6_b","1.0_CSx6_c","1.0_CSx6_d","0.9_CSx6_a","1.1_CSx6_a","1.0_ADNI_a")
plotData <- plotData %>% mutate(breakInfo = if_else(scanID %in% beforeBreakScans, "before", "after"))

###Clean up - remove atypical scans and NAs
regScans = c("0.8_CSx6_a","1.2_CSx6_a","1.0_CSx6_a","1.0_CSx6_b","1.0_CSx6_c","1.0_CSx6_d","0.9_CSx6_a","1.1_CSx6_a","1.0_ADNI_a", "0.8_CSx6_b","1.2_CSx6_b","1.0_CSx6_e","1.0_CSx6_f","1.0_CSx6_g","1.0_CSx6_h","0.9_CSx6_b","1.1_CSx6_b")
CS1mmScans <- regScans[c(3:6,12:15)]
CScrossResScans <-regScans[c(1:2,7:8,10:11,16:17)]
plotData <- plotData %>% filter(scanID %in% regScans)


#save out data for use here and other places
write_rds(plotData,"/path/to/data/ClusterLong/Data/ClusterCrossSectional_240803_allFSdataLong.rds")

## Use the line below if you just want to edit error script and you can skip all the data manipulation above
plotData <- read_rds("/path/to/data/ClusterLong/Data/ClusterCrossSectional_240803_allFSdataLong.rds")
  
###Save out baseline test retest errors for each paper region 

allTestRetest.IDs <- as_vector(read_csv(paste0("/path/to/data/ClusterLong/Data/allTestRetestID.csv")))
paperRegions<-read_csv("/path/to/data/ClusterLong/Data/paperRegions.csv")


plotData.allTestRetest <- plotData %>% filter(ID %in% allTestRetest.IDs, session < 3, dx != "MCInos" )

#### Baseline Errors for all paper regions for test retest subs
# ADNI 1 and CS 1 2 4 8

### CS errors
baselineTestRetestErrorData<-tibble()
for(numScans in c(1,2,4,8)){
  for(var in testRetestRegions){
    
    errorData <- plotData.allTestRetest %>% filter(region == var) %>% filter(scanID %in% CS1mmScans[1:numScans])
    errorData <- errorData %>% group_by(ID,dx, SessionWithinVisit,region,origData) %>% summarise(value = mean(value)) %>%
      ungroup() %>% mutate(SessionWithinVisit = str_to_title(SessionWithinVisit))
    errorData$brainVar<-errorData$region
    #errorData <- errorData %>% separate(c("Hemi", "region"),col = "region",sep = "-")  %>% distinct()  
    errorDataWide <- errorData %>% pivot_wider(names_from = SessionWithinVisit,values_from = value)
    errorDataWide <- errorDataWide %>% mutate(meanVal = (Test + Retest)/2)
    errorDataWide <- errorDataWide %>% mutate(error.raw = abs(Test - Retest), error.percent = abs(Test - Retest)/meanVal*100)
    errorDataWide$template <- paste0(numScans," CS")
    baselineTestRetestErrorData<-bind_rows(errorDataWide,baselineTestRetestErrorData)
  }
}

### ADNI errors
numScans=1
for(var in testRetestRegions){
  errorData <- plotData.allTestRetest %>% filter(region == var) %>% filter(scanID == "1.0_ADNI_a")
  errorData <- errorData %>% group_by(ID,dx, SessionWithinVisit,region,origData) %>% summarise(value = mean(value)) %>%
    ungroup() %>% mutate(SessionWithinVisit = str_to_title(SessionWithinVisit))
  errorData$brainVar<-errorData$region
  #errorData <- errorData %>% separate(c("Hemi", "region"),col = "region",sep = "-")  %>% distinct()  
  errorDataWide <- errorData %>% pivot_wider(names_from = SessionWithinVisit,values_from = value)
  errorDataWide <- errorDataWide %>% mutate(meanVal = (Test + Retest)/2)
  errorDataWide <- errorDataWide %>% mutate(error.raw = abs(Test - Retest), error.percent = abs(Test - Retest)/meanVal*100)
  errorDataWide$template <- paste0(numScans," ADNI")
  baselineTestRetestErrorData<-bind_rows(errorDataWide,baselineTestRetestErrorData)
}

baselineTestRetestErrorData.sum <- baselineTestRetestErrorData %>% group_by(brainVar,origData,template) %>% summarise(mean.error.raw = mean(error.raw),mean.error.percent = mean(error.percent),numIDs = n(),sd.error.raw = sd(error.raw),sd.error.percent = sd(error.percent))

write_rds(baselineTestRetestErrorData,"/path/to/data/ClusterLong/Data/baselineTestRetest.CrossSectional.participantErrors_240803.rds")

write_rds(baselineTestRetestErrorData.sum,"/path/to/data/ClusterLong/Data/baselineTestRetest.CrossSectional.ErrorSummary_240803.rds")



###################################
### Blood Biomarker Data

bloodData<-read_csv("/path/to/data/ClusterLong/Data/ADBloodBiomarkers_reformatted_240814.csv")

bloodData<-bloodData %>% mutate(BloodDrawDate = mdy(BloodDrawDate))

bloodData_mostRecent <- bloodData %>% arrange(MRI_Date,PIDN) %>% group_by(PIDN) %>% filter(BloodDrawDate==max(BloodDrawDate))

bloodData_mostRecent_TR <- bloodData_mostRecent %>% filter(PIDN %in% SixSeshIDs)

