######Application
rm(list=ls())

######1.NHANES data cleaning
library(Hmisc)

datahtn = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES\ Materials/NHANES/Htn/BPX_J.XPT')

###1.1 DEMOGRAPHICS ###
data18 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES\ Materials/NHANES/Demographics/DEMO_J.XPT')

data_LONG = merge(data18, datahtn, by = 'seqn', all.x = TRUE, all.y = TRUE)
data_SHORT = data_LONG[data_LONG$ridageyr>=18,]
data_SHORT = data_SHORT[data_SHORT$wtmec2yr != 0,]

data_SHORT$sys <- rowMeans(data_SHORT[,c('bpxsy1', 'bpxsy2','bpxsy3','bpxsy4')], na.rm=TRUE)
data_SHORT$dias <- rowMeans(data_SHORT[,c('bpxdi1', 'bpxdi2','bpxdi4','bpxdi4')], na.rm=TRUE)
data_SHORT$htn <- ifelse(data_SHORT$sys>=130 | data_SHORT$dias>=80, 1, 0)
data_SHORT = data_SHORT[!is.na(data_SHORT$htn),]
sum(data_SHORT$htn*data_SHORT$wtmec2yr, na.rm = T)/sum(data_SHORT$wtmec2yr) #0.3918549

data2 = sasxport.get('C:/Users/zgh98/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Demographics/DEMO_B.XPT')
data4 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Demographics/DEMO_C.XPT')
data6 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Demographics/DEMO_D.XPT')
data8 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Demographics/DEMO_E.XPT')
data10 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Demographics/DEMO_F.XPT')
data12 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Demographics/DEMO_G.XPT')
data14 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Demographics/DEMO_H.XPT')
data16 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Demographics/DEMO_I.XPT')
data18 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Demographics/DEMO_J.XPT')
ALL_NAMES = Reduce(intersect, list(names(data2), names(data4), names(data6), names(data8), names(data10), names(data12), names(data14), names(data16), names(data18)))
data2 = data2[,ALL_NAMES]
data4 = data4[,ALL_NAMES]
data6 = data6[,ALL_NAMES]
data8 = data8[,ALL_NAMES]
data10 = data10[,ALL_NAMES]
data12 = data12[,ALL_NAMES]
data14 = data14[,ALL_NAMES]
data16 = data16[,ALL_NAMES]
data18 = data18[,ALL_NAMES]
dataDEMOG = rbind(  data.frame(data18, Year = '2017-2018'),
                    data.frame(data16, Year = '2015-2016'), 
                    data.frame(data14, Year = '2013-2014'),
                    data.frame(data12, Year = '2011-2012'),
                    data.frame(data10, Year = '2009-2010'),
                    data.frame(data8, Year = '2007-2008'),
                    data.frame(data6, Year = '2005-2006'),
                    data.frame(data4, Year = '2003-2004'),
                    data.frame(data2, Year = '2001-2002'))

###1.2 BMI ###
data2 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/BMI/BMX_B.XPT')
data4 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/BMI/BMX_C.XPT')
data6 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/BMI/BMX_D.XPT')
data8 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/BMI/BMX_E.XPT')
data10 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/BMI/BMX_F.XPT')
data12 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/BMI/BMX_G.XPT')
data14 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/BMI/BMX_H.XPT')
data16 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/BMI/BMX_I.XPT')
data18 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/BMI/BMX_J.XPT')
ALL_NAMES = Reduce(intersect, list(names(data2), names(data4), names(data6), names(data8), names(data10), names(data12), names(data14), names(data16), names(data18)))
data2 = data2[,ALL_NAMES]
data4 = data4[,ALL_NAMES]
data6 = data6[,ALL_NAMES]
data8 = data8[,ALL_NAMES]
data10 = data10[,ALL_NAMES]
data12 = data12[,ALL_NAMES]
data14 = data14[,ALL_NAMES]
data16 = data16[,ALL_NAMES]
data18 = data18[,ALL_NAMES]
dataBMI = rbind(  data.frame(data18, Year = '2017-2018'),
                  data.frame(data16, Year = '2015-2016'), 
                  data.frame(data14, Year = '2013-2014'),
                  data.frame(data12, Year = '2011-2012'),
                  data.frame(data10, Year = '2009-2010'),
                  data.frame(data8, Year = '2007-2008'),
                  data.frame(data6, Year = '2005-2006'),
                  data.frame(data4, Year = '2003-2004'),
                  data.frame(data2, Year = '2001-2002'))

###1.3 SMOKING ###
data6 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Smoking/SMQ_D.XPT')
data8 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Smoking/SMQ_E.XPT')
data10 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Smoking/SMQ_F.XPT')
data12 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Smoking/SMQ_G.XPT')
data14 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Smoking/SMQ_H.XPT')
data16 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Smoking/SMQ_I.XPT')
data18 = sasxport.get('/Users/zhangguanghao/Dropbox/GuanghaoZhangResearchIdeas/Code/Application/NHANESdata/NHANES Materials/NHANES/Smoking/SMQ_J.XPT')
ALL_NAMES = Reduce(intersect, list(names(data6), names(data8), names(data10), names(data12), names(data14), names(data16), names(data18)))
data6 = data6[,ALL_NAMES]
data8 = data8[,ALL_NAMES]
data10 = data10[,ALL_NAMES]
data12 = data12[,ALL_NAMES]
data14 = data14[,ALL_NAMES]
data16 = data16[,ALL_NAMES]
data18 = data18[,ALL_NAMES]
dataSMOKE = rbind(  data.frame(data18, Year = '2017-2018'),
                    data.frame(data16, Year = '2015-2016'), 
                    data.frame(data14, Year = '2013-2014'),
                    data.frame(data12, Year = '2011-2012'),
                    data.frame(data10, Year = '2009-2010'),
                    data.frame(data8, Year = '2007-2008'),
                    data.frame(data6, Year = '2005-2006'))

###1.4 Merge Data ###
data_LONG = merge(dataDEMOG, subset(dataBMI, select = c(-Year)), by = 'seqn', all.x = TRUE, all.y = TRUE)
data_LONG = merge(data_LONG, subset(dataSMOKE, select = c(-Year)), by = 'seqn', all.x = TRUE, all.y = TRUE)

VARS_TO_SAVE = c('seqn', 'ridageyr', 'riagendr', 'wtint2yr','wtmec2yr', 'bmxbmi','smq040','smq020', 'Year', 'ridreth1')
data_SHORT = data_LONG[data_LONG$Year == '2017-2018' & data_LONG$ridageyr>=18,VARS_TO_SAVE]
data_SHORT = data_SHORT[data_SHORT$wtmec2yr != 0,]
data_SHORT$SAMP_NHANES = 1/data_SHORT$wtmec2yr  #these weights correspond to people with interviews and MEC example
data_SHORT$wtmec2yr = length(data_SHORT$wtmec2yr)*data_SHORT$wtmec2yr/sum(data_SHORT$wtmec2yr)

data_SHORT$riagendr = as.character(data_SHORT$riagendr)
data_SHORT$riagendr[data_SHORT$riagendr == 1] = 'Male'
data_SHORT$riagendr[data_SHORT$riagendr == 2] = 'Female'

data_SHORT$smq040 = as.character(data_SHORT$smq040)
data_SHORT$smq040[data_SHORT$smq040 %in% c('1', '2')] = 'Yes'
data_SHORT$smq040[data_SHORT$smq040 == '3'] = 'No'
data_SHORT$smq040[data_SHORT$smq040 == 9] = NA
data_SHORT$smq040[data_SHORT$smq040 == 7] = NA
data_SHORT$smq020 = as.character(data_SHORT$smq020)
data_SHORT$smq020[data_SHORT$smq020 == 1] = 'Yes'
data_SHORT$smq020[data_SHORT$smq020 == 2] = 'No'
data_SHORT$smq020[data_SHORT$smq020 == 9] = NA
data_SHORT$smq020[data_SHORT$smq020 == 7] = NA

data_SHORT$bmxbmi[data_SHORT$bmxbmi>86.2] = NA

data_SHORT$ridreth1 = as.character(data_SHORT$ridreth1)
data_SHORT$ridreth1[data_SHORT$ridreth1 %in% c('1', '2')] = 'Hispanic'
data_SHORT$ridreth1[data_SHORT$ridreth1 == '3'] = 'Non-Hispanic White'
data_SHORT$ridreth1[data_SHORT$ridreth1 == '4'] = 'Non-Hispanic Black'
data_SHORT$ridreth1[data_SHORT$ridreth1 == '5'] = 'Other'

#040 = do you now smoke cigarettes
#020 = have you smoked at least 100 cigarettes in life
table(data_SHORT$smq040,data_SHORT$smq020, useNA = 'ifany' )
data_SHORT$Smoking = rep(NA,length(data_SHORT$smq040))
data_SHORT$Smoking[data_SHORT$smq020 == 'No'] = 'Never'
data_SHORT$Smoking[data_SHORT$smq020 == 'Yes' & data_SHORT$smq040 == 'No'] = 'Former'
data_SHORT$Smoking[data_SHORT$smq020 == 'Yes' & data_SHORT$smq040 == 'Yes'] = 'Current'

#Age groupings
data_SHORT$NHANES_AGECAT = data.table::between(data_SHORT$ridageyr, 0,5)+
  2*data.table::between(data_SHORT$ridageyr, 6,11)+
  3*data.table::between(data_SHORT$ridageyr, 12,19)+
  4*data.table::between(data_SHORT$ridageyr, 20,39)+
  5*data.table::between(data_SHORT$ridageyr, 40,59)+
  6*data.table::between(data_SHORT$ridageyr, 60,150)

#BMI groupings
data_SHORT$NHANES_BMICAT = data.table::between(data_SHORT$bmxbmi, 0,18.499)+
  2*data.table::between(data_SHORT$bmxbmi, 18.5,24.999)+
  3*data.table::between(data_SHORT$bmxbmi, 25.0,29.999)+
  4*data.table::between(data_SHORT$bmxbmi, 30.0,120)

nh <- subset(data_SHORT, select = c(-smq040, -smq020, -Year, -wtint2yr, -wtmec2yr))
colnames(nh) <- c('id', 'age','sex','bmi','race','pr','smk','agec','bmic')
write.table(nh,"nhanes.txt",sep="\t",row.names=FALSE)
nh <- read.table("nhanes.txt", header = TRUE, sep = "\t", fill = TRUE)

nh$agec = data.table::between(nh$age, 18,39)+
  2*data.table::between(nh$age, 40,59)+
  3*data.table::between(nh$age, 60,150)

nh$racec[nh$race=='Non-Hispanic White'] <- 1
nh$racec[nh$race=='Non-Hispanic Black'] <- 2
nh$racec[nh$race=='Other'] <- 3
nh$racec[nh$race=='Hispanic'] <- 3

######2. MGI data
mgi <- read.table("mgi.txt", header = TRUE, sep = "\t", fill = TRUE)
mgi$agec = data.table::between(mgi$age, 18,39)+
  2*data.table::between(mgi$age, 40,59)+
  3*data.table::between(mgi$age, 60,150)

mgi$racec[mgi$race=='Non-Hispanic White'] <- 1
mgi$racec[mgi$race=='Non-Hispanic Black'] <- 2
mgi$racec[mgi$race=='Other'] <- 3
mgi$racec[mgi$race=='Hispanic'] <- 3

######3. Combine NHANES and MGI
TEMPDAT_ALL = data.frame(dataset = c(rep('NHANES', length(nh[,1])), rep('MGI', length(mgi[,1]))),
                         agec = c(nh$agec, mgi$agec),
                         sex = c(nh$sex, mgi$sex),
                         racec = c(nh$racec, mgi$racec),
                         smk = c(nh$smk, mgi$smk))

######4. Obtain lambda1
library(betareg)
selection_NHANES <- betareg(pr~factor(agec) + factor(racec) + factor(sex), data = nh,maxit = 10000)

mgiselect = glm(as.numeric(dataset == 'MGI')~factor(agec) + factor(racec) + factor(sex), 
                data = TEMPDAT_ALL, family = 'binomial')

p_Sext = predict(selection_NHANES, newdata = TEMPDAT_ALL[TEMPDAT_ALL$dataset == 'MGI',], type = 'response' )
p_MGI = predict(mgiselect, newdata = TEMPDAT_ALL[TEMPDAT_ALL$dataset == 'MGI',] , type = 'response')
SELECT_NHANES = p_Sext*(p_MGI/(1-p_MGI))

mgi$lambda1 <- SELECT_NHANES 
write.table(mgi,"mgi_cont.txt",sep="\t",row.names=FALSE)
write.table(mgi,"mgi_cat.txt",sep="\t",row.names=FALSE)
