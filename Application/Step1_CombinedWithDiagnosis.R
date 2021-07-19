mgi_demo <- read.table("HPI_4605_Demographics.txt", header = TRUE, sep = "\t", fill = TRUE)

mgi_smk <- read.table("HPI_4605_ClaritySocialHistory.txt", header = TRUE, sep = "\t", fill = TRUE)
mgi_smk <- mgi_smk[order(mgi_smk$Deid_ID, -abs(mgi_smk$DaysSinceBirth)), ] #get the latest smoking data
mgi_smk <- mgi_smk[!duplicated(mgi_smk$Deid_ID), ]

mgi_bmi <- read.table("HPI_4605_Anthropometrics.txt", header = TRUE, sep = "\t", fill = TRUE)
mgi_bmi <- mgi_bmi[order(mgi_bmi$Deid_ID, -abs(mgi_bmi$DaysSinceBirth)), ] #get the latest bmi data
mgi_bmi <- mgi_bmi[!duplicated(mgi_bmi$Deid_ID), ]

mgi = merge(mgi_demo, subset(mgi_smk, select = c(Deid_ID, SmokingStatus)), by = 'Deid_ID', all.x = TRUE, all.y = TRUE)
mgi = merge(mgi, subset(mgi_bmi, select = c(Deid_ID, BMI)), by = 'Deid_ID', all.x = TRUE, all.y = TRUE)

VARS_TO_SAVE = c('Deid_ID', 'Age', 'Sex', 'RaceName','EthnicityName', 'SmokingStatus','BMI')
mgi = mgi[mgi$Age>=18,VARS_TO_SAVE]

mgi$Sex = as.character(mgi$Sex)
mgi$Sex[mgi$Sex == 'M'] = 'Male'
mgi$Sex[mgi$Sex == 'F'] = 'Female'
mgi$Sex[mgi$Sex == 'U'] = NA

mgi$BMI[mgi$BMI>185.5] = NA

mgi$RaceName = as.character(mgi$RaceName)
mgi$RaceName[mgi$EthnicityName == 'Hispanic or Latino'] = 'Hispanic'
mgi$RaceName[mgi$RaceName == 'Caucasian' & mgi$EthnicityName != 'Hispanic or Latino'] = 'Non-Hispanic White'
mgi$RaceName[mgi$RaceName == 'African American' & mgi$EthnicityName != 'Hispanic or Latino'] = 'Non-Hispanic Black'
mgi$RaceName[mgi$RaceName != 'Hispanic' & mgi$RaceName != 'Non-Hispanic White' & mgi$RaceName != 'Non-Hispanic Black'] = 'Other'

mgi$SmokingStatus[mgi$SmokingStatus == 'Unknown'] = NA

mgi <- subset(mgi, select = c(-EthnicityName))

#Age groupings
mgi$agecat = data.table::between(mgi$Age, 0,5)+
  2*data.table::between(mgi$Age, 6,11)+
  3*data.table::between(mgi$Age, 12,19)+
  4*data.table::between(mgi$Age, 20,39)+
  5*data.table::between(mgi$Age, 40,59)+
  6*data.table::between(mgi$Age, 60,150)

#BMI groupings
mgi$bmicat = data.table::between(mgi$BMI, 0,18.499)+
  2*data.table::between(mgi$BMI, 18.5,24.999)+
  3*data.table::between(mgi$BMI, 25.0,29.999)+
  4*data.table::between(mgi$BMI, 30.0,120)

colnames(mgi) <- c('id','age','sex','race','smk','bmi','agec','bmic')

cc <- read.table("diagnosis.txt", header = TRUE, sep = "\t", fill = TRUE)
cc$y <- rep(1, length(cc$id))
cc <- subset(cc, select = c(id, y))

mgi.full <- read.table("mgi.txt", header = TRUE, sep = "\t", fill = TRUE)
mgi.full <- merge(cc, mgi.full, all.y = TRUE)
mgi.full$y[is.na(mgi.full$y)] <- 0

bd <- read.table("procedure.txt", header = TRUE, sep = "\t", fill = TRUE)
bd$iftest <- rep(1, length(bd$id))
bd <- subset(bd, select = c(id, iftest))

mgi.full <- merge(bd, mgi.full, all.y = TRUE)
mgi.full$iftest[is.na(mgi.full$iftest)] <- 0
mgi.full$y[mgi.full$iftest==0] <- -1 #we use -1 to represent NA

table(mgi.full$y)
table(mgi.full$iftest)

mgi <- mgi[which(mgi$y>-1),]
mgi <- mgi[which(mgi$bmic>0),]
mgi <- mgi[complete.cases(mgi), ]

write.table(mgi.full,"mgi.txt",sep="\t",row.names=FALSE)


#colonoscopy cpt code82270|82274|G0328|74261|74262|74263|45330|45331|45332|45333|45334|45335|45337|45338|45339|45340|45341|45342|45345|G0104|45.24|44388|44389|44390|44391|44392|44393|44394|44397|45355|45378|45379|45380|45381|45382|45384|45385|45386|45389|45390|45391|45392|45393|45398|G0105|G0121|45.22|45.23|45.25|45.42|45.43
#bp diagnosis: egrep -w 'I10|I11.0|I11.9|I12.0|I12.9|I13.0|I13.10|I13.11|I13.2|I15.0|I15.1|I15.2|I15.8|I15.9|401.0|401.1|401.9|402.00|402.01|402.10|402.11|402.90|402.91|403.00|403.01|403.10|403.11|403.90|403.91|404.00|404.01|404.02|404.03|404.10|404.11|404.12|404.13|404.90|404.91|404.92|404.93|405.01|405.09|405.11|405.19|405.91|405.99' HPI_4605_Diagnosis.txt > hpdg.txt
#bp cpt: egrep -w '93784|93786|93788|93790|99473|99474|90791|90792|90832|90834|90837|90839|90845|90880|92002|92004|92012|92014|96118|99201|99202|99203|99204|99205|99212|99213|99214|99281|99282|99283|99284|99285|99215|99304|99305|99306|99307|99308|99309|99310|99318|99324|99325|99326|99327|99328|99334|99335|99336|99337|99340|99341|99342|99343|99344|99345|99347|99348|99349|99350|D7140|D7210|G0101|G0402|G0438|G0439' HPI_4605_Procedures.txt > prhp.txt