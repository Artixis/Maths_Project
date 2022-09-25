library(tidyr)
library(MASS)
library(sf)
library(rgeos)
library(gmt)
library(lme4)
library(lmerTest)
library(ggplot2)

dt <- read.csv("32180DS0001_2001-21_Analysis.csv",skip=6)

dt <- dt[-1:-2,]
dt$SA2 <- dt$X.9

dt <- dt[,c("SA2","Area",paste0("X",2001:2021))]




dtlong <- dt %>% gather(year,pop,paste0("X",2001:2021))
dtlong$year <- as.numeric(gsub("X","",dtlong$year))
dtlong$pop <- as.numeric(dtlong$pop)
dtlong$Area <- as.numeric(dtlong$Area)
dtlong$density <- dtlong$pop/dtlong$Area

dtlong$yearsince2000 <- dtlong$year-2000 ##


syd <- st_read("/Users/russelthomson/Documents/WSU/nbn-trends/data/1270055001_sa2_2016_aust_shape/syd.shp")
syd <- syd[syd$GCC_NAME16=="Greater Sydney",]
centr <- st_centroid(syd)#, byid = TRUE)
centroid_lat <- unlist(lapply(centr$geometry,function(l) l[2]))
centroid_long <- unlist(lapply(centr$geometry,function(l) l[1]))

# Sydney 33.8688° S, 151.2093° E
# Parramatta 33.8148° S, 151.0017° E

syd$sydist <- geodist(-33.8688,151.2093,centroid_lat,centroid_long)
syd$parradist <- geodist(-33.8148,151.0017,centroid_lat,centroid_long)

dtlong$sydist <- syd$sydist[match(dtlong$SA2,syd$SA2_NAME16)]
dtlong$parradist <- syd$parradist[match(dtlong$SA2,syd$SA2_NAME16)]

dtlong$logsydist <- log(syd$sydist[match(dtlong$SA2,syd$SA2_NAME16)])
dtlong$logparradist <- log(syd$parradist[match(dtlong$SA2,syd$SA2_NAME16)])

# Drop SA2's without syddist
dtlong <- dtlong[!is.na(dtlong$sydist),]

# Drop SA2's with any pop of zero.
SA2s_with_zero <- unique(dtlong$SA2[dtlong$density==0])
dtlong <- dtlong[!dtlong$SA2 %in% SA2s_with_zero,]

# What are SA2's with smallest density?
head(dtlong[order(dtlong$density),])


bc <- boxcox(lm(density~yearsince2000*parradist+yearsince2000*sydist,data=dtlong[dtlong$density>0,]))
bc$x[which.max(bc$y)]
#Suggests sqrt transform

dtlong$trans_density <- sqrt(dtlong$density)

# Fixed Effect Model
fit_fixed <- lm(trans_density~yearsince2000*parradist+yearsince2000*sydist,data=dtlong)
summary(fit_fixed)

# Mixed Effect Model, Random Intercepts
fit_randint <- lmer(trans_density~yearsince2000*logparradist*logsydist+(1|SA2),data=dtlong)
summary(fit_randint)
# Random Intercept for each SA2  much bigger than residual variation
# this model better than fixed effect.

# Mixed Effect Model, Random Intercepts and Slopes
fit_randslope1 <- lmer(trans_density~yearsince2000*parradist*sydist+(yearsince2000|SA2),data=dtlong)
summary(fit_randslope1)
# intercept random effect is small. 

anova(fit_randint,fit_randslope1)
# Random slope model significantly better fit.

logsydist_1st_3rd_quartile <- quantile(dtlong$logsydist,c(0.25,.75))
logparradist_1st_3rd_quartile <- quantile(dtlong$logparradist,c(0.25,.75))


preds <- data.frame(yearsince2000=rep(1:21,4)
                    ,logsydist=rep(logsydist_1st_3rd_quartile,each=21)
                    ,logparradist=rep(logparradist_1st_3rd_quartile,each=42)
                    ,pred_trans_density=NA)
preds$sydist <- factor( round(exp(preds$logsydist)) )
preds$parradist <- factor( round(exp(preds$logparradist)) )
preds$year <- preds$yearsince2000+2000

preds$pred_trans_density <- predict(fit_randslope1,newdata=preds,re.form=NA)

# Back transform density
preds$pred_density <- preds$pred_trans_density^2 

# Plot Predictions
i <- ggplot(preds,aes(year,pred_density)) + expand_limits(y = 0)
i + geom_line(aes(colour=sydist,linetype=parradist))




# Test Covid Effect... (with sep. coefficients for 2020 and 2021)


# Plot COVID Effect
