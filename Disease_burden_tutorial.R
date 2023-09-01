rm(list=ls()) #run this to clear your environment 

#upload packages 
library(splines)
library(tidyverse)
library(zoo)

#set working directory - this is where we will save everything 
#setwd("C:/Users/hansencl/Desktop/ANISE 2023 activity") 


#read in the dataset we will use 
dat = read.csv("respiratory_mortality_SouthAfrica_2009_2018.csv") 

## check the dates variable to make sure R understand it as a date
class(dat$date) #reads as a character 
dat$date = as.Date(dat$date, format = '%m/%d/%Y') #convert to date 

#check again 
class(dat$date) #now it is a date 

#have a look at the data 
str(dat)
head(dat)
tail(dat)
View(dat)

dat = dat %>% 
  arrange(date) %>% 
  mutate(resp_rate = rollmean(resp_uc/population*100000,k=5,align="center",fill="extend"))#convert to a mortality rate and smooth



#make a figure to show respiratory mortality and influenza circulation - plot from slides 
plot1 = ggplot(data=dat)+
  theme_bw()+
  geom_line(aes(x=date, y=resp_rate,color="Underlying respiratory mortality rate"),linewidth=2)+
  geom_line(aes(x=date, y=perc_pos_flu*100, color="Percent or samples positive for influenza"))+
  labs(x=NULL, y="Mortality Rate and Viral Circulation")+
  scale_color_manual(name=NULL, values=c("aquamarine3","black"))+
  theme(legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15))
plot1 #have a look at the plot 
#save the plot 
ggsave(plot=plot1, "mortality and viral time series.png",height=5,width=13, units="in")



#Make a separate variable for flu and RSV for each season 
dat = dat %>% 
  arrange(date) %>% 
  mutate(season = year(date),
         flu = perc_pos_flu,
         rsv = perc_pos_rsv) %>% 
  pivot_wider(names_from = season, values_from=c("flu", "rsv"), values_fill=0)

#make a flu season matrix 
flu_matrix = as.matrix(dat[1:nrow(dat),c(grep("flu_",names(dat)))])
View(flu_matrix)



# Seasonal harmonic regression with covariates for flu and rsv --------
weeknum=seq(1,dim(dat)[1])
cos1= cos(2*pi/52.17*weeknum)
sin1= sin(2*pi/52.17*weeknum)

#single covariate for flu and single covariate for rsv
model_1 = glm(dat$resp_rate~ weeknum + weeknum^2 + cos1 + sin1  + dat$perc_pos_flu + dat$perc_pos_rsv, family=gaussian(link="identity"))
summary(model_1)



#try the same model with a covariate for each flu season 
model_2 = glm(dat$resp_rate~ weeknum + weeknum^2 + cos1 + sin1  + flu_matrix + dat$perc_pos_rsv, family=gaussian(link="identity"))
summary(model_2)



#try a coviariate for each flu subtype 
model_3 = glm(dat$resp_rate~ weeknum + weeknum^2 + cos1 + sin1  + dat$perc_pos_fluAh1n1 + dat$perc_pos_fluAh3n2 + dat$perc_pos_fluB + dat$perc_pos_rsv, family=gaussian(link="identity"))
summary(model_3)


#compare AIC
AIC(model_1)#single flu term
AIC(model_2)#term for each flu season
AIC(model_3)#term for each flu subtype



#Which is the best model? 
#Answer: The model with a separate covariate for each flu season 

#let's proceed using this model 
#create excess estimates table 
X = model.matrix(model_2)
coeff = model_2$coeff

## Setting all flu coeffs to zero for flu baseline
baseline.flu = model_2$coeff
baseline.flu[grep("flu",names(baseline.flu))] = 0

## Setting all RSV coeffs to zero for RSV baseline
baseline.rsv = model_2$coeff
baseline.rsv[grep("rsv",names(baseline.rsv))] = 0

### Estimating baseline, which is everything minus flu and RSV
seasonal.term = model_2$coeff
viral=c("flu","rsv")
seasonal.term[grepl(paste(viral,collapse="|"),names(seasonal.term))] = 0

data1 = as.data.frame(X)
data1$season = year(dat$date)
data1$date = dat$date
data1$population = dat$population
data1$flu = c((X %*% coeff)-(X %*% baseline.flu))
data1$rsv = c((X %*% coeff)-(X %*% baseline.rsv))
data1$baseline=c(X %*% seasonal.term)
data1$observed = dat$resp_rate 
data1$predicted = c(X%*%coeff)
pop = tapply(data1$population, data1$season, mean)
excess.flu.rate.estimated = tapply(data1$flu, data1$season, sum)
excess.rsv.rate.estimated = tapply(data1$rsv, data1$season, sum)
observed = tapply(data1$observed, data1$season, sum)
predicted = tapply(data1$predicted, data1$season, sum)
baseline=tapply(data1$baseline, data1$season, sum)

### Excess table
### Note here excess estimates are rates per 100,000 as this is what was modeled in glm	

out = as.data.frame(cbind(pop, 
                          excess.flu.rate.estimated, 
                          excess.rsv.rate.estimated, 
                          baseline,
                          observed,
                          predicted))
out$season = rownames(out)
View(out)

#why do you think RSV is negative? 
#Can you make the figure in plot 1 but with RSV instead of flu? 
#what do you notice about the RSV peaks and the mortality peaks? 

data1$baseline_withrsv = data1$baseline + data1$rsv
data1$flu.ex = data1$baseline_withrsv + data1$flu

plot2 <- ggplot(data=data1)+
  geom_area(aes(x=date, y=flu.ex, fill="Excess attributed to flu"))+
  geom_area(aes(x=date, y=baseline_withrsv, fill="Seasonal baseline"))+
  geom_line(aes(x=date, y=observed, color="Observed" ),linetype="dashed",size=1.5)+
  geom_line(aes(x=date, y=predicted, color="Modeled"))+
  geom_line(aes(x=date, y=baseline_withrsv),color="black")+
  scale_color_manual(name=NULL, values=c("hotpink","navy"))+
  scale_fill_manual(name=NULL, values=c("papayawhip","grey90"))+
  labs(x=NULL, y="Respiratory deaths per 100,000")+
  theme_bw()
plot2
ggsave(plot=plot2, "excess mortality rate - harmonic regression.png",height=4,width=13, units="in")



# Replace the linear and harmonic terms with a spline ---------------------
##### Splines provided through the "ns" function, which takes care of both seasonal trends and time trends
##### Spline degree has to be defined, suggest to start with 2 or 3 by default, and compare AIC
degree =3 
model.spline = glm(dat$resp_rate ~ flu_matrix + dat$perc_pos_rsv + ns(c(1:nrow(dat)), round(degree*(nrow(dat)/52.17))),
                     family=gaussian(link="identity"))
summary(model.spline)
AIC(model.spline)

#check AIC with other degrees 
#making the figures and estimates is the same as the harmonic regression example - give it a try  
#replace lines 96 - 153 with the name of your new model 



# Poisson example ---------------------------------------------------------
# for Poisson regression we will use the number of deaths instead of the mortality rate
# we will need to include an offset for the population size 
# we will need to exponentiate the results 
#put a rolling average on the count data 
dat$resp_smooth = round(rollmean(dat$resp_uc,k=5,align="center",fill="extend"))

model_poisson = glm(dat$resp_smooth~ weeknum + weeknum^2 + cos1 + sin1  + flu_matrix + dat$perc_pos_rsv + offset(log(dat$population)), family="poisson")
summary(model_poisson)
AIC(model_poisson)

#create excess estimates table 
X = model.matrix(model_poisson)
coeff = model_poisson$coeff

## Setting all flu coeffs to zero for flu baseline
baseline.flu = model_poisson$coeff
baseline.flu[grep("flu",names(baseline.flu))] = 0

## Setting all RSV coeffs to zero for RSV baseline
baseline.rsv = model_poisson$coeff
baseline.rsv[grep("rsv",names(baseline.rsv))] = 0

### Estimating baseline, which is everything minus flu and RSV
seasonal.term = model_poisson$coeff
viral=c("flu","rsv")
seasonal.term[grepl(paste(viral,collapse="|"),names(seasonal.term))] = 0

data1 = as.data.frame(X)
data1$season = year(dat$date)
data1$date = dat$date
data1$population = dat$population
data1$flu = exp(X %*% coeff + model_poisson$offset)-exp(X %*% baseline.flu+ model_poisson$offset)
data1$rsv = exp(X %*% coeff+ model_poisson$offset)-exp(X %*% baseline.rsv+ model_poisson$offset)
data1$baseline=exp(X %*% seasonal.term+ model_poisson$offset)
data1$observed = dat$resp_smooth
data1$predicted = exp(X%*%coeff+model_poisson$offset)
pop = tapply(data1$population, data1$season, mean)
excess.flu.rate.estimated = tapply(data1$flu, data1$season, sum)
excess.rsv.rate.estimated = tapply(data1$rsv, data1$season, sum)
observed = tapply(data1$observed, data1$season, sum)
predicted = tapply(data1$predicted, data1$season, sum)
baseline=tapply(data1$baseline, data1$season, sum)

### Excess table
#notice here the excess estimates are in number of deaths, not rates 
out = as.data.frame(cbind(pop, 
                          excess.flu.rate.estimated, 
                          excess.rsv.rate.estimated, 
                          baseline,
                          observed,
                          predicted))
out$season = rownames(out)

data1$baseline_withrsv = data1$baseline + data1$rsv
data1$flu.ex = data1$baseline_withrsv + data1$flu

#do you notice anything different about the figure? 
plot3 <- ggplot(data=data1)+
  geom_area(aes(x=date, y=flu.ex, fill="Excess attributed to flu"))+
  geom_area(aes(x=date, y=baseline_withrsv, fill="Seasonal baseline"))+
  geom_line(aes(x=date, y=observed, color="Observed" ),linetype="dashed",size=1.5)+
  geom_line(aes(x=date, y=predicted, color="Modeled"))+
  geom_line(aes(x=date, y=baseline_withrsv),color="black")+
  scale_color_manual(name=NULL, values=c("hotpink","navy"))+
  scale_fill_manual(name=NULL, values=c("papayawhip","grey90"))+
  labs(x=NULL, y="Respiratory deaths")+
  theme_bw()
plot3
ggsave(plot=plot2, "excess deaths - poisson regression.png",height=4,width=13, units="in")


# Challenge tasks  --------------------------------------------------------

#try repeating this activity with a different outcome (pi_uc, resp_cardio_uc, or natural)

#try a version of the model which has a covariate for each flu season AND a covariate for each RSV season 
  # which model has a lower AIC (single RSV term or multiple RSV terms)
