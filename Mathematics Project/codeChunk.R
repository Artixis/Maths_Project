
## ---- myrcode2
df_pop = read.csv("C:\\Users/Laura/OneDrive - Western Sydney University/Spring 2022/Mathematics Project/Sets Modified/Population_done.csv")
normData  = rnorm(n = 100, mean = 0, sd = 1)
par(mfrow=c(1,2))

hist(normData, main="Noraml Distribution", xlab="Data")
hist(df_pop$distance_from_syd, main="Our Data", xlab="Data")
