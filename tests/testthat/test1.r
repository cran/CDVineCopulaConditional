if(FALSE)#They work, just to check faster the package: FALSE
{
#Test 1

# Read data (data from \code{VineCopula} package)
data(dataset) 
plot(dataset$data)
data <- dataset$data[1:100,1:3]

# Define the variables Y and X. X are the conditioning variables and have to be positioned in the last columns of the \code{data.frame}
colnames(data) <- c("Y1","Y2","X3")

# Select a vine and fit the copula families, specifying that there are 2 conditioning variables
RVM <- CDVineCondFit(data,Nx=1)

# Set the values of the conditioning variables as those used for the calibration. Order them with respect to RVM$Matrix
cond1 <- data[,RVM$Matrix[1,1]]
condition <- cond1

# Simulate the variables
Sim <- CDVineCondSim(RVM,condition)

# Plot the simulated variables over the observed
Sim <- data.frame(Sim)
overplot(Sim,data)












#Test 2
#1 free variable
#C-Vine
data(dataset) 
data <- dataset$data[1:50,1:5]
colnames(data) <- c("Y1","X2","X3","X4","X5")
RVM <- CDVineCondFit(data,Nx=4,type="CVine")
RVM$Matrix
summary(RVM)
RVM$family
RVM$family[5,1]=1
RVM$par[5,1]=-0.9
RVM$par2[5,1]=0

N=500
dataRsim=data.frame(RVineSim(N,RVM))
#C-Vine
d=dim(RVM$Matrix)[1]
cond1 <- dataRsim[,RVM$Matrix[(d+1)-1,(d+1)-1]]
cond2 <- dataRsim[,RVM$Matrix[(d+1)-2,(d+1)-2]]
cond3 <- dataRsim[,RVM$Matrix[(d+1)-3,(d+1)-3]]
cond4 <- dataRsim[,RVM$Matrix[(d+1)-4,(d+1)-4]]
condition <- cbind(cond1,cond2,cond3,cond4)

dataCondSim <- data.frame(CDVineCondSim(RVM,N=N))
overplot(dataRsim,dataCondSim)

dataCondSim <- data.frame(CDVineCondSim(RVM,condition))
overplot(dataRsim,dataCondSim)


#Test 3
#DVine
data(dataset) 
data <- dataset$data[1:50,1:5]
colnames(data) <- c("Y1","X2","X3","X4","X5")
RVM <- CDVineCondFit(data,Nx=4,type="DVine")
RVM$Matrix
RVM$Matrix
summary(RVM)
RVM$family[4,1]=1
RVM$par[4,1]=-0.9
RVM$par2[4,1]=0

N=500
dataRsim=data.frame(RVineSim(N,RVM))
#D-Vine
cond1 <- dataRsim[,RVM$Matrix[1,1]]
cond2 <- dataRsim[,RVM$Matrix[2,2]]
cond3 <- dataRsim[,RVM$Matrix[3,3]]
cond4 <- dataRsim[,RVM$Matrix[4,4]]
condition <- cbind(cond1,cond2,cond3,cond4)

dataCondSim <- data.frame(CDVineCondSim(RVM,N=N))
overplot(dataRsim,dataCondSim)
plot(dataRsim)
plot(dataCondSim)

dataCondSim <- data.frame(CDVineCondSim(RVM,condition))
overplot(dataRsim,dataCondSim)












#Test 4
#2 free variables
#C-Vine
data(dataset) 
data <- dataset$data[1:50,1:5]
colnames(data) <- c("Y1","Y2","X3","X4","X5")
RVM <- CDVineCondFit(data,Nx=3,type="CVine")
RVM$Matrix
summary(RVM)
RVM$family
RVM$family[4,1]=1
RVM$par[4,1]=-0.9
RVM$par2[4,1]=0

N=500
dataRsim=data.frame(RVineSim(N,RVM))
#C-Vine
d=dim(RVM$Matrix)[1]
cond1 <- dataRsim[,RVM$Matrix[(d+1)-1,(d+1)-1]]
cond2 <- dataRsim[,RVM$Matrix[(d+1)-2,(d+1)-2]]
cond3 <- dataRsim[,RVM$Matrix[(d+1)-3,(d+1)-3]]
condition <- cbind(cond1,cond2,cond3)

dataCondSim <- data.frame(CDVineCondSim(RVM,N=N))
overplot(dataRsim,dataCondSim)

dataCondSim <- data.frame(CDVineCondSim(RVM,condition))
overplot(dataRsim,dataCondSim)


#Test 5
#DVine
data(dataset) 
data <- dataset$data[1:50,1:5]
colnames(data) <- c("Y1","Y2","X3","X4","X5")
RVM <- CDVineCondFit(data,Nx=3,type="DVine")
RVM$Matrix
RVM$Matrix
summary(RVM)
RVM$family
RVM$family[4,1]=1
RVM$par[4,1]=-0.9
RVM$par2[4,1]=0

N=1500
dataRsim=data.frame(RVineSim(N,RVM))
#D-Vine
cond1 <- dataRsim[,RVM$Matrix[1,1]]
cond2 <- dataRsim[,RVM$Matrix[2,2]]
cond3 <- dataRsim[,RVM$Matrix[3,3]]
condition <- cbind(cond1,cond2,cond3)

dataCondSim <- data.frame(CDVineCondSim(RVM,N=N))
overplot(dataRsim,dataCondSim)
plot(dataRsim)
plot(dataCondSim)

dataCondSim <- data.frame(CDVineCondSim(RVM,condition))
overplot(dataRsim,dataCondSim)



















#Test 6
data(dataset)
data <- dataset$data[1:20,c(1,2,3,4,5)]
colnames(data) <- c("Y1","Y2","Y3","Y4","Y5")
RVM <- CDVineCondFit(data,Nx=2,type="CVine")
RVM$Matrix
RVM$family
RVM$family=RVM$family/RVM$family
RVM$family[which(!is.finite(RVM$family))]=0
RVM$par=RVM$family/2
RVM$par2=RVM$par2*0

RVM$family[3,2]=17
RVM$par[3,2]=1.1
RVM$par2[3,2]=2.2
RVM$family[5,4]=33
RVM$par[5,4]=-7
RVM$par2[5,4]=0
RVM$family[5,3]=23
RVM$par[5,3]=-7
RVM$par2[5,3]=0
RVM$family[2,1]=214
RVM$par[2,1]=3
RVM$par2[2,1]=0.6
RVM$family[5,1]=104
RVM$par[5,1]=20
RVM$par2[5,1]=0.09

RVM$Matrix
RVM$family
RVM$par
RVM$par2

N=5000
dataRsim=RVineSim(N,RVM)
dataRsim=data.frame(dataRsim)
colnames(dataRsim) <- c("Y1","Y2","Y3","Y4","Y5")
RVM$Matrix
dataCondSim <- CDVineCondSim(RVM,N=N)
dataCondSim=data.frame(dataCondSim)
colnames(dataCondSim) <- c("Y1","Y2","Y3","Y4","Y5")
overplot(dataRsim,dataCondSim)
plot(dataRsim)
plot(dataCondSim)
hist(dataRsim[,1],breaks = 100)
hist(dataRsim[,2],breaks = 100)
hist(dataRsim[,3],breaks = 100)
hist(dataRsim[,4],breaks = 100)
hist(dataRsim[,5],breaks = 100)
hist(dataCondSim[,1],breaks = 100)
hist(dataCondSim[,2],breaks = 100)
hist(dataCondSim[,3],breaks = 100)
hist(dataCondSim[,4],breaks = 100)
hist(dataCondSim[,5],breaks = 100)

RVM2 <- CDVineCondFit(dataCondSim,Nx=2,type="CVine")
dataCondSim2 <- CDVineCondSim(RVM2,N=N)
dataCondSim2=data.frame(dataCondSim2)
colnames(dataCondSim2) <- c("Y1","Y2","Y3","Y4","Y5")
overplot(dataCondSim,dataCondSim2)















#Test 7
data(dataset) 
data <- dataset$data[1:20,c(1,2,3,4,5)]
colnames(data) <- c("Y1","Y2","Y3","Y4","Y5")
RVM <- CDVineCondFit(data,Nx=2,type="CVine")
RVM$Matrix
RVM$family
RVM$family=RVM$family/RVM$family
RVM$family[which(!is.finite(RVM$family))]=0
RVM$par=RVM$family/2
RVM$par2=RVM$par2*0

#RVM$family[3,2]=17
#RVM$par[3,2]=1.1
#RVM$par2[3,2]=2.2
#RVM$family[5,4]=33
#RVM$par[5,4]=-7
#RVM$par2[5,4]=0
#RVM$family[5,3]=23
#RVM$par[5,3]=-7
#RVM$par2[5,3]=0
RVM$family[3,1]=23
RVM$par[3,1]=-7
RVM$par2[3,1]=0
#RVM$family[2,1]=23
#RVM$par[2,1]=-7
#RVM$par2[2,1]=0
#RVM$family[2,1]=214
#RVM$par[2,1]=3
#RVM$par2[2,1]=0.6
#RVM$family[5,1]=104
#RVM$par[5,1]=20
#RVM$par2[5,1]=0.09

RVM$Matrix
RVM$family
RVM$par
RVM$par2

N=2500
dataRsim=RVineSim(N,RVM)
dataRsim=data.frame(dataRsim)
colnames(dataRsim) <- c("Y1","Y2","Y3","Y4","Y5")
RVM$Matrix
dataCondSim <- CDVineCondSim(RVM,N=N)
dataCondSim=data.frame(dataCondSim)
colnames(dataCondSim) <- c("Y1","Y2","Y3","Y4","Y5")
overplot(dataRsim,dataCondSim)
plot(dataRsim)
plot(dataCondSim)
hist(dataRsim[,1],breaks = 100)
hist(dataRsim[,2],breaks = 100)
hist(dataRsim[,3],breaks = 100)
hist(dataRsim[,4],breaks = 100)
hist(dataRsim[,5],breaks = 100)
hist(dataCondSim[,1],breaks = 100)
hist(dataCondSim[,2],breaks = 100)
hist(dataCondSim[,3],breaks = 100)
hist(dataCondSim[,4],breaks = 100)
hist(dataCondSim[,5],breaks = 100)













#Test 8
data(dataset) 
data <- dataset$data[1:20,c(1,2,3,4,5)]
colnames(data) <- c("Y1","Y2","Y3","Y4","Y5")
RVM <- CDVineCondFit(data,Nx=2,type="DVine")
RVM$Matrix
RVM$family
RVM$family=RVM$family/RVM$family
RVM$family[which(!is.finite(RVM$family))]=0
RVM$par=RVM$family/1.5
RVM$par2=RVM$par2*0

#RVM$family[3,2]=17
#RVM$par[3,2]=1.1
#RVM$par2[3,2]=2.2
#RVM$family[5,4]=23
#RVM$par[5,4]=-7
#RVM$par2[5,4]=0
#RVM$family[5,4]=33
#RVM$par[5,4]=-7
#RVM$par2[5,4]=0
#RVM$family[5,3]=23
#RVM$par[5,3]=-87
#RVM$par2[5,3]=0
RVM$family[5,1]=23
RVM$par[5,1]=-87
RVM$par2[5,1]=0
#RVM$family[2,1]=214
#RVM$par[2,1]=3
#RVM$par2[2,1]=0.6
#RVM$family[5,1]=104
#RVM$par[5,1]=20
#RVM$par2[5,1]=0.09

RVM$Matrix
RVM$family
RVM$par
RVM$par2

N=1000
dataRsim=RVineSim(N,RVM)
dataRsim=data.frame(dataRsim)
colnames(dataRsim) <- c("Y1","Y2","Y3","Y4","Y5")
RVM$Matrix
dataCondSim <- CDVineCondSim(RVM,N=N)
dataCondSim=data.frame(dataCondSim)
colnames(dataCondSim) <- c("Y1","Y2","Y3","Y4","Y5")
overplot(dataRsim,dataCondSim)
overplot(dataCondSim,dataRsim)
plot(dataRsim)
plot(dataCondSim)
hist(dataRsim[,1],breaks = 100)
hist(dataRsim[,2],breaks = 100)
hist(dataRsim[,3],breaks = 100)
hist(dataRsim[,4],breaks = 100)
hist(dataRsim[,5],breaks = 100)
hist(dataCondSim[,1],breaks = 100)
hist(dataCondSim[,2],breaks = 100)
hist(dataCondSim[,3],breaks = 100)
hist(dataCondSim[,4],breaks = 100)
hist(dataCondSim[,5],breaks = 100)
RVM$Matrix
RVM$family

#D-Vine
cond1 <- dataRsim[,RVM$Matrix[1,1]]
cond2 <- dataRsim[,RVM$Matrix[2,2]]
condition <- cbind(cond1,cond2)

dataCondSim <- data.frame(CDVineCondSim(RVM,condition))
overplot(dataRsim,dataCondSim)















#Test 9
data(dataset) 
data <- dataset$data[1:20,c(3,2,1)]
RVM <- CDVineCondFit(data,Nx=2,type="DVine")
RVM$Matrix
RVM$family
RVM$family=RVM$family/RVM$family
RVM$family[which(!is.finite(RVM$family))]=0
RVM$par=RVM$family/1.5
RVM$par2=RVM$par2*0

#RVM$family[3,2]=17
#RVM$par[3,2]=1.1
#RVM$par2[3,2]=2.2
#RVM$family[5,4]=23
#RVM$par[5,4]=-7
#RVM$par2[5,4]=0
#RVM$family[5,4]=33
#RVM$par[5,4]=-7
#RVM$par2[5,4]=0
#RVM$family[5,3]=23
#RVM$par[5,3]=-87
#RVM$par2[5,3]=0
RVM$family[3,1]=23
RVM$par[3,1]=-80
RVM$par2[3,1]=0
RVM$family[3,1]=214
RVM$par[3,1]=80
RVM$par2[3,1]=0.7
#RVM$family[2,1]=214
#RVM$par[2,1]=3
#RVM$par2[2,1]=0.6
#RVM$family[5,1]=104
#RVM$par[5,1]=20
#RVM$par2[5,1]=0.09

RVM$Matrix[1,1]=1
RVM$Matrix[2,2]=2
RVM$Matrix[3,3]=3
RVM$Matrix[2,1]=3
RVM$Matrix[3,2]=3
RVM$Matrix
RVM$family
RVM$par
RVM$par2

N=1000
dataRsim=RVineSim(N,RVM)
dataRsim=data.frame(dataRsim)
colnames(dataRsim) <- c("Y1","Y2","Y3")
RVM$Matrix
dataCondSim <- CDVineCondSim(RVM,N=N)
dataCondSim=data.frame(dataCondSim)
colnames(dataCondSim) <- c("Y1","Y2","Y3")
overplot(dataRsim,dataCondSim)
RVM$Matrix
RVM$family
RVM

overplot(dataCondSim,dataRsim)
plot(dataRsim)
plot(dataCondSim)
hist(dataRsim[,1],breaks = 100)
hist(dataRsim[,2],breaks = 100)
hist(dataRsim[,3],breaks = 100)
hist(dataCondSim[,1],breaks = 100)
hist(dataCondSim[,2],breaks = 100)
hist(dataCondSim[,3],breaks = 100)
RVM$Matrix
RVM$family
















    
    
    
    
    
    


    
    

    
    
    
    
        
#Test 10
data(dataset) 
data <- dataset$data[1:20,c(1,2,3,4,5)]
colnames(data) <- c("Y1","Y2","Y3","Y4","Y5")
RVM <- CDVineCondFit(data,Nx=2,type="DVine")
RVM$Matrix
RVM$family
RVM$family=RVM$family/RVM$family
RVM$family[which(!is.finite(RVM$family))]=0
RVM$par=RVM$family/1.5
RVM$par2=RVM$par2*0

#RVM$family[3,2]=17
#RVM$par[3,2]=1.1
#RVM$par2[3,2]=2.2
RVM$family[5,4]=23
RVM$par[5,4]=-7
RVM$par2[5,4]=0
#RVM$family[5,4]=33
#RVM$par[5,4]=-7
#RVM$par2[5,4]=0
#RVM$family[5,3]=23
#RVM$par[5,3]=-7
#RVM$par2[5,3]=0
#RVM$family[2,1]=214
#RVM$par[2,1]=3
#RVM$par2[2,1]=0.6
#RVM$family[5,1]=104
#RVM$par[5,1]=20
#RVM$par2[5,1]=0.09

RVM$Matrix
RVM$family
RVM$par
RVM$par2

N=1000
dataRsim=RVineSim(N,RVM)
dataRsim=data.frame(dataRsim)
colnames(dataRsim) <- c("Y1","Y2","Y3","Y4","Y5")
RVM$Matrix
dataCondSim <- CDVineCondSim(RVM,N=N)
dataCondSim=data.frame(dataCondSim)
colnames(dataCondSim) <- c("Y1","Y2","Y3","Y4","Y5")
overplot(dataRsim,dataCondSim)
overplot(dataCondSim,dataRsim)
plot(dataRsim)
plot(dataCondSim)
hist(dataRsim[,1],breaks = 100)
hist(dataRsim[,2],breaks = 100)
hist(dataRsim[,3],breaks = 100)
hist(dataRsim[,4],breaks = 100)
hist(dataRsim[,5],breaks = 100)
hist(dataCondSim[,1],breaks = 100)
hist(dataCondSim[,2],breaks = 100)
hist(dataCondSim[,3],breaks = 100)
hist(dataCondSim[,4],breaks = 100)
hist(dataCondSim[,5],breaks = 100)
RVM$Matrix
RVM$family

#D-Vine
cond1 <- dataRsim[,RVM$Matrix[1,1]]
cond2 <- dataRsim[,RVM$Matrix[2,2]]
condition <- cbind(cond1,cond2)

dataCondSim <- data.frame(CDVineCondSim(RVM,condition))
overplot(dataRsim,dataCondSim)







#Test 11
data(dataset) 
data <- dataset$data[1:20,c(1,2,3,4,5)]
colnames(data) <- c("Y1","Y2","Y3","Y4","Y5")
RVM <- CDVineCondFit(data,Nx=2,type="DVine")
RVM$Matrix
RVM$family
RVM$family=RVM$family/RVM$family
RVM$family[which(!is.finite(RVM$family))]=0
RVM$par=RVM$family/2
RVM$par2=RVM$par2*0

RVM$Matrix
RVM$family
RVM$par
RVM$par2

N=2500
dataRsim=RVineSim(N,RVM)
dataRsim=data.frame(dataRsim)
colnames(dataRsim) <- c("Y1","Y2","Y3","Y4","Y5")
RVM$Matrix
dataCondSim <- CDVineCondSim(RVM,N=N)
dataCondSim=data.frame(dataCondSim)
colnames(dataCondSim) <- c("Y1","Y2","Y3","Y4","Y5")
overplot(dataRsim,dataCondSim)
hist(dataRsim[,1],breaks = 100)
hist(dataRsim[,2],breaks = 100)
hist(dataRsim[,3],breaks = 100)
hist(dataRsim[,4],breaks = 100)
hist(dataRsim[,5],breaks = 100)
hist(dataCondSim[,1],breaks = 100)
hist(dataCondSim[,2],breaks = 100)
hist(dataCondSim[,3],breaks = 100)
hist(dataCondSim[,4],breaks = 100)
hist(dataCondSim[,5],breaks = 100)














#Test 12
data(dataset) 
data <- dataset$data[1:100,c(1,2,3,4,5)]
colnames(data) <- c("Y1","Y2","Y3","Y4","Y5")
RVM <- CDVineCondFit(data,Nx=2,type="DVine")
RVM$Matrix
RVM$family

RVM$Matrix
RVM$family
RVM$par
RVM$par2

N=8500
dataRsim=RVineSim(N,RVM)
dataRsim=data.frame(dataRsim)
colnames(dataRsim) <- c("Y1","Y2","Y3","Y4","Y5")
RVM$Matrix
dataCondSim <- CDVineCondSim(RVM,N=N)
dataCondSim=data.frame(dataCondSim)
colnames(dataCondSim) <- c("Y1","Y2","Y3","Y4","Y5")
overplot(dataRsim,dataCondSim)
overplot(dataCondSim,dataRsim)
par(mfrow=c(2,5))
hist(dataRsim[,1],breaks = 100)
hist(dataRsim[,2],breaks = 100)
hist(dataRsim[,3],breaks = 100)
hist(dataRsim[,4],breaks = 100)
hist(dataRsim[,5],breaks = 100)
hist(dataCondSim[,1],breaks = 100)
hist(dataCondSim[,2],breaks = 100)
hist(dataCondSim[,3],breaks = 100)
hist(dataCondSim[,4],breaks = 100)
hist(dataCondSim[,5],breaks = 100)
par(mfrow=c(1,1))




















#Test 13
data(dataset) 
data <- dataset$data[1:100,c(1,2)]
colnames(data) <- c("Y1","Y2")
RVM <- CDVineCondFit(data,Nx=2,type="DVine")
RVM$Matrix
RVM$family

RVM$Matrix
RVM$family[2,1]=23
RVM$par[2,1]=-5
RVM$par2[2,1]=0

N=850
dataRsim=RVineSim(N,RVM)
dataRsim=data.frame(dataRsim)
colnames(dataRsim) <- c("Y1","Y2")
RVM$Matrix
dataCondSim <- CDVineCondSim(RVM,N=N)
dataCondSim=data.frame(dataCondSim)
colnames(dataCondSim) <- c("Y1","Y2")
overplot(dataRsim,dataCondSim)

#D-Vine
cond1 <- dataRsim[,RVM$Matrix[1,1]]
condition <- cbind(cond1)
dataCondSim <- data.frame(CDVineCondSim(RVM,condition))
overplot(dataRsim,dataCondSim)

}



