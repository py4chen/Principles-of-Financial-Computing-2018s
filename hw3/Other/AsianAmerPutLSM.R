set.seed(0)

firstValueRow <-
function(x) {
	cumSumMat<-matrix(NA, nrow=dim(x)[1], ncol=dim(x)[2])
	for(i in 1:(dim(x)[1])) {
		cumSumMat[i,]<-cumsum(x[i,])
	}
	cumSumMat2<-cbind(matrix(0, nrow=dim(x)[1], ncol=1), cumSumMat[,-(dim(x)[2])])
	ResultMat<-matrix(NA, nrow=dim(x)[1], ncol=dim(x)[2])
	for(i in 1:dim(x)[2]) {
		ResultMat[,i]<-ifelse(cumSumMat2[,i]>0, 0, x[,i])
	}
	return(ResultMat)	
}

AsianAmerPutLSM <-
function (Spot=100, sigma=0.3, n=100000, m=100, Strike=100, r=0.05, dr=0.0, mT=1) {
	GBM<-matrix(NA, nrow=n, ncol=m)
	for(i in 1:n) {
		GBM[i,]<-Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
		# rnormal<-matrix(c(-1.3530204,-0.11307944, -1.21280941,  0.99634801,  0.2988387,0.69489436), nrow = 1)
		# GBM[i,]<-Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnormal)))
		}

	CFLp<-cbind(GBM[,1], matrix(NA, nrow=n, ncol=m-1))	


	for(i in 2:m) {
		CFLp[,i]<-rowMeans(GBM[,1:i])
		}

	CFL<-matrix(pmax(0,Strike-CFLp), nrow=n, ncol=m)

	X<-ifelse(CFL>0,GBM,NA)

	Xsh<-X[,-m]
	

	X2sh<-Xsh*Xsh
	Y1<-CFL*exp(-1*r*(mT/m))

# print(dim(Y1))

	Y2<-cbind((matrix(NA, nrow=n, ncol=m-1)), Y1[,m])
	

	# print(dim(Y2))

	CV<-matrix(NA, nrow=n, ncol=m-1)
	try(for(i in (m-1):1) {
		# test<-Xsh[,i]+X2sh[,i]
		# print(Xsh[,i])
		# quit()
		# print(Xsh[,i])
		# print(Y2[,i+1])


		reg1<-lm(Y2[,i+1]~Xsh[,i]+X2sh[,i])
		
		# print(reg1)


		CV[,i]<-(matrix(reg1$coefficients)[1,1])+((matrix(reg1$coefficients)[2,1])*Xsh[,i])+((matrix(reg1$coefficients)[3,1])*X2sh[,i])
		# print(CV[,i])
		# quit()

		CV[,i]<-(ifelse(is.na(CV[,i]),0,CV[,i]))
		Y2[,i]<-ifelse(CFL[,i]>CV[,i], Y1[,i], Y2[,i+1]*exp(-1*r*(mT/m)))
		}
	, silent = TRUE)
	CV<-ifelse(is.na(CV),0,CV)
	CVp<-cbind(CV, (matrix(0, nrow=n, ncol=1)))
	POF<-ifelse(CVp>CFL,0,CFL)
	FPOF<-firstValueRow(POF)
	dFPOF<-matrix(NA, nrow=n, ncol=m)
	for(i in 1:m) {
		dFPOF[,i]<-FPOF[,i]*exp(-1*mT/m*r*i)
		}
	PRICE<-mean(rowSums(dFPOF))
	std<-function(x) sqrt(var(x)/length(x))
    STD<-std(rowSums(dFPOF))
    print(rowSums(dFPOF))
	res<- list(price=(PRICE), std=STD, Spot, Strike, sigma, n, m, r, dr, mT)
	class(res)<-"AsianAmerPut"
	return(res)
}

ret <- AsianAmerPutLSM()
print(ret)
