.libPaths("~/RPackages")

libs <- c("nnls","SuperLearner","mvtnorm","survival","TH.data","MASS","Matrix",
          "foreach","ranger","lme4","geepack","doParallel","ggplot2",
          "multcomp","doBy","gam","future","data.table","optimx","clubSandwich")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}

load(file=paste0(cloudstor,"PhD/Paper 3 - Methods simulation/Simulation Data/data.RData"))

checkrange <- function(v) {
  v <- ifelse(v<0.001,0.001,v)
  v <- ifelse(v>0.999,ifelse(is.na(v),v,0.999),v)
}

compestimates <- function(data) {
  glmercont <- glmerControl(optimizer="optimx", calc.derivs=FALSE, optCtrl=list(method=c("Nelder-Mead","uobyqa","L-BFGS-B"), starttests=FALSE, kkt=FALSE))
  
  outmodel <- "y~a+la+l+ll+oa+ob+oc+obs"
  outmodel2 <- "y~a+la+ll+oa+ob+oc+obs"
  outmodelm <- "y~a+la+l+ll+oa+ob+oc"
  outmodel2m <- "y~a+la+ll+oa+ob+oc"
  
  propa0model <- "la~ll+oa+ob+oc+obs"
  propa1model <- "a~la+l+ll+oa+ob+oc+obs"
  propa0modelm <- "la~ll+oa+ob+oc"
  propa1modelm <- "a~la+l+ll+oa+ob+oc"
  
  create.Learner("SL.ranger", params = list(num.trees = 250))
  SLlib <- c("SL.glm","SL.glm.interaction","SL.gam","SL.ranger_1")
  
  D <- data
  D$obs <- factor(D$obs)
  I <- data$id
  Y <- as.vector(data$y)
  D1 <- data[,c(-1,-3)]
  
  # Correctly specified GLM G-computation
  GLlp <- glm(l~la+ll+oa+ob+oc,D,family=gaussian)
  GLlprop <- data.frame(cbind(l0=predict(GLlp,newdata=data.frame(D[,c(-5,-6)],la=0),type="response"),l1=predict(GLlp,newdata=data.frame(D[,c(-5,-6)],la=1),type="response")))
  
  GLout <- glm(formula=outmodel,D,family=gaussian)
  GLgcmpd <- data.frame(id=D$id,obs=D$obs,
                        glfy=predict(GLout,newdata=D[,-3],type="response"),
                        glfy00=(predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                        glfy10=(predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=1,a=0),type="response")*(1-GLlprop$l1)),
                        glfy01=(predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)),
                        glfy11=(predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=1),type="response")*GLlprop$l1) + (predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=1,a=1),type="response")*(1-GLlprop$l1)))
  
  # Incorrectly specified GLM G-computation
  GLoutm <- glm(formula=outmodelm,D,family=gaussian)
  GLgcmpdm <- data.frame(id=D$id,obs=D$obs,
                         glfy=predict(GLoutm,newdata=D[,-3],type="response"),
                         glfy00=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                         glfy10=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=1,a=0),type="response")*(1-GLlprop$l1)),
                         glfy01=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)),
                         glfy11=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=1),type="response")*GLlprop$l1) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=1,a=1),type="response")*(1-GLlprop$l1)))
  
  # Correctly specified Random-intercept G-computation
  RIlp <- glmer(l~la+ll+oa+ob+oc+(1|id),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIlprop <- data.frame(cbind(l0=predict(RIlp,newdata=data.frame(D[,c(-5,-6)],la=0),type="response"),l1=predict(RIlp,newdata=data.frame(D[,c(-5,-6)],la=1),type="response")))
  RIout <- lmer(formula=paste0(outmodel,"+(1|id)"),D)
  RIgcmpd <- data.frame(id=D$id,obs=D$obs,
                        mlfy=predict(RIout,newdata=D,type="response"),
                        mlfy00=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                        mlfy10=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                        mlfy01=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)),
                        mlfy11=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=1),type="response")*RIlprop$l1) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=1,a=1),type="response")*(1-RIlprop$l1)))
  
  # Correctly specified Random-intercept G-computation
  RIoutm <- lmer(formula=paste0(outmodelm,"+(1|id)"),D)
  RIgcmpdm <- data.frame(id=D$id,obs=D$obs,
                         mlfy=predict(RIoutm,newdata=D[,-3],type="response"),
                         mlfy00=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                         mlfy10=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                         mlfy01=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)),
                         mlfy11=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=1),type="response")*RIlprop$l0) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=1,a=1),type="response")*(1-RIlprop$l0)))
  
  # SuperLearner G-computation  
  SLlp <- SuperLearner(Y=as.vector(D$l),X=D[,c(-1,-3,-4,-6)],id=I,SL.library=SLlib,family=binomial)
  SLlprop <- data.frame(cbind(l0=as.vector(predict(SLlp,newdata=data.frame(D[,c(-3,-4,-5,-6)],la=0))$pred),l1=as.vector(predict(SLlp,newdata=data.frame(D[,c(-3,-4,-5,-6)],la=1))$pred)))
  SLout <- SuperLearner(Y=Y,X=D1,id=I,SL.library=SLlib,family=gaussian)
  SLgcmpd <- data.frame(id=D$id,obs=D$obs,
                        slfy=predict(SLout)$pred,
                        slfy00=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=0))$pred)*SLlprop$l0) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l0)),
                        slfy10=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=1,a=0))$pred)*SLlprop$l1) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l1)),
                        slfy01=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=1))$pred)*SLlprop$l0) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=1))$pred)*(1-SLlprop$l0)),
                        slfy11=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=1,a=1))$pred)*SLlprop$l1) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=1,a=1))$pred)*(1-SLlprop$l1)))
  
  # Correctly specified GLM IPTW
  GLexp0 <- glm(formula=propa0model,data=D,family=binomial)
  GLexp1 <- glm(formula=propa1model,data=D,family=binomial)
  GLexpp <- data.table(cbind(id=D$id,obs=D$obs,a=D$a,la=D$la,
                             propa=ifelse(D$a==1,checkrange(predict(GLexp1,type="response")),checkrange(1-predict(GLexp1,type="response"))),
                             propla=ifelse(D$la==1,checkrange(predict(GLexp0,type="response")),checkrange(1-predict(GLexp0,type="response")))))
  GLexpp$p <- GLexpp$propla*GLexpp$propa
  GLexpp$obs <- factor(GLexpp$obs,levels=c(1,2,3,4),labels=c("1","2","3","4"))
  GLexpp$GLwt <- 1/GLexpp$p
  
  GLexpp$Gla0 <- ifelse(GLexpp$la==1,1-GLexpp$propla,GLexpp$propla)
  GLexpp$Gla1 <- ifelse(GLexpp$la==1,GLexpp$propla,1-GLexpp$propla)
  GLexpp$Ga0 <- ifelse(GLexpp$a==1,1-GLexpp$propa,GLexpp$propa)
  GLexpp$Ga1 <- ifelse(GLexpp$a==1,GLexpp$propa,1-GLexpp$propa)
  
  GLiptw <- glm(y~a+la,data=merge(D,GLexpp[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLiptwv <- vcovCR(GLiptw, cluster=D$id, type = "CR3")
  
  # Incorrectly specified GLM IPTW
  GLexp0m <- glm(formula=propa0modelm,data=D,family=binomial)
  GLexp1m <- glm(formula=propa1modelm,data=D,family=binomial)
  GLexppm <- data.table(cbind(id=D$id,obs=D$obs,a=D$a,la=D$la,
                              propa=ifelse(D$a==1,checkrange(predict(GLexp1m,type="response")),checkrange(1-predict(GLexp1m,type="response"))),
                              propla=ifelse(D$la==1,checkrange(predict(GLexp0m,type="response")),checkrange(1-predict(GLexp0m,type="response")))))
  GLexppm$p <- GLexppm$propla*GLexppm$propa
  GLexppm$obs <- factor(GLexppm$obs,levels=c(1,2,3,4),labels=c("1","2","3","4"))
  GLexppm$GLwt <- 1/GLexppm$p
  
  GLexppm$Gla0 <- ifelse(GLexppm$la==1,1-GLexppm$propla,GLexppm$propla)
  GLexppm$Gla1 <- ifelse(GLexppm$la==1,GLexppm$propla,1-GLexppm$propla)
  GLexppm$Ga0 <- ifelse(GLexppm$a==1,1-GLexppm$propa,GLexppm$propa)
  GLexppm$Ga1 <- ifelse(GLexppm$a==1,GLexppm$propa,1-GLexppm$propa)
  
  GLiptwm <- glm(y~a+la,data=merge(D,GLexppm[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLiptwmv <- vcovCR(GLiptwm, cluster=D$id, type = "CR3")
  
  # DR GLM-IPTW
  # Correctly specified
  GLdriptw <- glm(outmodel2,data=merge(D,GLexpp[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdriptwv <- vcovCR(GLdriptw, cluster=D$id, type = "CR3")
  # Propensity incorrectly specified
  GLdriptwm1 <- glm(outmodel2,data=merge(D,GLexppm[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdriptwm1v <- vcovCR(GLdriptwm1, cluster=D$id, type = "CR3")
  # Outcome incorrectly specified
  GLdriptwm2 <- glm(outmodel2m,data=merge(D,GLexpp[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdriptwm2v <- vcovCR(GLdriptwm2, cluster=D$id, type = "CR3")
  # Both incorrectly specified
  GLdriptwm3 <- glm(outmodel2m,data=merge(D,GLexppm[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdriptwm3v <- vcovCR(GLdriptwm3, cluster=D$id, type = "CR3")
  
  # Correctly specified Random-intercept IPTW and DR-IPTW
  RIexp0 <- glmer(formula=paste0(propa0model,"+(1|id)"),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIexp1 <- glmer(formula=paste0(propa1model,"+(1|id)"),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIexpp <- data.table(cbind(id=D$id,obs=D$obs,a=D$a,la=D$la,
                             propa=ifelse(D$a==1,checkrange(predict(RIexp1,type="response")),checkrange(1-predict(RIexp1,type="response"))),
                             propla=ifelse(D$la==1,checkrange(predict(RIexp0,type="response")),checkrange(1-predict(RIexp0,type="response")))))
  RIexpp$p <- RIexpp$propla*RIexpp$propa
  RIexpp$obs <- factor(RIexpp$obs,levels=c(1,2,3,4),labels=c("1","2","3","4"))
  RIexpp$RIwt <- 1/RIexpp$p
  
  RIexpp$Gla0 <- ifelse(RIexpp$la==1,1-RIexpp$propla,RIexpp$propla)
  RIexpp$Gla1 <- ifelse(RIexpp$la==1,RIexpp$propla,1-RIexpp$propla)
  RIexpp$Ga0 <- ifelse(RIexpp$a==1,1-RIexpp$propa,RIexpp$propa)
  RIexpp$Ga1 <- ifelse(RIexpp$a==1,RIexpp$propa,1-RIexpp$propa)
  
  RIiptw <- lmer(y~a+la+(1|id),data=merge(D,RIexpp[,c(1,2,8)]),weight=RIwt)
  
  # Incorrectly specified Random-intercept IPTW and DR-IPTW
  RIexp0m <- glmer(formula=paste0(propa0modelm,"+(1|id)"),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIexp1m <- glmer(formula=paste0(propa1modelm,"+(1|id)"),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIexppm <- data.table(cbind(id=D$id,obs=D$obs,a=D$a,la=D$la,
                              propa=ifelse(D$a==1,checkrange(predict(RIexp1m,type="response")),checkrange(1-predict(RIexp1m,type="response"))),
                              propla=ifelse(D$la==1,checkrange(predict(RIexp0m,type="response")),checkrange(1-predict(RIexp0m,type="response")))))
  RIexppm$p <- RIexppm$propla*RIexppm$propa
  RIexppm$obs <- factor(RIexppm$obs,levels=c(1,2,3,4),labels=c("1","2","3","4"))
  RIexppm$RIwt <- 1/RIexppm$p
  
  RIexppm$Gla0 <- ifelse(RIexppm$la==1,1-RIexppm$propla,RIexppm$propla)
  RIexppm$Gla1 <- ifelse(RIexppm$la==1,RIexppm$propla,1-RIexppm$propla)
  RIexppm$Ga0 <- ifelse(RIexppm$a==1,1-RIexppm$propa,RIexppm$propa)
  RIexppm$Ga1 <- ifelse(RIexppm$a==1,RIexppm$propa,1-RIexppm$propa)
  
  RIiptwm <- lmer(y~a+la+(1|id),data=merge(D,RIexppm[,c(1,2,8)]),weight=RIwt)
  
  # DR Random-intercept-IPTW
  # Correctly specified
  RIdriptw <- lmer(paste0(outmodel2,"+(1|id)"),data=merge(D,RIexpp[,c(1,2,8)]),weight=RIwt)
  # Propensity incorrectly specified
  RIdriptwm1 <- lmer(paste0(outmodel2,"+(1|id)"),data=merge(D,RIexppm[,c(1,2,8)]),weight=RIwt)
  # Outcome incorrectly specified
  RIdriptwm2 <- lmer(paste0(outmodel2m,"+(1|id)"),data=merge(D,RIexpp[,c(1,2,8)]),weight=RIwt) 
  # Both incorrectly specified
  RIdriptwm3 <- lmer(paste0(outmodel2m,"+(1|id)"),data=merge(D,RIexppm[,c(1,2,8)]),weight=RIwt)
  
  # SuperLearner IPTW and DR-IPTW  
  SLexp0 <- SuperLearner(Y=as.vector(data$la),X=data[,c(-1,-3,-4,-5,-6)],id=I,SL.library=SLlib,family=binomial)
  SLexp1 <- SuperLearner(Y=as.vector(data$a),X=data[,c(-1,-3,-4)],id=I,SL.library=SLlib,family=binomial)
  SLexpp <- data.table(cbind(id=D$id,obs=D$obs,a=D$a,la=D$la,
                             propa=ifelse(D$a==1,checkrange(predict(SLexp1)$pred),checkrange(1-predict(SLexp1)$pred)),
                             propla=ifelse(D$la==1,checkrange(predict(SLexp0)$pred),checkrange(1-predict(SLexp0)$pred))))
  SLexpp$p <- SLexpp$propla*SLexpp$propa
  SLexpp$obs <- factor(SLexpp$obs,levels=c(1,2,3,4),labels=c("1","2","3","4"))
  SLexpp$SLwt <- 1/SLexpp$p
  
  SLexpp$Gla0 <- ifelse(SLexpp$la==1,1-SLexpp$propla,SLexpp$propla)
  SLexpp$Gla1 <- ifelse(SLexpp$la==1,SLexpp$propla,1-SLexpp$propla)
  SLexpp$Ga0 <- ifelse(SLexpp$a==1,1-SLexpp$propa,SLexpp$propa)
  SLexpp$Ga1 <- ifelse(SLexpp$a==1,SLexpp$propa,1-SLexpp$propa)
  
  SLiptw <- glm(y~a+la,data=merge(D,SLexpp[,c(1,2,8)]),weight=SLwt,family=gaussian)
  SLiptwv <- vcovCR(SLiptw, cluster=D$id, type = "CR3")
  SLdriptw <- glm(y~a+la+ll+oa+ob+oc,data=merge(D,SLexpp[,c(1,2,8)]),weight=SLwt,family=gaussian)
  SLdriptwv <- vcovCR(SLdriptw, cluster=D$id, type = "CR3")
  
  # TMLE clever covariate creation
  GLexpp$GLlah0 <- -(1-GLexpp$la)/GLexpp$Gla0
  GLexpp$GLlah1 <- GLexpp$la/GLexpp$Gla1
  GLexpp$GLah0 <- -(1-GLexpp$a)/GLexpp$Ga0
  GLexpp$GLah1 <- GLexpp$a/GLexpp$Ga1
  GLexpp$Hla <- GLexpp$GLlah1 + GLexpp$GLlah0
  GLexpp$Ha <- GLexpp$GLah1 + GLexpp$GLah0
  
  GLexppm$GLlah0 <- -(1-GLexppm$la)/GLexppm$Gla0
  GLexppm$GLlah1 <- GLexppm$la/GLexppm$Gla1
  GLexppm$GLah0 <- -(1-GLexppm$a)/GLexppm$Ga0
  GLexppm$GLah1 <- GLexppm$a/GLexppm$Ga1
  GLexppm$Hla <- GLexppm$GLlah1 + GLexppm$GLlah0
  GLexppm$Ha <- GLexppm$GLah1 + GLexppm$GLah0
  
  RIexpp$RIlah0 <- -(1-RIexpp$la)/RIexpp$Gla0
  RIexpp$RIlah1 <- RIexpp$la/RIexpp$Gla1
  RIexpp$RIah0 <- -(1-RIexpp$a)/RIexpp$Ga0
  RIexpp$RIah1 <- RIexpp$a/RIexpp$Ga1
  RIexpp$Hla <- RIexpp$RIlah1 + RIexpp$RIlah0
  RIexpp$Ha <- RIexpp$RIah1 + RIexpp$RIah0
  
  RIexppm$RIlah0 <- -(1-RIexppm$la)/RIexppm$Gla0
  RIexppm$RIlah1 <- RIexppm$la/RIexppm$Gla1
  RIexppm$RIah0 <- -(1-RIexppm$a)/RIexppm$Ga0
  RIexppm$RIah1 <- RIexppm$a/RIexppm$Ga1
  RIexppm$Hla <- RIexppm$RIlah1 + RIexppm$RIlah0
  RIexppm$Ha <- RIexppm$RIah1 + RIexppm$RIah0
  
  SLexpp$SLlah0 <- -(1-SLexpp$la)/SLexpp$Gla0
  SLexpp$SLlah1 <- SLexpp$la/SLexpp$Gla1
  SLexpp$SLah0 <- -(1-SLexpp$a)/SLexpp$Ga0
  SLexpp$SLah1 <- SLexpp$a/SLexpp$Ga1
  SLexpp$Hla <- SLexpp$SLlah1 + SLexpp$SLlah0
  SLexpp$Ha <- SLexpp$SLah1 + SLexpp$SLah0
  
  # GLM TMLE
  # Correctly specified
  DG1 <- merge(D[,c(1:3)],GLexpp[,c(1,2,13:18)])
  DG1 <- merge(DG1,GLgcmpd)
  GLtmle <- glm(data=DG1,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0 <- DG1$glfy00+GLtmle$coefficients["Hla"]*DG1$GLlah0
  GLfla1 <- DG1$glfy10+GLtmle$coefficients["Hla"]*DG1$GLlah1
  GLfa0 <- DG1$glfy00+GLtmle$coefficients["Ha"]*DG1$GLah0
  GLfa1 <- DG1$glfy01+GLtmle$coefficients["Ha"]*DG1$GLah1
  GLf00 <- DG1$glfy00+GLtmle$coefficients["Ha"]*DG1$GLah0+GLtmle$coefficients["Hla"]*DG1$GLlah1
  GLf11 <- DG1$glfy11+GLtmle$coefficients["Ha"]*DG1$GLah1+GLtmle$coefficients["Hla"]*DG1$GLlah1
  GLtmlela <- mean(GLfla1)-mean(GLfla0)
  GLtmlelaIC <- (DG1$GLlah1+DG1$GLlah0)*(DG1$y-DG1$glfy) + DG1$glfy10 - DG1$glfy00 - (GLfla1-GLfla0)
  GLtmlea <- mean(GLfa1)-mean(GLfa0)
  GLtmleaIC <- (DG1$GLah1+DG1$GLah0)*(DG1$y-DG1$glfy) + DG1$glfy01 - DG1$glfy00 - (GLfa1-GLfa0)
  GLtmlejtIC <- (DG1$GLlah1+DG1$GLlah0)*(DG1$y-DG1$glfy) + (DG1$GLah1+DG1$GLah0)*(DG1$y-DG1$glfy) + DG1$glfy11 - DG1$glfy00 - (GLf11-GLf00)
  # Propensity model incorrectly specified
  DG2 <- merge(D[,c(1:3)],GLexppm[,c(1,2,13:18)])
  DG2 <- merge(DG2,GLgcmpd)
  GLtmlem1 <- glm(data=DG2,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0m1 <- DG2$glfy00+GLtmlem1$coefficients["Hla"]*DG2$GLlah0
  GLfla1m1 <- DG2$glfy10+GLtmlem1$coefficients["Hla"]*DG2$GLlah1
  GLfa0m1 <- DG2$glfy00+GLtmlem1$coefficients["Ha"]*DG2$GLah0
  GLfa1m1 <- DG2$glfy01+GLtmlem1$coefficients["Ha"]*DG2$GLah1
  GLf00m1 <- DG2$glfy00+GLtmlem1$coefficients["Ha"]*DG2$GLah0+GLtmlem1$coefficients["Hla"]*DG2$GLlah0
  GLf11m1 <- DG2$glfy11+GLtmlem1$coefficients["Ha"]*DG2$GLah1+GLtmlem1$coefficients["Hla"]*DG2$GLlah1
  GLtmlelam1 <- mean(GLfla1m1)-mean(GLfla0m1)
  GLtmlelam1IC <- (DG2$GLlah1+DG2$GLlah0)*(DG2$y-DG2$glfy) + DG2$glfy10 - DG2$glfy00 - (GLfla1m1-GLfla0m1)
  GLtmleam1 <- mean(GLfa1m1)-mean(GLfa0m1)
  GLtmleam1IC <- (DG2$GLah1+DG2$GLah0)*(DG2$y-DG2$glfy) + DG2$glfy01 - DG2$glfy00 - (GLfa1m1-GLfa0m1)
  GLtmlejtm1IC <- (DG2$GLlah1+DG2$GLlah0)*(DG2$y-DG2$glfy) + (DG2$GLah1+DG2$GLah0)*(DG2$y-DG2$glfy) + DG2$glfy11 - DG2$glfy00 - (GLf11m1-GLf00m1)
  # Outcome model incorrectly specified
  DG3 <- merge(D[,c(1:3)],GLexpp[,c(1,2,13:18)])
  DG3 <- merge(DG3,GLgcmpdm)
  GLtmlem2 <- glm(data=DG3,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0m2 <- DG3$glfy00+GLtmlem2$coefficients["Hla"]*DG3$GLlah0
  GLfla1m2 <- DG3$glfy10+GLtmlem2$coefficients["Hla"]*DG3$GLlah1
  GLfa0m2 <- DG3$glfy00+GLtmlem2$coefficients["Ha"]*DG3$GLah0
  GLfa1m2 <- DG3$glfy01+GLtmlem2$coefficients["Ha"]*DG3$GLah1
  GLf00m2 <- DG3$glfy00+GLtmlem2$coefficients["Ha"]*DG3$GLah0+GLtmlem2$coefficients["Hla"]*DG3$GLlah0
  GLf11m2 <- DG3$glfy11+GLtmlem2$coefficients["Ha"]*DG3$GLah1+GLtmlem2$coefficients["Hla"]*DG3$GLlah1
  GLtmlelam2 <- mean(GLfla1m2)-mean(GLfla0m2)
  GLtmlelam2IC <- (DG3$GLlah1+DG3$GLlah0)*(DG3$y-DG3$glfy) + DG3$glfy10 - DG3$glfy00 - (GLfla1m2-GLfla0m2)
  GLtmleam2 <- mean(GLfa1m2)-mean(GLfa0m2)
  GLtmleam2IC <- (DG3$GLah1+DG3$GLah0)*(DG3$y-DG3$glfy) + DG3$glfy01 - DG3$glfy00 - (GLfa1m2-GLfa0m2)
  GLtmlejtm2IC <- (DG3$GLlah1+DG3$GLlah0)*(DG3$y-DG3$glfy) + (DG3$GLah1+DG3$GLah0)*(DG3$y-DG3$glfy) + DG3$glfy11 - DG3$glfy00 - (GLf11m2-GLf00m2)
  # Both incorrectly specified
  DG4 <- merge(D[,c(1:3)],GLexppm[,c(1,2,13:18)])
  DG4 <- merge(DG4,GLgcmpdm)
  GLtmlem3 <- glm(data=DG4,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0m3 <- DG4$glfy00+GLtmlem3$coefficients["Hla"]*DG4$GLlah0
  GLfla1m3 <- DG4$glfy10+GLtmlem3$coefficients["Hla"]*DG4$GLlah1
  GLfa0m3 <- DG4$glfy00+GLtmlem3$coefficients["Ha"]*DG4$GLah0
  GLfa1m3 <- DG4$glfy01+GLtmlem3$coefficients["Ha"]*DG4$GLah1
  GLf00m3 <- DG4$glfy00+GLtmlem3$coefficients["Ha"]*DG4$GLah0+GLtmlem3$coefficients["Hla"]*DG4$GLlah0
  GLf11m3 <- DG4$glfy11+GLtmlem3$coefficients["Ha"]*DG4$GLah1+GLtmlem3$coefficients["Hla"]*DG4$GLlah1
  GLtmlelam3 <- mean(GLfla1m3)-mean(GLfla0m3)
  GLtmlelam3IC <- (DG4$GLlah1+DG4$GLlah0)*(DG4$y-DG4$glfy) + DG4$glfy10 - DG4$glfy00 - (GLfla1m3-GLfla0m3)
  GLtmleam3 <- mean(GLfa1m3)-mean(GLfa0m3)
  GLtmleam3IC <- (DG4$GLah1+DG4$GLah0)*(DG4$y-DG4$glfy) + DG4$glfy01 - DG4$glfy00 - (GLfa1m3-GLfa0m3)
  GLtmlejtm3IC <- (DG4$GLlah1+DG4$GLlah0)*(DG4$y-DG4$glfy) + (DG4$GLah1+DG4$GLah0)*(DG4$y-DG4$glfy) + DG4$glfy11 - DG4$glfy00 - (GLf11m3-GLf00m3)
  
  # Random Intercept TMLE
  # Correctly specified
  DR1 <- merge(D[,c(1:3)],RIexpp[,c(1,2,13:18)])
  DR1 <- merge(DR1,RIgcmpd)
  RItmle <- glm(data=DR1,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0 <- DR1$mlfy00+RItmle$coefficients["Hla"]*DR1$RIlah0
  RIfla1 <- DR1$mlfy10+RItmle$coefficients["Hla"]*DR1$RIlah1
  RIfa0 <- DR1$mlfy00+RItmle$coefficients["Ha"]*DR1$RIah0
  RIfa1 <- DR1$mlfy01+RItmle$coefficients["Ha"]*DR1$RIah1
  RIf00 <- DR1$mlfy00+RItmle$coefficients["Ha"]*DR1$RIah0+RItmle$coefficients["Hla"]*DR1$RIlah0
  RIf11 <- DR1$mlfy11+RItmle$coefficients["Ha"]*DR1$RIah1+RItmle$coefficients["Hla"]*DR1$RIlah1
  RItmlela <- mean(RIfla1)-mean(RIfla0)
  RItmlelaIC <- (DR1$RIlah1+DR1$RIlah0)*(DR1$y-DR1$mlfy) + DR1$mlfy10 - DR1$mlfy00 - (RIfla1-RIfla0)
  RItmlea <- mean(RIfa1)-mean(RIfa0)
  RItmleaIC <- (DR1$RIah1+DR1$RIah0)*(DR1$y-DR1$mlfy) + DR1$mlfy01 - DR1$mlfy00 - (RIfa1-RIfa0)
  RItmlejtIC <- (DR1$RIlah1+DR1$RIlah0)*(DR1$y-DR1$mlfy) + (DR1$RIah1+DR1$RIah0)*(DR1$y-DR1$mlfy) + DR1$mlfy11 - DR1$mlfy00 - (RIf11-RIf00)
  # Propensity model incorrectly specified
  DR2 <- merge(D[,c(1:3)],RIexppm[,c(1,2,13:18)])
  DR2 <- merge(DR2,RIgcmpd)
  RItmlem1 <- glm(data=DR2,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0m1 <- DR2$mlfy00+RItmlem1$coefficients["Hla"]*DR2$RIlah0
  RIfla1m1 <- DR2$mlfy10+RItmlem1$coefficients["Hla"]*DR2$RIlah1
  RIfa0m1 <- DR2$mlfy00+RItmlem1$coefficients["Ha"]*DR2$RIah0
  RIfa1m1 <- DR2$mlfy01+RItmlem1$coefficients["Ha"]*DR2$RIah1
  RIf00m1 <- DR2$mlfy00+RItmlem1$coefficients["Ha"]*DR2$RIah0+RItmlem1$coefficients["Hla"]*DR2$RIlah0
  RIf11m1 <- DR2$mlfy11+RItmlem1$coefficients["Ha"]*DR2$RIah1+RItmlem1$coefficients["Hla"]*DR2$RIlah1
  RItmlelam1 <- mean(RIfla1m1)-mean(RIfla0m1)
  RItmlelam1IC <- (DR2$RIlah1+DR2$RIlah0)*(DR2$y-DR2$mlfy) + DR2$mlfy10 - DR2$mlfy00 - (RIfla1m1-RIfla0m1)
  RItmleam1 <- mean(RIfa1m1)-mean(RIfa0m1) 
  RItmleam1IC <- (DR2$RIah1+DR2$RIah0)*(DR2$y-DR2$mlfy) + DR2$mlfy01 - DR2$mlfy00 - (RIfa1m1-RIfa0m1)
  RItmlejtm1IC <- (DR2$RIlah1+DR2$RIlah0)*(DR2$y-DR2$mlfy) + (DR2$RIah1+DR2$RIah0)*(DR2$y-DR2$mlfy) + DR2$mlfy11 - DR2$mlfy00 - (RIf11m1-RIf00m1)
  # Outcome model incorrectly specified
  DR3 <- merge(D[,c(1:3)],RIexpp[,c(1,2,13:18)])
  DR3 <- merge(DR3,RIgcmpdm)
  RItmlem2 <- glm(data=DR3,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0m2 <- DR3$mlfy00+RItmlem2$coefficients["Hla"]*DR3$RIlah0
  RIfla1m2 <- DR3$mlfy10+RItmlem2$coefficients["Hla"]*DR3$RIlah1
  RIfa0m2 <- DR3$mlfy00+RItmlem2$coefficients["Ha"]*DR3$RIah0
  RIfa1m2 <- DR3$mlfy01+RItmlem2$coefficients["Ha"]*DR3$RIah1
  RIf00m2 <- DR3$mlfy00+RItmlem2$coefficients["Ha"]*DR3$RIah0+RItmlem2$coefficients["Hla"]*DR3$RIlah0
  RIf11m2 <- DR3$mlfy11+RItmlem2$coefficients["Ha"]*DR3$RIah1+RItmlem2$coefficients["Hla"]*DR3$RIlah1
  RItmlelam2 <- mean(RIfla1m2)-mean(RIfla0m2)
  RItmlelam2IC <- (DR3$RIlah1+DR3$RIlah0)*(DR3$y-DR3$mlfy) + DR3$mlfy10 - DR3$mlfy00 - (RIfla1m2-RIfla0m2)
  RItmleam2 <- mean(RIfa1m2)-mean(RIfa0m2)  
  RItmleam2IC <- (DR3$RIah1+DR3$RIah0)*(DR3$y-DR3$mlfy) + DR3$mlfy01 - DR3$mlfy00 - (RIfa1m2-RIfa0m2)
  RItmlejtm2IC <- (DR3$RIlah1+DR3$RIlah0)*(DR3$y-DR3$mlfy) + (DR3$RIah1+DR3$RIah0)*(DR3$y-DR3$mlfy) + DR3$mlfy11 - DR3$mlfy00 - (RIf11m2-RIf00m2)
  # Both incorrectly specified
  DR4 <- merge(D[,c(1:3)],RIexppm[,c(1,2,13:18)])
  DR4 <- merge(DR4,RIgcmpdm)
  RItmlem3 <- glm(data=DR4,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0m3 <- DR4$mlfy00+RItmlem3$coefficients["Hla"]*DR4$RIlah0
  RIfla1m3 <- DR4$mlfy10+RItmlem3$coefficients["Hla"]*DR4$RIlah1
  RIfa0m3 <- DR4$mlfy00+RItmlem3$coefficients["Ha"]*DR4$RIah0
  RIfa1m3 <- DR4$mlfy01+RItmlem3$coefficients["Ha"]*DR4$RIah1
  RIf00m3 <- DR4$mlfy00+RItmlem3$coefficients["Ha"]*DR4$RIah0+RItmlem3$coefficients["Hla"]*DR4$RIlah0
  RIf11m3 <- DR4$mlfy11+RItmlem3$coefficients["Ha"]*DR4$RIah1+RItmlem3$coefficients["Hla"]*DR4$RIlah1
  RItmlelam3 <- mean(RIfla1m3)-mean(RIfla0m3)
  RItmlelam3IC <- (DR4$RIlah1+DR4$RIlah0)*(DR4$y-DR4$mlfy) + DR4$mlfy10 - DR4$mlfy00 - (RIfla1m3-RIfla0m3)
  RItmleam3 <- mean(RIfa1m3)-mean(RIfa0m3) 
  RItmleam3IC <- (DR4$RIah1+DR4$RIah0)*(DR4$y-DR4$mlfy) + DR4$mlfy01 - DR4$mlfy00 - (RIfa1m3-RIfa0m3)
  RItmlejtm3IC <- (DR4$RIlah1+DR4$RIlah0)*(DR4$y-DR4$mlfy) + (DR4$RIah1+DR4$RIah0)*(DR4$y-DR4$mlfy) + DR4$mlfy11 - DR4$mlfy00 - (RIf11m3-RIf00m3)
  
  # SuperLearner TMLE    
  DS <- merge(D[,c(1:3)],SLexpp[,c(1,2,13:18)])
  DS <- merge(DS,SLgcmpd)
  SLtmle <- glm(data=DS,y~-1+Hla+Ha+offset(slfy),family="gaussian")
  SLfla0 <- DS$slfy00+SLtmle$coefficients["Hla"]*DS$SLlah0
  SLfla1 <- DS$slfy10+SLtmle$coefficients["Hla"]*DS$SLlah1
  SLfa0 <- DS$slfy00+SLtmle$coefficients["Ha"]*DS$SLah0
  SLfa1 <- DS$slfy01+SLtmle$coefficients["Ha"]*DS$SLah1
  SLf00 <- DS$slfy00+SLtmle$coefficients["Ha"]*DS$SLah0+SLtmle$coefficients["Hla"]*DS$SLlah0
  SLf11 <- DS$slfy11+SLtmle$coefficients["Ha"]*DS$SLah1+SLtmle$coefficients["Hla"]*DS$SLlah1
  SLtmlela <- mean(SLfla1)-mean(SLfla0)
  SLtmlelaIC <- (DS$SLlah1+DS$SLlah0)*(DS$y-DS$slfy) + DS$slfy10 - DS$slfy00 - (SLfla1-SLfla0)
  SLtmlea <- mean(SLfa1)-mean(SLfa0)
  SLtmleaIC <- (DS$SLah1+DS$SLah0)*(DS$y-DS$slfy) + DS$slfy01 - DS$slfy00 - (SLfa1-SLfa0)
  SLtmlejtIC <- (DS$SLlah1+DS$SLlah0)*(DS$y-DS$slfy) + (DS$SLah1+DS$SLah0)*(DS$y-DS$slfy) + DS$slfy11 - DS$slfy00 - (SLf11-SLf00)
  
  c(coef(summary(GLiptw))[2,1],sqrt(GLiptwv[2,2]),
    coef(summary(GLiptw))[3,1],sqrt(GLiptwv[3,3]),
    coef(summary(GLiptw))[2,1]+coef(summary(GLiptw))[3,1],sqrt(GLiptwv[2,2] + GLiptwv[3,3] + 2*GLiptwv[2,3]),
    coef(summary(GLiptwm))[2,1],sqrt(GLiptwmv[2,2]),
    coef(summary(GLiptwm))[3,1],sqrt(GLiptwmv[3,3]),
    coef(summary(GLiptwm))[2,1]+coef(summary(GLiptwm))[3,1],sqrt(GLiptwmv[2,2] + GLiptwmv[3,3] + 2*GLiptwmv[2,3]),
    coef(summary(RIiptw))[2,1],coef(summary(RIiptw))[2,2],
    coef(summary(RIiptw))[3,1],coef(summary(RIiptw))[3,2],
    coef(summary(RIiptw))[2,1]+coef(summary(RIiptw))[3,1],sqrt(vcov(RIiptw)[2,2] + vcov(RIiptw)[3,3] + 2*vcov(RIiptw)[2,3]),
    coef(summary(RIiptwm))[2,1],coef(summary(RIiptwm))[2,2],
    coef(summary(RIiptwm))[3,1],coef(summary(RIiptwm))[3,2],
    coef(summary(RIiptwm))[2,1]+coef(summary(RIiptwm))[3,1],sqrt(vcov(RIiptwm)[2,2] + vcov(RIiptwm)[3,3] + 2*vcov(RIiptwm)[2,3]),
    coef(summary(SLiptw))[2,1],coef(summary(SLiptw))[2,2],
    coef(summary(SLiptw))[3,1],coef(summary(SLiptw))[3,2],
    coef(summary(SLiptw))[2,1]+coef(summary(SLiptw))[3,1],sqrt(vcov(SLiptw)[2,2] + vcov(SLiptw)[3,3] + 2*vcov(SLiptw)[2,3]),
    coef(summary(GLdriptw))[2,1],sqrt(GLdriptwv[2,2]),
    coef(summary(GLdriptw))[3,1],sqrt(GLdriptwv[3,3]),
    coef(summary(GLdriptw))[2,1]+coef(summary(GLdriptw))[3,1],sqrt(GLdriptwv[2,2] + GLdriptwv[3,3] + 2*GLdriptwv[2,3]),
    coef(summary(GLdriptwm1))[2,1],sqrt(GLdriptwm1v[2,2]),
    coef(summary(GLdriptwm1))[3,1],sqrt(GLdriptwm1v[3,3]),
    coef(summary(GLdriptwm1))[2,1]+coef(summary(GLdriptwm1))[3,1],sqrt(GLdriptwm1v[2,2] + GLdriptwm1v[3,3] + 2*GLdriptwm1v[2,3]),
    coef(summary(GLdriptwm2))[2,1],sqrt(GLdriptwm2v[2,2]),
    coef(summary(GLdriptwm2))[3,1],sqrt(GLdriptwm2v[3,3]),
    coef(summary(GLdriptwm2))[2,1]+coef(summary(GLdriptwm2))[3,1],sqrt(GLdriptwm2v[2,2] + GLdriptwm2v[3,3] + 2*GLdriptwm2v[2,3]),
    coef(summary(GLdriptwm3))[2,1],sqrt(GLdriptwm3v[2,2]),
    coef(summary(GLdriptwm3))[3,1],sqrt(GLdriptwm3v[3,3]),
    coef(summary(GLdriptwm3))[2,1]+coef(summary(GLdriptwm3))[3,1],sqrt(GLdriptwm3v[2,2] + GLdriptwm3v[3,3] + 2*GLdriptwm3v[2,3]),
    coef(summary(RIdriptw))[2,1],coef(summary(RIdriptw))[2,2],
    coef(summary(RIdriptw))[3,1],coef(summary(RIdriptw))[3,2],
    coef(summary(RIdriptw))[2,1]+coef(summary(RIdriptw))[3,1],sqrt(vcov(RIdriptw)[2,2] + vcov(RIdriptw)[3,3] + 2*vcov(RIdriptw)[2,3]),
    coef(summary(RIdriptwm1))[2,1],coef(summary(RIdriptwm1))[2,2],
    coef(summary(RIdriptwm1))[3,1],coef(summary(RIdriptwm1))[3,2],
    coef(summary(RIdriptwm1))[2,1]+coef(summary(RIdriptwm1))[3,1],sqrt(vcov(RIdriptwm1)[2,2] + vcov(RIdriptwm1)[3,3] + 2*vcov(RIdriptwm1)[2,3]),
    coef(summary(RIdriptwm2))[2,1],coef(summary(RIdriptwm2))[2,2],
    coef(summary(RIdriptwm2))[3,1],coef(summary(RIdriptwm2))[3,2],
    coef(summary(RIdriptwm2))[2,1]+coef(summary(RIdriptwm2))[3,1],sqrt(vcov(RIdriptwm2)[2,2] + vcov(RIdriptwm2)[3,3] + 2*vcov(RIdriptwm2)[2,3]),
    coef(summary(RIdriptwm3))[2,1],coef(summary(RIdriptwm3))[2,2],
    coef(summary(RIdriptwm3))[3,1],coef(summary(RIdriptwm3))[3,2],
    coef(summary(RIdriptwm3))[2,1]+coef(summary(RIdriptwm3))[3,1],sqrt(vcov(RIdriptwm3)[2,2] + vcov(RIdriptwm3)[3,3] + 2*vcov(RIdriptwm3)[2,3]),
    coef(summary(SLdriptw))[2,1],coef(summary(SLdriptw))[2,2],
    coef(summary(SLdriptw))[3,1],coef(summary(SLdriptw))[3,2],
    coef(summary(SLdriptw))[2,1]+coef(summary(SLdriptw))[3,1],sqrt(vcov(SLdriptw)[2,2] + vcov(SLdriptw)[3,3] + 2*vcov(SLdriptw)[2,3]),
    GLtmlea,var(GLtmleaIC)/250,GLtmlela,var(GLtmlelaIC)/250,GLtmlea+GLtmlela,var(GLtmlejtIC)/250,
    GLtmleam1,var(GLtmleam1IC)/250,GLtmlelam1,var(GLtmlelam1IC)/250,GLtmleam1+GLtmlelam1,var(GLtmlejtm1IC)/250,
    GLtmleam2,var(GLtmleam2IC)/250,GLtmlelam2,var(GLtmlelam2IC)/250,GLtmleam2+GLtmlelam2,var(GLtmlejtm2IC)/250,
    GLtmleam3,var(GLtmleam3IC)/250,GLtmlelam3,var(GLtmlelam3IC)/250,GLtmleam3+GLtmlelam3,var(GLtmlejtm3IC)/250,
    RItmlea,var(RItmleaIC)/250,RItmlela,var(RItmlelaIC)/250,RItmlea+RItmlela,var(RItmlejtIC)/250,
    RItmleam1,var(RItmleam1IC)/250,RItmlelam1,var(RItmlelam1IC)/250,RItmleam1+RItmlelam1,var(RItmlejtm1IC)/250,
    RItmleam2,var(RItmleam2IC)/250,RItmlelam2,var(RItmlelam2IC)/250,RItmleam2+RItmlelam2,var(RItmlejtm2IC)/250,
    RItmleam3,var(RItmleam3IC)/250,RItmlelam3,var(RItmlelam3IC)/250,RItmleam3+RItmlelam3,var(RItmlejtm3IC)/250,
    SLtmlea,var(SLtmleaIC)/250,SLtmlela,var(SLtmlelaIC)/250,SLtmlea+SLtmlela,var(SLtmlejtIC)/250)
}

sim <- function(data) {
  outmodel <- "y~a+la+l+ll+oa+ob+oc+obs"
  outmodel2 <- "y~a+la+ll+oa+ob+oc+obs"
  outmodelm <- "y~a+la+l+ll+oa+ob+oc"
  outmodel2m <- "y~a+la+ll+oa+ob+oc"
  
  propa0model <- "la~ll+oa+ob+oc+obs"
  propa1model <- "a~la+l+ll+oa+ob+oc+obs"
  propa0modelm <- "la~ll+oa+ob+oc"
  propa1modelm <- "a~la+l+ll+oa+ob+oc"
  
  GLM <- glm(outmodel,data=data,family=gaussian)
  GLMm <- glm(outmodelm,data=data,family=gaussian)
  RI <- lmer(paste0(outmodel,"+(1|id)"),data=data)
  RIm <- lmer(paste0(outmodelm,"+(1|id)"),data=data)
  GEE <- geeglm(formula(outmodel),data=data,id=data$id,waves=data$obs,family=gaussian)
  GEEm <- geeglm(formula(outmodelm),data=data,id=data$id,waves=data$obs,family=gaussian)
  
  SL.est <- compestimates(data=data)
  
  results1 <- matrix(ncol = 36)
  
  results1[1]<-coef(summary(GLM))[2,1]
  results1[2]<-sqrt(vcovCR(GLM, cluster=data$id, type = "CR3")[2,2])
  results1[3]<-coef(summary(GLM))[3,1]
  results1[4]<-sqrt(vcovCR(GLM, cluster=data$id, type = "CR3")[3,3])
  results1[5]<-coef(summary(GLM))[2,1]+coef(summary(GLM))[3,1]
  results1[6]<-sqrt(vcovCR(GLM, cluster=data$id, type = "CR3")[2,2] + vcovCR(GLM, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(GLM, cluster=data$id, type = "CR3")[2,3])
  
  results1[7]<-coef(summary(GLMm))[2,1]
  results1[8]<-sqrt(vcovCR(GLMm, cluster=data$id, type = "CR3")[2,2])
  results1[9]<-coef(summary(GLMm))[3,1]
  results1[10]<-sqrt(vcovCR(GLMm, cluster=data$id, type = "CR3")[3,3])
  results1[11]<-coef(summary(GLMm))[2,1]+coef(summary(GLMm))[3,1]
  results1[12]<-sqrt(vcovCR(GLMm, cluster=data$id, type = "CR3")[2,2] + vcovCR(GLMm, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(GLMm, cluster=data$id, type = "CR3")[2,3])
  
  results1[13]<-coef(summary(RI))[2,1]
  results1[14]<-coef(summary(RI))[2,2]
  results1[15]<-coef(summary(RI))[3,1]
  results1[16]<-coef(summary(RI))[3,2]
  results1[17]<-coef(summary(RI))[2,1]+coef(summary(RI))[3,1]
  results1[18]<-sqrt(vcov(RI)[2,2] + vcov(RI)[3,3] + 2*vcov(RI)[2,3])
  
  results1[19]<-coef(summary(RIm))[2,1]
  results1[20]<-coef(summary(RIm))[2,2]
  results1[21]<-coef(summary(RIm))[3,1]
  results1[22]<-coef(summary(RIm))[3,2]
  results1[23]<-coef(summary(RIm))[2,1]+coef(summary(RIm))[3,1]
  results1[24]<-sqrt(vcov(RIm)[2,2] + vcov(RIm)[3,3] + 2*vcov(RIm)[2,3])
  
  results1[25]<-coef(summary(GEE))[2,1]
  results1[26]<-coef(summary(GEE))[2,2]
  results1[27]<-coef(summary(GEE))[3,1]
  results1[28]<-coef(summary(GEE))[3,2]
  results1[29]<-coef(summary(GEE))[2,1]+coef(summary(GEE))[3,1]
  results1[30]<-sqrt(summary(GEE)$cov.scaled[2,2] + summary(GEE)$cov.scaled[3,3] + 2*summary(GEE)$cov.scaled[2,3])
  
  results1[31]<-coef(summary(GEEm))[2,1]
  results1[32]<-coef(summary(GEEm))[2,2]
  results1[33]<-coef(summary(GEEm))[3,1]
  results1[34]<-coef(summary(GEEm))[3,2]
  results1[35]<-coef(summary(GEEm))[2,1]+coef(summary(GEEm))[3,1]
  results1[36]<-sqrt(summary(GEEm)$cov.scaled[2,2] + summary(GEEm)$cov.scaled[3,3] + 2*summary(GEEm)$cov.scaled[2,3])
  
  results <- c(results1,SL.est)
  results
}

col.names <- c("GLMaco","GLMase","GLMlaco","GLMlase","GLMjtco","GLMjtse",
               "RIaco","RIase","RIlaco","RIlase","RIjtco","RIjtse",
               "GEEaco","GEEase","GEElaco","GEElase","GEEjtco","GEEjtse",
               "GLMiptwaco","GLMiptwase","GLMiptwlaco","GLMiptwlase","GLMiptwjtco","GLMiptwjtse",
               "RIiptwaco","RIiptwase","RIiptwlaco","RIiptwlase","RIiptwjtco","RIiptwjtse",
               "SLiptwaco","SLiptwase","SLiptwlaco","SLiptwlase","SLiptwjtco","SLiptwjtse",
               "GLMdriptwaco","GLMdriptwase","GLMdriptwlaco","GLMdriptwlase","GLMdriptwjtco","GLMdriptwjtse",
               "RIdriptwaco","RIdriptwase","RIdriptwlaco","RIdriptwlase","RIdriptwjtco","RIdriptwjtse",
               "SLdriptwaco","SLdriptwase","SLdriptwlaco","SLdriptwlase","SLdriptwjtco","SLdriptwjtse",
               "GLMtmleaco","GLMtmlease","GLMtmlelaco","GLMtmlelase","GLMtmlejtco","GLMtmlejtse",
               "RItmleaco","RItmlease","RItmlelaco","RItmlelase","RItmlejtco","RItmlejtse",
               "SLtmleaco","SLtmlease","SLtmlelaco","SLtmlelase","SLtmlejtco","SLtmlejtse")

numcores <- future::availableCores()
cl <- makeCluster(numcores, type="FORK")
parallel::clusterEvalQ(cl, .libPaths("~/RPackages"))
parallel::clusterEvalQ(cl, library("nnls"))
parallel::clusterEvalQ(cl, library("SuperLearner"))
parallel::clusterEvalQ(cl, library("Matrix"))
parallel::clusterEvalQ(cl, library("foreach"))
parallel::clusterEvalQ(cl, library("ranger"))
parallel::clusterEvalQ(cl, library("lme4"))
parallel::clusterEvalQ(cl, library("geepack"))
parallel::clusterEvalQ(cl, library("ggplot2"))
parallel::clusterEvalQ(cl, library("mvtnorm"))
parallel::clusterEvalQ(cl, library("survival"))
parallel::clusterEvalQ(cl, library("TH.data"))
parallel::clusterEvalQ(cl, library("MASS"))
parallel::clusterEvalQ(cl, library("multcomp"))
parallel::clusterEvalQ(cl, library("doBy"))
parallel::clusterEvalQ(cl, library("splines"))
parallel::clusterEvalQ(cl, library("gam"))
parallel::clusterEvalQ(cl, library("future"))
parallel::clusterEvalQ(cl, library("stats"))
parallel::clusterEvalQ(cl, library("data.table"))
parallel::clusterEvalQ(cl, library("optimx"))
parallel::clusterEvalQ(cl, library("clubSandwich"))
parallel::clusterExport(cl, c("bootdata","checkrange","compestimates","sim"))
registerDoParallel(cl)

set.seed(269012,kind="L'Ecuyer-CMRG")
start1 <- Sys.time()
results1 <- data.frame(matrix(unlist(parLapply(cl,data[1:250], function (x) {sim(x)})),nrow=250,ncol=174,byrow=TRUE))
end1 <- Sys.time()
end1-start1
start1 <- Sys.time()
results2 <- data.frame(matrix(unlist(parLapply(cl,data[251:500], function (x) {sim(x)})),nrow=250,ncol=174,byrow=TRUE))
end2 <- Sys.time()
end2-start2

resultsB<-rbind(results1,results2)
save(resultsB,file="resultsB.RData")
resultsB1<-results[,c(1:6,13:18,25:30,37:42,49:54,61:66,67:72,91:96,115:120,121:126,145:150,169:174)]
colnames(resultsB1) <- col.names
save(resultsB1,file="resultsB1.RData")
resultsB2<-results[,c(1:6,13:18,25:30,43:48,55:60,61:66,73:78,97:102,115:120,127:132,151:156,169:174)]
colnames(resultsB2) <- col.names
save(resultsB2,file="resultsB2.RData")
resultsB3<-results[,c(7:12,19:24,31:36,37:42,49:54,61:66,79:84,103:108,115:120,133:138,157:162,169:174)]
colnames(resultsB3) <- col.names
save(resultsB3,file="resultsB3.RData")
resultsB4<-results[,c(7:12,19:24,31:36,43:48,55:60,61:66,85:90,109:114,115:120,139:144,163:168,169:174)]
colnames(resultsB4) <- col.names
save(resultsB4,file="resultsB4.RData")

