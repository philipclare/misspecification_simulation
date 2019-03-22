.libPaths("~/RPackages")

libs <- c("nnls","SuperLearner","mvtnorm","survival","TH.data","MASS","Matrix",
          "foreach","ranger","lme4","geepack","doParallel","ggplot2",
          "multcomp","doBy","gam","future","data.table","optimx")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing])
}

load(file="data.RData")

bootdata <- function(d) {
  index<-sample(unique(d$id), length(unique(d$id)), repl=TRUE)
  data<-lapply(index,function (x) {d[d$id==x,]})
  for (o in 1:length(unique(d$id))){
    data[[o]]$id <- o
  }
  data<-do.call(rbind.data.frame, data)
}

checkrange <- function(v) {
  v <- ifelse(v<0.001,0.001,v)
  v <- ifelse(v>0.999,ifelse(is.na(v),v,0.999),v)
}

compestimates <- function(data,outmodel,outmodel2,propa0model,propa1model) {
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
                        glfy10=(predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                        glfy01=(predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLgcmpla <- mean(GLgcmpd$glfy10)-mean(GLgcmpd$glfy00)
  GLgcmpa <- mean(GLgcmpd$glfy01)-mean(GLgcmpd$glfy00)
  
  # Incorrectly specified GLM G-computation
  GLoutm <- glm(formula=outmodelm,D,family=gaussian)
  GLgcmpdm <- data.frame(id=D$id,obs=D$obs,
                         glfy=predict(GLoutm,newdata=D[,-3],type="response"),
                         glfy00=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                         glfy10=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                         glfy01=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLgcmplam <- mean(GLgcmpdm$glfy10)-mean(GLgcmpdm$glfy00)
  GLgcmpam <- mean(GLgcmpdm$glfy01)-mean(GLgcmpdm$glfy00)
  
  # Correctly specified Random-intercept G-computation
  RIlp <- glmer(l~la+ll+oa+ob+oc+(1|id),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIlprop <- data.frame(cbind(l0=predict(RIlp,newdata=data.frame(D[,c(-5,-6)],la=0),type="response"),l1=predict(RIlp,newdata=data.frame(D[,c(-5,-6)],la=1),type="response")))
  RIout <- lmer(formula=paste0(outmodel,"+(1|id)"),D)
  RIgcmpd <- data.frame(id=D$id,obs=D$obs,
                        mlfy=predict(RIout,newdata=D[,-3],type="response"),
                        mlfy00=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                        mlfy10=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                        mlfy01=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIgcmpla <- mean(RIgcmpd$mlfy10)-mean(RIgcmpd$mlfy00)
  RIgcmpa <- mean(RIgcmpd$mlfy01)-mean(RIgcmpd$mlfy00)
  
  # Correctly specified Random-intercept G-computation
  RIoutm <- lmer(formula=paste0(outmodelm,"+(1|id)"),D)
  RIgcmpdm <- data.frame(id=D$id,obs=D$obs,
                         mlfy=predict(RIoutm,newdata=D[,-3],type="response"),
                         mlfy00=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                         mlfy10=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                         mlfy01=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIgcmplam <- mean(RIgcmpdm$mlfy10)-mean(RIgcmpdm$mlfy00)
  RIgcmpam <- mean(RIgcmpdm$mlfy01)-mean(RIgcmpdm$mlfy00)
  
  # SuperLearner G-computation  
  SLlp <- SuperLearner(Y=as.vector(D$l),X=D[,c(-1,-3,-4,-6)],id=I,SL.library=SLlib,family=binomial)
  SLlprop <- data.frame(cbind(l0=as.vector(predict(SLlp,newdata=data.frame(D[,c(-3,-4,-5,-6)],la=0))$pred),l1=as.vector(predict(SLlp,newdata=data.frame(D[,c(-3,-4,-5,-6)],la=1))$pred)))
  SLout <- SuperLearner(Y=Y,X=D1,id=I,SL.library=SLlib,family=gaussian)
  SLgcmpd <- data.frame(id=D$id,obs=D$obs,
                        slfy=predict(SLout)$pred,
                        slfy00=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=0))$pred)*SLlprop$l0) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l0)),
                        slfy10=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=1,a=0))$pred)*SLlprop$l1) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l1)),
                        slfy01=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=1))$pred)*SLlprop$l0) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=1))$pred)*(1-SLlprop$l0)))
  SLgcmpla <- mean(SLgcmpd$slfy10)-mean(SLgcmpd$slfy00)
  SLgcmpa <- mean(SLgcmpd$slfy01)-mean(SLgcmpd$slfy00)
  
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
  
  # DR GLM-IPTW
  # Correctly specified
  GLdriptw <- glm(outmodel2,data=merge(D,GLexpp[,c(1,2,8)]),weight=GLwt,family=gaussian) 
  # Propensity incorrectly specified
  GLdriptwm1 <- glm(outmodel2,data=merge(D,GLexppm[,c(1,2,8)]),weight=GLwt,family=gaussian) 
  # Outcome incorrectly specified
  GLdriptwm2 <- glm(outmodel2m,data=merge(D,GLexpp[,c(1,2,8)]),weight=GLwt,family=gaussian) 
  # Both incorrectly specified
  GLdriptwm3 <- glm(outmodel2m,data=merge(D,GLexppm[,c(1,2,8)]),weight=GLwt,family=gaussian) 
  
  # Correctly specified Random-intercept IPTW and DR-IPTW
  RIexp0 <- glmer(formula=paste0(propa0model,"+(1|id)"),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIexp1 <- glmer(formula=paste0(propa1model,"+(1|id)"),data=D,family="binomial",nAGQ=15)
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
  RIexp1m <- glmer(formula=paste0(propa1modelm,"+(1|id)"),data=D,family="binomial",nAGQ=15)
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
  SLdriptw <- glm(y~a+la+ll+oa+ob+oc,data=merge(D,SLexpp[,c(1,2,8)]),weight=SLwt,family=gaussian)
  
  # GLM DR-GCOMP
  # Correctly specified
  GLdrout <- glm(formula=outmodel,data=merge(D,GLexpp[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdrgcmpd <- data.frame(id=D$id,obs=D$obs,
                          glfy=predict(GLdrout,newdata=D[,-3],type="response"),
                          glfy00=(predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                          glfy10=(predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                          glfy01=(predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLdrgcmpla <- mean(GLdrgcmpd$glfy10)-mean(GLdrgcmpd$glfy00)
  GLdrgcmpa <- mean(GLdrgcmpd$glfy01)-mean(GLdrgcmpd$glfy00)
  
  # Propensity incorrectly specified
  GLdroutm1 <- glm(formula=outmodel,data=merge(D,GLexppm[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdrgcmpdm1 <- data.frame(id=D$id,obs=D$obs,
                            glfy=predict(GLdroutm1,newdata=D[,-3],type="response"),
                            glfy00=(predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                            glfy10=(predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                            glfy01=(predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLdrgcmplam1 <- mean(GLdrgcmpdm1$glfy10)-mean(GLdrgcmpdm1$glfy00)
  GLdrgcmpam1 <- mean(GLdrgcmpdm1$glfy01)-mean(GLdrgcmpdm1$glfy00)
  
  # Outcome incorrectly specified
  GLdroutm2 <- glm(formula=outmodelm,data=merge(D,GLexpp[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdrgcmpdm2 <- data.frame(id=D$id,obs=D$obs,
                            glfy=predict(GLdroutm2,newdata=D[,-3],type="response"),
                            glfy00=(predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                            glfy10=(predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                            glfy01=(predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLdrgcmplam2 <- mean(GLdrgcmpdm2$glfy10)-mean(GLdrgcmpdm2$glfy00)
  GLdrgcmpam2 <- mean(GLdrgcmpdm2$glfy01)-mean(GLdrgcmpdm2$glfy00)
  
  # Both incorrectly specified
  GLdroutm3 <- glm(formula=outmodelm,data=merge(D,GLexppm[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdrgcmpdm3 <- data.frame(id=D$id,obs=D$obs,
                            glfy=predict(GLdroutm3,newdata=D[,-3],type="response"),
                            glfy00=(predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                            glfy10=(predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                            glfy01=(predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLdrgcmplam3 <- mean(GLdrgcmpdm3$glfy10)-mean(GLdrgcmpdm3$glfy00)
  GLdrgcmpam3 <- mean(GLdrgcmpdm3$glfy01)-mean(GLdrgcmpdm3$glfy00)
  
  # Random-intercept DR-GCOMP
  # Correctly specified
  RIdrout <- lmer(formula=paste0(outmodel,"+(1|id)"),data=merge(D,RIexpp[,c(1,2,8)]),weight=RIwt)
  RIdrgcmpd <- data.frame(id=D$id,obs=D$obs,
                          mlfy=predict(RIdrout,newdata=D[,-3],type="response"),
                          mlfy00=(predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                          mlfy10=(predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                          mlfy01=(predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIdrgcmpla <- mean(RIdrgcmpd$mlfy10)-mean(RIdrgcmpd$mlfy00)
  RIdrgcmpa <- mean(RIdrgcmpd$mlfy01)-mean(RIdrgcmpd$mlfy00)
  # Propensity incorrectly specified
  RIdroutm1 <- lmer(formula=paste0(outmodel,"+(1|id)"),data=merge(D,RIexppm[,c(1,2,8)]),weight=RIwt)
  RIdrgcmpdm1 <- data.frame(id=D$id,obs=D$obs,
                            mlfy=predict(RIdroutm1,newdata=D[,-3],type="response"),
                            mlfy00=(predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                            mlfy10=(predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                            mlfy01=(predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIdrgcmplam1 <- mean(RIdrgcmpdm1$mlfy10)-mean(RIdrgcmpdm1$mlfy00)
  RIdrgcmpam1 <- mean(RIdrgcmpdm1$mlfy01)-mean(RIdrgcmpdm1$mlfy00)
  # Outcome incorrectly specified
  RIdroutm2 <- lmer(formula=paste0(outmodelm,"+(1|id)"),data=merge(D,RIexpp[,c(1,2,8)]),weight=RIwt)
  RIdrgcmpdm2 <- data.frame(id=D$id,obs=D$obs,
                            mlfy=predict(RIdroutm2,newdata=D[,-3],type="response"),
                            mlfy00=(predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                            mlfy10=(predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                            mlfy01=(predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIdrgcmplam2 <- mean(RIdrgcmpdm2$mlfy10)-mean(RIdrgcmpdm2$mlfy00)
  RIdrgcmpam2 <- mean(RIdrgcmpdm2$mlfy01)-mean(RIdrgcmpdm2$mlfy00)
  # Both incorrectly specified
  RIdroutm3 <- lmer(formula=paste0(outmodelm,"+(1|id)"),data=merge(D,RIexppm[,c(1,2,8)]),weight=RIwt)
  RIdrgcmpdm3 <- data.frame(id=D$id,obs=D$obs,
                            mlfy=predict(RIdroutm3,newdata=D[,-3],type="response"),
                            mlfy00=(predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                            mlfy10=(predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                            mlfy01=(predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIdrgcmplam3 <- mean(RIdrgcmpdm3$mlfy10)-mean(RIdrgcmpdm3$mlfy00)
  RIdrgcmpam3 <- mean(RIdrgcmpdm3$mlfy01)-mean(RIdrgcmpdm3$mlfy00)
  
  # SuperLearner DR-GOMP  
  SLdrout <- SuperLearner(Y=Y,X=D1,id=I,SL.library=SLlib,obsWeights=SLexpp$SLwt,family=gaussian)
  SLdrgcmpd <- data.frame(id=D$id,obs=D$obs,
                          slfy=predict(SLdrout)$pred,
                          slfy00=(as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=0))$pred)*SLlprop$l0) + (as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l0)),
                          slfy10=(as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=1,a=0))$pred)*SLlprop$l1) + (as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l1)),
                          slfy01=(as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=1))$pred)*SLlprop$l0) + (as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=1))$pred)*(1-SLlprop$l0)))
  SLdrgcmpla <- mean(SLdrgcmpd$slfy10)-mean(SLdrgcmpd$slfy00)
  SLdrgcmpa <- mean(SLdrgcmpd$slfy01)-mean(SLdrgcmpd$slfy00)
  
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
  GLtmlela <- mean(GLfla1)-mean(GLfla0)
  GLtmlea <- mean(GLfa1)-mean(GLfa0)  
  # Propensity model incorrectly specified
  DG2 <- merge(D[,c(1:3)],GLexppm[,c(1,2,13:18)])
  DG2 <- merge(DG2,GLgcmpd)
  GLtmlem1 <- glm(data=DG2,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0m1 <- DG2$glfy00+GLtmlem1$coefficients["Hla"]*DG2$GLlah0
  GLfla1m1 <- DG2$glfy10+GLtmlem1$coefficients["Hla"]*DG2$GLlah1
  GLfa0m1 <- DG2$glfy00+GLtmlem1$coefficients["Ha"]*DG2$GLah0
  GLfa1m1 <- DG2$glfy01+GLtmlem1$coefficients["Ha"]*DG2$GLah1
  GLtmlelam1 <- mean(GLfla1m1)-mean(GLfla0m1)
  GLtmleam1 <- mean(GLfa1m1)-mean(GLfa0m1)  
  # Outcome model incorrectly specified
  DG3 <- merge(D[,c(1:3)],GLexpp[,c(1,2,13:18)])
  DG3 <- merge(DG3,GLgcmpdm)
  GLtmlem2 <- glm(data=DG3,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0m2 <- DG3$glfy00+GLtmlem2$coefficients["Hla"]*DG3$GLlah0
  GLfla1m2 <- DG3$glfy10+GLtmlem2$coefficients["Hla"]*DG3$GLlah1
  GLfa0m2 <- DG3$glfy00+GLtmlem2$coefficients["Ha"]*DG3$GLah0
  GLfa1m2 <- DG3$glfy01+GLtmlem2$coefficients["Ha"]*DG3$GLah1
  GLtmlelam2 <- mean(GLfla1m2)-mean(GLfla0m2)
  GLtmleam2 <- mean(GLfa1m2)-mean(GLfa0m2)  
  # Both incorrectly specified
  DG4 <- merge(D[,c(1:3)],GLexppm[,c(1,2,13:18)])
  DG4 <- merge(DG4,GLgcmpdm)
  GLtmlem3 <- glm(data=DG4,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0m3 <- DG4$glfy00+GLtmlem3$coefficients["Hla"]*DG4$GLlah0
  GLfla1m3 <- DG4$glfy10+GLtmlem3$coefficients["Hla"]*DG4$GLlah1
  GLfa0m3 <- DG4$glfy00+GLtmlem3$coefficients["Ha"]*DG4$GLah0
  GLfa1m3 <- DG4$glfy01+GLtmlem3$coefficients["Ha"]*DG4$GLah1
  GLtmlelam3 <- mean(GLfla1m3)-mean(GLfla0m3)
  GLtmleam3 <- mean(GLfa1m3)-mean(GLfa0m3)  
  
  # Random Intercept TMLE
  # Correctly specified
  DR1 <- merge(D[,c(1:3)],RIexpp[,c(1,2,13:18)])
  DR1 <- merge(DR1,RIgcmpd)
  RItmle <- glm(data=DR1,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0 <- DR1$mlfy00+RItmle$coefficients["Hla"]*DR1$RIlah0
  RIfla1 <- DR1$mlfy10+RItmle$coefficients["Hla"]*DR1$RIlah1
  RIfa0 <- DR1$mlfy00+RItmle$coefficients["Ha"]*DR1$RIah0
  RIfa1 <- DR1$mlfy01+RItmle$coefficients["Ha"]*DR1$RIah1
  RItmlela <- mean(RIfla1)-mean(RIfla0)
  RItmlea <- mean(RIfa1)-mean(RIfa0)
  # Propensity model incorrectly specified
  DR2 <- merge(D[,c(1:3)],RIexppm[,c(1,2,13:18)])
  DR2 <- merge(DR2,RIgcmpd)
  RItmlem1 <- glm(data=DR2,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0m1 <- DR2$mlfy00+RItmlem1$coefficients["Hla"]*DR2$RIlah0
  RIfla1m1 <- DR2$mlfy10+RItmlem1$coefficients["Hla"]*DR2$RIlah1
  RIfa0m1 <- DR2$mlfy00+RItmlem1$coefficients["Ha"]*DR2$RIah0
  RIfa1m1 <- DR2$mlfy01+RItmlem1$coefficients["Ha"]*DR2$RIah1
  RItmlelam1 <- mean(RIfla1m1)-mean(RIfla0m1)
  RItmleam1 <- mean(RIfa1m1)-mean(RIfa0m1)  
  # Outcome model incorrectly specified
  DR3 <- merge(D[,c(1:3)],RIexpp[,c(1,2,13:18)])
  DR3 <- merge(DR3,RIgcmpdm)
  RItmlem2 <- glm(data=DR3,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0m2 <- DR3$mlfy00+RItmlem2$coefficients["Hla"]*DR3$RIlah0
  RIfla1m2 <- DR3$mlfy10+RItmlem2$coefficients["Hla"]*DR3$RIlah1
  RIfa0m2 <- DR3$mlfy00+RItmlem2$coefficients["Ha"]*DR3$RIah0
  RIfa1m2 <- DR3$mlfy01+RItmlem2$coefficients["Ha"]*DR3$RIah1
  RItmlelam2 <- mean(RIfla1m2)-mean(RIfla0m2)
  RItmleam2 <- mean(RIfa1m2)-mean(RIfa0m2)  
  # Both incorrectly specified
  DR4 <- merge(D[,c(1:3)],RIexppm[,c(1,2,13:18)])
  DR4 <- merge(DR4,RIgcmpdm)
  RItmlem3 <- glm(data=DR4,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0m3 <- DR4$mlfy00+RItmlem3$coefficients["Hla"]*DR4$RIlah0
  RIfla1m3 <- DR4$mlfy10+RItmlem3$coefficients["Hla"]*DR4$RIlah1
  RIfa0m3 <- DR4$mlfy00+RItmlem3$coefficients["Ha"]*DR4$RIah0
  RIfa1m3 <- DR4$mlfy01+RItmlem3$coefficients["Ha"]*DR4$RIah1
  RItmlelam3 <- mean(RIfla1m3)-mean(RIfla0m3)
  RItmleam3 <- mean(RIfa1m3)-mean(RIfa0m3)  
  
  # SuperLearner TMLE    
  DS <- merge(D[,c(1:3)],SLexpp[,c(1,2,13:18)])
  DS <- merge(DS,SLgcmpd)
  SLtmle <- glm(data=DS,y~-1+Hla+Ha+offset(slfy),family="gaussian")
  SLfla0 <- DS$slfy00+SLtmle$coefficients["Hla"]*DS$SLlah0
  SLfla1 <- DS$slfy10+SLtmle$coefficients["Hla"]*DS$SLlah1
  SLfa0 <- DS$slfy00+SLtmle$coefficients["Ha"]*DS$SLah0
  SLfa1 <- DS$slfy01+SLtmle$coefficients["Ha"]*DS$SLah1
  SLtmlela <- mean(SLfla1)-mean(SLfla0)
  SLtmlea <- mean(SLfa1)-mean(SLfa0)
  
  c(coef(summary(GLiptw))[2,1],coef(summary(GLiptw))[3,1],coef(summary(GLiptw))[2,1]+coef(summary(GLiptw))[3,1],
    coef(summary(GLiptwm))[2,1],coef(summary(GLiptwm))[3,1],coef(summary(GLiptwm))[2,1]+coef(summary(GLiptwm))[3,1],
    coef(summary(RIiptw))[2,1],coef(summary(RIiptw))[3,1],coef(summary(RIiptw))[2,1]+coef(summary(RIiptw))[3,1],
    coef(summary(RIiptwm))[2,1],coef(summary(RIiptwm))[3,1],coef(summary(RIiptwm))[2,1]+coef(summary(RIiptwm))[3,1],
    coef(summary(SLiptw))[2,1],coef(summary(SLiptw))[3,1],coef(summary(SLiptw))[2,1]+coef(summary(SLiptw))[3,1],
    coef(summary(GLdriptw))[2,1],coef(summary(GLdriptw))[3,1],coef(summary(GLdriptw))[2,1]+coef(summary(GLdriptw))[3,1],
    coef(summary(GLdriptwm1))[2,1],coef(summary(GLdriptwm1))[3,1],coef(summary(GLdriptwm1))[2,1]+coef(summary(GLdriptwm1))[3,1],
    coef(summary(GLdriptwm2))[2,1],coef(summary(GLdriptwm2))[3,1],coef(summary(GLdriptwm2))[2,1]+coef(summary(GLdriptwm2))[3,1],
    coef(summary(GLdriptwm3))[2,1],coef(summary(GLdriptwm3))[3,1],coef(summary(GLdriptwm3))[2,1]+coef(summary(GLdriptwm3))[3,1],
    coef(summary(RIdriptw))[2,1],coef(summary(RIdriptw))[3,1],coef(summary(RIdriptw))[2,1]+coef(summary(RIdriptw))[3,1],
    coef(summary(RIdriptwm1))[2,1],coef(summary(RIdriptwm1))[3,1],coef(summary(RIdriptwm1))[2,1]+coef(summary(RIdriptwm1))[3,1],
    coef(summary(RIdriptwm2))[2,1],coef(summary(RIdriptwm2))[3,1],coef(summary(RIdriptwm2))[2,1]+coef(summary(RIdriptwm2))[3,1],
    coef(summary(RIdriptwm3))[2,1],coef(summary(RIdriptwm3))[3,1],coef(summary(RIdriptwm3))[2,1]+coef(summary(RIdriptwm3))[3,1],
    coef(summary(SLdriptw))[2,1],coef(summary(SLdriptw))[3,1],coef(summary(SLdriptw))[2,1]+coef(summary(SLdriptw))[3,1],
    GLgcmpa,GLgcmpla,GLgcmpa+GLgcmpla,
    GLgcmpam,GLgcmplam,GLgcmpam+GLgcmplam,
    RIgcmpa,RIgcmpla,RIgcmpa+RIgcmpla,
    RIgcmpam,RIgcmplam,RIgcmpam+RIgcmplam,
    SLgcmpa,SLgcmpla,SLgcmpa+SLgcmpla,
    GLdrgcmpa,GLdrgcmpla,GLdrgcmpa+GLdrgcmpla,
    GLdrgcmpam1,GLdrgcmplam1,GLdrgcmpam1+GLdrgcmplam1,
    GLdrgcmpam2,GLdrgcmplam2,GLdrgcmpam2+GLdrgcmplam2,
    GLdrgcmpam3,GLdrgcmplam3,GLdrgcmpam3+GLdrgcmplam3,
    RIdrgcmpa,RIdrgcmpla,RIdrgcmpa+RIdrgcmpla,
    RIdrgcmpam1,RIdrgcmplam1,RIdrgcmpam1+RIdrgcmplam1,
    RIdrgcmpam2,RIdrgcmplam2,RIdrgcmpam2+RIdrgcmplam2,
    RIdrgcmpam3,RIdrgcmplam3,RIdrgcmpam3+RIdrgcmplam3,
    SLdrgcmpa,SLdrgcmpla,SLdrgcmpa+SLdrgcmpla,
    GLtmlea,GLtmlela,GLtmlea+GLtmlela,
    GLtmleam1,GLtmlelam1,GLtmleam1+GLtmlelam1,
    GLtmleam2,GLtmlelam2,GLtmleam2+GLtmlelam2,
    GLtmleam3,GLtmlelam3,GLtmleam3+GLtmlelam3,
    RItmlea,RItmlela,RItmlea+RItmlela,
    RItmleam1,RItmlelam1,RItmleam1+RItmlelam1,
    RItmleam2,RItmlelam2,RItmleam2+RItmlelam2,
    RItmleam3,RItmlelam3,RItmleam3+RItmlelam3,
    SLtmlea,SLtmlela,SLtmlea+SLtmlela,
    coef(summary(GLiptw))[2,2],coef(summary(GLiptw))[3,2],sqrt(vcov(GLiptw)[2,2] + vcov(GLiptw)[3,3] + 2*vcov(GLiptw)[2,3]),
    coef(summary(GLiptwm))[2,2],coef(summary(GLiptwm))[3,2],sqrt(vcov(GLiptwm)[2,2] + vcov(GLiptwm)[3,3] + 2*vcov(GLiptwm)[2,3]),
    coef(summary(RIiptw))[2,2],coef(summary(RIiptw))[3,2],sqrt(vcov(RIiptw)[2,2] + vcov(RIiptw)[3,3] + 2*vcov(RIiptw)[2,3]),
    coef(summary(RIiptwm))[2,2],coef(summary(RIiptwm))[3,2],sqrt(vcov(RIiptwm)[2,2] + vcov(RIiptwm)[3,3] + 2*vcov(RIiptwm)[2,3]),
    coef(summary(SLiptw))[2,2],coef(summary(SLiptw))[3,2],sqrt(vcov(SLiptw)[2,2] + vcov(SLiptw)[3,3] + 2*vcov(SLiptw)[2,3]),
    coef(summary(GLdriptw))[2,2],coef(summary(GLdriptw))[3,2],sqrt(vcov(GLdriptw)[2,2] + vcov(GLdriptw)[3,3] + 2*vcov(GLdriptw)[2,3]),
    coef(summary(GLdriptwm1))[2,2],coef(summary(GLdriptwm1))[3,2],sqrt(vcov(GLdriptwm1)[2,2] + vcov(GLdriptwm1)[3,3] + 2*vcov(GLdriptwm1)[2,3]),
    coef(summary(GLdriptwm2))[2,2],coef(summary(GLdriptwm2))[3,2],sqrt(vcov(GLdriptwm2)[2,2] + vcov(GLdriptwm2)[3,3] + 2*vcov(GLdriptwm2)[2,3]),
    coef(summary(GLdriptwm3))[2,2],coef(summary(GLdriptwm3))[3,2],sqrt(vcov(GLdriptwm3)[2,2] + vcov(GLdriptwm3)[3,3] + 2*vcov(GLdriptwm3)[2,3]),
    coef(summary(RIdriptw))[2,2],coef(summary(RIdriptw))[3,2],sqrt(vcov(RIdriptw)[2,2] + vcov(RIdriptw)[3,3] + 2*vcov(RIdriptw)[2,3]),
    coef(summary(RIdriptwm1))[2,2],coef(summary(RIdriptwm1))[3,2],sqrt(vcov(RIdriptwm1)[2,2] + vcov(RIdriptwm1)[3,3] + 2*vcov(RIdriptwm1)[2,3]),
    coef(summary(RIdriptwm2))[2,2],coef(summary(RIdriptwm2))[3,2],sqrt(vcov(RIdriptwm2)[2,2] + vcov(RIdriptwm2)[3,3] + 2*vcov(RIdriptwm2)[2,3]),
    coef(summary(RIdriptwm3))[2,2],coef(summary(RIdriptwm3))[3,2],sqrt(vcov(RIdriptwm3)[2,2] + vcov(RIdriptwm3)[3,3] + 2*vcov(RIdriptwm3)[2,3]),
    coef(summary(SLdriptw))[2,2],coef(summary(SLdriptw))[3,2],sqrt(vcov(SLdriptw)[2,2] + vcov(SLdriptw)[3,3] + 2*vcov(SLdriptw)[2,3]))
}

comperror <- function(data,outmodel,outmodel2,propa0model,propa1model) {
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
                        glfy10=(predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                        glfy01=(predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLgcmpla <- mean(GLgcmpd$glfy10)-mean(GLgcmpd$glfy00)
  GLgcmpa <- mean(GLgcmpd$glfy01)-mean(GLgcmpd$glfy00)
  
  # Incorrectly specified GLM G-computation
  GLoutm <- glm(formula=outmodelm,D,family=gaussian)
  GLgcmpdm <- data.frame(id=D$id,obs=D$obs,
                         glfy=predict(GLoutm,newdata=D[,-3],type="response"),
                         glfy00=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                         glfy10=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                         glfy01=(predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLgcmplam <- mean(GLgcmpdm$glfy10)-mean(GLgcmpdm$glfy00)
  GLgcmpam <- mean(GLgcmpdm$glfy01)-mean(GLgcmpdm$glfy00)
  
  # Correctly specified Random-intercept G-computation
  RIlp <- glmer(l~la+ll+oa+ob+oc+(1|id),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIlprop <- data.frame(cbind(l0=predict(RIlp,newdata=data.frame(D[,c(-5,-6)],la=0),type="response"),l1=predict(RIlp,newdata=data.frame(D[,c(-5,-6)],la=1),type="response")))
  RIout <- lmer(formula=paste0(outmodel,"+(1|id)"),D)
  RIgcmpd <- data.frame(id=D$id,obs=D$obs,
                        mlfy=predict(RIout,newdata=D[,-3],type="response"),
                        mlfy00=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                        mlfy10=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                        mlfy01=(predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIgcmpla <- mean(RIgcmpd$mlfy10)-mean(RIgcmpd$mlfy00)
  RIgcmpa <- mean(RIgcmpd$mlfy01)-mean(RIgcmpd$mlfy00)
  
  # Correctly specified Random-intercept G-computation
  RIoutm <- lmer(formula=paste0(outmodelm,"+(1|id)"),D)
  RIgcmpdm <- data.frame(id=D$id,obs=D$obs,
                         mlfy=predict(RIoutm,newdata=D[,-3],type="response"),
                         mlfy00=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                         mlfy10=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                         mlfy01=(predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIoutm,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIgcmplam <- mean(RIgcmpdm$mlfy10)-mean(RIgcmpdm$mlfy00)
  RIgcmpam <- mean(RIgcmpdm$mlfy01)-mean(RIgcmpdm$mlfy00)
  
  # SuperLearner G-computation  
  SLlp <- SuperLearner(Y=as.vector(D$l),X=D[,c(-1,-3,-4,-6)],id=I,SL.library=SLlib,family=binomial)
  SLlprop <- data.frame(cbind(l0=as.vector(predict(SLlp,newdata=data.frame(D[,c(-3,-4,-5,-6)],la=0))$pred),l1=as.vector(predict(SLlp,newdata=data.frame(D[,c(-5,-6)],la=1))$pred)))
  SLout <- SuperLearner(Y=Y,X=D1,id=I,SL.library=SLlib,family=gaussian)
  SLgcmpd <- data.frame(id=D$id,obs=D$obs,
                        slfy=predict(SLout)$pred,
                        slfy00=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=0))$pred)*SLlprop$l0) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l0)),
                        slfy10=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=1,a=0))$pred)*SLlprop$l1) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l1)),
                        slfy01=(as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=1))$pred)*SLlprop$l0) + (as.vector(predict(SLout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=1))$pred)*(1-SLlprop$l0)))
  SLgcmpla <- mean(SLgcmpd$slfy10)-mean(SLgcmpd$slfy00)
  SLgcmpa <- mean(SLgcmpd$slfy01)-mean(SLgcmpd$slfy00)
  
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
  
  # Correctly specified Random-intercept IPTW and DR-IPTW
  RIexp0 <- glmer(formula=paste0(propa0model,"+(1|id)"),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIexp1 <- glmer(formula=paste0(propa1model,"+(1|id)"),data=D,family="binomial",nAGQ=15)
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
  
  # Incorrectly specified Random-intercept IPTW and DR-IPTW
  RIexp0m <- glmer(formula=paste0(propa0modelm,"+(1|id)"),data=D,family="binomial",nAGQ=15,control=glmercont)
  RIexp1m <- glmer(formula=paste0(propa1modelm,"+(1|id)"),data=D,family="binomial",nAGQ=15)
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
  
  # GLM DR-GCOMP
  # Correctly specified
  GLdrout <- glm(formula=outmodel,data=merge(D,GLexpp[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdrgcmpd <- data.frame(id=D$id,obs=D$obs,
                          glfy=predict(GLdrout,newdata=D[,-3],type="response"),
                          glfy00=(predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                          glfy10=(predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                          glfy01=(predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLdrgcmpla <- mean(GLdrgcmpd$glfy10)-mean(GLdrgcmpd$glfy00)
  GLdrgcmpa <- mean(GLdrgcmpd$glfy01)-mean(GLdrgcmpd$glfy00)
  
  # Propensity incorrectly specified
  GLdroutm1 <- glm(formula=outmodel,data=merge(D,GLexppm[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdrgcmpdm1 <- data.frame(id=D$id,obs=D$obs,
                            glfy=predict(GLdroutm1,newdata=D[,-3],type="response"),
                            glfy00=(predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                            glfy10=(predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                            glfy01=(predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLdrgcmplam1 <- mean(GLdrgcmpdm1$glfy10)-mean(GLdrgcmpdm1$glfy00)
  GLdrgcmpam1 <- mean(GLdrgcmpdm1$glfy01)-mean(GLdrgcmpdm1$glfy00)
  
  # Outcome incorrectly specified
  GLdroutm2 <- glm(formula=outmodelm,data=merge(D,GLexpp[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdrgcmpdm2 <- data.frame(id=D$id,obs=D$obs,
                            glfy=predict(GLdroutm2,newdata=D[,-3],type="response"),
                            glfy00=(predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                            glfy10=(predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                            glfy01=(predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLdrgcmplam2 <- mean(GLdrgcmpdm2$glfy10)-mean(GLdrgcmpdm2$glfy00)
  GLdrgcmpam2 <- mean(GLdrgcmpdm2$glfy01)-mean(GLdrgcmpdm2$glfy00)
  
  # Both incorrectly specified
  GLdroutm3 <- glm(formula=outmodelm,data=merge(D,GLexppm[,c(1,2,8)]),weight=GLwt,family=gaussian)
  GLdrgcmpdm3 <- data.frame(id=D$id,obs=D$obs,
                            glfy=predict(GLdroutm3,newdata=D[,-3],type="response"),
                            glfy00=(predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*GLlprop$l0) + (predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l0)),
                            glfy10=(predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*GLlprop$l1) + (predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-GLlprop$l1)),
                            glfy01=(predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*GLlprop$l0) + (predict(GLdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-GLlprop$l0)))
  GLdrgcmplam3 <- mean(GLdrgcmpdm3$glfy10)-mean(GLdrgcmpdm3$glfy00)
  GLdrgcmpam3 <- mean(GLdrgcmpdm3$glfy01)-mean(GLdrgcmpdm3$glfy00)
  
  # Random-intercept DR-GCOMP
  # Correctly specified
  RIdrout <- lmer(formula=paste0(outmodel,"+(1|id)"),data=merge(D,RIexpp[,c(1,2,8)]),weight=RIwt)
  RIdrgcmpd <- data.frame(id=D$id,obs=D$obs,
                          mlfy=predict(RIdrout,newdata=D[,-3],type="response"),
                          mlfy00=(predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                          mlfy10=(predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                          mlfy01=(predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIdrgcmpla <- mean(RIdrgcmpd$mlfy10)-mean(RIdrgcmpd$mlfy00)
  RIdrgcmpa <- mean(RIdrgcmpd$mlfy01)-mean(RIdrgcmpd$mlfy00)
  # Propensity incorrectly specified
  RIdroutm1 <- lmer(formula=paste0(outmodel,"+(1|id)"),data=merge(D,RIexppm[,c(1,2,8)]),weight=RIwt)
  RIdrgcmpdm1 <- data.frame(id=D$id,obs=D$obs,
                            mlfy=predict(RIdroutm1,newdata=D[,-3],type="response"),
                            mlfy00=(predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                            mlfy10=(predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                            mlfy01=(predict(RIdroutm1,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIdrout,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIdrgcmplam1 <- mean(RIdrgcmpdm1$mlfy10)-mean(RIdrgcmpdm1$mlfy00)
  RIdrgcmpam1 <- mean(RIdrgcmpdm1$mlfy01)-mean(RIdrgcmpdm1$mlfy00)
  # Outcome incorrectly specified
  RIdroutm2 <- lmer(formula=paste0(outmodelm,"+(1|id)"),data=merge(D,RIexpp[,c(1,2,8)]),weight=RIwt)
  RIdrgcmpdm2 <- data.frame(id=D$id,obs=D$obs,
                            mlfy=predict(RIdroutm2,newdata=D[,-3],type="response"),
                            mlfy00=(predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                            mlfy10=(predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                            mlfy01=(predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIdroutm2,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIdrgcmplam2 <- mean(RIdrgcmpdm2$mlfy10)-mean(RIdrgcmpdm2$mlfy00)
  RIdrgcmpam2 <- mean(RIdrgcmpdm2$mlfy01)-mean(RIdrgcmpdm2$mlfy00)
  # Both incorrectly specified
  RIdroutm3 <- lmer(formula=paste0(outmodelm,"+(1|id)"),data=merge(D,RIexppm[,c(1,2,8)]),weight=RIwt)
  RIdrgcmpdm3 <- data.frame(id=D$id,obs=D$obs,
                            mlfy=predict(RIdroutm3,newdata=D[,-3],type="response"),
                            mlfy00=(predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=0),type="response")*RIlprop$l0) + (predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l0)),
                            mlfy10=(predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=1,a=0),type="response")*RIlprop$l1) + (predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=0),type="response")*(1-RIlprop$l1)),
                            mlfy01=(predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=1,la=0,a=1),type="response")*RIlprop$l0) + (predict(RIdroutm3,newdata=data.frame(D[,c(-4,-5,-6)],l=0,la=0,a=1),type="response")*(1-RIlprop$l0)))
  RIdrgcmplam3 <- mean(RIdrgcmpdm3$mlfy10)-mean(RIdrgcmpdm3$mlfy00)
  RIdrgcmpam3 <- mean(RIdrgcmpdm3$mlfy01)-mean(RIdrgcmpdm3$mlfy00)
  
  # SuperLearner DR-GOMP  
  SLdrout <- SuperLearner(Y=Y,X=D1,id=I,SL.library=SLlib,obsWeights=SLexpp$SLwt,family=gaussian)
  SLdrgcmpd <- data.frame(id=D$id,obs=D$obs,
                          slfy=predict(SLdrout)$pred,
                          slfy00=(as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=0))$pred)*SLlprop$l0) + (as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l0)),
                          slfy10=(as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=1,a=0))$pred)*SLlprop$l1) + (as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=0))$pred)*(1-SLlprop$l1)),
                          slfy01=(as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=1,la=0,a=1))$pred)*SLlprop$l0) + (as.vector(predict(SLdrout,newdata=data.frame(D1[,c(-2,-3,-4)],l=0,la=0,a=1))$pred)*(1-SLlprop$l0)))
  SLdrgcmpla <- mean(SLdrgcmpd$slfy10)-mean(SLdrgcmpd$slfy00)
  SLdrgcmpa <- mean(SLdrgcmpd$slfy01)-mean(SLdrgcmpd$slfy00)
  
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
  GLtmlela <- mean(GLfla1)-mean(GLfla0)
  GLtmlea <- mean(GLfa1)-mean(GLfa0)  
  # Propensity model incorrectly specified
  DG2 <- merge(D[,c(1:3)],GLexppm[,c(1,2,13:18)])
  DG2 <- merge(DG2,GLgcmpd)
  GLtmlem1 <- glm(data=DG2,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0m1 <- DG2$glfy00+GLtmlem1$coefficients["Hla"]*DG2$GLlah0
  GLfla1m1 <- DG2$glfy10+GLtmlem1$coefficients["Hla"]*DG2$GLlah1
  GLfa0m1 <- DG2$glfy00+GLtmlem1$coefficients["Ha"]*DG2$GLah0
  GLfa1m1 <- DG2$glfy01+GLtmlem1$coefficients["Ha"]*DG2$GLah1
  GLtmlelam1 <- mean(GLfla1m1)-mean(GLfla0m1)
  GLtmleam1 <- mean(GLfa1m1)-mean(GLfa0m1)  
  # Outcome model incorrectly specified
  DG3 <- merge(D[,c(1:3)],GLexpp[,c(1,2,13:18)])
  DG3 <- merge(DG3,GLgcmpdm)
  GLtmlem2 <- glm(data=DG3,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0m2 <- DG3$glfy00+GLtmlem2$coefficients["Hla"]*DG3$GLlah0
  GLfla1m2 <- DG3$glfy10+GLtmlem2$coefficients["Hla"]*DG3$GLlah1
  GLfa0m2 <- DG3$glfy00+GLtmlem2$coefficients["Ha"]*DG3$GLah0
  GLfa1m2 <- DG3$glfy01+GLtmlem2$coefficients["Ha"]*DG3$GLah1
  GLtmlelam2 <- mean(GLfla1m2)-mean(GLfla0m2)
  GLtmleam2 <- mean(GLfa1m2)-mean(GLfa0m2)  
  # Both incorrectly specified
  DG4 <- merge(D[,c(1:3)],GLexppm[,c(1,2,13:18)])
  DG4 <- merge(DG4,GLgcmpdm)
  GLtmlem3 <- glm(data=DG4,y~-1+Hla+Ha+offset(glfy),family="gaussian")
  GLfla0m3 <- DG4$glfy00+GLtmlem3$coefficients["Hla"]*DG4$GLlah0
  GLfla1m3 <- DG4$glfy10+GLtmlem3$coefficients["Hla"]*DG4$GLlah1
  GLfa0m3 <- DG4$glfy00+GLtmlem3$coefficients["Ha"]*DG4$GLah0
  GLfa1m3 <- DG4$glfy01+GLtmlem3$coefficients["Ha"]*DG4$GLah1
  GLtmlelam3 <- mean(GLfla1m3)-mean(GLfla0m3)
  GLtmleam3 <- mean(GLfa1m3)-mean(GLfa0m3)  
  
  # Random Intercept TMLE
  # Correctly specified
  DR1 <- merge(D[,c(1:3)],RIexpp[,c(1,2,13:18)])
  DR1 <- merge(DR1,RIgcmpd)
  RItmle <- glm(data=DR1,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0 <- DR1$mlfy00+RItmle$coefficients["Hla"]*DR1$RIlah0
  RIfla1 <- DR1$mlfy10+RItmle$coefficients["Hla"]*DR1$RIlah1
  RIfa0 <- DR1$mlfy00+RItmle$coefficients["Ha"]*DR1$RIah0
  RIfa1 <- DR1$mlfy01+RItmle$coefficients["Ha"]*DR1$RIah1
  RItmlela <- mean(RIfla1)-mean(RIfla0)
  RItmlea <- mean(RIfa1)-mean(RIfa0)
  # Propensity model incorrectly specified
  DR2 <- merge(D[,c(1:3)],RIexppm[,c(1,2,13:18)])
  DR2 <- merge(DR2,RIgcmpd)
  RItmlem1 <- glm(data=DR2,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0m1 <- DR2$mlfy00+RItmlem1$coefficients["Hla"]*DR2$RIlah0
  RIfla1m1 <- DR2$mlfy10+RItmlem1$coefficients["Hla"]*DR2$RIlah1
  RIfa0m1 <- DR2$mlfy00+RItmlem1$coefficients["Ha"]*DR2$RIah0
  RIfa1m1 <- DR2$mlfy01+RItmlem1$coefficients["Ha"]*DR2$RIah1
  RItmlelam1 <- mean(RIfla1m1)-mean(RIfla0m1)
  RItmleam1 <- mean(RIfa1m1)-mean(RIfa0m1)  
  # Outcome model incorrectly specified
  DR3 <- merge(D[,c(1:3)],RIexpp[,c(1,2,13:18)])
  DR3 <- merge(DR3,RIgcmpdm)
  RItmlem2 <- glm(data=DR3,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0m2 <- DR3$mlfy00+RItmlem2$coefficients["Hla"]*DR3$RIlah0
  RIfla1m2 <- DR3$mlfy10+RItmlem2$coefficients["Hla"]*DR3$RIlah1
  RIfa0m2 <- DR3$mlfy00+RItmlem2$coefficients["Ha"]*DR3$RIah0
  RIfa1m2 <- DR3$mlfy01+RItmlem2$coefficients["Ha"]*DR3$RIah1
  RItmlelam2 <- mean(RIfla1m2)-mean(RIfla0m2)
  RItmleam2 <- mean(RIfa1m2)-mean(RIfa0m2)  
  # Both incorrectly specified
  DR4 <- merge(D[,c(1:3)],RIexppm[,c(1,2,13:18)])
  DR4 <- merge(DR4,RIgcmpdm)
  RItmlem3 <- glm(data=DR4,y~-1+Hla+Ha+offset(mlfy),family="gaussian")
  RIfla0m3 <- DR4$mlfy00+RItmlem3$coefficients["Hla"]*DR4$RIlah0
  RIfla1m3 <- DR4$mlfy10+RItmlem3$coefficients["Hla"]*DR4$RIlah1
  RIfa0m3 <- DR4$mlfy00+RItmlem3$coefficients["Ha"]*DR4$RIah0
  RIfa1m3 <- DR4$mlfy01+RItmlem3$coefficients["Ha"]*DR4$RIah1
  RItmlelam3 <- mean(RIfla1m3)-mean(RIfla0m3)
  RItmleam3 <- mean(RIfa1m3)-mean(RIfa0m3)  
  
  # SuperLearner TMLE    
  DS <- merge(D[,c(1:3)],SLexpp[,c(1,2,13:18)])
  DS <- merge(DS,SLgcmpd)
  SLtmle <- glm(data=DS,y~-1+Hla+Ha+offset(slfy),family="gaussian")
  SLfla0 <- DS$slfy00+SLtmle$coefficients["Hla"]*DS$SLlah0
  SLfla1 <- DS$slfy10+SLtmle$coefficients["Hla"]*DS$SLlah1
  SLfa0 <- DS$slfy00+SLtmle$coefficients["Ha"]*DS$SLah0
  SLfa1 <- DS$slfy01+SLtmle$coefficients["Ha"]*DS$SLah1
  SLtmlela <- mean(SLfla1)-mean(SLfla0)
  SLtmlea <- mean(SLfa1)-mean(SLfa0)
  
  c(GLgcmpa,GLgcmpla,GLgcmpa+GLgcmpla,
    GLgcmpam,GLgcmplam,GLgcmpam+GLgcmplam,
    RIgcmpa,RIgcmpla,RIgcmpa+RIgcmpla,
    RIgcmpam,RIgcmplam,RIgcmpam+RIgcmplam,
    SLgcmpa,SLgcmpla,SLgcmpa+SLgcmpla,
    GLdrgcmpa,GLdrgcmpla,GLdrgcmpa+GLdrgcmpla,
    GLdrgcmpam1,GLdrgcmplam1,GLdrgcmpam1+GLdrgcmplam1,
    GLdrgcmpam2,GLdrgcmplam2,GLdrgcmpam2+GLdrgcmplam2,
    GLdrgcmpam3,GLdrgcmplam3,GLdrgcmpam3+GLdrgcmplam3,
    RIdrgcmpa,RIdrgcmpla,RIdrgcmpa+RIdrgcmpla,
    RIdrgcmpam1,RIdrgcmplam1,RIdrgcmpam1+RIdrgcmplam1,
    RIdrgcmpam2,RIdrgcmplam2,RIdrgcmpam2+RIdrgcmplam2,
    RIdrgcmpam3,RIdrgcmplam3,RIdrgcmpam3+RIdrgcmplam3,
    SLdrgcmpa,SLdrgcmpla,SLdrgcmpa+SLdrgcmpla,
    GLtmlea,GLtmlela,GLtmlea+GLtmlela,
    GLtmleam1,GLtmlelam1,GLtmleam1+GLtmlelam1,
    GLtmleam2,GLtmlelam2,GLtmleam2+GLtmlelam2,
    GLtmleam3,GLtmlelam3,GLtmleam3+GLtmlelam3,
    RItmlea,RItmlela,RItmlea+RItmlela,
    RItmleam1,RItmlelam1,RItmleam1+RItmlelam1,
    RItmleam2,RItmlelam2,RItmleam2+RItmlelam2,
    RItmleam3,RItmlelam3,RItmleam3+RItmlelam3,
    SLtmlea,SLtmlela,SLtmlea+SLtmlela)
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
  SL.bootdata <- foreach(1:500) %do% bootdata(d=data)
  SL.se <- as.data.frame(do.call(rbind,suppressWarnings(
    foreach(x=1:500, .errorhandling = 'remove') %do% comperror(data=SL.bootdata[[x]]))))
  
  results <- matrix(ncol = 258)
  
  results[1]<-coef(summary(GLM))[2,1]
  results[2]<-coef(summary(GLM))[2,2]
  results[3]<-coef(summary(GLM))[3,1]
  results[4]<-coef(summary(GLM))[3,2]
  results[5]<-coef(summary(GLM))[2,1]+coef(summary(GLM))[3,1]
  results[6]<-sqrt(vcov(GLM)[2,2] + vcov(GLM)[3,3] + 2*vcov(GLM)[2,3])
  
  results[7]<-coef(summary(GLMm))[2,1]
  results[8]<-coef(summary(GLMm))[2,2]
  results[9]<-coef(summary(GLMm))[3,1]
  results[10]<-coef(summary(GLMm))[3,2]
  results[11]<-coef(summary(GLMm))[2,1]+coef(summary(GLMm))[3,1]
  results[12]<-sqrt(vcov(GLMm)[2,2] + vcov(GLMm)[3,3] + 2*vcov(GLMm)[2,3])
  
  results[13]<-coef(summary(RI))[2,1]
  results[14]<-coef(summary(RI))[2,2]
  results[15]<-coef(summary(RI))[3,1]
  results[16]<-coef(summary(RI))[3,2]
  results[17]<-coef(summary(RI))[2,1]+coef(summary(RI))[3,1]
  results[18]<-sqrt(vcov(RI)[2,2] + vcov(RI)[3,3] + 2*vcov(RI)[2,3])
  
  results[19]<-coef(summary(RIm))[2,1]
  results[20]<-coef(summary(RIm))[2,2]
  results[21]<-coef(summary(RIm))[3,1]
  results[22]<-coef(summary(RIm))[3,2]
  results[23]<-coef(summary(RIm))[2,1]+coef(summary(RIm))[3,1]
  results[24]<-sqrt(vcov(RIm)[2,2] + vcov(RIm)[3,3] + 2*vcov(RIm)[2,3])
  
  results[25]<-coef(summary(GEE))[2,1]
  results[26]<-coef(summary(GEE))[2,2]
  results[27]<-coef(summary(GEE))[3,1]
  results[28]<-coef(summary(GEE))[3,2]
  results[29]<-coef(summary(GEE))[2,1]+coef(summary(GEE))[3,1]
  results[30]<-sqrt(summary(GEE)$cov.scaled[2,2] + summary(GEE)$cov.scaled[3,3] + 2*summary(GEE)$cov.scaled[2,3])
  
  results[31]<-coef(summary(GEEm))[2,1]
  results[32]<-coef(summary(GEEm))[2,2]
  results[33]<-coef(summary(GEEm))[3,1]
  results[34]<-coef(summary(GEEm))[3,2]
  results[35]<-coef(summary(GEEm))[2,1]+coef(summary(GEEm))[3,1]
  results[36]<-sqrt(summary(GEEm)$cov.scaled[2,2] + summary(GEEm)$cov.scaled[3,3] + 2*summary(GEEm)$cov.scaled[2,3])
  
  for (k in 1:42) {
    l<-k+111
    m<-(k*2)+35
    n<-m+1
    results[m]<-SL.est[k]
    results[n]<-SL.est[l]
  }
  
  for (k in 43:111) {
    l<-(k*2)+35
    m<-l+1
    n<-k-42
    results[l]<-SL.est[k]
    results[m]<-sd(SL.se[,n])
  }
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
               "GLMgcmpaco","GLMgcmpase","GLMgcmplaco","GLMgcmplase","GLMgcmpjtco","GLMgcmpjtse",
               "RIgcmpaco","RIgcmpase","RIgcmplaco","RIgcmplase","RIgcmpjtco","RIgcmpjtse",
               "SLgcmpaco","SLgcmpase","SLgcmplaco","SLgcmplase","SLgcmpjtco","SLgcmpjtse",
               "GLMdrgcmpaco","GLMdrgcmpase","GLMdrgcmplaco","GLMdrgcmplase","GLMdrgcmpjtco","GLMdrgcmpjtse",
               "RIdrgcmpaco","RIdrgcmpase","RIdrgcmplaco","RIdrgcmplase","RIdrgcmpjtco","RIdrgcmpjtse",
               "SLdrgcmpaco","SLdrgcmpase","SLdrgcmplaco","SLdrgcmplase","SLdrgcmpjtco","SLdrgcmpjtse",
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
parallel::clusterExport(cl, c("bootdata","checkrange","compestimates","comperror","sim"))
registerDoParallel(cl)

set.seed(269012,kind="L'Ecuyer-CMRG")
start1 <- Sys.time()
results1 <- data.frame(matrix(unlist(parLapply(cl,data[1:250], function (x) {sim(x)})),nrow=250,ncol=258,byrow=TRUE))
end1 <- Sys.time()
end1-start1
start1 <- Sys.time()
results2 <- data.frame(matrix(unlist(parLapply(cl,data[251:500], function (x) {sim(x)})),nrow=250,ncol=258,byrow=TRUE))
end2 <- Sys.time()
end2-start2

resultsA<-rbind(results1,results2)
save(resultsA,file="resultsA.RData")
resultsA1<-results[,c(1:6,13:18,25:30,37:42,49:54,61:66,67:72,91:96,115:120,121:126,133:138,145:150,151:156,175,180,199:204,205:210,229:234,253:258)]
colnames(resultsA1) <- col.names
save(resultsA1,file="resultsA1.RData")
resultsA2<-results[,c(1:6,13:18,25:30,43:48,55:60,61:66,73:78,97:102,115:120,121:126,133:138,145:150,157:162,181:186,199:204,211:216,235:240,253:258)]
colnames(resultsA2) <- col.names
save(resultsA2,file="resultsA2.RData")
resultsA3<-results[,c(7:12,19:24,31:36,37:42,49:54,61:66,79:84,103:108,115:120,127:132,139:144,145:150,163:168,187:192,199:204,217:222,241:246,253:258)]
colnames(resultsA3) <- col.names
save(resultsA3,file="resultsA3.RData")
resultsA4<-results[,c(7:12,19:24,31:36,43:48,55:60,61:66,85:90,109:114,115:120,127:132,139:144,145:150,169:174,193:198,199:204,223:228,247:252,253:258)]
colnames(resultsA4) <- col.names
save(resultsA4,file="resultsA4.RData")



