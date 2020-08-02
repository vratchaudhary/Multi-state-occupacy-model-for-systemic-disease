#############################################################################
#############################################################################
#The R code is provided as supplementary information for Chaudhary et al. 2020 (Journal of Applied Ecology, JAPPL-2020-00088). 
# A multi-state occupancy modeling framework for robust estimation of disease prevalence in multi-tissue disease systems
#This code reads the PrV detection history and associated covariate data (PRV_detection_data.csv; see S2.1 for a detailed description of the data), 
#and fits the general and constrained occupancy model described in this paper
############################################################################



set.seed(54321); epsln=1e-5; ofile='multistate_occ_modlsf.out';
NREPS=2

# Initial values
K=3 # K = Number of surveys
psi=c(.5,.4,.3,.6) # psi = occupancy probability /Prevalence
p=c(.2,.8,.3,.7) # p = detection history

#For our example of four tissues
psi=psi[1:4]; p=p[1:length(psi)]

S=length(psi)        #  number of sub-levels of occupancy (4)
SS=2^S               #  number of possible states, given sub-levels, including 0000
SSm1=SS-1            #  number of possible states, given sub-levels excluding 0000

#creating combination of truestates

trustate=NULL
for (i in 0:SSm1) trustate=rbind(trustate,as.integer(intToBits(i)[S:1]))
states=OBNG=apply(trustate,1,paste,collapse='') #OBNG represent oral, blood, nasal, and genital

#read data.
a= read.csv("PRV_detection_data.csv"); N=nrow(a)

head(a,10)
tail(a,10)

#covariates

sampcols= a[,6:17] #sample columns
covcols = a[,2:5]  #covariate columns

sampcols=c(6:17)
covcols=c(2:5)

##Creating detection history.
#survey 1
h1=cbind(a$sampling1.oral1,a$sampling1.blood1,a$sampling1.nasal1,a$sampling1.gen1)
#surevy 2
h2=cbind(a$sampling2.oral1,a$sampling2.blood1,a$sampling2.nasal1,a$sampling2.gen1)
#survey 3
h3=cbind(a$sampling3.oral1,a$sampling3.blood1,a$sampling3.nasal1,a$sampling3.gen1)

# combine samples...
# missing samples are the ones that were not collected or were not diagnosed
# and are designated as 2 here
h4=h40=cbind(h1,h2,h3)
h40[is.na(h4)]=0
h4[is.na(h4)]=2  #  h4=det.hist as matrix,

dhst=h4; dfrq=rep(1,nrow(dhst))
obsstate=array(as.integer(dhst),dim=c(nrow(dhst),S,K)); #  individual covariates
obsstate0=obsstate; obsstate0[obsstate==2]=0  #  obsstate=matrix of det. histories,
#obsstate0=same w/ missing replaced by 0's


#computes multinomial logit of A ( exp(A)/(1+sum(exp(A))))
mlogit<-function(a) {
  x=exp(a)/(1+sum(exp(a)))
  return(x)
}

#####Constrained model

#  computes the negative log-likelihood of multi-state occupancy model with independent sub-states
#parameters
# theta = psi and p
#desmat 1= design matrix for psi
#desmat 2 = design matrix for ps
like1<- function(theta,desmat,desmat2) {
  NN=nrow(obsstate0)
  loglike=0
  npar1=ncol(desmat) # psi
  npar2=ncol(desmat2) # p

  psi=p=rep(1,S) #initial values
  for (i in 1:NN) {  #   for each detection history/ animal, 510 in this case

    for (j in 1:S) {  # for each tissue

      psi[j] = plogis(sum(desmat[i+(j-1)*NN,] * theta[1:npar1])) # psi each tissue

      p[j] = plogis(sum(desmat2[i+(j-1)*NN,] * theta[npar1+1:npar2])) #p for each tissue
    }

    #calculating multi-state psi : 16 states = SS
    ms_psi=rep(1,SS); for (tsi in 2:SS) for (j in 1:S) ms_psi[tsi]= ms_psi[tsi]*c(1-psi[j],psi[j])[trustate[tsi,j]+1]

    #calculate ms psi and each of the 16 combinattion can be derived using the detection history

    ms_psi[1] = 1-sum(ms_psi[-1]); #prob of not being occuped 0000
    qp = cbind(1-p,p,1) #probability of detection

    # compute product of psi(tru-state) * p1(obs-state|tru-state)*p2(obs-state|tru-state)...

    prd=ms_psi
    for (srvy in 1:K) {    #for each survey, 3 survey all together
      for (tsi in 1:SS) { #  for each tru-state...   (tsi=trustate index -> 1=0000, 2=0001, 3=0010, 4=0011,... 16=1111)

# if obsstate is possible for given trustate-
        #i. don't compute product of p's if tru-state='0000'
        #ii for each of the 4 sample types (O,B,N,G) if tru-state for sample-type tiss == 1, then multiply
        #observed state for the sample-type.
        #So, the row of matrix qp is for the sample-type, and the column indicates the quantity to multiply:  (1-p) if not detected, p if detected, and 1 if missing data.

        #if  otherwise (obs-state not possible for tru-state e.g.,for true state 1011, obs state 1111 is not possible),set product for that tru-state to zero

      if (prod((trustate[tsi,]*obsstate0[i,,srvy])==obsstate0[i,,srvy])>0) {
          if (tsi>1)
            for (tiss in 1:S)

              if (trustate[tsi,tiss]>0) prd[tsi]=prd[tsi]*qp[tiss,obsstate[i,tiss,srvy]+1]
        } else
          prd[tsi]=0;
      }
    }
    loglike=loglike+dfrq[i]*log(sum(prd))  #  sum likelihood for det.history i
  }

  return(-loglike)
}




like2<- function(theta,desmat,desmat2) {  #  computes the negative log-likelihood of multi-state occupancy model with interaction among sub-states
  # ----------------------------------------------------------------
  ##they are dependent on 16 state
  # ----------------------------------------------------------------

  NN=nrow(obsstate0); loglike=0; npar1=ncol(desmat); npar2=ncol(desmat2)

  psi=rep(1,SSm1); p=rep(1,S)

  for (i in 1:NN) {              #   for each detection history...

    for (j in 1:SSm1) psi[j]=exp(sum(desmat[i+(j-1)*NN,] * theta[1:npar1]))
    for (j in 1:S)  p[j]=plogis(sum(desmat2[i+(j-1)*NN,] * theta[npar1+1:npar2]))

    ms_psi=psi/(1+sum(psi));
    ms_psi=c(1-sum(ms_psi),ms_psi);

    qp=cbind(1-p,p,1)

    prd=ms_psi                   #       compute product of psi(tru-state) * p1(obs-state|tru-state)*p2(obs-state|tru-state)...
    for (srvy in 1:K) {          #       for each survey...
      for (tsi in 1:SS) {        #          for each tru-state...   (tsi=trustate index -> 1=0000, 2=0001, 3=0010, 4=0011,... 16=1111)
        if (prod((trustate[tsi,]*obsstate0[i,,srvy])==obsstate0[i,,srvy])>0) {  #  if obsstate is possible for given trustate...
          if (tsi>1)            #              ... don't compute product of p's if tru-state='0000'
            for (tiss in 1:S)     #              for each of the 4 sample types (O,B,N,G)...
              #                  if tru-state for sample-tpe tiss == 1, then multiply product for
              #                  tru-state by either 1-p, p, or 1, depending on the
              #                  observed state for the sample-type.  So, the row
              #                  of matrix qp is for the sample-type, and the column
              #                  indicates the quantity to multiply: (1-p) if not detected,
              #                  p if detected, and 1 if missing data.

              if (trustate[tsi,tiss]>0) prd[tsi]=prd[tsi]*qp[tiss,obsstate[i,tiss,srvy]+1]
        } else                  #             otherwise (obs-state not possible for tru-state),
          prd[tsi]=0;        #                  set product for that tru-state to zero.
      }
    }
    loglike=loglike+dfrq[i]*log(sum(prd))  #  sum likelihood for det.history i
  }
  #cat(loglike,round(ms_psi,4),'  \n');
  return(-loglike)
}




print_psi_indep <- function(modl) {
  theta=modl$estimate; dm1=unique(modl$desmat); npar1=ncol(dm1);
  dm2=unique(modl$desmat2); npar2=ncol(dm2); npar=npar1+npar2
  psi=plogis(dm1 %*% theta[1:npar1]); p=plogis(dm2 %*% theta[-1:-npar1])

  # compute variances via delta method...
  grad=gradp=NULL;
  for (k in 1:npar) {
    theta[k]=theta[k]+epsln;                                 # increment k'th theta parameter
    psi2=plogis(dm1 %*% theta[1:npar1])
    p2=plogis(dm2 %*% theta[-1:-npar1])
    grad=cbind(grad,c(psi2-psi,p2-p)/epsln)                        # then, compute gradient as diff. in real parms / diff in theta's
    theta[k]=theta[k]-epsln                                  # reset theta back to orig. value
  }
  varcov.real=grad %*% modl$vcbeta %*% t(grad);   # delta-method: var-cov(real) = grad * var-cov(beta) * grad'
  se=sqrt(diag(varcov.real))
  x=cbind(psi,se[1:length(psi)]); colnames(x)=c('estimate','std.err');
  y=cbind(p,se[-1:-length(psi)]); colnames(y)=c('estimate','std.err');
  return(list(psi=x,p=y))
}


#  load snowfall package so we can run many replications of rand. init. value vectors at once...
#library(snowfall); sfInit(parallel=T,cpus=2); sfExportAll();




indep_mod <- function(psiform,pform) {                      #  independent tissues model (parms:psi1,psi2,psi3,psi4, p1,p2,p3,p4)
  mtx=a[,covcols]; mtx=data.frame(mtx[rep(row.names(mtx),S),])       #  create matrix of data so design matrix can be created using formula
  mtx$tissue=as.factor(rep(1:S,each=N));
  desmat=model.matrix(psiform,mtx); npar1=ncol(desmat);                    #  create design matrix for psi using psiform (~tissue+Sex)
  desmat2=model.matrix(pform,mtx); npar2=ncol(desmat2); npar=npar1+npar2   #  create design matrix for p using pform (~tissue)

  z= lapply(1:NREPS, function(ii) {                                       #  find min. neg. log-likelihood for a bunch of random initial

    ##################call like 2

    init_theta=theta=runif(npar)*6-3#; dyn.load('like2') ##for mac do so  instead of dll            #    value vectors.
    z<-nlm(like1,init_theta,desmat,desmat2,iterlim=500,steptol=1e-5,hessian=T)
    z$aic=z$minimum*2+2*npar; #z$tm=tm#; dyn.unload('like2')
    return(z)
  })
  jj=1; nlls=NULL                                                  #  find minimum neg. log-likelihood
  for (ii in 1:NREPS) {
    if (z[[ii]]$minimum<z[[jj]]$minimum) jj=ii
    nlls=c(nlls,z[[ii]]$minimum)                                  #     save neg. log-likelihood values in vector
  }
  zz=z[[jj]]; theta=zz$estimate
  zz$results=list(); izz=1                                         #  save model name and estimates in results list...
  zz$results[[1]]=paste0('===== independent tissues psi model === psi(',as.character(psiform)[-1],') p(',as.character(pform)[-1],') AIC=',zz$aic,collapse='')
  psi=p=rep(1,S); NN=nrow(obsstate0); psipmat=NULL
  for (i in 1:2) {              #   for each detection history...(only doing 1st two sites to save paper)
    for (j in 1:S) {
      psi[j]=plogis(sum(desmat[i+(j-1)*NN,] * theta[1:npar1]))      #  get psi's for each site
      p[j]=plogis(sum(desmat2[i+(j-1)*NN,] * theta[npar1+1:npar2])) #  get p's for each site
    }
    psipmat=cbind(psipmat,psi,p);           #  save psi's and p's in a matrix
  }

  ##################################33
  #interaction model 16 possible psis, 16 different beta parameters


  # compute variances via delta method
  i=which(rowSums(zz$hessian)!=0); hess=zz$hessian[i,i];           #  eliminate rows/cols w/ all zeros from hessian
  varcov.beta=matrix(0,npar,npar);
  try(varcov.beta[i,i]<-solve(hess))    #   then invert hessian to get var-cov matrix of beta's
  grad=NULL;
  for (k in 1:npar) {
    theta[k]=theta[k]+epsln; psipmat2=NULL                           # increment k'th theta parameter
    for (i in 1:2) {                                                 # then, compute each real parameter (only doing 1st sites to save paper)
      for (j in 1:S) {
        psi[j]=plogis(sum(desmat[i+(j-1)*NN,] * theta[1:npar1]))      #  get psi's for each site
        p[j]=plogis(sum(desmat2[i+(j-1)*NN,] * theta[npar1+1:npar2])) #  get p's for each site
      }
      psipmat2=cbind(psipmat2,psi,p);           #  save psi's and p's in a matrix
    }
    grad=cbind(grad,c(psipmat2-psipmat)/epsln)                        # then, compute gradient as diff. in real parms / diff in theta's
    theta[k]=theta[k]-epsln                                          # reset theta back to orig. value
  }
  varcov.real=grad %*% varcov.beta %*% t(grad);   # delta-method: var-cov(real) = grad * var-cov(beta) * grad'
  se=matrix(sqrt(diag(varcov.real)),nrow=S)
  psipmat2=cbind(psipmat[,1],se[,1],psipmat[,3],se[,3],psipmat[,2],se[,2],psipmat[,4],se[,4])
  rownames(psipmat2)=c('O','B','N','G');
  lbls=paste0(c(rep(c('psi','se'),2),rep(c('p','se'),2)),rep(1:2,each=2)); colnames(psipmat2)=lbls
  izz=izz+1; zz$results[[izz]]=psipmat2
  izz=izz+1; k=order(nlls); zz$results[[izz]]=nlls[k]
  names(zz$results)=c('modname','psi_p','negLL')
  rownames(desmat)=paste0('psi',rep(1:S,each=NN)); rownames(desmat2)=paste0('p',rep(1:S,each=NN))
  zz$desmat=desmat; zz$desmat2=desmat2; zz$vcbeta=varcov.beta
  return(zz)
}




########################################################################################

interaction_mod <- function(psiform,pform) {  #  interaction among tissues model (parms:psi0001,psi0010,psi0011,psi0100,..., p1,p2,p3,p4)
  mtx=a[,covcols]; mtx=data.frame(mtx[rep(row.names(mtx),SSm1),])       #  create matrix of data so design matrix can be created using formula
  mtx$state=as.factor(rep(1:SSm1,each=N));

  desmat<-model.matrix(psiform,mtx); npar1<-ncol(desmat);                    #  create design matrix for psi using psiform (~state+Sex)
  mtx2=a[,covcols]; mtx2=data.frame(mtx2[rep(row.names(mtx2),S),])
  mtx2$tissue=as.factor(rep(1:S,each=N));
  #   #  need a different data matrix, mtx2, since there are more psi parameters than p parameters
  desmat2=model.matrix(pform,mtx2); npar2=ncol(desmat2); npar=npar1+npar2    #  create design matrix for p using pform (~tissue)
  z=lapply(1:NREPS, function(ii) {                                         #  find min. neg. log-likelihood for a bunch of random initial
    init_theta=theta=runif(npar)*6-3 #; dyn.load('like3.so')                #    value vectors.
    emsg<-try(z<-nlm(like2,init_theta,desmat,desmat2,iterlim=500,steptol=1e-5,hessian=T))
    if (length(emsg)==6) {z$aic=z$minimum*2+2*npar;  } else z=list(aic=999999,minimum=999999)
   # dyn.unload('like2')
    return(z)
  })
  jj=1; nlls=NULL                                  #  find minimum neg. log-likelihood
  for (ii in 1:NREPS) {
    if (z[[ii]]$minimum<z[[jj]]$minimum) jj=ii
    nlls=c(nlls,z[[ii]]$minimum)                  #     save neg. log-likelihood values in vector
  }
  zz=z[[jj]];
  theta=zz$estimate
  zz$results=list(); izz=1                               #  save model name and estimates in results list...
  zz$results[[1]]=paste0('===== interaction psi model === psi(',as.character(psiform)[-1],') p(',as.character(pform)[-1],') AIC=',zz$aic,collapse='')
  psi=rep(1,SSm1); NN=nrow(obsstate0); psi1=p1=NULL
  for (i in 1:2) {              #   for each detection history...
    for (j in 1:S) p[j]=plogis(sum(desmat2[i+(j-1)*NN,] * theta[npar1+1:npar2])) # compute p's from model theta's
    for (j in 1:SSm1) psi[j]=exp(sum(desmat[i+(j-1)*NN,] * theta[1:npar1]))      # compute psi's from model theta's
    ms_psi=psi/(1+sum(psi)); ms_psi=c(1-sum(ms_psi),ms_psi);                     #   using mlogit transformation: exp(x)/(1+sum(exp(x)))
    psi1=c(psi1,ms_psi); p1=c(p1,p)
  }
  # compute variances via delta method
  i=which(rowSums(zz$hessian)!=0); hess=zz$hessian[i,i];           #  eliminate rows/cols w/ all zeros from hessian
  varcov.beta=matrix(0,npar,npar);
  k=try(varcov.beta[i,i]<-solve(hess))    #   then invert hessian to get var-cov matrix of beta's
  grad=NULL; psi2=p2=NULL
  for (k in 1:npar) {
    theta[k]=theta[k]+epsln; psi2=p2=NULL                            # increment k'th theta parameter
    for (i in 1:2) {                                                 # then, compute each real parameter (only doing 1st sites to save paper)
      for (j in 1:S) p[j]=plogis(sum(desmat2[i+(j-1)*NN,] * theta[npar1+1:npar2]))
      for (j in 1:SSm1) psi[j]=exp(sum(desmat[i+(j-1)*NN,] * theta[1:npar1]))
      ms_psi=psi/(1+sum(psi)); ms_psi=c(1-sum(ms_psi),ms_psi);
      psi2=c(psi2,ms_psi); p2=c(p2,p)
    }
    grad=cbind(grad,c(psi2-psi1,p2-p1)/epsln)                        # then, compute gradient as diff. in real parms / diff in theta's
    theta[k]=theta[k]-epsln                                          # reset theta back to orig. value
  }
  varcov.real=grad %*% varcov.beta %*% t(grad); se=round(sqrt(diag(varcov.real)),4)  # delta-method: var-cov(real) = grad * var-cov(beta) * grad'
  psi1=matrix(round(psi1,4),ncol=2)                                   # for printing psi, arrange as a matrix
  p1=matrix(round(p1,4),ncol=2)
  ii=1:(2*SS); psi.se=matrix(se[ii],ncol=2); p.se=matrix(se[-ii],ncol=2) #  compute se of psi and p, then save as a matrix
  psimat=cbind(OBNG,psi1[,1],psi.se[,1],psi1[,2],psi.se[,2]); colnames(psimat)=c('OBNG','site1','se1','site2','se2')
  izz=izz+1; zz$results[[izz]]=psimat
  pmat=cbind(c('O','B','N','G'),p1[,1],p.se[,1],p1[,2],p.se[,2]); colnames(pmat)=c('tissue','p1','se1','p2','se2')

  #   compute summed psi's and variances...
  sumpsi3=NULL
  for (k in 1:2) {                                                    #  for site's 1 and 2 (saving paper)
    sumpsi=rep(0,S); for (i in 1:S) sumpsi[i]=sum(psi1[trustate[,i]>0,k]);   #  calc sum of psi's for each tissue
    grad=matrix(0,S,SS)                                                      #
    for (i in 1:SS) {                                                        #  for each tru-state, 0001, 0010, 0011, 0100,...
      psi1[i,k]=psi1[i,k]+epsln                                             #    increment real parameter, psi
      sumpsi2=rep(0,S); for (j in 1:S) sumpsi2[j]=sum(psi1[trustate[,j]>0,k]);#  calc sum of psi's for each tissue using modified psi's
      grad[,i]=c(sumpsi2-sumpsi)/epsln                                      #    calc gradient
      psi1[i,k]=psi1[i,k]-epsln                                             #    reset psi's
    }
    varcov.sum=grad %*% varcov.real[1:SS,1:SS] %*% t(grad);         # delta method: var-cov(sum) = grad * var-cov(psi) * grad'
    se=sqrt(diag(varcov.sum))                                       #  std.err (se) = sqrt(variances)
    sumpsi3=cbind(sumpsi3,sumpsi,se)                                #  save sum of psi's and std.errors in matrix
  }
  rownames(sumpsi3)=c('O','B','N','G'); colnames(sumpsi3)=c('site1','se1','site2','se2')  # save sum of psi's in results
  izz=izz+1; zz$results[[izz]]=round(sumpsi3,4)
  izz=izz+1; zz$results[[izz]]=pmat
  izz=izz+1; k=order(nlls); zz$results[[izz]]=nlls[k]
  names(zz$results)=c('modname','ms_psi','sumpsi','p','negLL')
  return(zz)
}
m=list(); i=0 # m is all the models
i=i+1; m[[i]]=indep_mod(~tissue,~tissue); print(m[[i]]$results,quote=F)

i=i+1; m[[i]]=indep_mod(~tissue+Sex,~tissue); print(m[[i]]$results,quote=F)
i=i+1; m[[i]]=indep_mod(~tissue*Sex,~tissue); print(m[[i]]$results,quote=F)

i=i+1; m[[i]]=interaction_mod(~state,~tissue); print(m[[i]]$results,quote=F);
i=i+1; m[[i]]=interaction_mod(~state+Sex,~tissue); print(m[[i]]$results,quote=F);
i=i+1; m[[i]]=interaction_mod(~state*Sex,~tissue); print(m[[i]]$results,quote=F);

i=i+1; m[[i]]=indep_mod(~tissue+Age,~tissue); print(m[[i]]$results,quote=F)
i=i+1; m[[i]]=indep_mod(~tissue*Age,~tissue); print(m[[i]]$results,quote=F)

i=i+1; m[[i]]=interaction_mod(~state+Age,~tissue); print(m[[i]]$results,quote=F);
i=i+1; m[[i]]=interaction_mod(~state*Age,~tissue); print(m[[i]]$results,quote=F);

i=i+1; m[[i]]=indep_mod(~tissue+yr,~tissue); print(m[[i]]$results,quote=F)
i=i+1; m[[i]]=indep_mod(~tissue*yr,~tissue); print(m[[i]]$results,quote=F)

i=i+1; m[[i]]=interaction_mod(~state+yr,~tissue); print(m[[i]]$results,quote=F);
i=i+1; m[[i]]=interaction_mod(~state*yr,~tissue); print(m[[i]]$results,quote=F);

i=i+1; m[[i]]=indep_mod(~tissue+pdsi,~tissue); print(m[[i]]$results,quote=F)
i=i+1; m[[i]]=indep_mod(~tissue*pdsi,~tissue); print(m[[i]]$results,quote=F)

i=i+1; m[[i]]=interaction_mod(~state+pdsi,~tissue); print(m[[i]]$results,quote=F);
i=i+1; m[[i]]=interaction_mod(~state*pdsi,~tissue); print(m[[i]]$results,quote=F);

i=i+1; m[[i]]=indep_mod(~Sex,~tissue); print(m[[i]]$results,quote=F)
i=i+1; m[[i]]=indep_mod(~Age,~tissue); print(m[[i]]$results,quote=F)
i=i+1; m[[i]]=indep_mod(~yr,~tissue); print(m[[i]]$results,quote=F)
i=i+1; m[[i]]=indep_mod(~pdsi,~tissue); print(m[[i]]$results,quote=F)
i=i+1; m[[i]]=indep_mod(~1,~tissue); print(m[[i]]$results,quote=F)

nmods=i; aictbl=mnames=NULL
for (i in 1:nmods) { aictbl=c(aictbl,m[[i]]$aic); mnames=c(mnames,gsub(' AIC.+','',m[[i]]$results$modname))}
o=order(aictbl); minaic=min(aictbl); aictbl=cbind(aictbl,aictbl-minaic);
rownames(aictbl)=mnames; colnames(aictbl)=c('AIC','deltaAIC'); print(aictbl[o,])

x=print_psi_indep(m[[o[1]]]); cat('\nEstimates from top model:\n'); print(x)
