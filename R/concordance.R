concordance <-
function(db,test="Default",B=1000,alpha=0.05) {
  #Default è una funzione che non mi restituisce il pvalue#
  db2<-db*(db-1)
  sumcol<-colSums(db)
  sumrow<-rowSums(db)
  tot<-sum(db)
  vet<-list()
  
  #Calculate p mean
  
  pij<-matrix(,nrow=nrow(db),ncol=ncol(db))
  for (i in 1:length(sumrow)) {
    for (j in 1:length(sumcol)) {
      pij[i,j]<-db2[i,j]/(sumrow[i]*(sumrow[i]-1))
    }
  }
  pi<-rowSums(pij)
  p<-mean(pi) #p medio
  
  #Calculate pe 
  
  pj<-sumcol/tot
  pj2<-pj^2
  pe<-sum(pj2)
  fleiss.kappa<-(p-pe)/(1-pe) 
  s<-(ncol(db)*p-1)/(ncol(db)-1) 
  
  Default<-function(vet) {vet<-list(Fleiss=fleiss.kappa,S=s)
                          cat(paste("Inter-rater Agreement"),"\n")
                          cat(paste("Fleiss =",round(vet$Fleiss,3),"\n"))
                          cat(paste("     S =",round(vet$S,3),"\n"))
                          }
  
  
  Normal<-function(db) {
    stat.test<-s*nrow(db)*sqrt((ncol(db)-1)/  #modifica#
              (2*sum(1/(rowSums(db)*(rowSums(db)-1)))))
    bin<-stat.test>qnorm(1-alpha)
    pvalue<-1-pnorm(stat.test)
    vet<-list(Fleiss=fleiss.kappa,S=s,pvalue=pvalue)
    
    cat(paste("Inter-rater Agreement"),"\n")
    cat(paste("Fleiss =",round(vet$Fleiss,3),"\n"))
    cat(paste("     S =",round(vet$S,3),"\n"))
    cat(paste("pvalue =",round(vet$pvalue,5),"\n"))
  
    
  }
  
  MC<-function(db) {
    matrix.mc<-list()
    pij.mc<-list()
      for (l in 1:B) {
        
          matrix.mc[[l]]<-matrix(,nrow=nrow(db),ncol=ncol(db),byrow=T)
          pij.mc[[l]]<-matrix(,nrow=nrow(db),ncol=ncol(db),byrow=T)
      }
    
      for (k in 1:B) {
        for (j in 1:nrow(db))
          
        matrix.mc[[k]][j,]<-t(rmultinom(1,size=rowSums(db)[j],prob=rep(1/(ncol(db)),ncol(db))))
    }

      for(k in 1:B) {
        for (i in 1:nrow(db)) {
          for (j in 1:ncol(db)) {
      pij.mc[[k]][i,j]<-(matrix.mc[[k]][i,j])*((matrix.mc[[k]][i,j])-1)/(rowSums(db)[i]*(rowSums(db)[i]-1))
    }
  }
}   

    pi.mc<-list()
        for (k in 1:B) {
          pi.mc[[k]]<-matrix(,nrow=1,ncol=nrow(db))
}  

    for (k in 1:B) {
      for ( j in 1:nrow(db)) {
          pi.mc[[k]][,j]<-sum(pij.mc[[k]][j,])
  }
}

p.mc<-list()

    for (k in 1:B) {
      p.mc[[k]]<-mean(pi.mc[[k]]) 
}

s.mc<-c()

    for (k in 1:B) {
      s.mc[k]<-((p.mc[[k]]*ncol(db))-1)/(ncol(db)-1)
}

  crit.s.mc<-quantile(s.mc,0.95) #percentile di tutti le statistiche
  binary<-c()
    for (i in 1:length(s.mc)) {
      binary[i]<-(s.mc[i]>=s) #rivedere bene #modificato
}
pvalue<-sum(binary)/B

      vet<-list(Fleiss=fleiss.kappa,S=s,pvalue=pvalue)

      cat(paste("Inter-rater Agreement"),"\n")
      cat(paste("Fleiss =",round(vet$Fleiss,3),"\n"))
      cat(paste("     S =",round(vet$S,3),"\n"))
      cat(paste("pvalue =",round(vet$pvalue,5),"\n"))


} 
   

     
  
  
  #confidence intervals for fleiss if number of judges are the same
  
  if (length(unique(rowSums(db)))==1) {
    std.kappa<-sqrt(2/(tot*rowSums(db)[1]))
    ci.low<-fleiss.kappa-qnorm(0.975)*std.kappa
    ci.upp<-fleiss.kappa+qnorm(0.975)*std.kappa
      if (ci.upp>1 | ci.low<(-1/(rowSums(db)[1]-1))) {
        ci.low<-(1/(rowSums(db)[1]-1))
        ci.upp<-1
      }
      
    #result<-c(ci.low,fleiss.kappa,ci.upp)#
    #names(result)<-c("lower","fleiss","upper")#
    
    Chisq<-function(db) {
      
      stat.test<-nrow(db)*(ncol(db)-1)*((rowSums(db)-1)[1]*s+1)
      bin<-stat.test>qchisq(alpha,df=(nrow(db)*(ncol(db)-1)))
      
      pvalue<-1-pchisq(stat.test,df=(nrow(db)*(ncol(db)-1)))
      
      vet<-list(Lower=ci.low,Fleiss=fleiss.kappa,Upper=ci.upp,S=s,pvalue=pvalue)
      
      
      cat(paste("Inter-rater Agreement"),"\n")
      cat(paste(" Lower =",round(vet$Lower,3),"\n"))
      cat(paste("Fleiss =",round(vet$Fleiss,3),"\n"))
      cat(paste(" Upper =",round(vet$Upper,3),"\n"))
      cat(paste("     S =",round(vet$S,3),"\n"))
      cat(paste("pvalue =",round(vet$pvalue,5),"\n"))
      
      
      
    }
    
    
    Normal2<-function(db) {
      stat.test<-s*nrow(db)*sqrt((ncol(db)-1)/  #modifica#
                                   (2*sum(1/(rowSums(db)*(rowSums(db)-1)))))
      bin<-stat.test>qnorm(1-alpha)
      pvalue<-1-pnorm(stat.test)
    
      vet<-list(Lower=ci.low,Fleiss=fleiss.kappa,Upper=ci.upp,S=s,pvalue=pvalue)
      
      cat(paste("Inter-rater Agreement"),"\n")
      cat(paste(" Lower =",round(vet$Lower,3),"\n"))
      cat(paste("Fleiss =",round(vet$Fleiss,3),"\n"))
      cat(paste(" Upper =",round(vet$Upper,3),"\n"))
      cat(paste("     S =",round(vet$S,3),"\n"))
      cat(paste("pvalue =",round(vet$pvalue,5),"\n"))
      
      
    }
    
    
    
    
    
    MC2<-function(db) {
      matrix.mc<-list()
      pij.mc<-list()
      for (k in 1:B) {
        pij.mc[[k]]<-matrix(,nrow=nrow(db),ncol=ncol(db))
        matrix.mc[[k]]<-t(rmultinom(nrow(db),size=rowSums(db)[1],prob=rep(1/(ncol(db)),ncol(db))))
      }
      
      for(k in 1:B) {
        for (i in 1:nrow(db)) {
          for (j in 1:ncol(db)) {
            pij.mc[[k]][i,j]<-(matrix.mc[[k]][i,j])*((matrix.mc[[k]][i,j])-1)/(rowSums(db)[i]*(rowSums(db)[i]-1))
          }
        }
      }   
      
pi.mc<-list()
      for (k in 1:B) {
        pi.mc[[k]]<-matrix(,nrow=1,ncol=nrow(db))
      }  
      
      for (k in 1:B) {
        for ( j in 1:nrow(db)) {
          pi.mc[[k]][,j]<-sum(pij.mc[[k]][j,])
        }
      }
      
p.mc<-list()
      
      for (k in 1:B) {
        p.mc[[k]]<-mean(pi.mc[[k]]) 
      }
      
s.mc<-c()
      
      for (k in 1:B) {
        s.mc[k]<-((p.mc[[k]]*ncol(db))-1)/(ncol(db)-1)
      }
      
crit.s.mc<-quantile(s.mc,0.95) #percentile di tutti le statistiche
binary<-c()
      for (i in 1:length(s.mc)) {
        binary[i]<-(s.mc[i]>=s) #rivedere bene #modificato
      }
pvalue<-sum(binary)/B

    vet<-list(Lower=ci.low,Fleiss=fleiss.kappa,Upper=ci.upp,S=s,pvalue=pvalue)
    
    cat(paste("Inter-rater Agreement"),"\n")
    cat(paste(" Lower =",round(vet$Lower,3),"\n"))
    cat(paste("Fleiss =",round(vet$Fleiss,3),"\n"))
    cat(paste(" Upper =",round(vet$Upper,3),"\n"))
    cat(paste("     S =",round(vet$S,3),"\n"))
    cat(paste("pvalue =",round(vet$pvalue,5),"\n"))

  } 
    
    

    
  
    switch(test,
           Normal=Normal2(db),
           MC=MC2(db),
           Chisq=Chisq(db),
           Default=Default(vet)) 
}
    
 else { switch(test,
               Normal=Normal(db),
               MC=MC(db),
               Default=Default(vet)) }
  
  
}
