
setup <- function(k){
  
if(k=='A'){

snps = 200       #no of SNPs for X1
snpsc = 200         #No of SNPs for X2/X3
nobs = 100000
beta3 = 0
gamma31 = 0
}
  
  
  if(k=='B'){
    
    snps = 150       #no of SNPs for X1
    snpsc = 250         #No of SNPs for X2/X3
    nobs = 100000
    beta3 = 0
    gamma31 = 0
  }
  
  
  if(k=='C'){
    
    snps = 100       #no of SNPs for X1
    snpsc = 300        #No of SNPs for X2/X3
    nobs = 100000
    beta3 = 0
    gamma31 = 0
  }
  
  
  if(k=='D'){
    
    snps = 50       #no of SNPs for X1
    snpsc = 350         #No of SNPs for X2/X3
    nobs = 100000
    beta3 = 0
    gamma31 = 0
  }
  
  
  if(k=='E'){
    snps = 1       #no of SNPs for X1
    snpsc = 399         #No of SNPs for X2/X3               
    nobs = 100000
    beta3 = 0
    gamma31 = 0
  }
  
  if(k=='F'){
    
    snps = 200       #no of SNPs for X1
    snpsc = 200         #No of SNPs for X2/X3               
    nobs = 100000
    beta3 = 0.25
    gamma31 = 0.25
  }
  
  if(k=='G'){
    
    snps = 200       #no of SNPs for X1
    snpsc = 200         #No of SNPs for X2/X3           
    nobs = 100000
    beta3 = -0.25
    gamma31 = 0.25
  }
  
  if(k=='H'){
    
    snps = 200       #no of SNPs for X1
    snpsc = 200         #No of SNPs for X2/X3
    nobs = 300000
    beta3 = 0
    gamma31 = 0
  }
  
  return(c(snps, snpsc, nobs, beta3, gamma31))
}