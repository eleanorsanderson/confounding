
setup <- function(k){
  
if(k=='A'){

snps = 200       #no of SNPs for X1
snpsc = 200         #No of SNPs for X2/X3
nobs = 100000
}
  
  

  
  if(k=='B'){
    
    snps = 100       #no of SNPs for X1
    snpsc = 300        #No of SNPs for X2/X3
    nobs = 100000

  }
  
  
  if(k=='C'){
    
    snps = 50       #no of SNPs for X1
    snpsc = 350         #No of SNPs for X2/X3
    nobs = 100000

  }
  

  
  
  return(c(snps, snpsc, nobs))
}