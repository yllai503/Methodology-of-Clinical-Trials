library(MASS)
N=1000
subjid=100000+(1:N)
trt01pn=rep(c(0,1), each=(N/2))
trt01p=rep(c("placebo","active"), each=(N/2))

#simulate rd and missing
#5% rd
#5% missing
set.seed(999)
rd_miss_idx=sample(1:N, (N*0.1))
rd_idx=sample(rd_miss_idx, (N*0.05))
miss_idx=rd_miss_idx[!(rd_miss_idx %in% rd_idx)]

#construct a data frame with subjid, trt01pn, trt01p, rd(=0/1), miss(=0/1), base, aval, avisitn, last on-treatment visit
A=data.frame(subjid, trt01pn, trt01p, rd=rep(0, N), miss=rep(0, N), base=rep(NA, N), aval=rep(NA, N))
A[rd_idx, "rd"]=1
A[miss_idx, "miss"]=1

beta_i=c(0, -0.05, -0.1, -0.2, -0.25)
beta_t=c(0, -0.01, -0.05, -0.1, -0.2)
beta_mnar=c(0, 0, 0, 0, 0.25)
varmat=matrix(0.6, 5, 5)
diag(varmat)=1

A_long=NULL
for (i in 1:N){
  y_main=8.5+beta_t+beta_i*as.numeric(i>(N/2))+beta_mnar*as.numeric(i>(N/2))*as.numeric(A[i, "rd"]==1)
  y_i=mvrnorm(1, y_main, varmat)
  subj_full=do.call(rbind, replicate(5, A[i,], simplify=FALSE))
  subj_full$base=y_i[1]
  subj_full$aval=y_i
  subj_full$avisitn=c(0, 6, 12, 18, 26)
  subj_full$last_ontrt=26
  subj_full$last_ontrtval=y_i[5]
  if (A[i, "miss"]==1){
    #simulate last on-treatment visit to be either row 2, 3, 4
    miss_r=sample(c(2,3,4),1)
    #only visits till last on-treatment visit are kept to mimic clinical trial data
    subj_full=subj_full[1:miss_r, ]
    subj_full$last_ontrt=subj_full$avisitn[miss_r]
    subj_full$last_ontrtval=y_i[miss_r]
  }
  if (A[i, "rd"]==1){
    #simulate last on-treatment visit to be either row 2, 3, 4
    miss_r=sample(c(2,3,4),1)
    
    subj_full$last_ontrt=subj_full$avisitn[miss_r]
    subj_full$last_ontrtval=y_i[miss_r]
  }
  A_long=rbind(A_long, subj_full)
}
