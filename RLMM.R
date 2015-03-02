#rm(list=ls(all=TRUE)) 

library(MASS)
library(nlme)

my.psi<-function(u,c){
var.weights = rep(1, length(u))
sm<-median(abs(u/sqrt(var.weights)))/0.6745
w <- psi.huber(u/(sm * sqrt(var.weights)),c)
w*u
}



#Functions
REBLUP_SR<-function(y,x,group,var.weights = rep(1, nrow(x)),start.u=10, start.e=100,start.beta=c(rep(1,ncol(x))),tol=0.0001,maxit=100,k=1.345,m=20,K2=0.71,k_v=1.345)
{

assign("y",y,pos=1)
assign("x",x,pos=1)
assign("k",k,pos=1)
assign("m",m,pos=1)
assign("var.weights",var.weights,pos=1)
assign("group",group,pos=1)
assign("start.valuvu",start.valuvu,pos=1)
assign("start.valuve",start.valuve,pos=1)
assign("start.beta",start.beta,pos=1)
assign("tol",tol,pos=1)
assign("k_v",k_v,pos=1)
assign("maxit",maxit,pos=1)


n=length(y)
ni=table(group)
m=length(ni)
p=ncol(x)

z=matrix(0,n,m)
kk=0
for (j in 1:m){for(i in 1:ni[j]){
kk=kk+1
z[kk,j]=1}}

gstable=function(est){
sigma.v<-est[1]
sigma.e<-est[2]
V<-matrix(0,n,n)
V<-sigma.e*diag(1,n)+sigma.v*z%*%t(z)
U<-diag(diag(V),n,n)
svd.tmp.U=svd(U)
s<-matrix(0,2,1)
sqrt.U=t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d)))
sqrt.U.inv<-solve(t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d))))
Vi<-solve(V)
y.star1=sqrt.U.inv%*%y
x.star1=sqrt.U.inv%*%x
res.new<-y.star1-x.star1%*%betastim
res.robust=my.psi(res.new,k)

s[1,1]<-t(matrix(c(res.robust),n,1))%*%sqrt.U%*%Vi%*%z%*%t(z)%*%Vi%*%sqrt.U%*%matrix(c(res.robust),n,1)-K2*sum(diag(Vi%*%z%*%t(z)))
s[2,1]<-t(matrix(c(res.robust),n,1))%*%sqrt.U%*%Vi%*%diag(1,n,n)%*%Vi%*%sqrt.U%*%matrix(c(res.robust),n,1)-K2*sum(diag(Vi%*%diag(1,n,n)))


(s[1,1]^2)+(s[2,1]^2)

}

#STEP 1    
      

      beta.q<-matrix(start.beta,p,1)
      estsigma2u<-start.valuvu
      estsigma2e<-start.valuve
       
      diff.s<-1
      iter<-0
      while (abs(diff.s)>tol)
       {
       iter<-iter+1
#STEP 1
       diff.b<-1
       iter1<-0

       tmp=estsigma2e[iter]*diag(1,n)+estsigma2u[iter]*z%*%t(z)
       tmp.inv=solve(tmp)
       tmp.U=diag(diag(tmp),n,n)
       svd.tmp.U=svd(tmp.U)
       sqrt.ui=solve(t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d))))
       sqrt.u=t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d)))

       while (abs(diff.b)>tol)
       {iter1<-iter1+1
       
      
       y.star=sqrt.ui%*%y
       x.star=sqrt.ui%*%x
       
       res1<-(y.star-x.star%*%beta.q)
       w1<-diag(c(my.psi(res1,k)/res1),n,n)

       betastim=solve(t(x)%*%tmp.inv%*%sqrt.u%*%w1%*%sqrt.ui%*%x)%*%t(x)%*%tmp.inv%*%sqrt.u%*%w1%*%sqrt.ui%*%y

       

       diff.b<-sum((betastim-beta.q)^2)

       beta.q<-betastim

       if (iter1>maxit)
       {warning(paste("failed to converge in", maxit, "steps for beta"))
       break}

}
      
#STEP 2
        
     
opt.stable=optim(c(estsigma2u[iter],estsigma2e[iter]),gstable,method="Nelder-Mead")

estsigma2u[iter+1]<-opt.stable$par[1]
estsigma2e[iter+1]<-opt.stable$par[2]

diff.s<-sum((estsigma2u[iter+1]-estsigma2u[iter])^2+(estsigma2e[iter+1]-estsigma2e[iter])^2)
             
       if (iter>maxit)
       {warning(paste("failed to converge in", maxit, "steps"))
       break}
       }

R.tmp=estsigma2e[iter+1]*diag(1,n)
svd.R.tmp=svd(R.tmp)
sqrt.R.tmp.inv=solve(t(svd.R.tmp$v%*%(t(svd.R.tmp$u)*sqrt(svd.R.tmp$d))))

G.tmp=estsigma2u[iter+1]*diag(1,m)
svd.G.tmp=svd(G.tmp)
sqrt.G.tmp.inv=solve(t(svd.G.tmp$v%*%(t(svd.G.tmp$u)*sqrt(svd.G.tmp$d))))

tmp=estsigma2e[iter+1]*diag(1,n)+estsigma2u[iter+1]*z%*%t(z)
tmp.inv=solve(tmp)
tmp.U=diag(diag(tmp),n,n)
svd.tmp.U=svd(tmp.U)
sqrt.ui=solve(t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d))))
sqrt.u=t(svd.tmp.U$v%*%(t(svd.tmp.U$u)*sqrt(svd.tmp.U$d)))


vv.tmp<-G.tmp%*%t(z)%*%solve(R.tmp+z%*%G.tmp%*%t(z))%*%as.vector(y-x%*%beta.q)


diff.u<-1
iter2<-0
      while (abs(diff.u)>tol)
       {
       iter2<-iter2+1 
        v_robust=as.vector(vv.tmp)
        res1<-sqrt.R.tmp.inv%*%(y-x%*%beta.q-z%*%v_robust)
        res2<-sqrt.G.tmp.inv%*%v_robust
        w2<-diag(c(my.psi(res1,k_v)/res1),n,n)
        w3<-diag(c(my.psi(res2,k_v)/res2),m,m)
        A=t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv)%*%z+sqrt.G.tmp.inv%*%w3%*%sqrt.G.tmp.inv
        B<-t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv)%*%(y-x%*%beta.q)
        vv.tmp<-solve(A)%*%B
        diff.u<-sum(c((vv.tmp-v_robust)^2))
        if (iter2>maxit)
        {warning(paste("failed to converge in", maxit, "steps"))
        break}
 }
     
    rand.eff.robust=as.matrix(vv.tmp)

AA<-solve(t(x)%*%tmp.inv%*%sqrt.u%*%w1%*%sqrt.ui%*%x)%*%t(x)%*%tmp.inv%*%sqrt.u%*%w1%*%sqrt.ui
BB<-solve(t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv)%*%z+sqrt.G.tmp.inv%*%w3%*%sqrt.G.tmp.inv)%*%t(z)%*%(sqrt.R.tmp.inv)%*%w2%*%(sqrt.R.tmp.inv)

list(coefficients = beta.q, sigma2u=(estsigma2u[iter+1]),sigma2e=(estsigma2e[iter+1]),iterations=iter,randeff=rand.eff.robust,A.matrix=AA,B.matrix=BB,W1=w1,W2=w2,W3=w3)

}



