# Bayesian-Emulation-YT
The code for my project on Bayesian Emulation
# 1-D Emulator
xi=c(0.1,0.2,0.3,0.4,0.5)
xi=c(0,0.1,0.2,0.23,0.28,0.3,0.4,0.5)
f= function(n){exp(3.5*n)}
D=f(xi)
B=3.5
sig=1.5
theta=0.14
varf=sig^2
ED=c(rep(B,length(xi)))

covv= function(n){
v=numeric(length(xi))
	for(i in 1:length(v)){
		v[i]=varf*exp(-((n-xi[i])^2)/theta^2)}
return(v)}
m=matrix(,length(xi),length(xi))
for(i in 1:length(xi)){for(j in 1:length(xi)){
m[i,j]=varf*exp(-((xi[i]-xi[j])^2)/theta^2)}}
m
ed=function(n){g=B+(covv(n)%*%solve(m))%*%(D-ED)
return(as.numeric(g))}
vard=function(n){p=varf-(covv(n)%*%solve(m))%*%covv(n)
return(p)}

r=numeric(500)
for (i in 1:500){r[i]=ed(i/1000)}
y=seq(0.001,0.5,0.001)

rd1=function(n){ed(n)+3*sqrt(vard(n))}
rd2=function(n){ed(n)-3*sqrt(vard(n))}

s=numeric(500)
for (i in 1:500){s[i]=rd1(i/1000)}

d=numeric(500)
for (i in 1:500){d[i]=rd2(i/1000)}
pdf("3his1d.pdf", width=5, height=5)
plot(y,r,col="blue",type="l",lwd=2,
xlab="x",ylab="f(x)",cex.lab=1.5,cex.axis=1.5)
lines(y,s,col="red",lwd=2)
lines(y,d,col="red",lwd=2)
lines(y,thre,col="black",lwd=2)
for(i in 1:length(xi)){points(xi[i],ed(xi[i]),pch=19)}
dev.off()
z=2.5
im=function(x){imp=(ed(x)-z)^(2)/(vard(x)+0.05)
return(as.numeric(imp))}
thre=rep(2.5,500)
three=rep(3,500)
rd=numeric(500)
for (i in 1:500){rd[i]=sqrt(im(i/1000))}
pdf("4his1d.pdf", width=5, height=5)
plot(y,rd,col="blue",type='l',lwd=2,
xlab="x",ylab="Implausibility",cex.lab=1.5,cex.axis=1.5)
lines(y,three,col="black",lwd=2)
dev.off()

# 2-D Emulator
#Plotting Original
f=function(x){cos(x[1])+sin(x[2])}
#f=function(x){2*cos(x[1])-0.5*sin(2*x[2])}
xii=seq(0.1,10,0.1)
ma=matrix(,length(xii),length(xii))
for(i in 1:length(xii)){for(j in 1:length(xii)){
ma[i,j]=f(c(xii[i],xii[j]))}}
filled.contour(xii,xii,ma)

#Wave 1
m1=as.matrix(expand.grid(xi,xi))
mi=m1

#Wave 2
m2=matrix(c(8,8,8,2,5,8),3,2)
mi=rbind(m1,m2)

#Wave 3
m3=matrix(c(5,5,5,1,5,8),3,2)
mi=rbind(m1,m2,m3)

#Wave 4
m4=matrix(c(1,1,1,1,5,8),3,2)
mi=rbind(m1,m2,m3,m4)

#Parameters
xi=seq(0,10,10/3)
B=0
sig=2
varf=sig^2
theta=3
le=length(mi)/2
D=numeric(le)
for(i in 1:le){D[i]=f(mi[i,])}
ED=c(rep(B,le))
y=seq(0,10,0.1)
yi=seq(0,10,0.25)
z=1
matr=matrix(,41,41)
matri=matrix(,41,41)
mat=matrix(,41,41)
mtt=matrix(,41,41)
mtt2=matrix(,41,41)
mtt3=matrix(,41,41)
mtt4=matrix(,41,41)
sq=c(0,3,10,20,30,40,50,60)
sqq=c(-20,-15,-10,-5,-3,3,5,10,15,20)

#Functions
covv=function(x){
v=numeric(le)
for(i in 1:le){v[i]=varf*exp(-dist(rbind(x,mi[i,]))^2/theta^2)}
return(v)}

m=matrix(,le,le)
for(i in 1:le){for(j in 1:le){
m[i,j]=varf*exp(-((dist(rbind(mi[i,],mi[j,])))^2)/theta^2)}}

ed=function(n){g=B+(covv(n)%*%solve(m))%*%(D-ED)
return(as.numeric(g))}

vard=function(n){p=varf-(covv(n)%*%solve(m))%*%covv(n)
return(p)}

sd=function(n){l=sqrt(vard(n))
return(l)}

im=function(x){imp=sqrt((ed(x)-z)^(2)/(vard(x)+0.04^2))
return(as.numeric(imp))}

dia=function(x){cc=(ed(x)-f(x))/sd(x)
return(cc)}
ed(c(3,5))
#Matrix For Plots
for(i in 1:41){
	for(j in 1:41){
		matri[i,j]=sd(c(i/4.1,j/4.1))}}
for(i in 1:41){
	for(j in 1:41){
		mat[i,j]=ed(c(i/4.1,j/4.1))}}
for(i in 1:41){
	for(j in 1:41){
		matr[i,j]=im(c(i/4.1,j/4.1))}}
for(i in 1:41){
	for(j in 1:41){
		mtt[i,j]=dia(c(i/4.1,j/4.1))}}
length(which(mtt>3))
length(which(mtt< -3))

#Plots
pdf("W4Plots.pdf", width=8, height=7)
pdf("2DEX.pdf", width=6, height=5)
filled.contour(yi,yi,mat,#main = "Expectation",
    xlab = "x", ylab = "y")
dev.off()
pdf("2DSD.pdf", width=6, height=5)
filled.contour(yi,yi,matri,#main = "Standard Deviation",
    xlab = "x", ylab = "y",color.palette=terrain.colors,#cex.lab=2,
plot.axes = {axis(1);axis(2);points(mi,pch=19)})
dev.off()
filled.contour(yi,yi,matr,color.palette=heat.colors,
plot.title = title(main = "Implausibility",xlab = "X", ylab = "Y"),
levels=sq,key.axes=axis(4,sq))
filled.contour(yi,yi,mtt,color.palette=topo.colors,
plot.title = title(main = "Diagnostic",
xlab = "X", ylab = "Y"),levels=sqq,key.axes=axis(4,sqq))
dev.off()

pdf("New Implausibility.pdf", width=8, height=7)
filled.contour(yi,yi,matr,color.palette=heat.colors,
plot.title = title(main = "W1 Implausibility",xlab = "X", ylab = "Y"),
levels=sq,key.axes=axis(4,sq),
plot.axes = {axis(1);axis(2);points(m2,pch=19)})

filled.contour(yi,yi,mtt2,color.palette=heat.colors,
plot.title = title(main = "W2 Implausibility",xlab = "X", ylab = "Y"),
levels=sq,key.axes=axis(4,sq),
plot.axes = {axis(1);axis(2);points(m3,pch=19)})

filled.contour(yi,yi,mtt3,color.palette=heat.colors,
plot.title = title(main = "W3 Implausibility",xlab = "X", ylab = "Y"),
levels=sq,key.axes=axis(4,sq),
plot.axes = {axis(1);axis(2);points(m4,pch=19)})

filled.contour(yi,yi,mtt4,color.palette=heat.colors,
plot.title = title(main = "W4 Implausibility",xlab = "X", ylab = "Y"),
levels=sq,key.axes=axis(4,sq))

# 1-D Emulator with derivatives

#Initial function
f=function(x){0.5*x+cos(x)}
df=function(x){0.5-sin(x)}
sq=seq(0,20,0.01)

#Parameters
xi=c(0,5,10,15,20)
D=c(f(xi),df(xi))
B=5
dB=0
sig=4
theta=4
varf=sig^2
ED=c(rep(B,length(xi)),rep(dB,length(xi)))
st=varf/theta^2
st2=st/theta^2
le=length(xi)

#Functions
covv= function(n){
v=numeric(length(xi))
	for(i in 1:length(v)){
		v[i]=varf*exp(-((n-xi[i])^2)/theta^2)}
return(v)}

dcovv=function(n){
k=numeric(length(xi))
	for(i in 1:length(k)){
		k[i]=theta^(-2)*2*(n-xi[i])*varf*exp(-((n-xi[i])^2)/theta^2)}
return(k)}

cvr=function(x){b=c(covv(x),dcovv(x))
return(b)}

m=matrix(,2*le,2*le)
for(i in 1:le){for(j in 1:le){
m[i,j]=varf*exp(-((xi[i]-xi[j])^2)/theta^2)}
for(i in 1:le){for(j in (le+1):(2*le)){
m[i,j]=(theta^(-2)*-2*(xi[(j-le)]-xi[i]))*varf*exp(-((xi[i]-xi[(j-le)])^2)/theta^2)}}
for(i in (le+1):(2*le)){for(j in 1:le){
m[i,j]=(theta^(-2)*2*(xi[(j)]-xi[i-le]))*varf*exp(-((xi[i-le]-xi[(j)])^2)/theta^2)}}
for(i in (le+1):(2*le)){for(j in (le+1):(2*le)){
m[i,j]=(2*(st)-4*st2*(xi[(j-le)]-xi[i-le])^2)*exp(-((xi[i-le]-xi[j-le])^2)/theta^2)}}}
m

ed=function(n){g=B+(cvr(n)%*%solve(m))%*%(D-ED)
return(as.numeric(g))}
vard=function(n){p=varf-(cvr(n)%*%solve(m))%*%cvr(n)
return(p)}
rd1=function(n){ed(n)+3*sqrt(vard(n))}
rd2=function(n){ed(n)-3*sqrt(vard(n))}
av=numeric(2000)
for (i in 0:2000){av[i]=(vard(i/100))}
sum(av)/2000
s=numeric(2000)
for (i in 0:2000){s[i]=rd1(i/100)}

d=numeric(2000)
for (i in 0:2000){d[i]=rd2(i/100)}

r=numeric(2000)
for (i in 0:2000){r[i]=ed(i/100)}
y=seq(0.01,20,0.01)

pdf("DerPlot2.pdf", width=5, height=5)
plot(y,r,col="blue",type="l",lwd=2,xlab = "x", ylab = "f(x)",
cex.lab=1.5,cex.axis=1.5)
lines(y,s,col="red",lwd=2)
lines(y,d,col="red",lwd=2)
lines(sq,f(sq),col='green',type="l",lwd="2")
for(i in 1:length(xi)){points(xi[i],ed(xi[i]),pch=19,cex=1)}
dev.off()
#History Match
z=7
im=function(x){imp=(ed(x)-z)^(2)/(vard(x)+0.005)
return(as.numeric(sqrt(imp)))}
rd1=numeric(2000)
for (i in 1:2000){rd1[i]=im(i/100)}
plot(y,rd1,col="blue",type='l',lwd=2,
xlab="x",ylab="Implausibility (z=7)")
which(rd<3)
yyy=seq(11.8,16.37,0.01)

#Run 2 with derivatives
xi=c(0,5,10,13,15,20)

pdf("ImpDev3.pdf", width=8, height=7)
plot(y,rd,col="blue",type='l',lwd=2,
xlab="x",ylab="Implausibility (z=7)")

plot(yyy,rd[1180:1637],col="blue",type='l',lwd=2,
xlab="x",ylab="Implausibility (z=7)")

# Comparison with derivative emulator

#Initial function
f=function(x){0.5*x+cos(x)}
#f=function(x){0.5*x+cos(x+2*pi*x)}
sq=seq(0,20,0.01)
plot(sq,f(sq),type="l",lwd="2")

#Parameters
xi=seq(0,20,2.5)
#xi=seq(0,20,20/9)
#xi=seq(0,20,5)
length(xi)
D=f(xi)
B=5
sig=4
theta=4
varf=sig^2
ED=c(rep(B,length(xi)))
v=numeric(length(xi))

covv= function(n){
	for(i in 1:length(v)){
		v[i]=varf*exp(-((n-xi[i])^2)/theta^2)}
return(v)}
m=matrix(,length(xi),length(xi))
for(i in 1:length(xi)){for(j in 1:length(xi)){
m[i,j]=varf*exp(-((xi[i]-xi[j])^2)/theta^2)}}
m
ed=function(n){g=B+(covv(n)%*%solve(m))%*%(D-ED)
return(as.numeric(g))}
vard=function(n){p=varf-(covv(n)%*%solve(m))%*%covv(n)
return(p)}
rd1=function(n){ed(n)+3*sqrt(vard(n))}
rd2=function(n){ed(n)-3*sqrt(vard(n))}
av=numeric(2000)
for (i in 0:2000){av[i]=vard(i/100)}
sum(av)/2000
s=numeric(2000)
for (i in 0:2000){s[i]=rd1(i/100)}

d=numeric(2000)
for (i in 0:2000){d[i]=rd2(i/100)}
r=numeric(2000)
for (i in 1:2000){r[i]=ed(i/100)}
y=seq(0.01,20,0.01)
pdf("5Chap5.pdf", width=5, height=5)
plot(y,r,col="blue",type="l",lwd=2,xlab="x",ylab="f(x)")
lines(y,s,col="red",lwd=2)
lines(y,d,col="red",lwd=2)
lines(sq,f(sq),col='green',lwd="2")
for(i in 1:length(xi)){points(xi[i],ed(xi[i]),pch=19)}

# Emulating derivative of a function

#Initial function
f=function(x){0.5*x+cos(x)}
df=function(x){0.5-sin(x)}
sq=seq(0,20,0.01)
D
#Parameters
xi=seq(0,20,2.5)
D=c(df(xi))
B=0
dB=0
sig=7
theta=5
varf=sig^2
ED=c(rep(dB,length(xi)))
st=varf/theta^2
st2=st/theta^2
le=length(xi)
covv(4)
covv= function(x){
v=numeric(length(xi))
	for(i in 1:length(v)){
		v[i]=(2*(st)-4*st2*(xi[i]-x)^2)*exp(-((x-xi[i])^2)/theta^2)}
return(v)}
m=matrix(,length(xi),length(xi))
for(i in 1:length(xi)){for(j in 1:length(xi)){
m[i,j]=(2*(st)-4*(st2)*(xi[j]-xi[i])^2)*exp(-((xi[i]-xi[j])^2)/theta^2)}}
m[1,1]
ed=function(n){g=B+(covv(n)%*%solve(m))%*%(D-ED)
return(as.numeric(g))}

vard=function(n){p=2*varf/theta^2-(covv(n)%*%solve(m))%*%covv(n)
return(p)}
rd1=function(n){ed(n)+3*sqrt(vard(n))}
rd2=function(n){ed(n)-3*sqrt(vard(n))}

z=-0.49
im=function(x){imp=(ed(x)-z)^(2)/(vard(x)+0.005)
return(as.numeric(sqrt(imp)))}

rd=numeric(2000)
for (i in 1:2000){rd[i]=im(i/100)}

s=numeric(2000)
for (i in 0:2000){s[i]=rd1(i/100)}

d=numeric(2000)
for (i in 0:2000){d[i]=rd2(i/100)}

r=numeric(2000)
for (i in 0:2000){r[i]=ed(i/100)}
y=seq(0.01,20,0.01)
pdf("6Chap6.pdf", width=10, height=5)
plot(y,r,col="blue",type="l",lwd=2,xlab="x",ylab="f'(x)")
lines(y,s,col="red",lwd=2)
lines(y,d,col="red",lwd=2)
lines(sq,df(sq),col='green',type="l",lwd="2")
lines(y,rep(-0.49,length(y)))
dev.off()
pdf("7Chap6.pdf", width=10, height=5)
plot(y,rd,col="blue",type='l',lwd=2,
xlab="x",ylab="Implausibility (z=-0.49)")
lines(y,rep(3,length(y)))
dev.off()

# 2-D derivative emulator

f=function(x){2*cos(x[1])-0.5*sin(2*x[2])}
xii=seq(0.1,10,0.1)
ma=matrix(,length(xii),length(xii))
for(i in 1:length(xii)){for(j in 1:length(xii)){
ma[i,j]=f(c(xii[i],xii[j]))}}
pdf("3Chap7.pdf", width=6, height=5)
filled.contour(xii,xii,ma,xlab = expression('x'[1]), ylab = expression('x'[2]),
levels=sqqqq,key.axes=axis(4,sqqq)
)
dev.off()
d1f=function(x){-2*sin(x[1])}
d2f=function(x){cos(2*x[2])}
sqqqq=seq(-3.5,3,0.5)
sqqq=c(-3,-2,-1,0,1,2,3)
#Parameters
xi=seq(0,10,10/3)
m1=as.matrix(expand.grid(xi,xi))
mi=m1
le=length(mi)/2
D1=numeric(le)
for(i in 1:le){D1[i]=f(mi[i,])}
D2=numeric(le)
for(i in 1:le){D2[i]=d1f(mi[i,])}
D3=numeric(le)
for(i in 1:le){D3[i]=d2f(mi[i,])}
D=c(D1,D2,D3)
B=1
dB=0
sig=2
theta=4
th2=theta^(-2)
th4=theta^(-4)
varf=sig^2
ED=c(rep(B,le),rep(dB,2*le))

#Variance Matrix
#(1,1)
m=matrix(,3*le,3*le)
for(i in 1:le){for(j in 1:le){
m[i,j]=varf*exp(-((dist(rbind(mi[i,],mi[j,])))^2)/theta^2)}}
#(1,2)
for(i in 1:le){for(j in (le+1):(2*le)){
m[i,j]=theta^(-2)*2*(mi[i,1]-mi[j-le,1])*varf*
exp(-((dist(rbind(mi[i,],mi[j-le,])))^2)/theta^2)}}
#(1,3)
for(i in 1:le){for(j in (2*le+1):(3*le)){
m[i,j]=theta^(-2)*2*(mi[i,2]-mi[j-2*le,2])*varf*
exp(-((dist(rbind(mi[i,],mi[j-2*le,])))^2)/theta^2)}}
#(2,1)
for(i in (le+1):(2*le)){for(j in 1:le){
m[i,j]=theta^(-2)*-2*(mi[i-le,1]-mi[j,1])*varf*
exp(-((dist(rbind(mi[i-le,],mi[j,])))^2)/theta^2)}}
#(2,2)
for(i in (le+1):(2*le)){for(j in (le+1):(2*le)){
m[i,j]=varf*(2*th2-4*th4*(mi[(i-le),1]-mi[j-le,1])^2)*
exp(-((dist(rbind(mi[i-le,],mi[j-le,])))^2)/theta^2)}}
#(2,3)
for(i in (le+1):(2*le)){for(j in (2*le+1):(3*le)){
m[i,j]=varf*(-4*th4*(mi[(i-le),1]-mi[j-2*le,1])*(mi[(i-le),2]-mi[j-2*le,2]))*
exp(-((dist(rbind(mi[i-le,],mi[j-2*le,])))^2)/theta^2)}}
#(3,1)
for(i in (2*le+1):(3*le)){for(j in 1:le){
m[i,j]=theta^(-2)*-2*(mi[i-2*le,2]-mi[j,2])*varf*
exp(-((dist(rbind(mi[i-2*le,],mi[j,])))^2)/theta^2)}}
#(3,2)
for(i in (2*le+1):(3*le)){for(j in (le+1):(2*le)){
m[i,j]=varf*(-4*th4*(mi[(i-2*le),1]-mi[j-le,1])*(mi[(i-2*le),2]-mi[j-le,2]))*
exp(-((dist(rbind(mi[i-2*le,],mi[j-le,])))^2)/theta^2)}}
#(3,3)
for(i in (2*le+1):(3*le)){for(j in (2*le+1):(3*le)){
m[i,j]=varf*(2*th2-4*th4*(mi[(i-2*le),2]-mi[j-2*le,2])^2)*
exp(-((dist(rbind(mi[i-2*le,],mi[j-2*le,])))^2)/theta^2)}}

#Covariance Function
covv=function(x){
v=numeric(le)
for(i in 1:le){v[i]=varf*exp(-dist(rbind(x,mi[i,]))^2/theta^2)}
return(v)}

d1covv=function(x){
v1=numeric(le)
for(i in 1:le){v1[i]=theta^(-2)*2*(x[1]-mi[i,1])*varf*
exp(-((dist(rbind(x,mi[i,])))^2)/theta^2)}
return(v1)}

d2covv=function(x){
v2=numeric(le)
for(i in 1:le){v2[i]=theta^(-2)*2*(x[2]-mi[i,2])*varf*
exp(-((dist(rbind(x,mi[i,])))^2)/theta^2)}
return(v2)}

cvr=function(x){b=c(covv(x),d1covv(x),d2covv(x))
return(b)}

#Exp and Var
ed=function(n){g=B+(cvr(n)%*%solve(m))%*%(D-ED)
return(as.numeric(g))}
vard=function(n){p=varf-(cvr(n)%*%solve(m))%*%cvr(n)
return(p)}
sd=function(n){l=sqrt(vard(n))
return(l)}
rd1=function(n){ed(n)+3*sqrt(vard(n))}
rd2=function(n){ed(n)-3*sqrt(vard(n))}
mat=matrix(,41,41)
matri=matrix(,41,41)
matt=matrix(,41,41)
#Diagnostic
dia=function(x){cc=(ed(x)-f(x))/sd(x)
return(cc)}

#Matrix for plots
for(i in 1:41){
	for(j in 1:41){
		matri[i,j]=sd(c(i/4.1,j/4.1))}}
for(i in 1:41){
	for(j in 1:41){
		mat[i,j]=ed(c(i/4.1,j/4.1))}}
for(i in 1:41){
	for(j in 1:41){
		mtt[i,j]=dia(c(i/4.1,j/4.1))}}
yi=seq(0,10,0.25)

#Plots
pdf("1Chap7.pdf", width=6, height=5)
filled.contour(yi,yi,mat,
    xlab = expression('x'[1]), ylab = expression('x'[2]),
plot.axes = {axis(1);axis(2);points(mi,pch=19)})
dev.off()
pdf("2Chap7.pdf", width=6, height=5)
filled.contour(yi,yi,matri,
    xlab = expression('x'[1]), ylab = expression('x'[2])
,color.palette=terrain.colors,
plot.axes = {axis(1);axis(2);points(mi,pch=19)})
dev.off()
filled.contour(yi,yi,matri,
    xlab = expression('x'[1]), ylab = expression('x'[2])
,color.palette=topo.colors,
plot.axes = {axis(1);axis(2);points(mi,pch=19)}

