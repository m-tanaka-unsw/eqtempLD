## 
## Code for figures in the paper
## 
##   An equal tempered measure of linkage disequilibrium
##    by Mark M. Tanaka, UNSW Sydney
## 


## Multiple LD measures in a single function 
##  h = [ pAB, pAb, paB, pab ]  (proportions) 
ld <- function(h) {
    pA <- h[1]+h[2]
    pB <- h[1]+h[3]
    denom <- pA*(1-pA)*pB*(1-pB)
    D <- h[1]*h[4]-h[2]*h[3]
    r <- D/sqrt(denom) 
    if (denom==0) {
        cat("r squared denominator is zero\n")
        rsq <- 0      # to give it a value if denom=0
    } else { 
        rsq <- D^2/denom
    }
    if (D>0)
        Dmax <- min(pA*(1-pB), (1-pA)*pB)
    else 
        Dmax <- min(pA*pB, (1-pA)*(1-pB))
    Dprime <- D/Dmax
    zeta <- 2*(h[1] + h[4]) - 1
    list(D=D,Dp=Dprime,r=r,rsq=rsq,zeta=zeta,zeta.abs=abs(zeta))
}


## -----------------------------------------------
##   Changing frequencies by bits 
## -----------------------------------------------

## Make array of (relative) haplotype freqs, 
##  using vector of counts (absolute numbers)
## Rows are (a,b,c,d) + changes to frequencies 
makeH <- function(s, t="1a"){
    N <- sum(s)
    gap <- min(c(s, N-s))
    n <- c(gap:-gap) 
    H <- array(NA, c(2*gap+1, 4))
    if  (t=="1b") { 
        H[,1] <- s[1]+n
        H[,2] <- s[2]-n
        H[,3] <- s[3]
        H[,4] <- s[4]
    } else if (t=="1a") {
        H[,1] <- s[1]
        H[,2] <- s[2]
        H[,3] <- s[3]+n
        H[,4] <- s[4]-n
    } else { # assume it is "2" | preserve marginals via 2 changes
        H[,1] <- s[1]+n
        H[,2] <- s[2]-n
        H[,3] <- s[3]-n
        H[,4] <- s[4]+n
    }
    list(h=H/N, eps=n/N)
}

## Calculate LD for nearby frequencies 
## Input: frequency counts (absolute numbers) 
##   alter the frequency in a few ways
## Output: freq changes and measures of LD
LD_freq <- function(S1, type="1a") {
    mH <- makeH(S1, type)   # make array of rel freqs
    H <- mH$h
    eps <- mH$eps
    rr <- rsq <- DD <- Dp <- zeta <- zeta.abs <- c()
    for (i in 1:(dim(H)[1])) {
        s <- H[i,]
        LD <- ld(s)     # measure LD
        rr[i] <- LD$r
        rsq[i] <- LD$rsq
        DD[i] <- LD$D
        Dp[i] <- LD$Dp
        zeta[i] <- LD$zeta
        zeta.abs[i] <- LD$zeta.abs
    }
    list(eps=eps,rr=rr,rsq=rsq, zeta=zeta, azeta=zeta.abs
       , DD=DD, Dp=Dp) 
}

## Plot consequences on LD of altering frequencies
##   one/two bits at a time
change_freq <- function(cf, H=c(60, 10, 20, 25), suffix="cis"){
    ptcex <- 0.85
    onepan <- function(cf){ # plot one panel
        plot(cf$eps, cf$rsq
         , ylab="" #Linkage disequilibrium"
         , xlab="" #Change in marginal frequency (epsilon)"
         , ylim= range(cf$Dp,cf$rsq,cf$zeta,0)
         , type="n"   # don't do actual plot yet 
         , pch=19
         , las=1
           )
        abline(h=0, col="grey")
        points(cf$eps, cf$rsq, col=2, pch=23, bg=2, cex=ptcex)
        points(cf$eps, cf$Dp, col="grey55", pch=22, bg="grey55", cex=ptcex)
        points(cf$eps, cf$zeta, col=4, pch=21, bg=4, cex=ptcex)
        if (suffix=="trans") {
            points(cf$eps, abs(cf$Dp), col=1, pch=22, cex=ptcex)
            points(cf$eps, cf$azeta, col="blue3", pch=21, cex=ptcex)            
        }
    }
   ## Make the actual plots in a PDF
   fname <- paste("fig-ch-freq-LD-",suffix,".pdf", sep="")
   pdf(fname, width=5.0, height=2.5)
   par(mfrow=c(1,3), mar=c(2,2,1,1), oma=c(3,3,1,0.5)
     , cex=0.5)
   ## Panel 1
   cf <- LD_freq(H, "1a")  # generate the LD values from freqs   
   onepan(cf)
   mtext("A.  Change type 1a",3, cex=0.7, line=0.5, adj=-0.1)
   cf <- LD_freq(H, "1b")  # generate the LD values from freqs   
   onepan(cf)
   mtext("B.  Change type 1b",3, cex=0.7, line=0.5, adj=-0.1)
   if (suffix=="cis") { # smaller legend 
       legend("topleft", c("D'"
                         , expression(zeta)
                         , expression(r^2)) 
            , pch=c(22,21,23), pt.bg=c("grey55",4,2), col=c("grey55",4,2), cex=1)
   } else { # larger legend 
       legend("bottomright"
            , c("|D'|", "D'", 
                expression(paste("|",zeta,"|")),
                expression(zeta), expression(r^2)) 
            , pch=c(22,22,21,21,23)
            , col=c(1,"grey55","blue3",4,2)
            , pt.bg=c("white","grey55","white",4,2)
            , bg="white")
   }
   cf <- LD_freq(H, "2")  # generate the LD values from freqs   
   onepan(cf)
   mtext("C.  Change type 2",3, cex=0.7, line=0.5, adj=-0.1)
   mtext("Frequency change",1,line=1.5, outer=TRUE, cex=0.7)
   mtext("Linkage disequilibrium",2,line=1.5, outer=TRUE, cex=0.7)
   dev.off()
}

## Make some plots of LD and change freqs 
plot_change_freqs <- function(){
   change_freq(H=c(40, 10, 15, 35), suffix="cis")
   change_freq(H=c(30, 40, 10, 20), suffix="trans")
}


## -----------------------------------------------
##   Mixture of two subpopulations 
## -----------------------------------------------

## plot panel 
panel_two_subpops <- function(S1=c(5, 100, 5, 10),
                              S2=c(100, 10, 20, 80),
                              leg=1 # type of legend
                              ){
   s1 <- S1/sum(S1) # haplotype freqs in subpopulation 1
   s2 <- S2/sum(S2)
   rr <- rsq <- DD <- Dp <- zeta <- zeta.abs <- c()
   LD1 <- ld(s1)  # first subpop 
   LD2 <- ld(s2)  
   ms <- seq(0,1, length.out=101) # values of m parameter 
   for (i in 1:length(ms)) {
      sm <- s1*(1-ms[i]) +s2*ms[i]    # after migration or mixing 
      LDm <- ld(sm)
      rr[i] <- LDm$r
      rsq[i] <- LDm$rsq
      DD[i] <- LDm$D
      Dp[i] <- LDm$Dp 
      zeta[i] <- LDm$zeta
      zeta.abs[i] <- LDm$zeta.abs
   }
   ## Make plot for panel: 
   plot(ms,  zeta, t="l", col=4, lty=1, lwd=2
      , ylim=range(rsq,zeta,Dp,0)
      , xlim = c(0, 1.05)
      , xlab="Proportion from subpop 2, m"
      , ylab="Linkage disequilibrium"
      , las=1
        )
   abline(h=0, col="grey")
   lines(ms, Dp,  lty=1, col="grey30", lwd=1) 
   lines(ms, rsq, col=2, lwd=3)
   if (leg==1)  {# type 1: smaller legend
      legend("topleft"
           , c("D'", expression(zeta), expression(r^2))
           , lty=c(1,1,1), col=c("grey30",4,2), bg="white"  
           , lwd = c(1,2,3)
             )
   } else if (leg==2) { # larger legend with more series 
      lines(ms, abs(Dp), col=1, lty=2, lwd=2) 
      lines(ms, zeta.abs, lty=2, col="blue3", lwd=2.5)
      legend("bottomright"
           , c("|D'|", "D'",
               expression(paste("|",zeta,"|")),
               expression(zeta), expression(r^2))
           , lty=c(2,1, 2,1, 1), col=c(1,"grey30","blue3",4,2), bg="white"
           , lwd = c(2,1, 2.5,2, 3), seg.len=2
             )
   } else {
      ## no legend 
   }
}

## Make two files with plots for mixed populations 
plot_two_subpops <- function(){
   ## FIRST PDF 
   pdf("fig-mixed-pop-LD-cis-hump.pdf", width=6, height=3)
   op <- par(mfrow=c(1,2), mar=c(4,4,2,1)
           , cex=0.8) 
   panel_two_subpops(S1=c(55, 10, 25, 10), S2=c(50, 5, 5, 40), leg=1)
   mtext("A.", 3, adj=-0.05, line=0.5, cex=0.9)
   panel_two_subpops(S1=c(50, 10, 50, 10), S2=c(8.33, 16.67, 25, 50), leg=3) 
   mtext("B.", 3, adj=-0.05, line=0.5, cex=0.9)
   par(op)
   dev.off()
   ## -------------
   ## SECOND PDF 
   pdf("fig-mixed-pop-LD-trans.pdf", width=3, height=3) 
   op <- par(mar=c(4,4,0.5,0.5), cex=0.8)
   panel_two_subpops(S1=c(5, 70, 5, 20), S2=c(50, 5, 5, 40), leg=2)
   par(op) 
   dev.off()
}


## -----------------------------------------------
##  Squares and L-shapes to compare zeta and D 
## -----------------------------------------------

## Panel for unit square picture of zeta, D 
panel_square <- function(h){
   plot(c(1,1), t="n", xlim=c(0,1), ylim=c(0,1)
      , xaxs="i", yaxs="i"  # plot exactly to xlim,ylim 
      , xaxt="n", xlab=""
      , yaxt="n", ylab=""
        )
   ## label the axes, top and left side
   ## top 
    mtext(expression(p[AB]), 3, adj=0.13, line=0.3)
    mtext(expression(p[ab]), 3, adj=0.47, line=0.3)
    mtext(expression(p[Ab]), 3, adj=0.77, line=0.3)
    mtext(expression(p[aB]), 3, adj=0.97, line=0.3)
    ## left 
    mtext(expression(p[AB]), 2, padj=-6, las=1, line=0.7)
    mtext(expression(p[ab]), 2, padj=-0.5, las=1, line=0.7)
    mtext(expression(p[Ab]), 2, padj=5, las=1, line=0.7)
    mtext(expression(p[aB]), 2, padj=8, las=1, line=0.7)
    ## x boundaries going left to right
    bx1 <- h[1]
    bx2 <- h[1]+h[4]
    bx3 <- h[1]+h[4]+h[2]
    ## y boundaries going top to bottom
    by1 <- h[2]+h[3]+h[4]
    by2 <- h[2]+h[3]
    by3 <- h[3]
    polygon(c(0,0,bx2,bx2),
            c(by2,1,1,by2)
          , col="blue"
          , border=NA)
    polygon(c(bx2,bx2,1,1),
            c(0,by2,by2,0)
          , col="coral"
          , border=NA)
    ## Now the D rectangles
    polygon(c(0,0,bx1,bx1),
            c(by2,by1,by1,by2)
          , col="darkblue", border=NA)
    polygon(c(bx1,bx1,bx2,bx2),
            c(by1,1,1,by1)
          , col="darkblue", border=NA)
    polygon(c(bx2,bx2,bx3,bx3),
            c(0,by3,by3,0)
          , col="tomato3"
          , border=NA)
    polygon(c(bx3,bx3,1,1),
            c(by3,by2,by2,by3)
          , col="tomato3", border=NA)
    abline(h=c(by3, by2, by1)
         , col="grey")
    abline(v=c(bx1,bx2,bx3)
         , col="grey")
    }

## Panel for unit square elbows (L shapes)
panel_elbow <- function(h){
   plot(c(1,1), t="n", xlim=c(0,1), ylim=c(0,1)
      , xaxs="i", yaxs="i"  # plot exactly to xlim,ylim 
      , xaxt="n", xlab=""
      , yaxt="n", ylab=""
        )
   ## label the axes, top and left side
   ## top 
   mtext(expression(p[AB]), 3, adj=0.13, line=0.3)
   mtext(expression(p[ab]), 3, adj=0.47, line=0.3)
   mtext(expression(p[Ab]), 3, adj=0.77, line=0.3)
   mtext(expression(p[aB]), 3, adj=0.97, line=0.3)
   ## left 
   mtext(expression(p[AB]), 4, padj=-6, las=1, line=0.6)
   mtext(expression(p[ab]), 4, padj=-0.5, las=1, line=0.6)
   mtext(expression(p[Ab]), 4, padj=5, las=1, line=0.6)
   mtext(expression(p[aB]), 4, padj=8, las=1, line=0.6)
   ## x boundaries going left to right
   bx1 <- h[1]
   bx2 <- h[1]+h[4]
   bx3 <- h[1]+h[4]+h[2]
   ## y boundaries going top to bottom
   by1 <- h[2]+h[3]+h[4]
   by2 <- h[2]+h[3]
   by3 <- h[3]
   ## 
   polygon(c(0,0,bx1,bx1),
           c(by2,1,1,by2)
         , col="darkgreen", border=NA)
   polygon(c(bx1,bx1,bx2,bx2),
           c(by2,by1,by1,by2)
         , col="darkgreen", border=NA)
   polygon(c(bx2,bx2,bx3,bx3),
           c(0,by2,by2,0)
         , col="maroon", border=NA)
   polygon(c(bx3,bx3,1,1),
           c(0,by3,by3,0)
         , col="maroon", border=NA)
   abline(h=c(by3, by2, by1)
        , col="grey")
   abline(v=c(bx1,bx2,bx3)
        , col="grey")
}

## Put the panels together for the figure
plot_squares <- function(h=c(0.3,0.15,0.2,0.35)){
    pdf("fig-D-squares-elbows.pdf", height=2.7, width=5.4)
    op <- par(mfrow=c(1,2), mar=c(0.9,0.9,3,2), oma=c(0,3,1.1,3)
            , cex=0.4)
    panel_square(h)
    mtext("A.", 3, adj=-0.0, line=2.1, cex=1.1)
    panel_elbow(h)
    mtext("B.", 3, adj=-0.0, line=2.1, cex=1.1)
    par(op)
    dev.off()
}



## -----------------------------------------------
## DYNAMICS: Deterministic selective sweep model
## -----------------------------------------------
## Simulate dynamics of linkage disequilibrium measures
##  under a selective sweep 

## Set the parameters of the model 
setp <- function(
                 s=0.05,     # selective effect on a in homozygote
                 hs=0.025,   # selective effect on a in heterozygote
                 c=0.001,    # recomb fraction
                 D0=10^-4,   # initial D 
                 t.end=400,  # end of process
                 ict=2,      # fixed init freqs
                 p10=0.99,   # initial A freq 
                 p20=0.1     # initial B freq
                 ){
   list(s=s, hs=hs, c=c, D0=D0, t.end=t.end
      , ict=ict, p10=p10, p20=p20)
}

## Set initial conditions 
setic <- function(p){
   with(p,{
       if (ict==1) { # random initial frequencies
           p1 <- runif(1,0,1)
           p2 <- runif(1,0,1)
       } else {
           p1 <- p10
           p2 <- p20
       }
       x1 <- p1*p2 +D0
       x2 <- p1*(1-p2) -D0
       x3 <- (1-p1)*p2 -D0
       x4 <- (1-p1)*(1-p2) +D0
       if (x1<0 || x2<0 || x3<0 || x4<0 
           || x1>1 || x2>1 || x3>1 || x4>1)
           stop("Initial frequencies out of bounds") 
       c(x1, x2, x3, x4) 
   })
}


## Deterministic model of sweep
##  (see also Thomson 1977 and Maynard Smith and Haigh 1974) 
## 
##  Two loci (2 alleles each). 
## Haplotype frequencies: 
##  x1 : AB
##  x2 : Ab
##  x3 : aB
##  x4 : ab
## Marginal frequencies: 
##     p1 = x1+x2 = freq(A) 
##     p2 = x1+x3 = freq(B)
## Fitnesses:
##    A/A    1
##    A/a    1+hs
##    a/a    1+s
##  B and b are neutral 
det.sweep.model <- function(p){
   v <- 1+p$hs
   w <- 1+p$s
   ic <- setic(p) # initial D set here
   x1 <- ic[1]
   x2 <- ic[2]
   x3 <- ic[3]
   x4 <- ic[4]
   pop <- as.matrix(t(ic))
   ## record LD: D, r^2, zeta 
   lds <- ld(ic) # get initial LD values 
   LD <- as.matrix(t(c(lds$D, lds$rsq, lds$Dp, lds$zeta)))
   for (i in 1:p$t.end) {
      D <- x1*x4 - x2*x3
      x1.next <-   x1*(x1+x2) + v*x1*(x3+x4) - v*p$c*D
      x2.next <-   x2*(x1+x2) + v*x2*(x3+x4) + v*p$c*D
      x3.next <- v*x3*(x1+x2) + w*x3*(x3+x4) + v*p$c*D
      x4.next <- v*x4*(x1+x2) + w*x4*(x3+x4) - v*p$c*D
      Wbar <- (x1+x2)^2 +2*v*(x1+x2)*(x3+x4) + w*(x3+x4)^2
      ##      Wbar.alt <- x1.next +x2.next +x3.next +x4.next
      x1 <- x1.next/Wbar
      x2 <- x2.next/Wbar
      x3 <- x3.next/Wbar
      x4 <- x4.next/Wbar
      freq <- c(x1,x2,x3,x4)
      pop <- rbind(pop, unlist(freq))
      lds <- ld(freq) 
      LD <- rbind(LD, c(lds$D, lds$rsq, lds$Dp, lds$zeta)) 
   }
   list(pop=pop, LD=LD)
}


## Plot variables from numerical solution of the process
plot_dynamics <- function(p){
    y <- det.sweep.model(p) 
    pdf("fig-sweep-dynamics.pdf", width=6.4,height=2.3)
    op <- par(mfrow=c(1,3), mar=c(2, 4.5, 1.1, 0.6),
              oma=c(2.1, 0.2, 0.2, 0), cex=0.65) 
    time <- c(0:p$t.end) 

    ## First panel: frequencies of the four haplotypes 
    hcol <- c("chocolate","darkorchid",
              "darkolivegreen","dodgerblue4")
    matplot(time, y$pop, t="l"
          , ylab="" #Haplotype frequencies"
          , xlab="" #Time (generations)"
          , lwd=1.5
          , col=hcol
          , yaxt="n"
            )
    axis(2, las=2)
    mtext("A.", adj=-0.05, line=0.1, cex=0.8) 
    title(ylab="Haplotype frequencies", line=3.3)
    ## reverse order: 
    legend("right", c("ab","aB","Ab","AB"),
           lwd=1.5, col=rev(hcol), lty=c(4:1)
         , box.lwd=0.5)

    ## Second panel: D and r^2 
    matplot(time, y$LD[,1:2], t="l"
       , ylab="" 
       , xlab="" 
       , lwd=c(1,2), lty=1, col=c(1,2)
       , yaxt='n'
         )
    mtext("B.", adj=-0.05, line=0.1, cex=0.8) 
    title(ylab=expression(paste("D and ",r^2)), line=3.2)
    axis(2, at=c(0, 0.001, 0.002), las=2, cex.axis=0.85)
    legend("topright", c("D",expression(r^2))
         , lwd=c(1,2), col=c(1,2), lty=c(1,1)
         , box.lwd=0.5)
    mtext("Time (generations)", 1, cex=0.8, line=3)

   ## Third panel: zeta and D' 
    plot(time, y$LD[,4], t="l"
       , ylab="" 
       , xlab="" 
       , col=4, lwd=2     
       , yaxt="n")
    axis(2, las=2)
    lines(time, y$LD[,3], col="gray45", lty=1, lwd=1)
    legend("bottomright", c(expression(zeta),"D'")
         , lwd=c(2,1), col=c(4,"gray45"), lty=c(1,1)
         , box.lwd=0.5
           )
    mtext("C.", adj=-0.05, line=0.1, cex=0.7) 
    title(ylab=expression(paste(zeta," and D'")), line=2.7)

    par(op)
    dev.off()
}

## Plot zeta vs recomb frac at t.end for different
##  initial freq of b/B
plot_zeta_recomb <- function(t.end=1500){
    pts <- 40   # points to sample for c 
    cc <- 10^seq(log10(0.5),-4, length.out=pts)
    ss <- c(0.02, 0.05, 0.1)
    hss <- ss/2  # assume h=1/2 
    len <- length(ss)
    ## matrix: recomb x s   
    LDcsh <- matrix(NA, nrow=pts, ncol=len)
    for (i in 1:length(cc)) {
       for (j in 1:len) {
          y <- det.sweep.model(setp(c=cc[i], 
                                    s=ss[j], 
                                    hs=hss[j],
                                    t.end=t.end)) 
          LDcsh[i,j] <- y$LD[t.end,4]
       }
    }
    pdf("fig-zeta-recomb.pdf", height=2.5, width=2.5)
    op <- par(mar=c(3.5, 4.5, 0.2, 0.2) , cex=0.75)
    matplot(cc,LDcsh[,1:len], log="x", t="l", 
            col=c("dodgerblue","dodgerblue3","dodgerblue4"),
            xaxt="n", 
            lwd=c(1, 1.5, 2), 
            ylab="", 
            xlab="", 
            las=2
           )
    title(xlab="Recombination fraction", line=2.5)
    title(ylab=expression(paste(zeta," after sweep")), line=3.4)
    sfsmisc::eaxis(1, n.axp=1)
    ## ADD LEGEND 
    legend("topright", c("s=0.02","s=0.05","s=0.1"),
           lwd=c(1, 1.5, 2),             
           col=c("dodgerblue","dodgerblue3","dodgerblue4"),
           lty=c(1:3), bty="n")
    par(op)
    dev.off()
}



## -----------------------------------------------
##   Generate all the figures 
## -----------------------------------------------
make_figures <- function(){
   plot_change_freqs()
   plot_two_subpops()
   plot_squares()
   plot_dynamics(setp())
   plot_zeta_recomb()
}
