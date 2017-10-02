inv_logit3 <- function(m,mstar,a,b,c) {
  z = (m-mstar)*cos(atan(b))-a*sin(atan(b))
  #zstar = mstar*cos(atan(b))-a*sin(atan(b))
  #1/(1+ (z/mstar)^-c)
  1/(1+exp(-(c*z)))
}
sc <- function(x) x/max(x,na.rm=T)
lw <- function(w) (w/0.01)^(1/3)
wl <- function(l) 0.01*l^3

O2_plotfun <- function(etas,Tops){
  grid <- expand.grid(eta = etas,
                      Topt = Tops)
  
  O2fun <- function(x) {
    ts = with(x,((30-seq(1,30,l=1000))/(30-Topt))^eta*exp(-eta*(30-seq(1,30,l=1000))/(30-Topt)))
    ts = ts/max(ts)
  }
  
  
  
  require(dplyr)
  require(purrr)
  require(ggplot2)
  g <- grid %>% 
    by_row(~O2fun(.x), .collate='rows') %>%
    group_by(.row) %>%
    mutate(Temperature=seq(1,30,l=1000)) 
  
  ggplot(g) + 
    geom_line(aes(x = Temperature, 
                  y = .out,
                  linetype=factor(eta))) + 
    theme_bw() + 
    facet_grid(Topt~.) +
    scale_linetype_discrete(bquote(eta)) + 
    ylab(paste(expression(O_2, supply))) # + 
  #scale_color_discrete(bquote(T[opt]), h.start = 180)
  
}

plot_data_temp <- function(v){
  
  lm = seq(10^v$min_m,10^v$max_m,l=v$lm)
  
  O2_tcor <- O2_fact(v$temp,5)
  f=O2_supply(10*O2_tcor,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=v$temp,delta=v$lO,omega=v$shape,P50=v$P50,n=v$n)
  #browser()
  max_tau <- eval_tau_max_temp(f=f,
                               temp=v$temp,
                               temp_ref=v$temp_ref,
                               Ea=v$Ea,
                               omega = v$omega,
                               gamma=v$gamma,
                               delta=v$delta,
                               phi=v$phi,
                               h=v$h,
                               beta=v$beta,
                               k=v$k,
                               p=v$p,
                               q=v$q,
                               n=v$n,
                               m=v$m)
  tau <- eval_tau_eq_temp(temp=v$temp,
                          temp_ref=v$temp_ref,
                          Ea=v$Ea,
                          gamma=v$gamma,
                          delta=v$delta,
                          phi=v$phi,
                          h=v$h,
                          beta=v$beta,
                          k=v$k,
                          p=v$p,
                          q=v$q,
                          n=v$n,
                          m=v$m,
                          M=v$M,
                          v=v$v)
  
  tau_max = apply(cbind(tau,max_tau),1,min)
  tau_max[tau_max<0] <- 0
  tau_max[tau_max>1] <- 1
  
  model_frame <- model_out(tau_max = tau_max,
                           temp=v$temp,
                           temp_ref=v$temp_ref,
                           Ea=v$Ea,
                           gamma=v$gamma,
                           delta=v$delta,
                           phi=v$phi,
                           h=v$h,
                           beta=v$beta,
                           k=v$k,
                           p=v$p,
                           q=v$q,
                           n=v$n,
                           m=v$m)
  
  
  #browser()
  scope <- f*v$m^v$n-model_frame$Metabolism*v$omega
  scope[scope<0.0001] <- 0
  #browser()
  
  
  
  bind_cols(model_frame,
            data_frame(
              Temperature=v$temp,
              Realised = tau_max,
              Limit = sapply(sapply(max_tau,max,0),min,1),
              M = (tau_max*v$v+v$M)*v$m^v$nu,
              Optimum = sapply(sapply(tau,max,0),min,1),
              Scope=ifelse(scope<0.0001,0,scope),
              `Max Metabolism` = f*v$m^v$n,
              `Active Metabolism` = model_frame$Metabolism*v$omega,
              `Std Metabolism` = model_frame$Std*v$omega,
              Viable = as.numeric(tau_max>0.0001 & model_frame[['C for growth']]>0.0001))
  )
  
}

plot_data_growth_tO2 <- function(v){
  
  n_int <- length(v$temp)
  lm = seq(10^v$min_m,10^v$max_m,l=v$lm)
  O2_range <- seq(1,10,l=n_int)
  
  #### Step 1 get mstar=mo+slope*age for all Temp
  #browser()
  
  mout <- model_out_growth(temp=v$temp,
                           temp_ref=v$temp_ref,
                           l=v$lm,
                           Ea=v$Ea,
                           gamma=v$gamma,
                           delta=v$delta,
                           phi=v$phi,
                           h=v$h,
                           beta=v$beta,
                           k=v$k,
                           p=v$p,
                           q=v$q,
                           n=v$n,
                           tmax = v$tmax,
                           slope=v$slope,
                           tr=v$tr,
                           v=v,
                           lm=lw(lm),
                           dt = v$dt)
  
  winfs <- mout$winfs
  
  ### Step 2 calc mat for selected response 
  
  #browser()
  
  mouts <- model_out_growth_check(temp=winfs$Temperature,
                                  temp_ref=v$temp_ref,
                                  Ea=v$Ea,
                                  gamma= seq(v$gamma*(1-v$c),v$gamma*(1+v$c),l=n_int),
                                  delta=v$delta,
                                  phi=v$phi,
                                  h=v$h,
                                  beta=v$beta,
                                  k=v$k,
                                  p=v$p,
                                  q=v$q,
                                  n=v$n,
                                  mstar = lw(lm[winfs$opt[which.min(abs(winfs$Temperature-v$temp_ref))]]),
                                  tmax = v$tmax,
                                  slope=v$slope,
                                  tr=v$tr,
                                  v=v,
                                  lm=lw(lm),
                                  dt = v$dt)
  
  list(winfs=winfs,
       growth=mout$growth,
       g_growth=mouts$g_growth,
       t_growth=mouts$t_growth)
  
}


# Activity model

eval_tau_eq <- function(gamma=50,
                        delta=2,
                        phi=10,
                        h=30,
                        beta=0.75,
                        k=2,
                        p=0.8,
                        q=0.9,
                        n=0.8,
                        m=100,
                        M=0.2,
                        v=1){
  
  
  tc=1
  -(M*delta*h*k*m^(n - p + q)*tc - h*k*m^(n - p + q)*tc*v - sqrt(-((beta - 1)*h*k*m^(n - 2*p + 3*q) + h*k*m^(n - 2*p + 3*q)*phi)*tc^2*v^2 - ((beta - 1)*delta*gamma*k*m^(n - p + 2*q) + delta*gamma*k*m^(n - p + 2*q)*phi)*M^2*tc + (((beta - 1)*delta*h*k*m^(n - 2*p + 3*q) + delta*h*k*m^(n - 2*p + 3*q)*phi)*tc^2 + (gamma*h*m^(-p + 3*q)*phi^2 + (beta - 1)*gamma*k*m^(n - p + 2*q) + (beta^2 - 2*beta + 1)*gamma*h*m^(-p + 3*q) + (2*(beta - 1)*gamma*h*m^(-p + 3*q) + gamma*k*m^(n - p + 2*q))*phi)*tc)*M*v)*h)/(M*delta*gamma*k*m^n - ((beta - 1)*gamma*h*m^q + gamma*h*m^q*phi + gamma*k*m^n)*v)
  
  
}

eval_tau_eq_temp <- function(Ea,
                             temp,
                             temp_ref=15,
                             gamma=50,
                             delta=2,
                             phi=10,
                             h=30,
                             beta=0.75,
                             k=2,
                             p=0.8,
                             q=0.9,
                             n=0.8,
                             m=100,
                             M=0.2,
                             v=1){
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  -(M*delta*h*k*m^(n - p + q)*tc - h*k*m^(n - p + q)*tc*v - sqrt(-((beta - 1)*h*k*m^(n - 2*p + 3*q) + h*k*m^(n - 2*p + 3*q)*phi)*tc^2*v^2 - ((beta - 1)*delta*gamma*k*m^(n - p + 2*q) + delta*gamma*k*m^(n - p + 2*q)*phi)*M^2*tc + (((beta - 1)*delta*h*k*m^(n - 2*p + 3*q) + delta*h*k*m^(n - 2*p + 3*q)*phi)*tc^2 + (gamma*h*m^(-p + 3*q)*phi^2 + (beta - 1)*gamma*k*m^(n - p + 2*q) + (beta^2 - 2*beta + 1)*gamma*h*m^(-p + 3*q) + (2*(beta - 1)*gamma*h*m^(-p + 3*q) + gamma*k*m^(n - p + 2*q))*phi)*tc)*M*v)*h)/(M*delta*gamma*k*m^n - ((beta - 1)*gamma*h*m^q + gamma*h*m^q*phi + gamma*k*m^n)*v)
  
}




model_out <- function(tau_max,
                      temp,
                      temp_ref=15,
                      Ea,
                      r=0.2,
                      gamma=50,
                      delta=2,
                      phi=0.15,
                      h=30,
                      beta=0.2,
                      k=2,
                      p=0.8,
                      q=0.9,
                      n=0.8,
                      m=100){
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  f <- tau_max*gamma*m^p/(tau_max*gamma*m^p+tc*h*m^q)
  inp <- (1-phi-beta)*f*tc*h*m^q
  out <- (1+tau_max*delta)*k*tc*m^n
  e <-  inp -out
  efficiency <- e/f
  efficiency[efficiency<0] <- 0
  predation_rate <- f/phi
  
  met = beta*f*tc*h*m^q+ ((1+tau_max*delta)*k*tc*m^n)
  
  #browser()
  
  data_frame(`Feeding level`= f, 
             Consumption = inp,
             `C used for Metabolism` = out,
             `C for growth` = e, 
             Efficiency = efficiency, 
             `Predation rate`=predation_rate, 
             Metabolism = met,
             Std = k*tc*m^n)
}

model_out_growth <- function(temp,
                             temp_ref=15,
                             l,
                             lm,
                             Ea,
                             gamma=50,
                             delta=2,
                             phi=0.15,
                             h=30,
                             beta=0.2,
                             k=2,
                             p=0.8,
                             q=0.9,
                             n=0.8,
                             tmax=10,
                             slope=0.05,
                             tr = 1,
                             v=NULL,
                             dt=100){
  
  
  temps = length(temp)
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  ts <- seq(0,tmax,l=dt)
  
  
  s <- array(0, c(temps,l,dt))
  allocs <- array(0, c(temps,l,dt))
  
  R0 <- array(0, c(temps,l,dt))
  surv <- array(0, c(temps,l,dt))
  
  s[,,1] <- min(wl(lm))
  
  
  dts <- (tmax/(dt-1))
  
  for(t in 2:dt) {
    tm1 <- get_taus(v,1,10,temp,s[,,t-1])
    f <- tm1*gamma*s[,,t-1]^p/(tm1*gamma*s[,,t-1]^p+tc*h*s[,,t-1]^q)
    inp <- (1-phi-beta)*f*tc*h*s[,,t-1]^q
    out <- (1+ tm1*delta)*k*tc*s[,,t-1]^n 
    Es  <- inp-out
    allocs[,,t] <- pmax(allocs[,,t-1],t(apply(lw(s[,,t-1]),1,inv_logit3,lm,ts[t],slope,tr)))
    
    mt <- (tm1*v$v+v$M)*s[,,t-1]^v$nu
    surv[,,t] <- surv[,,t-1] + dts*mt
    R0[,,t] <- R0[,,t-1] + dts*allocs[,,t]*Es*exp(-surv[,,t])
    
    s[,,t] <- s[,,t-1]+dts*(1-allocs[,,t])*Es
  }
  ls <- lw(s)
  opt <- apply(R0[,,dt],1,function(x) ifelse(any(!is.nan(x)),which.max(x),NA))
  
  #browser()
  
  s=t(sapply(1:temps,function(x) s[x,opt[x],]))
  allocs=t(sapply(1:temps,function(x) allocs[x,opt[x],]))
  R0s=t(sapply(1:temps,function(x) R0[x,opt[x],]))
  
  lss <- reshape2::melt(s)
  colnames(lss) <- c('Temperature','t','size')
  lss$Temperature <- temp[lss$Temperature]
  lss$t <- ts[lss$t]
  lss$opt <- opt[lss$Temperature]
  
  alloc <- reshape2::melt(allocs)
  colnames(alloc) <- c('Temperature','t','allocs')
  alloc$Temperature <- temp[alloc$Temperature]
  alloc$t <- ts[alloc$t]
  alloc$opt <- opt[alloc$Temperature]
  
  R0s <- reshape2::melt(R0s)
  colnames(R0s) <- c('Temperature','t','R0')
  R0s$Temperature <- temp[R0s$Temperature]
  R0s$opt <- opt[R0s$Temperature]
  R0s$t <- ts[R0s$t]
  
  #browser()
  
  growth <- inner_join(inner_join(lss,alloc),R0s) %>% arrange(t,Temperature)
  
  
  winfs <- growth %>% group_by(Temperature) %>%
    summarise(l=size[ifelse(any(abs(allocs-0.5)<0.05),which.min(abs(allocs-0.5)),NA)-1],
              t=t[ifelse(any(abs(allocs-0.5)<0.05),which.min(abs(allocs-0.5)),NA)-1],
              opt=unique(opt))
  
  G <- rep(NA,length(unique(winfs$Temperature)))
  
  for(t in 2:(length(winfs$Temperature)-1)) {
    tau = winfs$Temperature[t]
    this.l <- winfs$l[t]
    if (is.na(this.l)) next
    tPM <- R0[t,opt[t],dt]
    lPM <- R0[t-1,opt[t],dt]
    nPM <- R0[t+1,opt[t],dt]
    
    sl <- abs((nPM-lPM)/(winfs$Temperature[t+1]-winfs$Temperature[t-1]))
    
    G[t] <- 0.004*this.l*sl/tPM
  }       
  
  
  sG <- sign(G)
  winfs$G <- G/max(winfs$t,1)
  winfs$L <- 10*(lw(winfs$l + winfs$G)-lw(winfs$l))
  
  
  
  list(winfs=winfs, growth=growth)
  
}


model_out_growth_check <- function(temp,
                                   temp_ref=15,
                                   l,
                                   lm,
                                   Ea,
                                   gamma=50,
                                   delta=2,
                                   phi=0.15,
                                   h=30,
                                   beta=0.2,
                                   k=2,
                                   p=0.8,
                                   q=0.9,
                                   n=0.8,
                                   mstar=1000,
                                   tmax=10,
                                   slope=0.05,
                                   tr = 1,
                                   v=NULL,
                                   dt=100){
  
  
  temps = length(temp)
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  ts <- seq(0,tmax,l=dt)
  
  s <- array(0, c(temps,temps,dt))
  allocs <- array(0, c(temps,temps,dt))
  
  s[,,1] <- min(wl(lm))
  
  dts <- (tmax/(dt-1))
  
  for(t in 2:dt) {
    tm1 <- sapply(1:temps,function(g) {
      w <- v
      w$gamma <- gamma[g]
      get_taus(w,1,10,temp,s[,g,t-1])
    })
    f <- tm1*gamma*s[,,t-1]^p/(tm1*gamma*s[,,t-1]^p+tc*h*s[,,t-1]^q)
    inp <- (1-phi-beta)*f*tc*h*s[,,t-1]^q
    out <- (1+ tm1*delta)*k*tc*s[,,t-1]^n 
    Es=inp-out
    allocs[,,t] <- pmax(allocs[,,t-1],t(apply(lw(s[,,t-1]),1,inv_logit3,mstar,ts[t],slope,tr)))
    
    s[,,t] <- s[,,t-1]+dts*(1-allocs[,,t])*Es
  }
  ls <- lw(s)
  
  #browser()
  
  ref = which.min(abs(temp-temp_ref))
  
  gls=t(sapply(1:temps,function(x) ls[ref,x,]))
  gallocs=t(sapply(1:temps,function(x) allocs[ref,x,]))
  
  tls=t(sapply(1:temps,function(x) ls[x,ref,]))
  tallocs=t(sapply(1:temps,function(x) allocs[x,ref,]))
  
  lss <- reshape2::melt(tls)
  colnames(lss) <- c('Temperature','t','t_length')
  lss$Temperature <- temp[lss$Temperature]
  lss$Gamma <- gamma[ref]
  lss$t <- ts[lss$t]
  
  alloc <- reshape2::melt(tallocs)
  colnames(alloc) <- c('Temperature','t','allocs')
  alloc$Temperature <- temp[alloc$Temperature]
  alloc$Gamma <- gamma[ref]
  alloc$t <- ts[alloc$t]
  
  
  lgs <- reshape2::melt(gls)
  colnames(lgs) <- c('Gamma','t','g_length')
  lgs$Temperature <- temp[ref]
  lgs$Gamma <- gamma[lgs$Gamma]
  lgs$t <- ts[lgs$t]
  
  galloc <- reshape2::melt(gallocs)
  colnames(galloc) <- c('Gamma','t','allocs')
  galloc$Temperature <- temp[ref]
  galloc$Gamma <- gamma[galloc$Gamma]
  galloc$t <- ts[galloc$t]
  
  #browser()
  
  
  list(t_growth = inner_join(lss,alloc) %>% arrange(t,Temperature),
       g_growth = inner_join(lgs,galloc) %>% arrange(t,Temperature))
  
}

O2_supply <- function(O2 = 1:100,O2crit=20,P50 = 40, Tmax=30,Topt=15,T,omega=1.870,delta=1038,n=0.8){
  
  level <- delta*((Tmax-T)/(Tmax-Topt))^omega*exp(-omega*(Tmax-T)/(Tmax-Topt))/exp(-omega)
  365*24*level*(1-exp(-(O2-O2crit)/(-(P50-O2crit)/log(0.5))))/1000
  
}

# plot(O2_supply(level=200,P50 = 10),t='l',xlab='Disolved O2',ylab='02 supply')

eval_tau_max_o2 <- function(f=O2_supply(),
                            omega = 0.4,
                            gamma=50,
                            delta=2,
                            h=30,
                            phi=0.15,
                            beta=0.25,
                            k=2,
                            p=0.8,
                            q=0.9,
                            n=0.8,
                            m=100){
  
  
  tau_max<-(m^(-n-p)*(sqrt((gamma*(-f)*m^(n+p)+gamma*k*omega*m^(n+p)+delta*k*h*omega*m^(n+q)+beta*gamma*h*omega*m^(p+q))^2-4*gamma*delta*k*omega*m^(n+p)*(k*h*omega*m^(n+q)-f*h*m^(n+q)))+gamma*f*m^(n+p)-gamma*k*omega*m^(n+p)-delta*k*h*omega*m^(n+q)-beta*gamma*h*omega*m^(p+q)))/(2*gamma*delta*k*omega)
  
  tau_max
  
}

O2_fact <- function(temp,Tref=15){
  
  exp(-0.01851*(temp-Tref))
  
}

eval_tau_max_temp <- function(f=O2_supply(),
                              Ea = 0.52,
                              temp=seq(5,10,l=100),
                              temp_ref=15,
                              omega = 0.4,
                              gamma=50,
                              delta=2,
                              h=30,
                              phi=0.15,
                              beta=0.25,
                              k=2,
                              p=0.8,
                              q=0.9,
                              n=0.8,
                              m=100){
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  tau_max<-(m^(-n-p)*(sqrt((gamma*(-f)*m^(n+p)+gamma*tc*k*omega*m^(n+p)+delta*k*tc^2*h*omega*m^(n+q)+beta*tc*gamma*h*omega*m^(p+q))^2-4*gamma*delta*tc*k*omega*m^(n+p)*(k*h*tc^2*omega*m^(n+q)-f*h*tc*m^(n+q)))+gamma*f*m^(n+p)-gamma*k*omega*tc*m^(n+p)-delta*k*h*tc^2*omega*m^(n+q)-beta*gamma*h*tc*omega*m^(p+q)))/(2*gamma*delta*tc*k*omega)
  
  tau_max
  
}

get_taus <- function(v,tau_uc,O2_in,temp_in,m=10^seq(0,6,l=1000)){
  #browser()
  
  O2_tcor <- O2_fact(temp_in,5)
  O2 = O2_supply(O2=10*O2_tcor,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=temp_in,delta=v$lO,omega=v$shape,P50=v$P50,n=v$n)
  max_tau <- eval_tau_max_temp(f=O2,
                               temp=temp_in,
                               Ea=v$Ea,
                               omega = v$omega,
                               gamma=v$gamma,
                               delta=v$delta,
                               phi=v$phi,
                               h=v$h,
                               beta=v$beta,
                               k=v$k,
                               p=v$p,
                               q=v$q,
                               n=v$n,
                               m=m)
  tau <- eval_tau_eq_temp(temp=temp_in,
                          Ea=v$Ea,
                          gamma=v$gamma,
                          delta=v$delta,
                          phi=v$phi,
                          h=v$h,
                          beta=v$beta,
                          k=v$k,
                          p=v$p,
                          q=v$q,
                          n=v$n,
                          m=m,
                          M=v$M,
                          v=v$v)
  
  tau_max = pmin(tau,max_tau)
  tau_max[tau_max<0] <- 0
  tau_max[tau_max>1] <- 1
  
  tau_max
  
}
