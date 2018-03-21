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
  f=O2_supply(10*O2_tcor,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=v$temp,delta=v$lO,omega=v$shape,P50=v$P50)
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
                               o=v$o,
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
                          o=v$o,
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
                           o=v$o,
                           m=v$m)
  
  
  #browser()
  scope <- f*v$m-model_frame$std*v$omega
  scope[scope<0.0001] <- 0
  #browser()
  
  
  #browser()
  bind_cols(model_frame,
            data_frame(
              Temperature=v$temp,
              Realised = tau_max,
              Limit = sapply(sapply(max_tau,max,0),min,1),
              M = (tau_max*v$v+v$M)*v$m^v$nu,
              Optimum = sapply(sapply(tau,max,0),min,1),
              Scope=ifelse(scope<0.0001,0,scope),
              `Max` = f*v$m,
              `Active` = model_frame$Met*v$omega,
              `Std` = model_frame$std*v$omega,
              Viable = as.numeric(tau_max>0.0001 & model_frame[['C for growth']]>0.0001))
  )
  
}

plot_data_growth_tO2 <- function(v,mstar=NULL){
  
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
                           o=v$o,
                           tmax = v$tmax,
                           slope=v$slope,
                           tr=v$tr,
                           v=v,
                           lm=lw(lm),
                           dt = v$dt)
  
  winfs <- mout$winfs
  
  ### Step 2 calc mat for selected response 
  
  #browser()
  
  mstar = ifelse(!is.null(mstar),mstar, lw(lm[winfs$opt[which.min(abs(winfs$Temperature-v$temp_ref))]]))
  
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
                                  o=v$o,
                                  mstar = mstar,
                                  tmax = v$tmax,
                                  slope=v$slope,
                                  tr=v$tr,
                                  v=v,
                                  lm=lw(lm),
                                  dt = v$dt)
  
  list(winfs=winfs,
       growth=mout$growth,
       R0=mout$R0s,
       mstar=mstar,
       g_growth=mouts$g_growth,
       t_growth=mouts$t_growth)
  
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
                             o=1,
                             m=100,
                             M=0.2,
                             v=1){
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  -(M*delta*h*k*m^(o + p + q)*tc - h*k*m^(n + p + q)*tc*v - sqrt(-((beta - 1)*h*k*m^(n + 2*p + 3*q) + h*k*m^(n + 2*p + 3*q)*phi)*tc^2*v^2 - ((beta - 1)*delta*gamma*k*m^(o + 3*p + 2*q) + delta*gamma*k*m^(o + 3*p + 2*q)*phi)*M^2*tc + (((beta - 1)*delta*h*k*m^(o + 2*p + 3*q) + delta*h*k*m^(o + 2*p + 3*q)*phi)*tc^2 + (gamma*h*m^(3*p + 3*q)*phi^2 + (beta - 1)*gamma*k*m^(n + 3*p + 2*q) + (beta^2 - 2*beta + 1)*gamma*h*m^(3*p + 3*q) + (2*(beta - 1)*gamma*h*m^(3*p + 3*q) + gamma*k*m^(n + 3*p + 2*q))*phi)*tc)*M*v)*h)/(M*delta*gamma*k*m^(o + 2*p) - ((beta - 1)*gamma*h*m^(2*p + q) + gamma*h*m^(2*p + q)*phi + gamma*k*m^(n + 2*p))*v)
  
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
                              o=1,
                              m=100){
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  -1/2*(delta*h*k*m^(o + q)*omega*tc^2 - f*gamma*m^(p + 1) + (beta*gamma*h*m^(p + q) + gamma*k*m^(n + p))*omega*tc - sqrt(delta^2*h^2*k^2*m^(2*o + 2*q)*omega^2*tc^4 + 2*(beta*delta*gamma*h^2*k*m^(o + p + 2*q) - delta*gamma*h*k^2*m^(n + o + p + q))*omega^2*tc^3 + f^2*gamma^2*m^(2*p + 2) - 2*(beta*f*gamma^2*h*m^(2*p + q + 1) + f*gamma^2*k*m^(n + 2*p + 1))*omega*tc + (2*delta*f*gamma*h*k*m^(o + p + q + 1)*omega + (beta^2*gamma^2*h^2*m^(2*p + 2*q) + 2*beta*gamma^2*h*k*m^(n + 2*p + q) + gamma^2*k^2*m^(2*n + 2*p))*omega^2)*tc^2))*m^(-o - p)/(delta*gamma*k*omega*tc)

  }

get_taus <- function(v,tau_uc,O2_in,temp_in,m=10^seq(0,6,l=1000)){
  #browser()
  
  O2_tcor <- O2_fact(temp_in,5)
  O2 = O2_supply(O2=10*O2_tcor,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=temp_in,delta=v$lO,omega=v$shape,P50=v$P50)
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
                               o=v$o,
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
                          o=v$o,
                          m=m,
                          M=v$M,
                          v=v$v)
  
  tau_max = pmin(tau,max_tau)
  tau_max[tau_max<0] <- 0
  tau_max[tau_max>1] <- 1
  
  tau_max
  
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
                      o=v$o,
                      m=100){
  
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  f <- tau_max*gamma*m^p/(tau_max*gamma*m^p+tc*h*m^q)
  inp <- (1-phi-beta)*f*tc*h*m^q
  out <- k*tc*m^n + tau_max*delta*k*tc*m^o
  e <-  inp -out
  efficiency <- e/f
  efficiency[efficiency<0] <- 0
  predation_rate <- f/phi
  
  met = beta*f*tc*h*m^q+ out
  #browser()
  
  data_frame(`Feeding level`= f, 
             Consumption = inp,
             `C used for Metabolism` = out,
             `C for growth` = e, 
             Efficiency = efficiency, 
             `Predation rate`=predation_rate, 
             Met = met,
             std = k*tc*m^n)
}


get_model_out <- function(v,tau){
  pars <- model_out(tau,
            temp = 15,
            temp_ref=15,
            Ea = v$Ea,
            r=v$r,
            gamma=v$gamma,
            delta=v$delta,
            phi=v$phi,
            h=v$h,
            beta=v$beta,
            k=v$k,
            p=v$p,
            q=v$q,
            n=v$n,
            o=v$o,
            m=10)
  pars$M <- (tau*v$v+v$M)*v$m^v$nu
  pars$tau <- tau
  pars
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
                             o=1,
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
  mu <- array(0, c(temps,l,dt))
  
  dts <- (tmax/(dt-1))
  s[,,1] <- min(wl(lm))
  
  
  for(t in 2:dt) {
    tm1 <- get_taus(v,1,10,temp,s[,,t-1])
    Es <- (1-phi-beta)*(tm1*gamma*s[,,t-1]^p/(tm1*gamma*s[,,t-1]^p+tc*h*s[,,t-1]^q))*tc*h*s[,,t-1]^q -k*tc*s[,,t-1]^n-tm1*delta*k*tc*s[,,t-1]^o
    allocs[,,t] <- pmax(allocs[,,t-1],t(apply(lw(s[,,t-1]),1,inv_logit3,lm,ts[t],slope,tr)))
    
    
    s[,,t] <- s[,,t-1]+dts*(1-allocs[,,t])*Es
    mu[,,t] <- mu[,,t-1] + dts*(tm1*v$v+v$M)*s[,,t]^v$nu
    R0[,,t] <- R0[,,t-1] + dts*allocs[,,t]*Es*exp(-mu[,,t])
    
  }
  ls <- lw(s)
  opt <- apply(R0[,,dt],1,function(x) ifelse(any(!is.nan(x)),which.max(x),NA))
  
  #browser()
  
  s=t(sapply(1:temps,function(x) s[x,opt[x],]))
  allocs=t(sapply(1:temps,function(x) allocs[x,opt[x],]))
  R0s=t(sapply(1:temps,function(x) R0[x,opt[x],]))
  
  lss <- reshape2::melt(s)
  colnames(lss) <- c('Temperature','t','size')
  lss$opt <- opt[lss$Temperature]
  lss$Temperature <- temp[lss$Temperature]
  lss$t <- ts[lss$t]
  
  alloc <- reshape2::melt(allocs)
  colnames(alloc) <- c('Temperature','t','allocs')
  alloc$t <- ts[alloc$t]
  alloc$opt <- opt[alloc$Temperature]
  alloc$Temperature <- temp[alloc$Temperature]
  
  R0s <- reshape2::melt(R0s)
  colnames(R0s) <- c('Temperature','t','R0')
  R0s$opt <- opt[R0s$Temperature]
  R0s$Temperature <- temp[R0s$Temperature]
  R0s$t <- ts[R0s$t]
  
  #browser()
  
  growth <- inner_join(inner_join(lss,alloc),R0s) %>% arrange(t,Temperature)
  
  
  winfs <- growth %>% group_by(Temperature) %>%
    summarise(l=size[ifelse(any(abs(allocs-0.5)<0.05),which.min(abs(allocs-0.5)),NA)-1],
              t=t[ifelse(any(abs(allocs-0.5)<0.05),which.min(abs(allocs-0.5)),NA)-1],
              opt=unique(opt))
  
  G <- rep(NA,length(unique(winfs$Temperature)))
  #browser()
  for(t in 2:(length(winfs$Temperature)-1)) {
    tau = winfs$Temperature[t]
    this.l <- winfs$l[t]
    if (is.na(this.l)) next
    tPM <- R0[t,opt[t],dt]
    lPM <- R0[t-1,opt[t],dt]
    nPM <- R0[t+1,opt[t],dt]
    
    sl <- abs((nPM-lPM)/(winfs$Temperature[t+1]-winfs$Temperature[t-1]))
    
    G[t] <- 0.04*this.l*(sl/tPM)
  }       
  
  
  sG <- sign(G)
  winfs$G <- G/winfs$t
  winfs$L <- 10*(lw(winfs$l + winfs$G)-lw(winfs$l))
  
  
  
  list(winfs=winfs, growth=growth, R0s  = R0[,which.min(abs(winfs$Temperature-v$temp_ref)),dt])
  
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
                                   o=1,
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
    Es <- (1-phi-beta)*(tm1*gamma*s[,,t-1]^p/(tm1*gamma*s[,,t-1]^p+tc*h*s[,,t-1]^q))*tc*h*s[,,t-1]^q -k*tc*s[,,t-1]^n-tm1*delta*k*tc*s[,,t-1]^o
   
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

O2_supply <- function(O2 = 1:100,O2crit=20,P50 = 40, Tmax=30,Topt=15,T,omega=1.870,delta=1038){
  
  level <- delta*((Tmax-T)/(Tmax-Topt))^omega*exp(-omega*(Tmax-T)/(Tmax-Topt))/exp(-omega)
  365*24*level*(1-exp(-(O2-O2crit)/(-(P50-O2crit)/log(0.5))))/1000
  
}

O2_fact <- function(temp,Tref=15){
  
  exp(-0.01851*(temp-Tref))
  
}

