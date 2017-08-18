inv_logit3 <- function(m,mstar,a,b,c) {
  z = (m-mstar)*cos(atan(b))-a*sin(atan(b))
  #zstar = mstar*cos(atan(b))-a*sin(atan(b))
  #1/(1+ (z/mstar)^-c)
  1/(1+exp(-(c*z)))
}

lw <- function(w) (w/0.01)^(1/3)


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
  
  lm = 10^seq(-2,5,l=v$lm)
  
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

plot_data_growth <- function(v){
  
  tau_uc = eval_tau_eq(gamma=v$gamma,
                       delta=v$delta,
                       phi=v$phi,
                       h=v$h,
                       beta=v$beta,
                       k=v$k,
                       p=v$p,
                       q=v$q,
                       n=v$n,
                       m=10^seq(1,5,l=1000),
                       M=v$M,
                       v=v$v)
  
  tau_uc[tau_uc<0] <- 0
  tau_uc[tau_uc>1] <- 1
  
  taus <- get_taus(v,tau_uc,10,v$temp_ref)
  
  model_out_par(tau_max = taus$tau_max,
                tau = taus$tau,
                tau_uc = tau_uc,
                tau_o2 = taus$tau_o2,
                M=v$M,
                temp=v$temp_ref,
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
                m = 10^seq(1,5,l=1000)
  )
}

plot_data_growth_tO2 <- function(v,lm = 10^seq(-2,6,l=1000),n_int=100){
  
  #browser()
  O2_range <- seq(1,10,l=n_int)
  winfs <- data.frame(matrix(NA,n_int,3))
  winfs[,1] <- O2_range
  winfs[,2] <- seq(min(v$temp),max(v$temp),l=n_int)
  
  taus_ref <- get_taus(v,1,10,v$temp_ref,m=lm)
  
  mout_ref <- model_out_par(tau_max = taus_ref$tau_max,
                            M=v$M,
                            v=v$v,
                            nu=v$nu,
                            temp=v$temp_ref,
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
                            n=ifelse(v$n==v$q,v$n+0.1,v$n),
                            m = lm)
  
  mouts <- parallel::mclapply(1:n_int, function(i){#parallel::mclapply(1:n_int, function(i){
    taus <- get_taus(v,1,winfs[i,1],winfs[i,2],m=lm)
    
    mout <- model_out_par(tau_max = taus$tau_max,
                          M=v$M,
                          v=v$v,
                          nu=v$nu,
                          temp=winfs[i,2],
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
                          m = lm)
    #browser()
    mout$gout <- model_out_growth(temp=winfs[i,2],
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
                                  mstar = mout_ref$winf,
                                  tmax = v$tmax,
                                  slope=v$slope,
                                  tr=v$tr,
                                  v=v)
    
    # Increment the progress bar, and update the detail text.
    
    mout
    
  },mc.cores=4)
  
  
  winfs[,3] <- do.call(c,lapply(mouts, function(x) x[['winf']]))
  
  colnames(winfs) <- c('O2','Temperature','Temp')
  
  dPM <- bind_rows(lapply(mouts, function(x) x[['dPm']]))
  
  alloc <- bind_rows(lapply(mouts, function(x) x[['gout']]))
  
  aa <- alloc %>% 
    group_by(Temperature) %>% 
    summarise(ns=sum(is.na(ls)),
              na = sum(is.na(allocs)),
              ni = sum(!is.finite(allocs))) %>% 
    filter(na>1 | ns >0 | ni > 1)
  
  alloc <- alloc %>% 
    filter(!Temperature %in% aa$Temperature)
  dPM <- dPM %>% 
    filter(!Temperature %in% aa$Temperature)
  winfs <- winfs %>% 
    filter(!Temperature %in% aa$Temperature)
  
  dpj <- dPM %>% inner_join(winfs)
  G <- rep(NA,length(winfs$Temperature))
  for(t in 2:(length(winfs$Temperature)-1)) {
    tau = winfs$Temperature[t]
    this.mass <- unique(dpj$Temp[dpj$Temperature == tau])
    if (is.na(this.mass)) next
    tPM <- dpj$Pm[dpj$Temperature == tau & dpj$Mass == this.mass]
    lPM <- dpj$Pm[dpj$Temperature == winfs$Temperature[t-1] & dpj$Mass == this.mass]
    nPM <- dpj$Pm[dpj$Temperature == winfs$Temperature[t+1] & dpj$Mass == this.mass]
    
    sl <- (nPM-lPM)/(winfs$Temperature[t+1]-winfs$Temperature[t-1])
    
    G[t] <- 0.04*this.mass*sl/tPM
  }       
  
  #browser()
  norm <- alloc %>%
    group_by(Temperature) %>% 
    summarise(ts=t[ifelse(any(abs(allocs-0.5)<0.1),which.min(abs(allocs-0.5)),NA)-1],
              m=ls[ifelse(any(abs(allocs-0.5)<0.1),which.min(abs(allocs-0.5)),NA)-1]) 
  
  sG <- sign(G)
  winfs$G <- G/norm$ts
  winfs$L <- 10*(lw(winfs$Temp + winfs$G)-lw(winfs$Temp))
  list(winfs=winfs, dPM=dPM,alloc=alloc)
  
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

get_dPdm <- function(m=NULL,tau,
                     tc,
                     nu=-0.25,
                     v=0.2,
                     M=0.2,
                     gamma=50,
                     delta=2,
                     phi=0.15,
                     h=30,
                     beta=0.2,
                     k=2,
                     p=0.8,
                     q=0.9,
                     n=0.8) ((beta + phi - 1)*gamma*h*m^p*m^q*tau*tc/(gamma*m^p*tau + h*m^q*tc) + (delta*tau + 1)*k*m^n*tc)*m^(nu - 1)*nu/((tau*v + M)*m^(2*nu)) - ((beta + phi - 1)*gamma*h*m^(p - 1)*m^q*p*tau*tc/(gamma*m^p*tau + h*m^q*tc) + (beta + phi - 1)*gamma*h*m^p*m^(q - 1)*q*tau*tc/(gamma*m^p*tau + h*m^q*tc) - (gamma*m^(p - 1)*p*tau + h*m^(q - 1)*q*tc)*(beta + phi - 1)*gamma*h*m^p*m^q*tau*tc/(gamma*m^p*tau + h*m^q*tc)^2 + (delta*tau + 1)*k*m^(n - 1)*n*tc)/((tau*v + M)*m^nu)

model_out_par <- function(tau_max,
                          temp,
                          temp_ref=15,
                          Ea,
                          M=0.2,
                          v=1,
                          nu=-0.25,
                          gamma=50,
                          delta=2,
                          phi=0.15,
                          h=30,
                          beta=0.2,
                          k=2,
                          p=0.8,
                          q=0.9,
                          n=0.8,
                          m=10^seq(0,6,l=1000)){
  
  #browser()
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  
  f <- tau_max*gamma*m^p/(tau_max*gamma*m^p+tc*h*m^q)
  inp <- (1-phi-beta)*f*tc*h*m^q
  out <- (1+tau_max*delta)*k*tc*m^n 
  E=inp/out
  
  #browser()
  mm <- get_dPdm(m=m,tau=tau_max,tc=tc,v=v,M=M,nu=nu,gamma=gamma,delta=delta,
                 phi=phi,h=h,beta=beta,k=k,p=p,q=q,n=n)
  mtemp <- (tau_max*v+M)*m^nu
  
  test <- abs(mm-1)
  cond <- order(test[test<1])[1]
  winf <- ifelse(!is.na(cond), m[test<1][cond],NA)
  
  list(winf=winf, 
       dPm = data_frame(Temperature=temp,
                        Mass=m,
                        dPm=mm,
                        Pm=(inp-out)/mtemp))
}

model_out_growth <- function(temp,
                             temp_ref=15,
                             l,
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
                             v=NULL){
  
  #browser()
  tc = exp(Ea*((temp+(288.2-temp_ref))-288.2)/(8.6173324*10^(-5)*(temp+(288.2-temp_ref))*288.2))
  
  ts <- seq(0,tmax,l=l)
  s <- vector()
  allocs <- vector()
  s[1] <- 0.01
  
  
  for(t in 2:l) {
    tm1 <- get_taus(v,1,10,temp,s[t-1])$tau_max
    f <- tm1*gamma*s[t-1]^p/(tm1*gamma*s[t-1]^p+tc*h*s[t-1]^q)
    inp <- (1-phi-beta)*f*tc*h*s[t-1]^q
    out <- (1+ tm1*delta)*k*tc*s[t-1]^n 
    Es=inp-out
    allocs[t] = min(max(allocs[t-1],inv_logit3(lw(s[t-1]),lw(mstar),ts[t],slope,tr),na.rm=T),1)
    
    s[t] = s[t-1]+(tmax/l)*(1-allocs[t])*Es
  }
  ls = lw(s)
  data_frame(Temperature=temp,
             Growth=gamma,
             ls=ls,
             allocs=allocs,
             t=ts)
  
}


model_out_growth_check <- function(temp,
                                   l,
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
                                   nu,
                                   M,
                                   v,
                                   tmax=10,
                                   dat=NULL){
  
  #browser()
  tc = 1
  
  ts <- seq(1,tmax,l=l)
  s <- vector()
  allocs <- vector()
  dpMs <- vector()
  dpMs[1] <- Inf
  
  s[1] <- 0.01
  for(t in 2:l) {
    tm1 <- get_taus(dat,1,10,temp,s[t-1])$tau_max
    f <- tm1*gamma*s[t-1]^p/(tm1*gamma*s[t-1]^p+tc*h*s[t-1]^q)
    inp <- (1-phi-beta)*f*tc*h*s[t-1]^q
    out <- (1+ tm1*delta)*k*tc*s[t-1]^n 
    Es=inp-out
    dpMs[t] <- get_dPdm(m=s[t-1],tau=tm1,tc=tc,v=v,M=M,nu=nu,gamma=gamma,delta=delta,
                        phi=phi,h=h,beta=beta,k=k,p=p,q=q,n=n)
    allocs[t] = ifelse(abs(dpMs[t]-1)>=abs(dpMs[t-1]-1),1,0)
    
    s[t] = s[t-1]+(tmax/l)*(1-allocs[t])*Es
  }
  ls = lw(s)
  #browser()
  data_frame(Growth=gamma,
             ls=ls,
             allocs=allocs,
             t=ts)
  
  #  (Eq*h*m^(q - 1)*nu - Eq*h*m^(q - 1)*q - (h*m^(q - 1)*nu - h*m^(q - 1)*q)*E - (delta*k*m^(n - 1)*n - delta*k*m^(n - 1)*nu)*tau + (delta*k*m^(n - 1)*n - delta*k*m^(n - 1)*nu)*tauq)*M/(k*m^(n - 1)*n - k*m^(n - 1)*nu + (h*m^(q - 1)*nu - h*m^(q - 1)*q)*E + (delta*k*m^(n - 1)*n - delta*k*m^(n - 1)*nu)*tau)
  
  #  -(Eq*h*m^(q - 1)*nu*phi + (Eq*beta - Eq)*h*m^(q - 1)*nu - ((beta - 1)*h*m^(q - 1)*nu + h*m^(q - 1)*nu*phi - ((beta - 1)*h*m^(q - 1) + h*m^(q - 1)*phi)*q)*E - (Eq*h*m^(q - 1)*phi + (Eq*beta - Eq)*h*m^(q - 1))*q + (delta*k*m^(n - 1)*n - delta*k*m^(n - 1)*nu)*tau - (delta*k*m^(n - 1)*n - delta*k*m^(n - 1)*nu)*tauq)*M/(k*m^(n - 1)*n - k*m^(n - 1)*nu - ((beta - 1)*h*m^(q - 1)*nu + h*m^(q - 1)*nu*phi - ((beta - 1)*h*m^(q - 1) + h*m^(q - 1)*phi)*q)*E + (delta*k*m^(n - 1)*n - delta*k*m^(n - 1)*nu)*tau)
  
  #  -(((beta - 1)*h*m^(q - 1)*nu + h*m^(q - 1)*nu*phi - ((beta - 1)*h*m^(q - 1) + h*m^(q - 1)*phi)*q)*E*tauq*v - (Eq*h*m^(q - 1)*nu*phi + (Eq*beta - Eq)*h*m^(q - 1)*nu - ((beta - 1)*h*m^(q - 1)*nu + h*m^(q - 1)*nu*phi - ((beta - 1)*h*m^(q - 1) + h*m^(q - 1)*phi)*q)*E - (Eq*h*m^(q - 1)*phi + (Eq*beta - Eq)*h*m^(q - 1))*q + (delta*k*m^(n - 1)*n - delta*k*m^(n - 1)*nu)*tau - (delta*k*m^(n - 1)*n - delta*k*m^(n - 1)*nu)*tauq)*M - ((Eq*h*m^(q - 1)*nu*phi - k*m^(n - 1)*n + ((Eq*beta - Eq)*h*m^(q - 1) + k*m^(n - 1))*nu - (Eq*h*m^(q - 1)*phi + (Eq*beta - Eq)*h*m^(q - 1))*q)*tau + (k*m^(n - 1)*n - k*m^(n - 1)*nu)*tauq)*v)/(((beta - 1)*h*m^(q - 1)*nu + h*m^(q - 1)*nu*phi - ((beta - 1)*h*m^(q - 1) + h*m^(q - 1)*phi)*q)*E*tauq - (k*m^(n - 1)*n - k*m^(n - 1)*nu + (delta*k*m^(n - 1)*n - delta*k*m^(n - 1)*nu)*tau)*tauq)
  
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
  O2 = O2_supply(O2=O2_in,Topt=v$Topt,O2crit=v$O2crit,Tmax=v$temp[length(v$temp)],T=v$Topt,delta=v$lO,omega=v$shape,P50=v$P50,n=v$n)
  tau_o2 = eval_tau_max_o2(f=O2,omega = v$omega,
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
  
  tau_o2 <- apply(cbind(tau_uc,tau_o2),1,min)
  
  tau_o2[tau_o2<0] <- 0
  tau_o2[tau_o2>1] <- 1
  
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
  
  tau_max = apply(cbind(tau,max_tau),1,min)
  tau_max[tau_max<0] <- 0
  tau_max[tau_max>1] <- 1
  
  list(tau_max=tau_max, tau_o2=tau_o2, tau=tau)
  
}
