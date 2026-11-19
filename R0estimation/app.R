library(shiny)
library(SimInf)
library(RColorBrewer)
library(EpiEstim)

# Estimating R0: simulate an SEIR outbreak and estimate R0 from it by the exponential
# growth rate (linear regression on log cases), the final size and EpiEstim, for one
# epidemic or the distribution of estimates over many. The simulation code follows the
# SEIRmodel app (SimInf for the stochastic model, a base-R rk4 for the deterministic
# one) so it also runs in the browser via shinylive

# Classic 4th-order Runge-Kutta with a fixed internal step h. Same interface and
# output layout (time, then one column per variable) as deSolve::lsoda, but base R
# only, so the app also runs in the browser via shinylive
rk4 <- function(y, times, func, parms, h = 0.02)
{
    deriv <- function(t, y) func(t, y, parms)[[1]]
    out <- matrix(NA, length(times), length(y) + 1)
    out[1, ] <- c(times[1], y)
    for(k in seq_along(times)[-1])
    {
        n  <- ceiling((times[k] - times[k-1])/h)
        dt <- (times[k] - times[k-1])/n
        t  <- times[k-1]
        for(j in 1:n)
        {
            k1 <- deriv(t, y)
            k2 <- deriv(t + dt/2, y + dt/2*k1)
            k3 <- deriv(t + dt/2, y + dt/2*k2)
            k4 <- deriv(t + dt, y + dt*k3)
            y  <- y + dt/6*(k1 + 2*k2 + 2*k3 + k4)
            t  <- t + dt
        }
        out[k, ] <- c(times[k], y)
    }
    out
}

mass_actionD <- function(beta,a,g,S0,E0,I0,R0,time_vec)
{
    SI.dyn <- function(t,var,par) {
        S <- var[1]
        E <- var[2]
        I <- var[3]
        R <- var[4]

        N <- S + E + I + R

        beta <- par[1]
        a    <- par[2]
        g    <- par[3]

        dS <- - beta*S*I/N
        dE <-   beta*S*I/N - a*E
        dI <-   a*E - g*I
        dR <-   g*I
        return(list(c(dS,dE,dI,dR)))
    }

    SI.par  <- c(beta,a,g)
    SI.init <- c(S0,E0,I0,R0)

    oot<-rk4(SI.init,time_vec,SI.dyn,SI.par)

    return(oot)
}

# Run n stochastic replicates of SimInf's built-in SEIR model, which has the transitions
#   S -> beta*S*I/(S+E+I+R) -> E
#   E -> epsilon*E -> I
#   I -> gamma*I -> R
# It is precompiled, so unlike mparse() it needs no C compiler (and runs in the browser).
# SimInf draws its random numbers from R's generator, so set.seed() makes runs repeatable.
# Returns a list of time x replicate matrices S, E, I and R, starting at t=0
sim_seir <- function(n, beta, te, ti, init, t_run)
{
    initial_state <- data.frame(S=rep(init[1],n),E=rep(init[2],n),
                                I=rep(init[3],n),R=rep(init[4],n))
    # SimInf's first recorded point is already one time unit into the run, so record
    # from t=1 and prepend the initial state at t=0 ourselves
    model <- SEIR(u0 = initial_state, tspan = seq(from = 1, to = t_run, by = 1),
                  beta = as.double(beta), epsilon = 1.0/te, gamma = 1.0/ti)
    traj <- trajectory(run(model = model))
    # Sort by replicate then time, so each compartment reshapes to one column per replicate
    traj <- traj[order(traj$node, traj$time),]
    lapply(c(S='S',E='E',I='I',R='R'), function(v)
        rbind(rep(init[match(v,c('S','E','I','R'))],n), matrix(traj[[v]], ncol=n)))
}

# As sim_seir, but if major_only is set, re-run any replicate whose outbreak infects
# fewer than thresh*N individuals (a minor outbreak that fades out early), up to
# max_tries times. Returns the matrices plus how many replicates never took off
sim_major <- function(n, beta, te, ti, init, t_run, major_only, thresh, max_tries = 50)
{
    out <- sim_seir(n, beta, te, ti, init, t_run)
    N   <- sum(init)
    is_minor <- function(sims) sum(init[1:3]) - sims$S[nrow(sims$S),] < thresh*N
    redo  <- if(major_only) which(is_minor(out)) else integer(0)
    tries <- 1
    while(length(redo) > 0 && tries < max_tries)
    {
        new <- sim_seir(length(redo), beta, te, ti, init, t_run)
        for(v in names(out)) out[[v]][,redo] <- new[[v]]
        redo  <- redo[is_minor(new)]
        tries <- tries + 1
    }
    out$n_minor <- if(major_only) length(redo) else sum(is_minor(out))
    out
}

# The three data streams we can estimate R0 from, worked out from the compartments
# (columns are replicates). New cases are counted at symptom onset, i.e. on entering I:
#   incidence  = new onsets per day = daily change in I+R (only I flows into R)
#   prevalence = number currently infectious, I
#   cumulative = onsets since t=0 = (I+R) - (I+R at t=0)
series_names <- c(incidence  = 'Daily incidence (new onsets)',
                  prevalence = 'Prevalence (infectious, I)',
                  cumulative = 'Cumulative cases (onsets)')

data_stream <- function(I, R, stream)
{
    IR <- as.matrix(I + R)
    switch(stream,
           incidence  = rbind(NA, diff(IR)),
           prevalence = as.matrix(I),
           cumulative = sweep(IR, 2, IR[1,]))
}

# Colours for S, E, I and R as in the SEIR app; each data stream takes the colour of
# the compartment it is counted from (onsets are E -> I, prevalence is I, and
# cumulative onsets have mostly moved on to R)
seir_cols   <- brewer.pal(4,'Set1')[c(1,4,2,3)]
stream_cols <- c(incidence = seir_cols[2], prevalence = seir_cols[3], cumulative = seir_cols[4])

# Transparent version of a colour, for drawing many replicates
alpha_col <- function(col, alpha) {
    rgb_val <- col2rgb(col)/255
    rgb(rgb_val[1], rgb_val[2], rgb_val[3], alpha)
}

# ---- Estimating R0 ------------------------------------------------------------
# All the functions below work on one replicate: vectors indexed by time 0, 1, 2, ...
# (so the value at time t is element t+1)

# Default rule for the exponential growth phase, applied to cumulative onsets so that
# every data stream (and method) uses the same days:
#   start = first day with at least c_min cumulative cases (skips the noisy, chance-
#           driven start of the outbreak)
#   end   = first day cumulative cases reach p_max of the population (so susceptibles
#           have fallen by at most about p_max and growth hasn't slowed much), or the
#           peak of the 7-day average incidence if that comes first (and use_peak is set)
# Returns c(start, end) in days, or NULL if the outbreak never reaches c_min cases
growth_window <- function(cum, N, c_min = 10, p_max = 0.05, use_peak = TRUE)
{
    start <- which(cum >= c_min)[1] - 1
    if(is.na(start)) return(NULL)
    end <- which(cum >= p_max*N)[1] - 1
    if(is.na(end)) end <- length(cum) - 1
    if(use_peak)
    {
        inc  <- c(NA, diff(cum))
        avg7 <- stats::filter(inc, rep(1/7, 7), sides = 2)
        if(any(!is.na(avg7))) end <- min(end, which.max(avg7) - 1)
    }
    c(start, max(start, end))
}

# Exponential growth: y(t) = y0*exp(r*t), so log(y) = log(y0) + r*t is a straight line.
# Fit it by least squares over the window (days with y = 0 have no log, so are left out)
# and turn the growth rate r into R0. For an SEIR model with exponentially distributed
# latent (mean TE) and infectious (mean TI) periods
#   R0 = (1 + r*TE)*(1 + r*TI)
# which is increasing in r, so the confidence limits for r map onto limits for R0
# The formula only holds while both factors are positive: a growth rate below
# -1/max(TE, TI) (e.g. the lower CI limit of a very uncertain fit) corresponds to R0 = 0
r_to_R0 <- function(r, te, ti) ifelse(r <= -1/max(te, ti), 0, (1 + r*te)*(1 + r*ti))

fit_growth <- function(y, win, te, ti)
{
    if(is.null(win)) return(NULL)
    t    <- win[1]:win[2]
    yw   <- y[t + 1]
    keep <- !is.na(yw) & yw > 0
    if(sum(keep) < 3) return(NULL)
    fit <- lm(log(yw[keep]) ~ t[keep])
    r   <- unname(coef(fit)[2])
    ci  <- unname(confint(fit)[2,])
    list(fit = fit, r = r, r_lo = ci[1], r_hi = ci[2],
         R0 = r_to_R0(r, te, ti), R0_lo = r_to_R0(ci[1], te, ti), R0_hi = r_to_R0(ci[2], te, ti),
         n_used = sum(keep), n_dropped = sum(!keep))
}

# Final size: R0 = -log(1 - z)/z, where z is the proportion of the population ever infected
# (final size / N, with final size = S + E + I at the start - S at the end, so the initial
# infections are counted)
final_size_R0 <- function(final_size, N)
{
    z <- final_size/N
    if(is.na(z) || z <= 0 || z >= 1) return(NA)
    -log(1 - z)/z
}

# Generation interval (time from being infected to infecting someone else) for this
# model: an exponential latent period (mean TE) followed by a time within an
# exponential infectious period (mean TI), so it's the sum of two exponentials. EpiEstim
# needs it on whole days, with no mass on day 0. Each day k takes the probability
# within a day either side, weighted towards k (the usual way of discretising, which
# keeps the mean); the small day-0 mass (infecting someone within a day of being
# infected) is moved to day 1. That biases R slightly upwards when TE and TI are both
# short (about 1% for Biggles). G is the integral of the cumulative distribution
# function, so each weight is a second difference of G
gen_interval <- function(te, ti, tol = 1e-6)
{
    a <- 1/te; g <- 1/ti
    G <- function(x)
    {
        x <- pmax(x, 0)
        if(abs(a - g) < 1e-8*max(a, g))   # equal rates: a Gamma(2, a) distribution
            x - (2*(1 - exp(-a*x))/a - x*exp(-a*x))
        else
            x - (a*(1 - exp(-g*x))/g - g*(1 - exp(-a*x))/a)/(a - g)
    }
    # Go out until the cumulative probability is within tol of 1
    K <- 1
    while(1 - (G(K + 1) - G(K)) > tol && K < 2000) K <- K + 1
    k <- 0:K
    w <- G(k + 1) - 2*G(k) + G(k - 1)
    w <- c(0, w[2] + w[1], w[-(1:2)])
    w/sum(w)
}

# EpiEstim (Cori et al. 2013) estimates the reproduction number from incidence and
# the generation interval using the renewal equation. The cases at t=0 are the ones
# that started the outbreak, so they're entered as imported. If t_start/t_end are
# given, R is estimated over those windows (in days); by default over weekly sliding
# windows. Returns EpiEstim's table of estimates, with times in days
epiestim_R <- function(inc, I0, si, t_start = NULL, t_end = NULL)
{
    n <- length(inc)
    # (rounded, as EpiEstim needs whole numbers of cases; simulated counts already are)
    incid <- data.frame(local = c(0, round(inc[-1])), imported = c(round(I0), rep(0, n - 1)))
    if(is.null(t_start))
    {
        t_start <- 1:(n - 7)     # in days; EpiEstim windows can't start on the first day
        t_end   <- t_start + 6
    }
    config <- make_config(list(si_distr = si, t_start = t_start + 1, t_end = t_end + 1))
    est <- suppressWarnings(suppressMessages(
        estimate_R(incid, method = 'non_parametric_si', config = config)))
    out <- est$R
    out$t_start <- out$t_start - 1
    out$t_end   <- out$t_end - 1
    out
}

# Renewal-equation likelihood with susceptible depletion, fitted to the whole epidemic
# (no growth window needed). New onsets on day t are Poisson, or negative binomial
# (allowing more variation), with mean
#   mu_t = R0 * s_t * Lambda_t
#   Lambda_t = sum over k of w_k * onsets_{t-k}, the infection pressure from earlier
#              cases (w is the generation interval, as for EpiEstim)
#   s_t      = fraction of the population susceptible when the day-t cases were infected
# Cases were infected a latent period before onset, and by then more people had been
# infected than had shown symptoms, so s_t is not just N minus cumulative onsets. With
# exponential latent periods, on average onsets per day = E/TE, so the number infected
# but not yet showing symptoms is about TE * daily onsets, giving
#   S(u) = S0 - (cumulative onsets by u) - TE * (onsets around u)
# which is then averaged over the latent period before day t. (Using cumulative onsets
# alone underestimates R0, by about 3% for Biggles and 10% when R0 = 3.) Daily data
# can't keep up with very fast epidemics, so it underestimates when growth is very
# fast (R0 = 5 with Biggles' periods gives about 4.2).
# R0 is fitted by maximum likelihood, with a 95% profile likelihood interval.
# Returns the estimate, interval, overdispersion k (negative binomial) and fitted mean
fit_renewal <- function(inc, init, N, te, gi, family = 'negbin')
{
    n  <- length(inc)
    y  <- c(init[3], inc[-1])        # day 0: the initially infectious cases
    o  <- c(0, inc[-1])              # onsets (none at day 0)
    w  <- gi[-1]                     # generation interval for lags 1, 2, ...
    Lam <- numeric(n)
    for(k in seq_len(min(length(w), n - 1)))
        Lam[(k + 1):n] <- Lam[(k + 1):n] + w[k]*y[1:(n - k)]
    # Susceptibles on each day, then averaged over the latent period before each onset
    o_next <- c(o[-1], 0)
    S_day  <- pmax(init[1] - cumsum(o) - te*(o + o_next)/2, 0)
    d  <- 0:ceiling(30*te)
    fL <- pexp(d + 0.5, 1/te) - pexp(pmax(d - 0.5, 0), 1/te)
    s  <- sapply(1:n, function(t) { dd <- d[d < t]; sum(fL[dd + 1]*S_day[t - dd])/sum(fL[dd + 1]) })/N
    x  <- s*Lam
    keep <- (1:n) > 1 & x > 0
    yy <- y[keep]
    xx <- x[keep]
    if(sum(yy) == 0) return(NULL)

    # Poisson: the likelihood is maximised at R0 = sum(onsets)/sum(s*Lambda)
    Rp <- sum(yy)/sum(xx)
    if(family == 'poisson')
    {
        nll <- function(R) -sum(dpois(yy, R*xx, log = TRUE))
        R0  <- Rp
        k   <- NA
        prof <- nll
        best <- nll(R0)
    }else{
        nll2 <- function(lR, lk) -sum(dnbinom(yy, mu = exp(lR)*xx, size = exp(lk), log = TRUE))
        o2   <- optim(c(log(Rp), 0), function(p) nll2(p[1], p[2]))
        R0   <- exp(o2$par[1])
        k    <- exp(o2$par[2])
        prof <- function(R) optimize(function(lk) nll2(log(R), lk), c(-5, 12))$objective
        best <- o2$value
    }
    # 95% interval: the values of R0 whose (profile) log-likelihood is within
    # qchisq(0.95, 1)/2 of the maximum
    f  <- function(R) 2*(prof(R) - best) - qchisq(0.95, 1)
    lo <- tryCatch(uniroot(f, c(R0/10, R0))$root, error = function(e) NA)
    hi <- tryCatch(uniroot(f, c(R0, R0*10))$root, error = function(e) NA)
    list(R0 = R0, lo = lo, hi = hi, k = k, mu = R0*x, n_days = sum(keep))
}

# All the estimates for one replicate, as a one-row data frame (NA where a method
# fails, e.g. a minor outbreak that never reaches the start of the window).
# win is the growth window, gi the generation interval, family the renewal
# likelihood's distribution; epiestim = FALSE skips EpiEstim (the slowest)
estimate_replicate <- function(S, I, R, pars, win, gi, family = 'negbin', epiestim = TRUE)
{
    out <- data.frame(start = NA, end = NA)
    if(!is.null(win)) out[1,] <- win
    for(s in names(series_names))
    {
        g <- fit_growth(data_stream(I, R, s)[,1], win, pars$te, pars$ti)
        out[[paste0(s, '_R0')]] <- if(is.null(g)) NA else g$R0
        out[[paste0(s, '_lo')]] <- if(is.null(g)) NA else g$R0_lo
        out[[paste0(s, '_hi')]] <- if(is.null(g)) NA else g$R0_hi
    }
    out$final_size <- sum(pars$init[1:3]) - S[length(S)]
    out$final_size_R0 <- final_size_R0(out$final_size, pars$N)
    inc <- data_stream(I, R, 'incidence')[,1]
    out$epiestim_R0 <- NA
    out$epiestim_lo <- NA
    out$epiestim_hi <- NA
    if(epiestim && !is.null(win))
    {
        e <- epiestim_R(inc, pars$init[3], gi, max(1, win[1]), max(1, win[2]))
        out$epiestim_R0 <- e$`Mean(R)`
        out$epiestim_lo <- e$`Quantile.0.025(R)`
        out$epiestim_hi <- e$`Quantile.0.975(R)`
    }
    rn <- tryCatch(fit_renewal(inc, pars$init, pars$N, pars$te, gi, family), error = function(e) NULL)
    out$renewal_R0 <- if(is.null(rn)) NA else rn$R0
    out$renewal_lo <- if(is.null(rn)) NA else rn$lo
    out$renewal_hi <- if(is.null(rn)) NA else rn$hi
    out
}

# Methods as reported in the tables and plots: label and column prefix
methods <- c('Growth rate: daily incidence' = 'incidence',
             'Growth rate: prevalence'      = 'prevalence',
             'Growth rate: cumulative'      = 'cumulative',
             'Final size'                   = 'final_size',
             'EpiEstim (growth window)'     = 'epiestim',
             'Renewal likelihood'           = 'renewal')
method_cols <- c(stream_cols, final_size = 'grey30', epiestim = brewer.pal(8, 'Set1')[5],
                 renewal = '#E7298A')

# Line styles shared by the plots: models drawn for reference are solid grey, the
# regression fit is black (its CI dotted), the renewal model dashed magenta and
# EpiEstim orange (its growth-window estimate a darker orange)
det_col   <- 'grey60'
epiestim_dark <- '#A63603'

# In the shinylive (browser) version a plain click on a download link bypasses the
# service worker that runs the app, so the browser gets "File wasn't available on
# site". Fetch the file from the app instead (fetch() does go through the service
# worker) and save it from a blob. Delegated, so it also covers the re-rendered button
# The Update button turns orange and pulses while there are changes it hasn't applied
update_js <- tags$head(
    tags$style(HTML("
        #update { font-weight: bold; margin-top: 4px; }
        #new_seed { margin-top: 6px; }
        #update.pending { background-color: #e67e22; border-color: #d35400;
                          animation: update-pulse 1.2s infinite; }
        @keyframes update-pulse {
            0%   { box-shadow: 0 0 0 0 rgba(230, 126, 34, 0.7); }
            70%  { box-shadow: 0 0 0 8px rgba(230, 126, 34, 0); }
            100% { box-shadow: 0 0 0 0 rgba(230, 126, 34, 0); } }
        #update_hint { display: none; color: #d35400; font-size: 12px; margin-top: 3px; }
        #update_hint.pending { display: block; }")),
    tags$script(HTML("
        Shiny.addCustomMessageHandler('update_pending', function(pending) {
            ['update', 'update_hint'].forEach(function(id) {
                var el = document.getElementById(id);
                if (el) el.classList.toggle('pending', pending);
            });
        });")))

download_js <- tags$script(HTML("
document.addEventListener('click', async function(e) {
    const link = e.target.closest('a.shiny-download-link');
    if(!link) return;
    e.preventDefault();
    const resp = await fetch(link.href);
    if(!resp.ok) return;
    const m = (resp.headers.get('Content-Disposition') || '').match(/filename=\"([^\"]+)\"/);
    const a = document.createElement('a');
    a.href = URL.createObjectURL(await resp.blob());
    a.download = m ? m[1] : 'download.csv';
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(function() { URL.revokeObjectURL(a.href); }, 1000);
});
"))

# R0 with the 0 as a subscript: R0_html for the page and tables, R0_sub for text that
# can't hold HTML (the browser tab title); plots use plotmath, R[0]
R0_tag  <- HTML('R<sub>0</sub>')   # one string, so no space creeps in before the 0
R0_html <- 'R<sub>0</sub>'
R0_sub  <- 'R\u2080'

# Default parameters: the fictional disease Biggles in a boarding school of 1,300
# (R0 = 2, mean latent period 4 days, mean infectious period 1 day)
default_seed <- 2027

ui <- fluidPage(

    # Application title
    titlePanel(tagList('Estimating ', R0_tag), windowTitle = paste('Estimating', R0_sub)),
    download_js,
    update_js,

    # Sidebar: model inputs, simulation options, what to plot and the growth window
    sidebarLayout(
        sidebarPanel(
            fluidRow(column(width=4,
                            numericInput('R0',
                                         R0_tag,
                                         min = 0,
                                         max = 20,
                                         value = 2,
                                         step = 0.1),
                            numericInput('te',
                                         'TE (latent period)',
                                         min=0.1,
                                         max=30,
                                         value=4,
                                         step=0.1),
                            numericInput('ti',
                                         'TI (infectious period)',
                                         min=0.1,
                                         max=30,
                                         value=1,
                                         step=0.1),
                            uiOutput('beta_info')),
                     column(width=4,
                            numericInput('S0',
                                         'Initial S',
                                         value = 1299,
                                         min = 0,
                                         max = 1000000),
                            numericInput('E0',
                                         'Initial E',
                                         value = 0,
                                         min = 0,
                                         max = 1000000),
                            numericInput('I0',
                                         'Initial I',
                                         value = 1,
                                         min = 0,
                                         max = 1000000),
                            numericInput('Rinit',
                                         'Initial R',
                                         value = 0,
                                         min = 0,
                                         max = 1000000)),
                     column(width=4,
                            numericInput('t_run',
                                         'Time',
                                         value = 100,
                                         min = 1,
                                         max = 1000),
                            numericInput('reps',
                                         'Replicates',
                                         value = 1,
                                         min = 1,
                                         max = 500),
                            numericInput('seed',
                                         'Random seed',
                                         value = default_seed,
                                         min = 1,
                                         max = 100000),
                            actionButton('update', 'Update', icon = icon('play'),
                                         class = 'btn-primary', width = '100%'),
                            tags$div(id = 'update_hint', 'Press Update to apply your changes'),
                            actionButton('new_seed',
                                         'New epidemic', width = '100%')),width=6),
            fluidRow(column(width=6,
                            checkboxInput('major_only',
                                          'Major outbreaks only',
                                          value = TRUE),
                            numericInput('major_thresh',
                                         'Major outbreak: final size at least (% of N)',
                                         value = 10,
                                         min = 0,
                                         max = 100)),
                     column(width=6,
                            radioButtons('stream',
                                         'Data',
                                         choiceNames = lapply(names(series_names), function(s){
                                             span(series_names[[s]], style=paste0('color:',stream_cols[[s]],';font-weight:bold'))}),
                                         choiceValues = names(series_names),
                                         selected = 'incidence'),
                            checkboxInput('log_y', 'Log scale', value = FALSE),
                            checkboxInput('show_det', 'Show deterministic', value = TRUE),
                            checkboxInput('show_cum', 'Plot cumulative infections (final size)', value = FALSE))),
            tags$hr(),
            tags$h4('Exponential growth window', style = 'font-size:16px'),
            fluidRow(column(width=6,
                            numericInput('c_min',
                                         'Start: cumulative cases at least',
                                         value = 10,
                                         min = 1),
                            numericInput('p_max',
                                         'End: cumulative cases reach (% of N)',
                                         value = 5,
                                         min = 0.1,
                                         max = 100),
                            checkboxInput('use_peak',
                                          'or end at peak incidence, if earlier',
                                          value = TRUE)),
                     column(width=6,
                            checkboxInput('auto_window',
                                          'Use this rule',
                                          value = TRUE),
                            conditionalPanel('!input.auto_window',
                                             sliderInput('window',
                                                         'Window (days)',
                                                         min = 0,
                                                         max = 100,
                                                         value = c(0, 30),
                                                         step = 1)),
                            helpText('The window is shaded on the plots. With more than one',
                                     'replicate the rule is applied to each one.'))),
            tags$h4('Renewal likelihood', style = 'font-size:16px'),
            radioButtons('renewal_family',
                         'Daily onsets vary as',
                         choices = c('Negative binomial' = 'negbin', 'Poisson' = 'poisson'),
                         selected = 'negbin',
                         inline = TRUE),

            downloadButton("downloadData", "Download Simulated Data")
        ),
        mainPanel(
            tabsetPanel(id = 'tabs',
                # Replicate 1: data, fitted growth and the estimates of R0
                tabPanel('Epidemic',
                         plotOutput("tsPlot",hover = "plot_hover"),
                         textOutput("info"),
                         conditionalPanel('input.show_cum',
                                          plotOutput("cumPlot", height = 320, hover = "cum_hover"),
                                          textOutput("cum_info")),
                         tags$hr(),
                         tags$h4('Estimates of ', R0_tag),
                         tableOutput("estimates"),
                         tableOutput("growth"),
                         uiOutput("warnings"),
                         helpText(tags$b('Growth rate:'), 'a straight line log(y) = a + r t is fitted to the',
                                  'chosen data over the window; for this model', R0_tag, '= (1 + r TE)(1 + r TI).',
                                  tags$br(), tags$b('Final size:'), R0_tag, '= −ln(1 − z)/z, where z is the proportion of the',
                                  'population ever infected, including the initial cases (final size / N).',
                                  tags$br(), tags$b('EpiEstim:'), 'R over the whole window from daily incidence and',
                                  'the generation interval (see the R(t) tab).',
                                  tags$br(), tags$b('Renewal likelihood:'), 'onsets on day t vary around', R0_tag, '× s(t) ×',
                                  'Σ w(k) onsets(t − k), where w is the generation interval and s(t) the fraction',
                                  'susceptible when those cases were infected.', R0_tag, 'is fitted by maximum likelihood to',
                                  'the whole epidemic, so no window is needed.'),
                         tags$h4('Outbreak'),
                         tableOutput("summary")),
                # EpiEstim's R(t) over weekly windows, and the generation interval it uses
                tabPanel('R(t)',
                         plotOutput("rtPlot"),
                         helpText('EpiEstim estimates of R over sliding 7-day windows, plotted at the end of',
                                  'each window, with 95% credible intervals. Windows with fewer than 12 cases',
                                  'are left out. The grey line is', R0_tag, '× S(t)/N, how many people each case',
                                  'would infect given the susceptibles left.'),
                         plotOutput("siPlot", height = 250)),
                # The same estimates for every replicate
                tabPanel('Many simulations',
                         uiOutput('many_note'),
                         checkboxInput('many_epiestim', 'Include EpiEstim (slower)', value = TRUE),
                         plotOutput("manyPlot"),
                         tableOutput("manyTable"),
                         plotOutput("fsHist", height = 300),
                         downloadButton("downloadEstimates", "Download Estimates")))
        )
    )
)


server <- function(input, output, session) {

# Parameters are only read when Update (or New epidemic) is pressed, so a half-typed
# value doesn't trigger a run
rv <- reactiveValues(sims = NULL, det = NULL, pars = NULL)

run_model <- function(seed)
{
    init  <- c(input$S0, input$E0, input$I0, input$Rinit)
    beta  <- input$R0/input$ti
    t_run <- input$t_run
    reps  <- max(1, round(input$reps))

    if(any(is.na(c(init, beta, input$te, input$ti, t_run, seed, reps))) || sum(init) <= 0 ||
       input$te <= 0 || input$ti <= 0 || t_run < 1)
    {
        showNotification('Check the inputs: all values must be filled in, the periods must be positive and the population non-empty.', type = 'error')
        return()
    }

    withProgress(message = 'Simulating', value = 0.5, {
        set.seed(seed)
        rv$sims <- sim_major(reps, beta, input$te, input$ti, init, t_run,
                             input$major_only, input$major_thresh/100)
        rv$det  <- mass_actionD(beta,1.0/input$te,1.0/input$ti,init[1],init[2],init[3],init[4],seq(0,t_run))
    })
    rv$pars <- list(R0 = input$R0, beta = beta, te = input$te, ti = input$ti, init = init,
                    N = sum(init), t_run = t_run, reps = reps, seed = seed,
                    major_only = input$major_only, major_thresh = input$major_thresh)
}

# Run once at start-up and whenever Update is pressed
observeEvent(input$update, run_model(input$seed), ignoreNULL = FALSE)

# A new random seed: show it (so the epidemic can be recreated) and run with it
observeEvent(input$new_seed, {
    seed <- as.integer((as.numeric(Sys.time())*1000) %% 100000) + 1L
    updateNumericInput(session, 'seed', value = seed)
    run_model(seed)
})

# Flag the Update button while the inputs differ from those of the epidemic shown
observe({
    p <- rv$pars
    req(p)
    now <- list(input$R0, input$te, input$ti, input$S0, input$E0, input$I0, input$Rinit, input$t_run,
                max(1, round(input$reps)), input$seed, input$major_only, input$major_thresh)
    old <- list(p$R0, p$te, p$ti, p$init[1], p$init[2], p$init[3], p$init[4], p$t_run,
                p$reps, p$seed, p$major_only, p$major_thresh)
    session$sendCustomMessage('update_pending', !isTRUE(all.equal(now, old, check.attributes = FALSE)))
})

output$beta_info <- renderUI({
    if(is.na(input$R0) || is.na(input$ti) || input$ti <= 0) return(NULL)
    HTML(paste0('Beta = ', R0_html, '/TI = ', signif(input$R0/input$ti, 3)))
})

# ---- Replicate 1 ----

# The chosen data stream for all replicates (time x replicate) and the deterministic model
stream_data <- reactive({
    req(rv$sims)
    list(time = 0:rv$pars$t_run,
         stoch = data_stream(rv$sims$I, rv$sims$R, input$stream),
         det   = data_stream(rv$det[,4], rv$det[,5], input$stream)[,1])
})

# The window rule's settings, or NULL if one is missing or out of range
rule_settings <- reactive({
    if(any(is.na(c(input$c_min, input$p_max))) || input$c_min < 0 || input$p_max <= 0) return(NULL)
    list(c_min = input$c_min, p_max = input$p_max/100, use_peak = input$use_peak)
})

window_rule <- function(I, R, rule)
{
    if(is.null(rule)) return(NULL)
    growth_window(data_stream(I, R, 'cumulative')[,1], rv$pars$N, rule$c_min, rule$p_max, rule$use_peak)
}

rule_win <- reactive({
    req(rv$sims)
    window_rule(rv$sims$I[,1], rv$sims$R[,1], rule_settings())
})

# Keep the manual slider in step with the rule (and the run length) while the rule is
# in use, so switching to manual starts from the rule's window
observe({
    req(rv$pars)
    if(isTRUE(input$auto_window))
    {
        w <- rule_win()
        if(is.null(w)) w <- c(0, rv$pars$t_run)
        updateSliderInput(session, 'window', min = 0, max = rv$pars$t_run, value = w)
    }
})

# The growth window used for replicate 1: the rule's, or the one set by hand
win <- reactive({
    if(isTRUE(input$auto_window)) rule_win() else input$window
})

growth_fit <- reactive({
    req(rv$sims)
    fit_growth(stream_data()$stoch[,1], win(), rv$pars$te, rv$pars$ti)
})

gen_int <- reactive({
    req(rv$pars)
    gen_interval(rv$pars$te, rv$pars$ti)
})

incidence1 <- reactive({
    req(rv$sims)
    data_stream(rv$sims$I[,1], rv$sims$R[,1], 'incidence')[,1]
})

# EpiEstim over the growth window (the estimate of R0) and over weekly windows (R(t)).
# NULL if it fails, e.g. too few days
epiestim_window <- reactive({
    w <- win()
    if(is.null(w)) return(NULL)
    tryCatch(epiestim_R(incidence1(), rv$pars$init[3], gen_int(), max(1, w[1]), max(1, w[2])),
             error = function(e) NULL)
})

# Renewal likelihood for replicate 1, over the whole epidemic
renewal_fit <- reactive({
    req(rv$sims)
    tryCatch(fit_renewal(incidence1(), rv$pars$init, rv$pars$N, rv$pars$te, gen_int(), input$renewal_family),
             error = function(e) NULL)
})

epiestim_sliding <- reactive({
    if(rv$pars$t_run < 9) return(NULL)
    tryCatch(epiestim_R(incidence1(), rv$pars$init[3], gen_int()), error = function(e) NULL)
})

output$tsPlot <- renderPlot({
    d <- stream_data()
    col <- stream_cols[[input$stream]]
    y1  <- d$stoch[,1]
    n_reps <- ncol(d$stoch)
    w <- win()
    g <- growth_fit()

    # Zeros can't be shown on a log scale
    if(input$log_y)
    {
        d$stoch[d$stoch <= 0] <- NA
        d$det[d$det <= 0] <- NA
        y1[y1 <= 0] <- NA
    }
    yvals <- c(d$stoch, if(input$show_det) d$det)
    if(all(is.na(yvals))) yvals <- 1

    ylim <- range(yvals,na.rm=TRUE)
    # On the log scale start at half a case, so the deterministic curve's early
    # fractions of a case don't squash the data
    if(input$log_y) ylim[1] <- max(ylim[1], 0.5)
    plot(range(d$time),ylim,type='n',xlab='Time (days)',yaxt='n',
         ylab=series_names[[input$stream]],log=if(input$log_y) 'y' else '')
    axis(2, at=axTicks(2), labels=format(axTicks(2), scientific=FALSE, drop0trailing=TRUE, trim=TRUE))
    # Shade the growth window
    if(!is.null(w))
    {
        usr <- par('usr')
        ybox <- if(input$log_y) 10^usr[3:4] else usr[3:4]
        rect(w[1] - 0.5, ybox[1], w[2] + 0.5, ybox[2], col = 'grey92', border = NA)
    }
    grid()
    # Other replicates faintly, the first (the one the estimates use) on top
    if(n_reps > 1)
        matlines(d$time, d$stoch[,-1,drop=FALSE], col=alpha_col(col,0.15), lty=1)
    if(input$show_det)
        lines(d$time, d$det, col=det_col, lwd=2.5)
    if(input$stream == 'incidence')
        points(d$time, y1, pch=19, col=col)
    else
        lines(d$time, y1, col=col, lwd=3)
    # Fitted exponential growth, with its 95% confidence band
    if(!is.null(g))
    {
        tt <- seq(w[1], w[2], length.out = 100)
        b  <- coef(g$fit)
        V  <- vcov(g$fit)
        se <- sqrt(V[1,1] + 2*tt*V[1,2] + tt^2*V[2,2])
        q  <- qt(0.975, df.residual(g$fit))
        lines(tt, exp(b[1] + b[2]*tt), lwd = 2.5)
        lines(tt, exp(b[1] + b[2]*tt - q*se), lty = 3, lwd = 1.5)
        lines(tt, exp(b[1] + b[2]*tt + q*se), lty = 3, lwd = 1.5)
    }
    # Renewal model's expected onsets (each day given the cases before it), for incidence
    rn <- if(input$stream == 'incidence') renewal_fit() else NULL
    if(!is.null(rn))
    {
        mu <- rn$mu
        if(input$log_y) mu[mu <= 0] <- NA
        lines(d$time[-1], mu[-1], col = method_cols[['renewal']], lwd = 2.5, lty = 2)
    }

    leg_text <- c(if(n_reps > 1) 'Replicate 1' else 'Simulated', if(n_reps > 1) 'Other replicates',
                  if(input$show_det) 'Deterministic', if(!is.null(g)) c('Fit', '95% CI'), if(!is.null(rn)) 'Renewal model')
    leg_col  <- c(col, if(n_reps > 1) alpha_col(col,0.4), if(input$show_det) det_col, if(!is.null(g)) c('black', 'black'),
                  if(!is.null(rn)) method_cols[['renewal']])
    leg_lty  <- c(if(input$stream == 'incidence') NA else 1, if(n_reps > 1) 1, if(input$show_det) 1, if(!is.null(g)) c(1, 3),
                  if(!is.null(rn)) 2)
    leg_pch  <- c(if(input$stream == 'incidence') 19 else NA, if(n_reps > 1) NA, if(input$show_det) NA, if(!is.null(g)) c(NA, NA),
                  if(!is.null(rn)) NA)
    leg_lwd  <- c(3, if(n_reps > 1) 1, if(input$show_det) 2.5, if(!is.null(g)) c(2.5, 1.5), if(!is.null(rn)) 2.5)
    # Trailing spaces keep each entry clear of the next one's line
    legend('top',inset=c(0,-0.12),xpd=TRUE,horiz=TRUE,bty='n',text.width=NA,legend=paste0(leg_text,'    '),
           col=leg_col,lty=leg_lty,pch=leg_pch,lwd=leg_lwd)
})

# Everyone ever infected (S, E and I at the start minus S now), which ends at the final
# size. Other replicates are drawn faintly, with their final sizes as a sideways
# histogram to the right of the time axis
output$cumPlot <- renderPlot({
    req(rv$sims, isTRUE(input$show_cum))
    p    <- rv$pars
    t    <- 0:p$t_run
    n0   <- sum(p$init[1:3])
    cum  <- n0 - as.matrix(rv$sims$S)
    det  <- n0 - rv$det[,2]
    fs   <- cum[nrow(cum),]
    n_reps <- ncol(cum)
    col  <- seir_cols[4]

    # Room on the right for the histogram of final sizes
    x_hist <- p$t_run*c(1.04, 1.22)
    par(mar = c(4, 4, 2.5, 4.5))
    plot(c(0, if(n_reps > 1) x_hist[2] else p$t_run), c(0, p$N), type = 'n', xaxt = 'n',
         xlab = 'Time (days)', ylab = 'Ever infected')
    axis(1, at = pretty(t)[pretty(t) <= p$t_run])
    z_at <- seq(0, 1, 0.2)
    axis(4, at = z_at*p$N, labels = z_at)
    mtext('Proportion infected, z', side = 4, line = 2.5)
    grid()
    if(n_reps > 1)
        matlines(t, cum[,-1,drop=FALSE], col = alpha_col(col, 0.15), lty = 1)
    if(input$show_det)
        lines(t, det, col = det_col, lwd = 2.5)
    lines(t, cum[,1], col = col, lwd = 3)
    # Final size of replicate 1, labelled so it can be read off
    abline(h = fs[1], lty = 2, col = col)
    text(0, fs[1], sprintf('Final size %d (z = %.3f)', round(fs[1]), fs[1]/p$N),
         adj = c(0, if(fs[1] > 0.85*p$N) 1.4 else -0.5), col = col, font = 2)
    if(n_reps > 1)
    {
        h <- hist(fs, breaks = seq(0, p$N, length.out = 27), plot = FALSE)
        abline(v = mean(c(p$t_run, x_hist[1])), col = 'grey70')
        w <- diff(x_hist)*h$counts/max(h$counts)
        rect(x_hist[1], h$breaks[-length(h$breaks)], x_hist[1] + w, h$breaks[-1],
             col = alpha_col(col, 0.6), border = 'white')
        text(mean(x_hist), par('usr')[4], paste0('Final sizes (', n_reps, ' runs)'),
             xpd = TRUE, pos = 3, cex = 0.85)
    }
    leg_text <- c(if(n_reps > 1) 'Replicate 1' else 'Simulated', if(n_reps > 1) 'Other replicates',
                  if(input$show_det) 'Deterministic', 'Final size')
    leg_col  <- c(col, if(n_reps > 1) alpha_col(col, 0.4), if(input$show_det) det_col, col)
    leg_lty  <- c(1, if(n_reps > 1) 1, if(input$show_det) 1, 2)
    leg_lwd  <- c(3, if(n_reps > 1) 1, if(input$show_det) 2.5, 1)
    legend('top', inset = c(0, -0.14), xpd = TRUE, horiz = TRUE, bty = 'n', text.width = NA,
           legend = paste0(leg_text, '    '), col = leg_col, lty = leg_lty, lwd = leg_lwd)
})

output$cum_info <- renderText({
    e <- input$cum_hover
    if(is.null(e) || is.null(rv$pars) || e$x > rv$pars$t_run) return("")
    paste0("Time = ", round(e$x, 1), ", Ever infected = ", round(e$y), " (z = ",
           formatC(e$y/rv$pars$N, digits = 3, format = 'f'), ")")
})

output$info <- renderText({

    # Read-out of the cursor position; blank when the cursor is off the plot
    e <- input$plot_hover
    if(is.null(e)) return("")
    paste0("Time = ", round(e$x, 1), ", Number = ", round(e$y, 1))

})

fmt2 <- function(x) ifelse(is.na(x), '–', formatC(x, digits = 2, format = 'f'))
fmt_ci <- function(lo, hi) ifelse(is.na(lo), '–', paste0(fmt2(lo), ' – ', fmt2(hi)))

output$estimates <- renderTable({
    req(rv$sims)
    p <- rv$pars
    g <- growth_fit()
    e <- epiestim_window()
    S <- rv$sims$S[,1]
    final_size <- sum(p$init[1:3]) - S[length(S)]
    z  <- final_size/p$N
    fs <- final_size_R0(final_size, p$N)
    rn <- renewal_fit()
    data.frame(Method = c('Growth rate (regression)', 'Final size', 'EpiEstim (growth window)',
                          paste0('Renewal likelihood (', if(input$renewal_family == 'poisson') 'Poisson' else 'negative binomial', ')')),
               Data = c(series_names[[input$stream]], paste0('Proportion infected z = ', formatC(z, digits = 3, format = 'f')), 'Daily incidence',
                        'Daily incidence, whole epidemic'),
               `R<sub>0</sub> estimate` = c(fmt2(if(is.null(g)) NA else g$R0), fmt2(fs),
                                 fmt2(if(is.null(e)) NA else e$`Mean(R)`), fmt2(if(is.null(rn)) NA else rn$R0)),
               `95% interval` = c(if(is.null(g)) '–' else fmt_ci(g$R0_lo, g$R0_hi), 'none from one epidemic',
                                  if(is.null(e)) '–' else fmt_ci(e$`Quantile.0.025(R)`, e$`Quantile.0.975(R)`),
                                  if(is.null(rn)) '–' else fmt_ci(rn$lo, rn$hi)),
               `True R<sub>0</sub>` = format(p$R0),
               check.names = FALSE)
}, striped = TRUE, sanitize.text.function = function(x) x)

output$growth <- renderTable({
    g <- growth_fit()
    w <- win()
    if(is.null(w)) return(NULL)
    data.frame(Quantity = c('Growth window', 'Days fitted', 'Growth rate r (per day)', 'Doubling time (days)'),
               Value = c(paste0('day ', w[1], ' to ', w[2]),
                         if(is.null(g)) '–' else paste0(g$n_used, if(g$n_dropped > 0) paste0(' (', g$n_dropped, ' with no cases left out)')),
                         if(is.null(g)) '–' else paste0(formatC(g$r, digits = 3, format = 'f'), ' (95% CI ',
                                                             formatC(g$r_lo, digits = 3, format = 'f'), ' – ',
                                                             formatC(g$r_hi, digits = 3, format = 'f'), ')'),
                         if(is.null(g) || g$r <= 0) '–' else formatC(log(2)/g$r, digits = 1, format = 'f')))
}, striped = TRUE, sanitize.text.function = function(x) x)

# Summary of replicate 1 (the epidemic the estimates are based on)
output$summary <- renderTable({
    req(rv$sims)
    p   <- rv$pars
    S   <- rv$sims$S[,1]
    inc <- incidence1()
    final_size <- sum(p$init[1:3]) - S[length(S)]
    data.frame(Quantity = c('Beta', 'Population size N', 'Final size (ever infected)',
                            'Proportion infected z', 'Peak daily incidence', 'Day of peak incidence', 'Random seed'),
               Value = c(format(signif(p$beta,3)), format(p$N), format(final_size),
                         formatC(final_size/p$N, digits = 3, format = 'f'), format(max(inc, na.rm=TRUE)),
                         format(which.max(inc) - 1), format(p$seed)))
}, striped = TRUE, sanitize.text.function = function(x) x)

# Warn about anything that makes the estimates unreliable
output$warnings <- renderUI({
    req(rv$sims)
    p   <- rv$pars
    msgs <- character(0)
    still_going <- rv$sims$E[nrow(rv$sims$E),] + rv$sims$I[nrow(rv$sims$I),] > 0
    if(still_going[1])
        msgs <- c(msgs, 'The epidemic had not finished by the end of the run, so the final size estimate is too low: increase Time.')
    else if(any(still_going))
        msgs <- c(msgs, paste(sum(still_going), 'of', p$reps, 'replicates had not finished by the end of the run.'))
    if(p$major_only && rv$sims$n_minor > 0)
        msgs <- c(msgs, paste0(rv$sims$n_minor, ' of ', p$reps, ' replicates were still minor outbreaks (final size below ',
                               p$major_thresh, '% of N) after 50 attempts. Is ', R0_html, ' above 1, and is Time long enough?'))
    w <- win()
    g <- growth_fit()
    if(is.null(w))
        msgs <- c(msgs, 'The outbreak never reached the number of cases that starts the growth window.')
    else if(is.null(g))
        msgs <- c(msgs, 'Too few days with cases in the growth window to fit a line: widen the window.')
    else if(g$n_used < 5)
        msgs <- c(msgs, 'Only a few days in the growth window, so the growth rate is very uncertain.')
    if(length(msgs) == 0) return(NULL)
    tags$div(style='color:#b35900', lapply(msgs, function(m) tags$p(HTML(m))))
})

# R(t) from EpiEstim, with the growth window and the R0 estimate from it
output$rtPlot <- renderPlot({
    e <- epiestim_sliding()
    validate(need(!is.null(e), 'EpiEstim could not estimate R(t) for this epidemic.'))
    p   <- rv$pars
    inc <- incidence1()
    # Leave out windows with too few cases to say much
    n_cases <- sapply(seq_len(nrow(e)), function(i) sum(inc[(e$t_start[i]:e$t_end[i]) + 1]))
    e <- e[n_cases >= 12,]
    validate(need(nrow(e) > 0, 'Too few cases to estimate R(t).'))
    t  <- e$t_end
    lo <- e$`Quantile.0.025(R)`
    hi <- e$`Quantile.0.975(R)`
    w  <- win()
    ew <- epiestim_window()
    S  <- rv$sims$S[,1]

    plot(range(t), c(0, max(c(hi, 1.2*p$R0), na.rm = TRUE)), type = 'n', xlab = 'Time (days)', ylab = 'R(t)')
    if(!is.null(w))
    {
        usr <- par('usr')
        rect(w[1] - 0.5, usr[3], w[2] + 0.5, usr[4], col = 'grey92', border = NA)
    }
    grid()
    polygon(c(t, rev(t)), c(lo, rev(hi)), col = alpha_col(method_cols[['epiestim']], 0.3), border = NA)
    lines(t, e$`Mean(R)`, col = method_cols[['epiestim']], lwd = 3)
    lines(0:p$t_run, p$R0*S/p$N, col = det_col, lwd = 2.5)
    abline(h = 1, col = 'grey50')
    if(!is.null(ew))
    {
        segments(w[1], ew$`Mean(R)`, w[2], ew$`Mean(R)`, lwd = 3, col = epiestim_dark)
        rect(w[1], ew$`Quantile.0.025(R)`, w[2], ew$`Quantile.0.975(R)`, border = epiestim_dark, lty = 2, lwd = 1.5)
    }
    legend('top', inset = c(0, -0.12), xpd = TRUE, horiz = TRUE, bty = 'n', text.width = NA,
           legend = c(expression('R(t), 7-day windows    '), expression(R[0] %*% S(t)/N * '    '),
                      expression('Growth window estimate    ')),
           col = c(method_cols[['epiestim']], det_col, epiestim_dark), lty = c(1, 1, 1), lwd = c(3, 2.5, 3))
})

# The generation interval EpiEstim is given
output$siPlot <- renderPlot({
    si <- gen_int()
    k  <- seq_along(si) - 1
    show <- k <= max(which(cumsum(si) < 0.995))
    par(mar = c(4, 4, 2, 1))
    plot(k[show], si[show], type = 'h', lwd = 6, col = 'grey40', xlab = 'Days from infection',
         ylab = 'Probability', main = paste0('Generation interval (latent + infectious period), mean ',
                                             formatC(sum(k*si), digits = 1, format = 'f'), ' days'),
         font.main = 1, cex.main = 1)
    grid()
})

# ---- Many simulations ----

# The same rule and estimators applied to every replicate
many_est <- reactive({
    req(rv$sims, rv$pars$reps > 1)
    rule <- rule_settings()
    gi   <- gen_int()
    ee   <- isTRUE(input$many_epiestim)
    fam  <- input$renewal_family
    p    <- rv$pars
    n    <- p$reps
    withProgress(message = 'Estimating', value = 0, {
        est <- lapply(1:n, function(i) {
            if(i %% 10 == 0) incProgress(10/n)
            I <- rv$sims$I[,i]
            R <- rv$sims$R[,i]
            tryCatch(estimate_replicate(rv$sims$S[,i], I, R, p, window_rule(I, R, rule), gi, fam, ee),
                     error = function(e) estimate_replicate(rv$sims$S[,i], I, R, p, NULL, gi, fam, FALSE))
        })
    })
    cbind(replicate = 1:n, do.call(rbind, est))
})

output$many_note <- renderUI({
    if(is.null(rv$pars) || rv$pars$reps > 1) return(NULL)
    helpText('Set Replicates to more than 1 and press Update to see how the estimates vary from one simulated epidemic to the next.')
})

# Estimates for each method in a list, leaving out methods with none (e.g. EpiEstim switched off)
many_values <- function(est)
{
    v <- lapply(methods, function(m) { x <- est[[paste0(m, '_R0')]]; x[!is.na(x)] })
    v[sapply(v, length) > 0]
}

output$manyPlot <- renderPlot({
    est <- many_est()
    v <- many_values(est)
    validate(need(length(v) > 0, 'No estimates: did any outbreak reach the growth window?'))
    v <- rev(v)  # first method at the top
    cols <- method_cols[methods[names(v)]]
    par(mar = c(4, 13, 3.5, 1))
    boxplot(v, horizontal = TRUE, las = 1, outline = FALSE, col = alpha_col('grey', 0.3), border = 'grey40',
            xlab = expression('Estimated ' * R[0]), ylim = range(c(unlist(v), rv$pars$R0)))
    title(main = paste('Estimates from', rv$pars$reps, 'simulated epidemics'), font.main = 1, cex.main = 1, line = 1.5)
    grid(ny = NA)
    stripchart(v, method = 'jitter', jitter = 0.2, pch = 19, cex = 0.6,
               col = sapply(cols, alpha_col, 0.4), add = TRUE)
    abline(v = rv$pars$R0, col = 'red', lwd = 2, lty = 2)
    mtext(expression('True ' * R[0]), side = 3, line = 0.1, at = rv$pars$R0, col = 'red', cex = 0.9)
})

output$manyTable <- renderTable({
    est <- many_est()
    R0  <- rv$pars$R0
    rows <- lapply(names(methods), function(lab) {
        m <- methods[[lab]]
        x <- est[[paste0(m, '_R0')]]
        if(all(is.na(x))) return(NULL)
        lo <- est[[paste0(m, '_lo')]]
        hi <- est[[paste0(m, '_hi')]]
        cover <- if(is.null(lo)) '–' else
            sprintf('%.0f%%', 100*mean(lo <= R0 & R0 <= hi, na.rm = TRUE))
        data.frame(Method = lab, Estimates = sum(!is.na(x)),
                   Mean = fmt2(mean(x, na.rm = TRUE)), Median = fmt2(median(x, na.rm = TRUE)),
                   `Middle 95%` = fmt_ci(quantile(x, 0.025, na.rm = TRUE), quantile(x, 0.975, na.rm = TRUE)),
                   `Intervals containing true R<sub>0</sub>` = cover,
                   check.names = FALSE)
    })
    do.call(rbind, rows)
}, striped = TRUE, sanitize.text.function = function(x) x)

output$fsHist <- renderPlot({
    est <- many_est()
    hist(est$final_size, breaks = 30, col = 'black', border = 'white', xlab = 'Final size (ever infected)',
         main = paste0('Final size distribution', if(rv$pars$major_only) ' (major outbreaks only)'),
         font.main = 1, cex.main = 1)
    grid()
})

# All replicates in long format, with the three data streams alongside the compartments
output$downloadData <- downloadHandler(
    filename = function() {
        paste('R0-simulation-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
        s <- rv$sims
        n_reps <- ncol(s$S)
        out <- data.frame(replicate = rep(1:n_reps, each = nrow(s$S)),
                          time = rep(0:(nrow(s$S)-1), n_reps),
                          S = c(s$S), E = c(s$E), I = c(s$I), R = c(s$R),
                          incidence  = c(data_stream(s$I, s$R, 'incidence')),
                          cumulative = c(data_stream(s$I, s$R, 'cumulative')))
        write.csv(out, con, row.names=FALSE)
    }
)

# One row per replicate: growth window, each method's estimate and interval
output$downloadEstimates <- downloadHandler(
    filename = function() {
        paste('R0-estimates-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
        write.csv(many_est(), con, row.names=FALSE)
    }
)

}

# Run the application
shinyApp(ui = ui, server = server)
