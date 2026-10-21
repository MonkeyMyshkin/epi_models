
if(!requireNamespace("igraph", quietly=TRUE)) install.packages("igraph")

library(shiny)

library(igraph)

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

mass_actionD <- function(beta,k,g,N0,I0,time_vec)
{
SI.dyn <- function(t,var,par) {
    S <- var[1]
    I <- var[2]
    N <- S + I

    beta <- par[1]
    k    <- par[2]
    g    <- par[3]

    dS <- g*I - k*beta*S*I/N
    dI <- k*beta*S*I/N - g*I

    return(list(c(dS,dI)))
}

SI.par<-c(beta,k,g)
SI.init<-c(N0-I0,I0)

oot<-rk4(SI.init,time_vec,SI.dyn,SI.par)

return(oot)
}

gillespie_stepSIS <- function(net, beta, gamma, tsim, intro)
{
    t = 0

    infected = V(net)[status=='I']

    # Only edges between an infected and a susceptible person can transmit
    SIpairs  = E(net)[infected %--% V(net)[status=='S']]
    rate = gamma*length(infected)+beta*length(SIpairs)

    while(t < tsim)
    {
        # If disease has gone extinct, reintroduce at rate intro
        if(length(infected) == 0)
        {
            # With no reintroduction nothing more can happen in this step
            if(intro == 0){break;}
            delta_t = rexp(1,intro)
            if((t + delta_t) > tsim){break;}
            susceptible = V(net)[status=='S']
            V(net)[susceptible[sample.int(length(susceptible),1)]]$status = 'I'
            t = t + delta_t
        }else
        {
            delta_t = rexp(1,rate)
            if((t + delta_t) > tsim){break;}

            if(runif(1,0,rate) < gamma*length(infected))
            {
                # Recovery
                V(net)[infected[sample.int(length(infected),1)]]$status = 'S'
            }else
            {
                # Infection: the susceptible end of a random S-I edge
                edge = SIpairs[sample.int(length(SIpairs),1)]
                v = head_of(net,edge)
                if(is.element(v,infected))
                {
                    v = tail_of(net,edge)
                }
                V(net)[v]$status = 'I'
            }
            t = t + delta_t
        }

        infected = V(net)[status=='I']
        SIpairs  = E(net)[infected %--% V(net)[status=='S']]
        rate = gamma*length(infected)+beta*length(SIpairs)
    }

    return(net)
}


ui <- fluidPage(

    # Application title
    titlePanel("SIS Model"),

    sidebarLayout(
        sidebarPanel(
            selectInput('network',
                        'Contact Network',
                        choices=c('Pathology 1B Class' = 'pathology1B',
                                  'High School' = 'hschool',
                                  'High School 2' = 'hschool2',
                                  'Random' = 'mass_action')),
            uiOutput('degreeUI'),
            uiOutput('degree_msg'),
            numericInput('I0',
                         'Initial Infected',
                         value = 5,
                         min = 1,
                         max = 100),
            uiOutput('I0_msg'),
            numericInput("beta",
                        "Transmission Rate",
                        min = 0,
                        max = 2,
                        value = 0.04),
            uiOutput('beta_msg'),
            numericInput('ti',
                        'Infectious Period',
                        min=0.5,
                        max=365,
                        value=10),
            uiOutput('ti_msg'),
            numericInput('intro',
                         'Reintroduction Rate',
                         min=0,
                         max=2,
                         value=0.5),
            uiOutput('intro_msg'),
            actionButton('simulate',
                         'Stop/Go!'),
            actionButton('step',
                         "Step"),
            actionButton('reset',
                         "Reset"),
            uiOutput('pending_msg'),
            plotOutput("timeseries")
            ),
        # Network plot
        mainPanel(
           # Fill the panel width and the window height below the title (at
           # least 400px); igraph keeps the network's aspect ratio inside this area
           plotOutput("netPlot",width='100%',height='calc(max(400px, 100vh - 100px))')
        )
    )
)


edgelist <- read.table('SchoolAE.csv')

nodes    <- read.table('SchoolAN.csv')
nodes    <- data.frame(id=nodes[,1],status='S',stringsAsFactors = FALSE)

edgelist2 <- read.table('SchoolBE.csv')
# Remove edges with missing endpoints
edgelist2 <- edgelist2[complete.cases(edgelist2$from,edgelist2$to),]

nodes2    <- read.table('SchoolBN.csv')
nodes2    <- data.frame(id=nodes2[,1],status='S',stringsAsFactors = FALSE)
# Remove nodes with missing ids
nodes2    <- nodes2[!is.na(nodes2$id),]

edgelist1B <- read.csv('Path1BE.csv',header=TRUE,stringsAsFactors = FALSE)[,2:3]

nodes1B = data.frame(id=sort(unique(c(edgelist1B$V1,edgelist1B$V2))),status='S',stringsAsFactors = FALSE)

# Empirical networks, built once
# simplify() merges repeated contacts into a single edge and removes self-loops
empirical_nets <- list(pathology1B = simplify(graph_from_data_frame(d=edgelist1B,vertices=nodes1B,directed=F)),
                       hschool     = simplify(graph_from_data_frame(d=edgelist,vertices=nodes,directed=F)),
                       hschool2    = simplify(graph_from_data_frame(d=edgelist2,vertices=nodes2,directed=F)))

# Random network: same size as the Pathology 1B class, with N*k/2 edges placed
# uniformly at random so the average degree is k
random_net <- function(k)
{
    N <- vcount(empirical_nets$pathology1B)
    sample_gnm(N,round(N*k/2))
}

# TRUE if x is a single number between lo and hi (and a whole number if whole=TRUE)
valid_input <- function(x,lo,hi,whole=FALSE)
{
    !is.null(x) && length(x)==1 && !is.na(x) && x >= lo && x <= hi && (!whole || x == round(x))
}

# Red message shown under a box holding an invalid value
range_msg <- function(ok,lo,hi,whole=FALSE)
{
    if(ok){return(NULL)}
    div(style='color:#c00; margin-top:-10px; margin-bottom:10px;',
        sprintf('Enter a %s between %s and %s. Using the last valid value.',
                if(whole) 'whole number' else 'value',lo,hi))
}


server <- function(input, output, session) {

    # beta, ti, intro and I0 hold the parameter values in use, which are only
    # updated from the boxes on Reset or a network change; n0 is the number
    # actually infected at the start of the current run; N is the network size;
    # network is the currently selected network
    rv <- reactiveValues(net=empirical_nets$pathology1B,
                         coord_fix=NULL,
                         sim_snap=list(),time=0,
                         toggle=TRUE,
                         beta=0.04,ti=10,intro=0.5,I0=5,n0=5,
                         N=vcount(empirical_nets$pathology1B),network=NULL)

    compute_layout <- function()
    {
        net <- isolate(rv$net)
        rv$coord_fix=as.matrix(layout_nicely(net))
    }

    # Copy valid box values into rv; invalid ones keep the value in use
    apply_params <- function()
    {
        isolate({
            if(valid_input(input$I0,1,rv$N,whole=TRUE)){rv$I0 = input$I0}
            if(valid_input(input$beta,0,2)){rv$beta = input$beta}
            if(valid_input(input$ti,0.5,365)){rv$ti = input$ti}
            if(valid_input(input$intro,0,2)){rv$intro = input$intro}
        })
    }

    # Clear the network, infect I0 random people and restart the clock
    reset_sim <- function()
    {
        rv$n0 = min(rv$I0,vcount(rv$net))
        V(rv$net)$status = 'S'
        i0 <- sample(1:length(V(rv$net)),rv$n0)
        V(rv$net)$status[i0] = 'I'
        rv$time = 0
        rv$sim_snap = list()
    }

    # Well-mixed solution up to tmax. The solution so far is cached and only
    # the new time steps are integrated, so each plot update stays cheap; the
    # cache is rebuilt when any of the model parameters change
    det_cache <- NULL
    det_key <- NULL
    well_mixed <- function(beta,k,g,N,n0,tmax)
    {
        key <- c(beta,k,g,N,n0)
        if(!identical(key,det_key)){
            det_cache <<- mass_actionD(beta,k,g,N,n0,0)
            det_key <<- key
        }
        last <- det_cache[nrow(det_cache),]
        if(tmax > last[1]){
            ext <- mass_actionD(beta,k,g,N,last[3],seq(last[1],tmax))
            det_cache <<- rbind(det_cache,ext[-1,,drop=FALSE])
        }
        det_cache[det_cache[,1] <= tmax,,drop=FALSE]
    }

    # toggle is TRUE when stopped; Stop/Go flips it
    observeEvent(input$simulate,{
        rv$toggle = !rv$toggle
        })

    # Single observer that steps while running, pausing 100ms between steps.
    # It has a lower priority than the plots, so after each step it reruns
    # (because rv$time changed) only once they have been redrawn, and starts
    # the pause then. Starting the pause before redrawing would leave the next
    # step already due whenever a step plus redraw takes over 100ms, and in the
    # browser (shinylive) the app then never gets to handle clicks or inputs
    last_time <- -1
    observe({
        if(rv$toggle){return()}
        if(rv$time == last_time){
            # Woken by the pause ending: take a step
            isolate({
                rv$net <- gillespie_stepSIS(rv$net,rv$beta,1.0/rv$ti,1.0,rv$intro)
                rv$time=rv$time+1
                rv$sim_snap[[rv$time]] = V(rv$net)$status=='I'
            })
        }else{
            # A step has just been taken and drawn (or Go pressed): pause
            last_time <<- rv$time
            invalidateLater(100)
        }
    }, priority = -100)

    observeEvent(input$step,{
        rv$net <- gillespie_stepSIS(rv$net,rv$beta,1.0/rv$ti,1.0,rv$intro)
        rv$time=rv$time+1
        rv$sim_snap[[rv$time]] = V(rv$net)$status=='I'
    })

    observeEvent(input$reset,{
        apply_params()
        reset_sim()
        compute_layout()
    })

    # Degree box: editable for the random network, fixed at the average degree otherwise
    output$degreeUI <- renderUI({
        if(input$network == 'mass_action'){
            k <- isolate(input$degree)
            if(is.null(k)){k <- 6}
            numericInput('degree',
                         'Network Degree',
                         value = k,
                         min = 2,
                         max = 12,
                         step = 1)
        }else{
            box <- numericInput('degree_fixed',
                                'Network Degree (average)',
                                value = round(mean(degree(empirical_nets[[input$network]])),2))
            htmltools::tagQuery(box)$find('input')$addAttrs(disabled=NA)$allTags()
        }
    })

    observeEvent(c(input$network,input$degree),{
        if(input$network == 'mass_action'){
            k <- input$degree
            # Wait for the degree box to appear, and ignore values outside 2-12
            # (e.g. while typing), keeping the current network
            if(is.null(k) || is.na(k) || k < 2 || k > 12){return()}
            rv$net = random_net(k)
        }else{
            rv$net = empirical_nets[[input$network]]
        }
        # A new network starts a new run, so pending parameter changes apply too.
        # Initial Infected goes back to the default of 5 on every network switch
        # (but not when only the random network's degree changes), and can be
        # at most the network size
        rv$N = vcount(rv$net)
        apply_params()
        if(!identical(input$network,rv$network)){
            rv$network = input$network
            rv$I0 = 5
            updateNumericInput(session,'I0',value=5,max=rv$N)
        }
        reset_sim()
        compute_layout()
    })

    # Parameter boxes don't affect the current run: valid changes are pending
    # until Reset (or a network change) applies them
    output$pending_msg <- renderUI({
        changed <- function(x,lo,hi,current,whole=FALSE){valid_input(x,lo,hi,whole) && x != current}
        if(changed(input$I0,1,rv$N,rv$I0,whole=TRUE) ||
           changed(input$beta,0,2,rv$beta) ||
           changed(input$ti,0.5,365,rv$ti) ||
           changed(input$intro,0,2,rv$intro)){
            div(style='color:#06c; margin-top:8px;',
                'Parameter changes pending: press Reset to apply them.')
        }
    })

    output$I0_msg   <- renderUI(range_msg(valid_input(input$I0,1,rv$N,whole=TRUE),1,rv$N,whole=TRUE))
    output$beta_msg <- renderUI(range_msg(valid_input(input$beta,0,2),0,2))
    output$ti_msg   <- renderUI(range_msg(valid_input(input$ti,0.5,365),0.5,365))
    output$intro_msg <- renderUI(range_msg(valid_input(input$intro,0,2),0,2))
    output$degree_msg <- renderUI({
        if(input$network == 'mass_action' && !is.null(input$degree)){
            range_msg(valid_input(input$degree,2,12),2,12)
        }
    })

    output$netPlot <- renderPlot({
        # Drawn with vectorised base graphics rather than plot.igraph, which is
        # several times slower and limits the Stop/Go frame rate. Same picture:
        # grey by default, red for infected people and I-I edges, blue for S-I edges
        infected <- V(rv$net)$status=='I'
        el <- as_edgelist(rv$net,names=FALSE)
        n_inf <- infected[el[,1]] + infected[el[,2]]
        edge_col <- c('grey','blue','red')[n_inf+1]
        xy <- norm_coords(rv$coord_fix)
        # Node radius is a fraction of the plot width, so it already scales with
        # the plot; shrink it for larger networks to limit overlap.
        # Edge and outline widths are in pixels, so scale them with the plot
        # (1 at 800px, never below 0.5)
        radius <- 4*sqrt(300/vcount(rv$net))/200
        plot_px <- min(session$clientData$output_netPlot_width,session$clientData$output_netPlot_height)
        lwd <- max(0.5,plot_px/800)
        plot(0,0,type='n',xlab='',ylab='',xlim=c(-1,1),ylim=c(-1,1),axes=FALSE,frame.plot=FALSE,asp=1)
        segments(xy[el[,1],1],xy[el[,1],2],xy[el[,2],1],xy[el[,2],2],col=edge_col,lwd=lwd)
        symbols(xy[,1],xy[,2],circles=rep(radius,nrow(xy)),inches=FALSE,add=TRUE,
                bg=ifelse(infected,'red','grey'),fg='black',lwd=lwd)
        })

    output$timeseries <- renderPlot({
        if(length(rv$sim_snap) > 0){
        N <- length(V(rv$net))
        k <- mean(degree(rv$net))
        sim_prev <- 100*sapply(rv$sim_snap,sum)/N
        # Well-mixed model with the same R0 = beta * k * infectious period
        R0 <- rv$beta*k*rv$ti
        detSol <- well_mixed(rv$beta,k,1.0/rv$ti,N,rv$n0,rv$time)
        det_prev <- 100*detSol[,3]/N
        plot(seq(1,rv$time,1),sim_prev,xlab='Time',ylab='% Prevalence',xlim=c(0,rv$time),ylim=c(0,max(c(sim_prev,det_prev))),pch=19,
             main=sprintf('R0 = %.2f (mean degree %.2f)',R0,k))
        lines(detSol[,1],det_prev,col='red',lwd=2)
        legend('bottomright',c('Network simulation','Well-mixed model'),pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c('black','red'),bty='n')

        }else{plot(0,0,xlab='Time',ylab='% Prevalence',pch='',xlim=c(0,1),ylim=c(0,100))}
        })
    }


# Run the application
shinyApp(ui = ui, server = server)
