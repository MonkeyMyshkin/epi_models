library(shiny)
library(SimInf)
library(RColorBrewer)

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

mass_actionD <- function(beta,g,S0,I0,R0,time_vec) 
{
    SI.dyn <- function(t,var,par) {
        S <- var[1]
        I <- var[2]
        R <- var[3]
        
        N <- S + I + R
        
        beta <- par[1]
        g    <- par[2]
        
        dS <- - beta*S*I/N
        dI <- beta*S*I/N - g*I
        dR <- g*I 
        return(list(c(dS,dI,dR)))
    }
    
    SI.par<-c(beta,g)
    SI.init<-c(S0,I0,R0)
    
    oot<-rk4(SI.init,time_vec,SI.dyn,SI.par)
    
    return(oot)
}

# Colours for S, I and R, shared by the checkbox labels and the plot
sir_cols <- brewer.pal(3,'Set1')

ui <- fluidPage(
    
    # Application title
    titlePanel("SIR Model"),
    
# Sidebar: final size histogram and model inputs
sidebarLayout(
    sidebarPanel(
        fluidRow(plotOutput("histogram"),width=6),
        fluidRow(column(width=4,numericInput("beta",
                    'Beta',
                    min = 0,
                    max = 5,
                    value = 1.5,
                    step = 0.05),
                    numericInput('ti',
                    'TI',
                    min=0.5,
                    max=30,
                    value=1,
                    step=0.1)),
        column(width=4,
        numericInput('S0',
                     'Initial S',
                     value = 99,
                     min = 0,
                     max = 1000000),
        numericInput('I0',
                     'Initial I',
                     value = 1,
                     min = 1,
                     max = 1000000),
        numericInput('R0',
                     'Initial R',
                     value = 0,
                     min = 0,
                     max = 1000000)),
        column(width=4,
        numericInput('t_run',
                     'Time',
                     value = 30,
                     min = 1,
                     max = 1000),
        numericInput('reps',
                     'Replicates',
                     value = 0,
                     min = 0,
                     max = 100),
        checkboxGroupInput('plot_these',
                           "Plot",
                           choiceNames = lapply(1:3, function(i){
                               span(c('S','I','R')[i], style=paste0('color:',sir_cols[i],';font-weight:bold'))}),
                           choiceValues = c('S','I','R'),
                           selected = c('I')),
        actionButton('update',
                     'Update')),width=6),

        uiOutput("downloadUI")
    ),
    # Time series of S, I and R (deterministic and stochastic)
    mainPanel(
        plotOutput("tsPlot",hover = "plot_hover"),textOutput("info")
    )
)
)


server <- function(input, output) {

cb_safepal <- col2rgb(brewer.pal(4,'Set1'))/255    
    
rv <- reactiveValues(detSol = mass_actionD(1.5,1.0/1.0,99,1,0,seq(0,30)), stoch_out=NULL)

detSol <- observeEvent(input$update,{
    rv$detSol = mass_actionD(input$beta,1.0/input$ti,input$S0,input$I0,input$R0,seq(0,input$t_run)) 
    if(input$reps>0)
    {
    initial_state  	<- data.frame(S=rep(input$S0,input$reps),I=rep(input$I0,input$reps),R=rep(input$R0,input$reps))
    # SimInf's first recorded point is already one time unit into the run, so record
    # from t=1 and prepend the initial state at t=0 ourselves
    tspan         		<- seq(from = 1, to=input$t_run, by = 1)

    # SimInf's built-in SIR model has the transitions
    #   S -> beta*S*I/(S+I+R) -> I
    #   I -> gamma*I -> R
    # It is precompiled, so unlike mparse() it needs no C compiler (and runs in the browser)
    model  <- SIR(u0 = initial_state, tspan = tspan, beta = as.double(input$beta), gamma = 1.0/input$ti)
    traj <- trajectory(run(model = model))
    traj0 <- data.frame(node=1:input$reps, time=0, initial_state)
    traj <- rbind(traj0, traj)
    rv$stoch_out <- traj[order(traj$node, traj$time),]
    }else{rv$stoch_out <- NULL}
})

# Download the stochastic sims if there are any, otherwise the deterministic solution
output$downloadUI <- renderUI({
    downloadButton("downloadData",
                   if(is.null(rv$stoch_out)) "Download Deterministic Solution" else "Download Stochastic Sims")
})

output$downloadData <- downloadHandler(
    filename = function() {
        paste(if(is.null(rv$stoch_out)) 'deterministic-' else 'stochastic-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
        if(is.null(rv$stoch_out))
        {
            det_out <- as.data.frame(rv$detSol[,1:4])
            names(det_out) <- c('time','S','I','R')
            write.csv(det_out, con, row.names=FALSE)
        }else{
            write.csv(rv$stoch_out, con, row.names=FALSE)
        }
    }
)

output$tsPlot <- renderPlot({
    
    if(!is.null(input$plot_these))
    {
    cols = numeric()
    if(is.element('S',input$plot_these)){cols = c(cols,2)}
    if(is.element('I',input$plot_these)){cols = c(cols,3)}
    if(is.element('R',input$plot_these)){cols = c(cols,4)}
    
    # Determinstic solution
    plot(0,0,pch=' ',xlim=range(rv$detSol[,1]),ylim=range(c(rv$detSol[,cols],rv$stoch_out[,1+cols])),xlab='Time',ylab='')
    if(is.element('S',input$plot_these)){lines(rv$detSol[,1],rv$detSol[,2],
                                               col=rgb(cb_safepal[1,1],cb_safepal[2,1],cb_safepal[3,1]),lwd=4)}
    if(is.element('I',input$plot_these)){lines(rv$detSol[,1],rv$detSol[,3],
                                               col=rgb(cb_safepal[1,2],cb_safepal[2,2],cb_safepal[3,2]),lwd=4)}
    if(is.element('R',input$plot_these)){lines(rv$detSol[,1],rv$detSol[,4],
                                               col=rgb(cb_safepal[1,3],cb_safepal[2,3],cb_safepal[3,3]),lwd=4)}
    grid()   
    # Add stochastic replicates. stoch_out is sorted by replicate then time, so each
    # compartment reshapes to a matrix with one column per replicate, drawn in one go
    if(!is.null(rv$stoch_out))
    {
    n_reps <- max(rv$stoch_out$node)
    times  <- rv$stoch_out$time[rv$stoch_out$node==1]
    stoch_cols <- c(rgb(cb_safepal[1,1],cb_safepal[2,1],cb_safepal[3,1],0.2),
                   rgb(cb_safepal[1,2],cb_safepal[2,2],cb_safepal[3,2],0.2),
                   rgb(cb_safepal[1,3],cb_safepal[2,3],cb_safepal[3,3],0.2))
    for(c in cols)
    {
        matlines(times,matrix(rv$stoch_out[[1+c]],ncol=n_reps),col=stoch_cols[c-1],type='l',lty=1)
    }
    }

    # Legend: one entry per compartment shown, plus line styles if stochastic runs are shown
    leg_text <- c('S','I','R')[cols-1]
    leg_col  <- sir_cols[cols-1]
    leg_lwd  <- rep(4,length(cols))
    if(!is.null(rv$stoch_out))
    {
        leg_text <- c(leg_text,'Deterministic','Stochastic run')
        leg_col  <- c(leg_col,'black',rgb(0,0,0,0.3))
        leg_lwd  <- c(leg_lwd,4,1)
    }
    legend('top',inset=c(0,-0.12),xpd=TRUE,horiz=TRUE,bty='n',text.width=NA,legend=leg_text,col=leg_col,lwd=leg_lwd)

    }
})

output$histogram <- renderPlot({
    if(!is.null(rv$stoch_out))
    {
    if(max(rv$stoch_out$node)>1)
    {
     # Final size: R at the end minus R at the start, one column per replicate
     R <- matrix(rv$stoch_out$R,ncol=max(rv$stoch_out$node))
     hist(R[nrow(R),] - R[1,],
              col='black',xlab='Cases',main='Final Size Distribution')  
    grid()
    }
    }
    
    
})

output$info <- renderText({
    
    # Read-out of the cursor position; blank when the cursor is off the plot
    e <- input$plot_hover
    if(is.null(e)) return("")
    paste0("Time = ", round(e$x, 1), ", Number = ", round(e$y, 1))

})





}

# Run the application 
shinyApp(ui = ui, server = server)
