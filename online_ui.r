library(shiny)
library(flexdashboard)
library(plotrix)
library(matrixStats)

pred_id = c(208100, 208575, 209675, 327375, 498400, 
            509800, 510100, 544200, 544450, 820775)
panel_dim = c(4,1)
inset_dim = c(0,-.18)

ui <- navbarPage( title = "Internal Bleeding Investigation",
                  header = selectInput("subject",
                                       label = "Choose a patient: ",
                                       choices = pred_id,
                                       selected = 327375),
                  tabPanel(title = "EKG Analysis",
                           # Sidebar layout with input and output definitions ----
                           column(6,
                                  fluidRow(splitLayout(cellWidths = c("100%"), plotOutput("vitals", height = "25%")))
                           ),
                           column(6,
                                  fluidRow(splitLayout(cellWidths = c("100%"), plotOutput("probs", height = "25%")))
                           )
                  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    output$vitals <- renderPlot({
        
        # Look for new data every 5 seconds
        invalidateLater(5000,session)
        load('Model_out/info_list.rda')
        
        par(mfrow = c(2,1), mar = c(2, 2.1, 2.1, 1), bg='black', fg='green')
        
        i = input$subject
        ind_i = which(pred_id == i)
        
        indices_i = (info_list$use_data[,'EID']==i)
        n_i = sum(indices_i)
        
        t_grid = round(info_list$use_data[indices_i, 'time'] / 60, digits = 3)
        t_grid_bar = 1:length(t_grid)
        rbc_times_bar = which(info_list$use_data[info_list$use_data[,'EID']==i, 'RBC_ordered'] != 0)
        
        rbc_admin_times_bar = which(info_list$use_data[info_list$use_data[,'EID']==i, 'RBC_admin'] != 0) 
        
        rbc_times = t_grid[rbc_times_bar]
        rbc_admin_times = t_grid[rbc_admin_times_bar]
        
        curr_post_prob = rbind(colMeans(info_list$B_chain[, indices_i] == 1),
                               colMeans(info_list$B_chain[, indices_i] == 2),
                               colMeans(info_list$B_chain[, indices_i] == 3),
                               colMeans(info_list$B_chain[, indices_i] == 4),
                               colMeans(info_list$B_chain[, indices_i] == 5))
        prev_post_prob = info_list$prev_state_avg[2:6, info_list$prev_state_avg[1,] == i]
        
        if(ncol(prev_post_prob) == 48) {
            if(sum(info_list$ind_list[[ind_i]] %in% info_list$prev_ind_list[[ind_i]]) != length(info_list$ind_list[[ind_i]])) {
                prev_post_prob = prev_post_prob[,-1]
            }
        }
        curr_post_prob[,1:ncol(prev_post_prob)] = 0.5 * (curr_post_prob[,1:ncol(prev_post_prob)] + prev_post_prob)  
        
        pb = barplot(curr_post_prob, 
                     col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                     xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
        
        
        # Heart Rate and MAP double plot ---------------------------------------
        if(mean(info_list$use_data[indices_i, 'clinic_rule']) != 0) {
            title_name = paste0('HR & MAP: ', i, ', RBC = ', 
                                info_list$RBC_rule_vec[ind_i], ', clinic = ', 
                                mean(info_list$use_data[indices_i, 'clinic_rule']))
        } else {
            title_name = paste0('HR & MAP: ', i, ', RBC = ', 
                                info_list$RBC_rule_vec[ind_i])
        }
        
        hr_upper = colQuantiles( info_list$Hr_chain[, indices_i, drop=F], probs=.975)
        hr_lower = colQuantiles( info_list$Hr_chain[, indices_i, drop=F], probs=.025)
        bp_upper = colQuantiles( info_list$Map_chain[, indices_i, drop=F], probs=.975)
        bp_lower = colQuantiles( info_list$Map_chain[, indices_i, drop=F], probs=.025)
        
        hr_map_ylim = c(min(hr_lower, bp_lower), max(hr_upper, bp_upper))
        
        # Make a new plot to add the background color
    
        plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
             xlab='time', ylab=NA, xaxt='n', col.main='green',
             col.axis='green')
        
        plotCI( x = pb, y=colMeans(info_list$Hr_chain[, indices_i, drop=F]), 
                ui=hr_upper, li=hr_lower,
                main=title_name,xlab='time', ylab=NA, xaxt='n', col.main='green',
                col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
                xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim, add =T) 
        plotCI( x = pb, y=colMeans(info_list$Map_chain[, indices_i, drop=F]), 
                ui=bp_upper, li=bp_lower, main=title_name,
                xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
                col = 'orange',
                xlim = range(pb) + c(-0.5,0.5), add = T) 
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'HR', 'MAP'), pch=15, pt.cex=1.5, 
                col=c( 'aquamarine', 'orange'))
        
        grid(nx = NA, ny = NULL, col = "white")
        axis( side=1, at=pb, col.axis='green', labels=t_grid)
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        
        # Hemoglobin and Lactate double plot -----------------------------------
        if(mean(info_list$use_data[indices_i, 'clinic_rule']) != 0) {
            title_name = paste0('Hemo. & Lact: ', i, ', RBC = ', 
                                info_list$RBC_rule_vec[ind_i], ', clinic = ',
                                mean(info_list$use_data[indices_i, 'clinic_rule']))
        } else {
            title_name = paste0('Hemo. & Lact: ', i, ', RBC = ',
                                info_list$RBC_rule_vec[ind_i])
        }
        
        hc_upper = colQuantiles( info_list$Hc_chain[, indices_i, drop=F], probs=.975)
        hc_lower = colQuantiles( info_list$Hc_chain[, indices_i, drop=F], probs=.025)
        la_upper = colQuantiles( info_list$La_chain[, indices_i, drop=F], probs=.975)
        la_lower = colQuantiles( info_list$La_chain[, indices_i, drop=F], probs=.025)
        
        hr_map_ylim = c(min(hc_lower, la_lower), max(hc_upper, la_upper))
        
        plot(NULL, xlim=range(pb) + c(-0.5,0.5), ylim=hr_map_ylim, main=title_name,
             xlab='time', ylab=NA, xaxt='n', col.main='green',
             col.axis='green')
        
        plotCI(x = pb, y = colMeans(info_list$Hc_chain[, indices_i, drop=F]), ui=hc_upper, li=hc_lower,
               main=title_name,
               xlab='time', ylab=NA, xaxt='n', col.main='green',
               col.axis='green', pch=20, cex=1, sfrac=.0025, col = 'aquamarine',
               xlim = range(pb) + c(-0.5,0.5), ylim = hr_map_ylim, add = T) 
        plotCI( x = pb, y=colMeans(info_list$La_chain[, indices_i, drop=F]), ui=la_upper, li=la_lower,
                main=title_name,
                xlab='time', ylab=NA, xaxt='n', pch=20, cex=1, sfrac=.0025,
                col = 'orange',
                xlim = range(pb) + c(-0.5,0.5), add = T) 
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'hemo', 'lactate'), pch=15, pt.cex=1.5, 
                col=c( 'aquamarine', 'orange'))
        # grid( nx=20, NULL, col='white')
        grid(nx = NA, ny = NULL, col = "white")
        axis( side=1, at=pb, col.axis='green', labels=t_grid)
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
    }, height = 420, width = 700)
    
    output$probs <- renderPlot({
        
        # Look for new data every 5 seconds
        invalidateLater(5000,session)
        load('Model_out/info_list.rda')
        
        par(mfrow = c(2,1), mar = c(2, 2.1, 2.1, 1), bg='black', fg='green')
        i = input$subject
        ind_i = which(pred_id == i)
        
        indices_i = (info_list$use_data[,'EID']==i)
        n_i = sum(indices_i)
        
        t_grid = round(info_list$use_data[indices_i, 'time'] / 60, digits = 3)
        t_grid_bar = 1:length(t_grid)
        rbc_times_bar = which(info_list$use_data[info_list$use_data[,'EID']==i, 'RBC_ordered'] != 0)
        
        rbc_admin_times_bar = which(info_list$use_data[info_list$use_data[,'EID']==i, 'RBC_admin'] != 0) 
        
        rbc_times = t_grid[rbc_times_bar]
        rbc_admin_times = t_grid[rbc_admin_times_bar]
        
        curr_post_prob = rbind(colMeans(info_list$B_chain[, indices_i] == 1),
                               colMeans(info_list$B_chain[, indices_i] == 2),
                               colMeans(info_list$B_chain[, indices_i] == 3),
                               colMeans(info_list$B_chain[, indices_i] == 4),
                               colMeans(info_list$B_chain[, indices_i] == 5))
        prev_post_prob = info_list$prev_state_avg[2:6, info_list$prev_state_avg[1,] == i]
        
        if(ncol(prev_post_prob) == 48) {
            if(sum(info_list$ind_list[[ind_i]] %in% info_list$prev_ind_list[[ind_i]]) != length(info_list$ind_list[[ind_i]])) {
                prev_post_prob = prev_post_prob[,-1]
            }
        }
        curr_post_prob[,1:ncol(prev_post_prob)] = 0.5 * (curr_post_prob[,1:ncol(prev_post_prob)] + prev_post_prob)  
        
        pb = barplot(curr_post_prob, 
                     col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                     xlab='time', space=0, col.main='green', border=NA, axes = F, plot = F) 
        
        
        # BAR PLOTS --------------------------------------------------------------
        barplot(curr_post_prob, 
                col=c( 'dodgerblue', 'firebrick1', 'yellow2', 'green', 'darkgray'), 
                xlab='time', space=0, col.main='green', border=NA,
                xlim=range(pb) + c(-0.5,0.5)) 
        # grid( nx=20, NULL, col='white')
        grid(nx = NA, ny = NULL, col = "white")
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'Baseline', 'State 2', 'State 3', 'State 4', 'State 5'), 
                pch=15, pt.cex=1.5, 
                col=c( 'dodgerblue', 'firebrick1', 'yellow2','green', 'darkgray'))
        legend( 'topleft', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c( 'RBC order', 'RBC admin'), pch=15, pt.cex=1.5, 
                col=c( 'darkorchid1', 'aquamarine'))				
        axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
        axis( side=2, at=0:1, col.axis='green')
        
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        
        # Cumulative PLOTS ---------------------------------------------------------
        cumulative_post_prob = matrix(nrow = 2, ncol = n_i)
        ind = 1
        win_length = 0
        c = 0.257
        
        indices_i_new = which(indices_i == T)
        for(w in 1:length(indices_i_new)) {
            start_index = indices_i_new[1]
            end_index = indices_i_new[w] 
            if(w - win_length > 0) start_index = indices_i_new[w - win_length]
            
            y_or_n_2 = apply(info_list$B_chain[, start_index:end_index, drop=F],
                             1, function(x) (2 %in% x))
            prob_2 = mean(y_or_n_2)
            
            cumulative_post_prob[, ind] = c(prob_2, 1-prob_2)
            ind = ind + 1
        }
        
        barplot( cumulative_post_prob,
                 col=c('darkred', 'black'),
                 main=paste0('cumulative prob.'), xlab='time', space=0, 
                 col.main='green', border=NA,
                 xlim=range(pb) + c(-0.5,0.5))
        # grid( nx=20, NULL, col='white')
        grid(nx = NA, ny = NULL, col = "white")
        legend( 'topright', inset=inset_dim, xpd=T, horiz=T, bty='n', x.intersp=.75,
                legend=c('', 'State 2'), pch=15, pt.cex=1.5,
                col=c('black', 'darkred'))
        axis( side=1, at=t_grid_bar-0.5, col.axis='green', labels = t_grid)
        axis( side=2, at=0:1, col.axis='green')
        
        abline(v = rbc_times_bar-0.5, col = 'darkorchid1', lwd = 1)
        abline(v = rbc_admin_times_bar-0.5, col = 'aquamarine', lwd = 1)
        abline(h = c, col = 'yellow', lwd = 2)
        
    }, height = 420, width = 700)
}

# Run the application 
shinyApp(ui = ui, server = server)