for(run in 1:10){
  
  #Required packages
  
  #library(plotly)
  
  
  #Clean environment
  
  
  setwd("~/Documents/Cluster/")
  
  
  ##Set fixed parameters
  
  ## Density of animals (ind/km2) at carring capacity
  ## from Distance Sampling (1/16 to fit the cells in the simulation matrix) based on data from 7 localities
  
  abund.small <- 40/16
  abund.medium <- 1.7/16
  abund.monkeys <- 36/16
  abund.large <- 0.3/16
  
  ## Camaritagua Indigenous Reserve
  ### Load populations and carrying capacity (relative productivity of habitats)
  load(file = "Initial-Abundance-Camaritagua.Rdata")
  
  
  
  template.mtx <- matrix(0, nrow= nrow(Carr.capacity), ncol=ncol(Carr.capacity))
  cell.nums <- matrix(1:length(template.mtx), nrow(template.mtx), ncol(template.mtx),  byrow=T)
  #Load required functions
  
  source("~/Documents/Cluster//Hunting_Model_Functions_2-5.R")
  
  rodent.correction <- matrix(0, nrow(template.mtx), ncol(template.mtx), byrow = TRUE)
  
  for(rod in 1:length(near[,1])){
    rodent.correction[near[rod,1], near[rod,2]] <- (abund.small*0.75)
  }
  
  Carr.capacity.small <- Carr.capacity.small + rodent.correction
  
  
  
  ##Simulates hunting events
  people <- 0
  daily <- 1.93e-2
  village.location <- as.matrix(c(54,37))
  close <- 0.6
  far <- 0.4
  consumed <- 0
  Abundance.yearly <- NULL
  timeline.small <- NULL
  timeline.medium <- NULL
  timeline.large <- NULL
  timeline.monkeys <- NULL
  
  events <- NULL
  
  for(month in 1:720){
    #month <- period
    #print(c("month", month))
    if(month==120){
      people <- 20
    }
    if(month==480){
      daily <- (daily*2)
    }
    
    budget <- people*(30*daily)
    effort <- 0
    #consumed <- 0
    effort_limit <- (people/8)*10000
    
    #print(c("saved from last month", consumed))
    kills <- NULL
    while(consumed<budget){
      if(effort > effort_limit){
        break
      }
      
      trip.type <- sample(c(close, far), 1, prob = c(close, far))
      if (trip.type==close){
        start.cell <- village.location
        kill <- hunt.event.close(start.cell, 1500)
      }
      else{
        start.cell.random <- far.away[sample(1:nrow(far.away)),,drop=FALSE]
        start.cell.nu <- sample(1:length(start.cell.random[,2]), 1)
        start.cell <- start.cell.random[start.cell.nu,]
        start.cell <- as.matrix(c(start.cell[1], start.cell[2]))
        #print(c("s_cell", start.cell))
        distance <- sample(c(4000, 5000, 6000, 8000), 1, replace = T, prob = (c(0.32, 0.32, 0.32, 0.04)))
        kill <- hunt.event(start.cell, distance)
      }
      
      
      
      if (kill[4]!= "none"){
        #print(c("today's dinner is", kill[4]))
        #print(c(as.integer(kill[2]), as.integer(kill[3])))
        #dd <- c(as.integer(kill[2]), as.integer(kill[3]))
        #print(Curr.abundance.All[[kill[4]]][as.integer(kill[2]), as.integer(kill[3])])
        Curr.abundance.All[[kill[4]]][as.integer(kill[2]), as.integer(kill[3])] <- (Curr.abundance.All[[kill[4]]][as.integer(kill[2]), as.integer(kill[3])])-1
        #print(Curr.abundance.All[[kill[4]]][as.integer(kill[2]), as.integer(kill[3])])
      }
      kills <- rbind(kills, kill)
      consumed <-  consumed + as.numeric(kill[5])
    }
    
    events <- rbind(events, kills)
    
    if(month%%12 == 0){
      
      Curr.abundance.All$small <- popdynamics(Curr.abundance.All$small, Carr.capacity.small, 0.2, 2000)
      Curr.abundance.All$medium <- popdynamics(Curr.abundance.All$medium, Carr.capacity.medium, 0.5, 4000)
      Curr.abundance.All$large <- popdynamics(Curr.abundance.All$large, Carr.capacity.large, 0.5, 8000)
      Curr.abundance.All$monkeys <- popdynamics(Curr.abundance.All$monkeys, Carr.capacity.monkeys, 0.2, 2000)
      
      
      
      people <- people + round(people*0.05)
      #print(c("people", people))
    }
    if(is.element(-1, unique(Curr.abundance.All$small))){
      print("menos 1")
      break
    }
    
    timeline.small[[month]] <- Curr.abundance.All$small
    timeline.medium[[month]] <- Curr.abundance.All$medium
    timeline.large[[month]] <- Curr.abundance.All$large
    timeline.monkeys[[month]] <- Curr.abundance.All$monkeys
    Curr.abundance.All <-  movement.month(Curr.abundance.All, Carr.capacity.small, Carr.capacity.medium,
                                          Carr.capacity.large, Carr.capacity.monkeys)
    
    #Abundance.yearly[[month]] <- Curr.abundance.All
    consumed <- consumed-budget
    if(consumed<0){
      consumed=0
    }
    
    
  }
  
  run_name <- as.character(paste("~/Documents/Cluster/Results/Cam_Hg_", run, ".Rdata", sep = "")) 
  save(events, timeline.small, timeline.monkeys, timeline.medium, timeline.large, file = run_name)
  rm(list=ls())
  
  
  
}
