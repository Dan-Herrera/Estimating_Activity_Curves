rm(list=ls())
source("./code/activityFunctions.R") #source functions (rename this eventually)

data <- UtahData() #load data

#format observed activity density datasets for all species
cottontail <- format_observedPDF(df = data,
                                 sp.col = "species",
                                 sp.name = "Mountain Cottontail",
                                 time.start = "18:00",
                                 time.stop = "12:00",
                                 time.col = "obs.time")
fox <- format_observedPDF(df = data,
                          sp.col = "species",
                          sp.name = "Red Fox",
                          time.start = "18:00",
                          time.stop = "12:00",
                          time.col = "obs.time")
coyote <- format_observedPDF(df = data,
                             sp.col = "species",
                             sp.name = "Coyote",
                             time.start = "18:00",
                             time.stop = "12:00",
                             time.col = "obs.time")

#format potential ideal curves
ideal.list <- list()

#left-beta
ideal.list[[1]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "beta",
                                   param1 = 2,
                                   param2 = 3)

#left-kumaraswamy
ideal.list[[2]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "kumaraswamy",
                                   param1 = 2,
                                   param2 = 5)

#normal
ideal.list[[3]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "normal",
                                   param1 = 0.4,
                                   param2 = 1.5) #WARNING MESSAGE IS EXPECTED

#triangular
ideal.list[[4]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 1,
                                   param3 = 0.2)

#uniform
ideal.list[[5]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "uniform")

#bimodal
ideal.list[[6]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "bimodal",
                                   param1 = 0.2,
                                   param2 = 0.6,
                                   param3 = 0.8)

#format predator list
pred.list <- list()
pred.list[[1]] <- fox
#pred.list[[1]] <- coyote #alternate bewteen species (in the future, it would be nice to have both in the same list, but not possible now)

results.beta <- estimate_activity(idealPDF = ideal.list[[1]],
                                  disturbance = pred.list,
                                  observedPDF = cottontail,
                                  confInt = TRUE,
                                  nIter = 1000,
                                  AIC = TRUE)

results.kumaraswamy <- estimate_activity(idealPDF = ideal.list[[2]],
                                         disturbance = pred.list,
                                         observedPDF = cottontail,
                                         confInt = TRUE,
                                         nIter = 1000,
                                         AIC = TRUE)

results.normal <- estimate_activity(idealPDF = ideal.list[[3]],
                                    disturbance = pred.list,
                                    observedPDF = cottontail,
                                    confInt = TRUE,
                                    nIter = 1000,
                                    AIC = TRUE)

results.triangular <- estimate_activity(idealPDF = ideal.list[[4]],
                                        disturbance = pred.list,
                                        observedPDF = cottontail,
                                        confInt = TRUE,
                                        nIter = 1000,
                                        save.boots = TRUE,
                                        AIC = TRUE)

results.uniform <- estimate_activity(idealPDF = ideal.list[[5]],
                                     disturbance = pred.list,
                                     observedPDF = cottontail,
                                     confInt = TRUE,
                                     nIter = 1000,
                                     AIC = TRUE)

results.bimodal <- estimate_activity(idealPDF = ideal.list[[6]],
                                     disturbance = pred.list,
                                     observedPDF = cottontail,
                                     confInt = TRUE,
                                     nIter = 1000,
                                     AIC = TRUE)

sens.list <- list()
sens.list[[1]] <- 0.10

pred.prediction <- predict_activity(idealPDF = ideal.list[[4]],
                                           disturbance = pred.list,
                                           sensitivity = sens.list)

################## testing variation within a single distribution

pred.list <- list()
pred.list[[1]] <- fox

#format potential ideal curves
ideal.list <- list()

#position 1
ideal.list[[1]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 0.8,
                                   param3 = 0.1)

#position 2
ideal.list[[2]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 0.8,
                                   param3 = 0.2)

#position 3
ideal.list[[3]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 0.8,
                                   param3 = 0.3)


#position 4
ideal.list[[4]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 0.9,
                                   param3 = 0.1)

#position 5
ideal.list[[5]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 0.9,
                                   param3 = 0.2)

#position 6
ideal.list[[6]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 0.9,
                                   param3 = 0.3)

#position 7
ideal.list[[7]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 1,
                                   param3 = 0.1)

#position 8
ideal.list[[8]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 1,
                                   param3 = 0.2)

#position 9
ideal.list[[9]] <- format_idealPDF(time.start = "18:00",
                                   time.stop = "12:00",
                                   observedPDF = cottontail,
                                   dist = "triangular",
                                   param1 = 0,
                                   param2 = 1,
                                   param3 = 0.3)

#estimate s for each competing distribution
results.p1 <- estimate_activity(idealPDF = ideal.list[[1]],
                                disturbance = pred.list,
                                observedPDF = cottontail,
                                confInt = TRUE,
                                nIter = 1000)

results.p2 <- estimate_activity(idealPDF = ideal.list[[2]],
                                disturbance = pred.list,
                                observedPDF = cottontail,
                                confInt = TRUE,
                                nIter = 1000)

results.p3 <- estimate_activity(idealPDF = ideal.list[[3]],
                                disturbance = pred.list,
                                observedPDF = cottontail,
                                confInt = TRUE,
                                nIter = 1000)

results.p4 <- estimate_activity(idealPDF = ideal.list[[4]],
                                disturbance = pred.list,
                                observedPDF = cottontail,
                                confInt = TRUE,
                                nIter = 1000)

results.p5 <- estimate_activity(idealPDF = ideal.list[[5]],
                                disturbance = pred.list,
                                observedPDF = cottontail,
                                confInt = TRUE,
                                nIter = 1000)

results.p6 <- estimate_activity(idealPDF = ideal.list[[6]],
                                disturbance = pred.list,
                                observedPDF = cottontail,
                                confInt = TRUE,
                                nIter = 1000)

results.p7 <- estimate_activity(idealPDF = ideal.list[[7]],
                                disturbance = pred.list,
                                observedPDF = cottontail,
                                confInt = TRUE,
                                nIter = 1000)

results.p8 <- estimate_activity(idealPDF = ideal.list[[8]],
                                disturbance = pred.list,
                                observedPDF = cottontail,
                                confInt = TRUE,
                                nIter = 1000)

results.p9 <- estimate_activity(idealPDF = ideal.list[[9]],
                                disturbance = pred.list,
                                observedPDF = cottontail,
                                confInt = TRUE,
                                nIter = 1000)

#combine results
p.results <- rbind(results.p1,
                   results.p2,
                   results.p3,
                   results.p4,
                   results.p5,
                   results.p6,
                   results.p7,
                   results.p8,
                   results.p9)

#add information about the distribution positon
p.results$position <- c(1,2,3,4,5,6,7,8,9)

#rearrange according to neg log liklihood
sort(p.results$negLL)
