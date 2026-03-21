library(NicheMapR) # load the NicheMapR package
micro <- micro_terra(ystart = 2000,
                     yfinish = 2000,
                     clearsky = 0,
                     cap = 1,
                     snowmodel = 0,
                     elevation = 270,
                     runmoist = 0,
                     scenario = 0)
metout <- as.data.frame(micro$metout)
soil <- as.data.frame(micro$soil)
write.csv(metout, file = 'c:/git/BiophysicalGrids.jl/test/data/micro_terra/metout_monthly_terra.csv')
write.csv(soil, file = 'c:/git/BiophysicalGrids.jl/test/data/micro_terra/soil_monthly_terra.csv')
