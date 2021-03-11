library(raster)
dir = "/home/liorgr/PycharmProjects/land_mammal_biomass"
setwd(dir)
habitat_raster = raster::raster("iucn_habitatclassification_composite_lvl2_ver001.tif")
# 100  101  102  103  104  105  106  107  108  109  201  202  300   500  502  503  
# 504  505  506  507  510  511  513  514  
# # 0 
# 
# 
# 
##

maskvalue = c(0, 301,  302,  303,  304,  305,  306,  307 , 308,
              400, 401, 402, 403, 404, 405, 406, 407, 
              600, 801, 802, 803, 140, 1402, 1403, 1405)


plot(habitat_raster)
