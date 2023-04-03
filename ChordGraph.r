library(circlize)
library(readxl)
ChordChartSampleChemoDownQval01 <- read_excel("E:/qvalues/chemokine_network_qval0.1_pval0.05.xlsx", sheet=2)
View(ChordChartSampleChemoDownQval01)
chordDiagram(ChordChartSampleChemoDownQval01, directional = 1,
             transparency = 0.5, grid.col=c("11 - Ccl5"="#9999FF","2 - Ccl2"="#009900","2 - Ccl4"="#009900",
                                         "2 - Ccl5"="#009900","2 - Ccl8"="#009900","2 - Ccl12"="#009900",
                                         "3 - Ccl3"="#66FFFF","3 - Ccl4"="#66FFFF","4 - Ccl3"="#994C00",
                                         "5 - Ccl7"= "#FF6666","6 - Ccl5"="#999900","9 - Ccl3"="#FF9933",
                                         "9 - Ccl4"="#FF9933","11 - Ccr2"="#9999FF","3 - Ccr5"="#66FFFF",
                                         "5 - Ccr5"= "#FF6666"))

ChordChartSampleCytoDownQval01 <- read_excel("E:/qvalues/cytokine_network_qval0.1_pval0.05.xlsx", sheet=2)
View(ChordChartSampleCytoDownQval01)
chordDiagram(ChordChartSampleCytoDownQval01, directional = 1,
             transparency = 0.5,grid.col=c("1 - Il1a"="#4C9900","1 - Il1b"="#4C9900","2 - Inhba"="#009900",
                                           "2 - Ifnb1"="#009900","2 - Tnf"="#009900","2 - Il1a"="#009900",
                                           "2 - Il1rn"="#009900","2 - Tnfsf13b"="#009900",
                                           "5 - Gdf11" = "#FF6666","5 - Il1b" = "#FF6666",
                                           "5 - Il1rn" = "#FF6666","6 - Il1b" = "#999900",
                                           "7 - Il1b" = "#BFBFBF","7 - Il1rn" = "#BFBFBF","9 - Tnfsf11" = "#FF9933",
                                           "9 - Cd40lg" = "#FF9933","9 - Tnf" = "#FF9933","9 - Il1b" = "#FF9933",
                                           "1 - Tgfbr1"="#4C9900","2 - Acvr2a"="#009900",
                                           "2 - Tnfrsf11a"="#009900","3 - Ifnar1" = "#66FFFF","3 - Tnfrsf1b" = "#66FFFF",
                                           "3 - Tgfbr1" = "#66FFFF","3 - Cd40" = "#66FFFF","3 - Tnfrsf11a" = "#66FFFF",
                                           "5 - Il1r2"= "#FF6666","5 - Acvr2a"= "#FF6666","5 - Tnfrsf13c"= "#FF6666",
                                           "5 - Tgfbr1"= "#FF6666","9 - Il1r1"= "#FF9933","9 - Ifnar1"= "#FF9933"),
             order=c("1 - Il1a","1 - Il1b","2 - Inhba",
                     "2 - Ifnb1","2 - Tnf","2 - Il1a",
                     "2 - Il1rn","2 - Tnfsf13b",
                     "5 - Il1b",
                     "5 - Il1rn","5 - Gdf11","6 - Il1b",
                     "7 - Il1b","7 - Il1rn","9 - Tnfsf11",
                     "9 - Cd40lg","9 - Tnf","9 - Il1b",
                     "1 - Tgfbr1","2 - Acvr2a",
                     "2 - Tnfrsf11a","3 - Ifnar1","3 - Tnfrsf1b",
                     "3 - Tgfbr1","3 - Cd40","3 - Tnfrsf11a",
                     "5 - Il1r2","5 - Acvr2a","5 - Tnfrsf13c",
                     "5 - Tgfbr1","9 - Il1r1","9 - Ifnar1"))            