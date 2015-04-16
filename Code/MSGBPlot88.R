# Script to Create MSGB Plot

# Load Libraries and Set Parameters ---------------------------------------
library(png)
library(grid)

outputPath = "Output"
codePath = "Code"
dataPath = "Data"

# Functions ---------------------------------------------------------------
BubblePlot = function(DoseInfo, SNPInfo, ProxyInfo, GeneMatch, GeneInfo, FID, Woman, Man, ShowGenes = TRUE, labelSize = .6){
    
  ### Running Requred Functions
  Labels = function(Genes, a, r, length, offset = .2){
    r = r + offset
    degrees = seq(0, 360, length.out = (length+1))[-(length+1)]
    for (i in 1:length){
      if (degrees[i] > 90 & degrees[i] < 270){
        text((r*cos(a))[i], (r*sin(a))[i], Genes[i], cex = labelSize, srt = degrees[i]+180, font = 1)
      }
      else{
        text((r*cos(a))[i], (r*sin(a))[i], Genes[i], cex = labelSize, srt = degrees[i], font = 1)
      }
    }
  }
  Distribution = function(distribution, MSrange){
    dens = density(distribution)
    dens$x = ((dens$x - dens$x[1]) / diff(range(dens$x))) * diff(MSrange) + MSrange[1]
    x = dens$x
    y = dens$y
    yMod = sqrt(y)
    yPercentiles = seq(0.05, 0.95, length.out = 17)
    yPerc = cumsum(yMod) / sum(yMod)
    xPoints = sapply(yPercentiles, function(yP_i) x[which(yPerc-yP_i>0)[1]])
    xPoints = jitter(xPoints, amount=0.02)
    return(xPoints)
  }
  rescale = function(x, y, a, b){
    min = min(x,y)
    x  = x - min
    y = y -  min
    max = max(x,y)
    x = x/max
    y = y/max
    x = a + (b-a)*x
    y = a + (b-a)*y
    return(list(x,y))
  }
  PlotPeople = function(image, x = 8, y=0, xlim = c(-.5,.5), ylim = c(-.5,.5), height = .07, color = "#444444BB"){
    x = grconvertX(x, from="user", to="ndc")
    y = grconvertY(y, from="user", to="ndc")
    image = as.raster(image)
    image[image!="#00000000"] = color
    grid.raster(x = x,y = y,image = image, height=height)
    par(new = T)
  }
  
  ### Checking to see if FID is valid
  N = which(DoseInfo$FID == FID)
  if (length(N) == 0){
    stop("Invalid FID")  
  }
  
  ### Setting up Colors to Use
  ColorSet1 = c("#D95F02","#E6AB02","#E7298A","#66A61E","#7570B3","#377EB8", "#A6761D", "#666666") #Functions
  ColorSet2 = c("#FD8D3C55","#F768A155","#225EA855","#74C47655","#810F7C55") #Tissue Labels
  ColorSet3 = c("#FD8D3C","#F768A1","#225EA8","#74C476","#810F7C") #Tissues
  ColorSet4 = c("#99D8C9BB", "#66C2A4BB") #Controls
  ColorSet5 = c("#FA9FB5BB", "#FF8888BB") #Cases
  ColorSet6 = "#444444" #Distribution Info
  
  ### Get Case and Control Data
  Control = DoseInfo$Score[DoseInfo$Pheno == "HC"]
  Case = DoseInfo$Score[DoseInfo$Pheno == "MS"]
  
  ### Get rsIDs found in DoseInfo and SNPInfo
  start = grep("rs", colnames(DoseInfo))[1]
  end = tail(grep("rs", colnames(DoseInfo)),1)
  length = length(start:end)
  SNPnames = colnames(DoseInfo)[start:end]
  index = SNPInfo$rsIDforJoseph %in% SNPnames
  SNPnames = SNPInfo$rsIDforJoseph[index]
  
  ### Get Genes, Functions, Tissues, OR, RAF from SNPInfo 
  Genes = sapply(strsplit(sapply(strsplit(gsub("^\\s+","",as.character(SNPInfo$Gene))," "), "[", 1), "\\("), "[", 1)[index]
  Tissues = as.factor(SNPInfo$Tissue)[index]
  Tissues = c("Brain", "B-Cell", "T-Cell", "B and T", "Other")
  # Tissues = sample(Tissues, length, replace = TRUE)  
  Tissues = as.character(GeneInfo$Tissues)
  Functions = as.factor(SNPInfo$Function)[index]
  FuncColors = ColorSet1[as.numeric(Functions)]
  OR = SNPInfo$OR[index]
  RAF = SNPInfo$RAF[index]  
  
  ### Getting Dose for each SNP (Original + Proxy)
  common = rep(0,nrow(SNPInfo))
  for (i in 1:nrow(SNPInfo)){
    cut = ProxyInfo[ProxyInfo$rsIDforJoseph == unique(ProxyInfo$rsIDforJoseph)[i],]
    index = cut$Proxy %in% colnames(GeneMatch)
    if(sum(index) == 0){
      common[i] = NA    
    }
    else{
      common[i] = as.character(cut$Proxy[which(cut$Order == min(cut$Order[index]))])
    }
  }
  Dose = round(as.numeric(unlist(GeneMatch[N, which(colnames(GeneMatch) %in% common)])))  
  
  ### Getting Info Unique to Individual from DoseInfo
  Sex = DoseInfo$SEX[N]
  MSGB = round(DoseInfo$Score[N],1)
  Pheno = DoseInfo$Pheno[N]
  
  ### Put Info in Data Frame and Sort by Tissue Type and by Increasing OR
  Data = data.frame(SNPnames, Genes, Tissues, Functions, FuncColors, OR, RAF, Dose)
  Data$Genes = as.character(Genes)
  empty = Data[1,]
  empty$OR = NA
  empty$Genes = ""
  empty$Dose = NA
  d1 = Data[which(Data$Tissues == "Brain"),]
  d1 = d1[order(d1$OR),]
  d1[nrow(d1)+1, ] = empty
  d2 = Data[which(Data$Tissues == "B-Cell"),]
  d2 = d2[order(d2$OR),]
  d2[nrow(d2)+1, ] = empty
  d3 = Data[which(Data$Tissues == "B and T"),]
  d3 = d3[order(d3$OR),]
  d3[nrow(d3)+1, ] = empty
  d4 = Data[which(Data$Tissues == "T-Cell"),]
  d4 = d4[order(d4$OR),]
  d4[nrow(d4)+1, ] = empty
  d5 = Data[which(Data$Tissues == "Other"),]
  d5 = d5[order(d5$OR),]
  d5[nrow(d5)+1, ] = empty
  Data = rbind(d1,d2,d3,d4,d5)
  
  length = nrow(Data)
  ### Create Vectors of Angles that will be used throughout plotting
  aShort = seq(0, 2*pi, length.out = length + 1)[-(length+1)]  
  aLong = seq(0, 2*pi, length.out = length + 1)
  
  ### Create Shaded Donut Shaped Plot, each shaded portion corresponds to different Tissue types
  small = data.frame("cos.a." = 1.25*cos(aLong), "sin.a." = 1.25*sin(aLong))
  large = data.frame("cos.a." = 1.85*cos(aLong), "sin.a." = 1.85*sin(aLong))
  n1 = nrow(d1)
  n2 = nrow(d2) + n1
  n3 = nrow(d3) + n2
  n4 = nrow(d4) + n3
  n5 = nrow(d5) + n4
  par(xpd = T)
  plot(0,0,xaxt = 'n', xlab = '', ylab = '',yaxt = 'n', bty = 'n', 
       xlim = c(-1.7, 1.7), ylim = c(-1.7, 1.7), col = "white")
  polygon(rbind(small[1:n1,],large[n1:1,]), col = ColorSet2[1], border = ColorSet2[1])
  polygon(rbind(small[n1:n2,],large[n2:n1,]), col = ColorSet2[2], border = ColorSet2[2])
  polygon(rbind(small[n2:n3,],large[n3:n2,]), col = ColorSet2[3], border = ColorSet2[3])
  polygon(rbind(small[n3:n4,],large[n4:n3,]), col = ColorSet2[4], border = ColorSet2[4])
  polygon(rbind(small[n4:n5,],large[n5:n4,]), col = ColorSet2[5], border = ColorSet2[5])    
  polygon(rbind(small[n5:(n5+1),],large[(n5+1):n5,]), col = ColorSet2[1], border = ColorSet2[1])    
  
  ### Create XY Coordinates for Circle with different Radii (further away = lower RAF)
  r = as.numeric(as.factor(round(Data$RAF,1)))
  r = seq(1.8, 1.3, length = length(unique(r)))[r]
  xvec = r*cos(aShort)
  yvec = r*sin(aShort)
  
  ### Draw solid/dotted lines around cirle corresponding to different RAF levels
  n = 1
  r1 = sort(unique(r))
  for (i in 1:length(r1)){
    n = 3 - n
    xvec1 = r1[i]*cos(aLong)
    yvec1 = r1[i]*sin(aLong)
    lines(xvec1,yvec1, col = "#F0F0F0", lty = n)
  }
  
  ### RAF Legend
  RAFLeg = seq(-1.3, -1.8, length = length(unique(r))) + .025
  text(0,RAFLeg[2],"90%", cex = .6, col = "white")
  text(0,RAFLeg[4],"70%", cex = .6, col = "white")
  text(0,RAFLeg[6],"50%", cex = .6, col = "white")
  text(0,RAFLeg[8],"30%", cex = .6, col = "white")
  text(0,RAFLeg[10],"10%", cex = .6, col = "white") 
  
  ### Set up cex, pch, lwd for each SNP, larger cex = greater OR, pch = "x" if Individual has 0 SNPs for that SNP, "circle" otherwise
  cex = seq(1.4,3.2,length.out = length(levels(as.factor(round(Data$OR,1)))))[as.numeric(as.factor(round(Data$OR,1)))]
  cex[Data$Dose == 0] = .8
  pch = rep(21,length)
  pch[Data$Dose == 0] = 4
  lwd = rep(3, length)
  lwd[Data$Dose == 0] = 2
  
  ### Plot SNPs around circle
  points(xvec, yvec, col = as.character(Data$FuncColors), cex = cex, pch = pch, lwd = lwd,
         xaxt = 'n', xlab = '', ylab = '',yaxt = 'n', bty = 'n', bg = "white",
         xlim = c(-1.7, 1.7), ylim = c(-1.7, 1.7)) 
  
  ### Plot a second ring inside SNPs when Individual has two of that SNP
  index = which(Data$Dose != 2)  
  points(xvec[-index], yvec[-index], col = as.character(Data$FuncColors)[-index],
         cex = cex[-index]-1.1, pch = 21, lwd = 3,
         xaxt = 'n', xlab = '', ylab = '',yaxt = 'n', bty = 'n', bg = "white",
         xlim = c(-1.7, 1.7), ylim = c(-1.7, 1.7)) 
  
  ### Add Gene Names to Plot if ShowGenes = TRUE  
  if(ShowGenes){
    Labels(as.character(Data$Genes), aShort, r, length, offset = .2)    
  }
  
  ### Tissue Legend
  t1 = ceiling(nrow(d1)/2)
  t2 = ceiling(nrow(d2)/2) + n1
  t3 = ceiling(nrow(d3)/2) + n2
  t4 = ceiling(nrow(d4)/2) + n3
  t5 = ceiling(nrow(d5)/2) + n4
  text(2.15*cos(aLong[t1]),2.15*sin(aLong[t1]), "Brain", col = ColorSet3[1], cex = 1)
  text(2.15*cos(aLong[t2]),2.15*sin(aLong[t2]), "B-Cell", col = ColorSet3[2], cex = 1)
  text(2.2*cos(aLong[t3]),2.2*sin(aLong[t3]), "B & T", col = ColorSet3[3], cex = 1)
  text(2.2*cos(aLong[t4]),2.2*sin(aLong[t4]), "T-Cell", col = ColorSet3[4], cex = 1)
  text(2.2*cos(aLong[t5]),2.2*sin(aLong[t5]), "Other", col = ColorSet3[5], cex = 1)
  
  ### MSGB Score
  text(0,.8, paste("Your MSGB Score:", as.character(MSGB)), cex = .8, font = 1, col = "dark grey")
  
  ### FID
  text(1.8,-2.2, paste("FID:", as.character(FID)), cex = .8, font = 1, col = "dark grey")
  
  ### Functions Legend
  legend(-2.2, 2.2, levels(Functions), pch = 21, col = ColorSet1, pt.bg = "white", pt.lwd = 2, 
         pt.cex = 1, bty = "n", cex = labelSize, y.intersp = 1.1)  
  #.5
  ### OR Legend
  legend(1.6, 2.25, legend=c("Low OR",NA,NA,NA,"High OR"), pch=21, pt.cex = c(.6,.95,1.3,1.65, 2),cex = labelSize, bty = "n", 
         pt.bg = 'light grey', y.intersp = c(1,1.2,1.3,1.4,1.5), x.intersp = 1.35)

  ### Dose Legend
  legend(-2.2, -1.9, legend=c("0 SNPs","1 SNP","2 SNPs"), pch=c(4,21,21), pt.cex = c(.8,2,2),cex = labelSize, bty = "n", 
         pt.bg = c("black","white", "white"), y.intersp = 2, pt.lwd = c(2,2,2), x.intersp = 1.35)
  legend(-2.2, -1.9, legend=c(NA,NA,NA), pch=c(NA,NA,21), pt.cex = c(NA,NA,1.2),cex = labelSize, bty = "n", 
         pt.bg = c(NA,NA, "white"), y.intersp = 2, pt.lwd = c(NA,NA,2), x.intersp = 1.35)
  
  ### Rescale Distribution of Control and Case Data to fit between -1.1 and 1.1
  ControlDist = Distribution(Control,MSrange = c(min(Control), max(Control)))
  CaseDist = Distribution(Case, MSrange = c(min(Case), max(Case)))
  DistList = rescale(ControlDist, CaseDist, -1.1, 1.1)
  xControl = DistList[[1]]
  xCase = DistList[[2]]
  
  ### Make Vectors of Colors for Controls and Cases
  ControlColor = sample(ColorSet4, length(xControl), replace = T)
  CaseColor = sample(ColorSet5, length(xCase), replace = T)
  
  ### Plot Men and Women according to Control Distribution
  MaleFem = sample(c(1,2), length(xControl), replace = T)
  PlotPeople(Man, x = xControl[which(MaleFem == 1)][which(ControlColor == ColorSet4[1])], y = .25, color = ColorSet4[1])
  PlotPeople(Man, x = xControl[which(MaleFem == 1)][which(ControlColor == ColorSet4[2])], y = .25, color = ColorSet4[2])
  PlotPeople(Woman, x = xControl[which(MaleFem == 2)][which(ControlColor == ColorSet4[1])], y = .25, color = ColorSet4[1])
  PlotPeople(Woman, x = xControl[which(MaleFem == 2)][which(ControlColor == ColorSet4[2])], y = .25, color = ColorSet4[2])
  
  ### Plot Men and Women according to Case Distribution
  MaleFem = sample(c(1,2), length(xCase), replace = T)
  PlotPeople(Man, x = xCase[which(MaleFem == 1)][which(CaseColor == ColorSet5[1])], y = -.35, color = ColorSet5[1])
  PlotPeople(Man, x = xCase[which(MaleFem == 1)][which(CaseColor == ColorSet5[2])], y = -.35, color = ColorSet5[2])
  PlotPeople(Woman, x = xCase[which(MaleFem == 2)][which(CaseColor == ColorSet5[1])], y = -.35, color = ColorSet5[1])
  PlotPeople(Woman, x = xCase[which(MaleFem == 2)][which(CaseColor == ColorSet5[2])], y = -.35, color = ColorSet5[2])
  
  ### Plot Individual in Case and Control Data
  Min = min(ControlDist,CaseDist)
  Max = max(ControlDist,CaseDist)
  Score = (MSGB-Min)/(Max-Min)*(1.1+1.1)+-1.1  
  if(Sex == 1){
    PlotPeople(Man, x = Score, y=c(-.35,0.25))
  }
  else{
    PlotPeople(Woman, x = Score, y=c(-.35,0.25))
  }
  
  ### Add Info to Middle Distribution
  arrows(.45, 0, 1, 0, lwd = 2.5, col = ColorSet6, length = .18)
  segments(-.8, 0, -.45, 0, lwd = 2.5, col = ColorSet6)
  text(0, 0, "Increasing Genetic Risk for MS", cex = .6, col = ColorSet6)
  text(0, .5, "Healthy Controls", col = ColorSet6, cex = 1)
  text(0, -.6, "Patients with Multiple Sclerosis", col = ColorSet6, cex = 1)
  par(new = F)  
}
pngPlots = function(FID, outputPath){
  filePath = file.path(outputPath, paste0("MSGBPlot", FID, ".png"))
  png(filePath, width = 800, height = 800 , res = 100)
  BubblePlot(DoseInfo, SNPInfo, ProxyInfo, GeneMatch, GeneInfo, FID, Woman, Man, ShowGenes = TRUE, labelSize = .6)
  dev.off()
}

# Read in Data ------------------------------------------------------------
SNPInfo = read.table(file.path(dataPath, "SNPInfo.txt"), header = T)
DoseInfo = read.table(file.path(dataPath, "DoseInfo.txt"), header = T)
ProxyInfo = read.table(file.path(dataPath, "ProxyInfo.txt"), header = T)
GeneMatch = read.table(file.path(dataPath, "GeneMatch.txt"), header = T)
GeneInfo = read.csv(file.path(dataPath, "GeneInfo.csv"), header = T)
Woman = readPNG(file.path(dataPath, "Woman.png"))
Man = readPNG(file.path(dataPath, "Man.png"))

### MSGB Plot for FID = 1
BubblePlot(DoseInfo, SNPInfo, ProxyInfo, GeneMatch, GeneInfo, 1, Woman, Man, ShowGenes = TRUE)

### Create a MSGB Plot for each Patient
FID = DoseInfo$FID
sapply(FID, pngPlots, outputPath = outputPath)
