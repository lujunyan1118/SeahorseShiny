library(shiny)
library(Biobase)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Rtsne)
library(genefilter)
library(DESeq2)
library(piano)
library(DT)
library(limma)

load("shinyData.RData")

#Pre-process for correlation test between Seahorse and RNAseq
overPat <- intersect(colData(ddsNorm)$PatID,colnames(seaMain))
ddsSea <- ddsNorm[,colData(ddsNorm)$PatID %in% overPat]
seaRNA <- seaMain[,overPat]

#funtion to perform t test
tTest <- function(geneList, seaList, eVar = TRUE) {
  geneList <- geneList[!is.na(geneList)]
  seaList <- seaList[names(geneList)]
  tt <- t.test(seaList~geneList, var.equal = eVar)
  tt$p.value
}

function(input,output){
  
  #Reactive object for the t-test panel
  tResult.raw <- reactive({
    
    #Get the useful seahorse measurement
    seaTest <- seaMain[input$seleSea1,]
    seaTest <- seaTest[!is.na(seaTest)]
    
    if (input$treatIGHV == TRUE) {
      ighv <- as.vector(t(geneMain[,"IGHV.Uppsala.U.M",drop=FALSE]))
      names(ighv) <- rownames(geneMain)
      ighv <- as.factor(ighv[!is.na(ighv)])
      seaTest <- seaTest[intersect(names(seaTest),names(ighv))]
      ighv <- ighv[names(seaTest)]
      if (input$methodIGHV == "Remove"){
        seaTest <- removeBatchEffect(t(as.matrix(seaTest)),batch = ighv)[1,]  #remove ighv effect
      } else if (input$methodIGHV == "Unmutated only") {
        seaTest <- seaTest[names(ighv[ighv == 0])]
      } else if (input$methodIGHV == "Mutated only") {
        seaTest <- seaTest[names(ighv[ighv == 1])]
      }
    }
    
    #Get useful genetic backgroudn matrix
    geneTest <- geneMain[names(seaTest),]
    geneTest <- geneTest[,colSums(geneTest == 0, na.rm = TRUE) >=3 & colSums(geneTest == 1,na.rm = TRUE) >=3] #at least 3 unmutated value and 3 mutated samples
    geneTest$Methylation_Cluster <- NULL 
    geneTest <- data.frame(t(geneTest))
    
    #perform t-test
    if (input$tMethod == "Equal variance") eVar <- TRUE else eVar <- FALSE
    pRes <- apply(geneTest, 1, function(x) tTest(x,seaTest,eVar))
    pRes <- sort(pRes)
    pTabRes <- data.frame(cbind(pRes,p.adjust(pRes,method = "BH")))
    colnames(pTabRes) <- c("p.raw","p.adj.BH")
    list(tTab = pTabRes, seaVec = seaTest)
  })
  
    tResult <- reactive({
      rawTab <- tResult.raw()$tTab
      if (input$filterP == TRUE) {
        newTab <- rawTab[rawTab$p.raw <= input$rawPcut & rawTab$p.adj.BH <= input$adjPcut,]
      } else newTab <- rawTab
      
      list(tTab = newTab, seaVec = tResult.raw()$seaVec)
    })

  selectTab <- reactive({
    if (input$screenSet == "main"){
      list(viabTab = viabMain, tabCor = mainCor, geneTab = geneMain, seaTab = seaMain)
    } else {
      list(viabTab = viabEMBL, tabCor = emblCor, geneTab = geneEMBL, seaTab = seaEMBL)
    }
  })
  
  corFiltered <- reactive({
    if (input$filterCor) {
    selectTab()$tabCor[abs(selectTab()$tabCor$coef) >= input$coefCut & selectTab()$tabCor$pval <= input$pCut,]}
    else selectTab()$tabCor 
    })
  
  drugID <- reactive(selectTab()$tabCor[selectTab()$tabCor$name == input$drugName & selectTab()$tabCor$conc == input$drugConc,]$drugID[1])

  corRNA.filtered <- reactive({
    origin <- ddsCor[[input$corMethod]][[input$seleSea]]
    if (input$filterRNA == TRUE) {
      filtered <- origin[(abs(origin$coef) >= input$coefCut1 & origin$p.raw <= input$pCut1), ]
    } else filtered <- origin
    filtered
  })
  
  output$downloadTable <- downloadHandler(
    filename = function() { paste('CorTable', '.csv', sep='') },
    content = function(file) {
      write.csv2(corRNA.filtered(), file)
    }
  )
  
  resGSEA <- eventReactive(input$goGSEA, {
    #parameters for GSEA
    nPerm <- input$permNum
    
    #readin geneset database
    inGMT <- loadGSC(input$setBox,type="gmt")
    
    #proceccing correlation data
    corTab <- corRNA.filtered()
    corTab <- corTab[order(corTab$p.raw),]
    corTab <- corTab[!duplicated(corTab$symbol),] #remove genes with the same name, since the genes are already ordered by p value, the genes with low pvalue will be kept.
    rownames(corTab) <- corTab$symbol
    corTab <- corTab[order(corTab$coef,decreasing = TRUE),] #re-rank by coefficient 
    
    #GSEA analysis
    myCoef <- corTab[,"coef",drop=FALSE]
    res <- runGSA(geneLevelStats = myCoef,geneSetStat = "gsea",adjMethod = "fdr", gsc=inGMT, signifMethod = 'geneSampling', nPerm = nPerm)
  })
    
  gseaList <- reactive({
    if (input$ifGSEA == TRUE) {
    corGene <- corRNA.filtered()
    setName <- input$tab1_row_last_clicked
    geneList <- loadGSC(input$setBox,type="gmt")$gsc[[setName]]
    geneTab <- corGene[corGene$symbol %in% geneList,]
    } 
  })
  
  plotTable <- reactive({
    meta <- selectTab()$seaTab[input$measureName,]
    viab <- exprs(selectTab()$viabTab)[drugID(),]
    #stopifnot(all(names(viab) == names(meta)) & all(names(viab) == colnames(selectTab()$geneTab)))
    as.data.frame(cbind(meta,viab,selectTab()$geneTab))
  })
  
  #prepare a table for scatter plot in panel 3
  scatterTable <- reactive({
    if (input$dataSet1 == "Seahorse") {
      seaMain <- seaMain[order(rowMeans(seaMain)),]
      allTab <- cbind(t(seaMain),geneMain)
      allTab <- melt(data = allTab)
    }
    if (input$dataSet1 == "Main screen") {
      if (input$patSet == "All patients") {
        viabSub <- viabAll
        Diagnosis <- pData(viabSub)$Diagnosis
        #allTab <- cbind(t(exprs(viabAll)),geneAll,Diagnosis)
      } else {
        viabSub <- viabAll[,pData(viabAll)$Diagnosis == "CLL"]
        Diagnosis <- pData(viabSub)$Diagnosis
      }
      if (input$drugSet == "High variant drugs") {
        sds <- rowSds(exprs(viabSub))
        viabSub <- viabSub[sds > quantile(sds,1-input$varCut/100),]
        viabSub <- viabSub[order(rowMedians(exprs(viabSub))),]
        allTab <- cbind(t(exprs(viabSub)),geneAll[colnames(viabSub),],Diagnosis)
      } else if (input$drugSet == "Drug set 1 (small)") {
        drugs <- c("duvelisib","idelalisib","MK-2206","everolimus","ibrutinib","spebrutinib","selumetinib","encorafenib") #a list of drugs, representing 6 major pathways: PI3K, AKT, mTOR, BTK, MEK, BRAF
        listID <- rownames(viabSub[(fData(viabSub)$name %in% drugs & fData(viabSub)$subtype %in% c(2,3,4,5)),]) #Only using the lowest two concentrations
        viabSub <- viabSub[listID,]
        viabSub <- viabSub[order(rowMedians(exprs(viabSub))),]
        allTab <- cbind(t(exprs(viabSub)),geneAll[colnames(viabSub),],Diagnosis)
      } else if (input$drugSet == "Drug set 2 (large)") {
        targetedDrugNames <- c("ibrutinib", "idelalisib",  "PRT062607 HCl", "duvelisib", 
                               "spebrutinib", "selumetinib", "MK-2206",  "everolimus", "encorafenib")
        targetedDrugs <- rownames(fData(viabAll)[fData(viabAll)$name %in% targetedDrugNames & fData(viabAll)$subtype %in% c(4,5),])
        chemoDrugNames <- c("fludarabine", "doxorubicine",  "nutlin-3")
        chemoDrugs <- rownames(fData(viabAll)[fData(viabAll)$name %in% chemoDrugNames & fData(viabAll)$subtype %in% c(3,4,5),])
        badDrugs <- c(bortezomib = "D_008", `NSC 74859` = "D_025")
        candDrugs <- rownames(viabAll)[fData(viabAll)$type=="viab" & !(fData(viabAll)$id %in% badDrugs) & fData(viabAll)$subtype %in% paste(2:5)]
        thresh <- list(effectVal = 0.7, effectNum = 4, viab = 0.6, maxval = 1.1)
        overallMean  <- rowMeans(exprs(viabAll)[candDrugs, ])
        nthStrongest <- apply(exprs(lpdAll)[candDrugs, ], 1, function(x) sort(x)[thresh$effectNum])
        eligibleDrugs <- candDrugs[ overallMean >= thresh$viab & nthStrongest <= thresh$effectVal ] %>%
          union(targetedDrugs) %>% union(chemoDrugs)
        viabSub <- viabSub[eligibleDrugs,]
        viabSub <- viabSub[order(rowMedians(exprs(viabSub))),]
        allTab <- cbind(t(exprs(viabSub)),geneAll[colnames(viabSub),],Diagnosis)
      } else {
        viabSub <- viabSub[order(rowMedians(exprs(viabSub))),]
        allTab <- cbind(t(exprs(viabSub)),geneAll[colnames(viabSub),],Diagnosis)
        }
      allTab <- melt(data = allTab)
      allTab$variable <- as.character(allTab$variable)
      #annotation drug names
      allTab$variable <- sprintf("%s(%s)",fData(viabSub)[allTab$variable,]$name,fData(viabSub)[allTab$variable,]$subtype)
      }
      allTab
  })
  
  tsneResult <- reactive({
    if (input$dataSet == "Seahorse"){
      if (input$clusterMethod == "ifTSNE") {
        tsne <- Rtsne(distSea, perplexity = input$perplexity, theta = input$theta, max_iter = input$maxIter, is_distance = TRUE)
        tsne <- tsne$Y
        rownames(tsne) <- labels(distSea)
        colnames(tsne) <- c("x","y")
        outTab <- as.data.frame(cbind(tsne,geneMain[rownames(tsne),]))}
      if (input$clusterMethod == "ifPCA") {
        pca <- prcomp(na.omit(t(seaMain)))
        varExp <- (pca$sdev^2 / sum(pca$sdev^2))*100
        pca <- data.frame(pca$x)
        names(varExp) <- colnames(pca)
        outTab <- list( data = as.data.frame(cbind(pca,geneMain[rownames(pca),])), percent = varExp)
      }
      if (input$clusterMethod == "ifKmeans") {
        pca <- prcomp(na.omit(t(seaMain)))
        pca <- pca$x[,1:2]
        cluster <- as.factor(kmeans(na.omit(t(seaMain)),centers = input$clusterNum)$cluster)
        outTab <- as.data.frame(cbind(pca,cluster,geneMain[rownames(pca),]))
      }
    } else if (input$dataSet %in% c("Main screen","EMBL screen")) {
        if (input$dataSet == "Main screen"){
          viabTab <- log(exprs(viabAll))
          viabTab <- viabTab - rowMedians(viabTab, na.rm = TRUE) # center
          viabTab <- t(apply(viabTab,1,function(x) x/mad(x,na.rm = TRUE)))
          Diagnosis <- pData(viabAll)$Diagnosis
        } else if (input$dataSet == "EMBL screen") {
          thresh <- list(effectVal = 0.8, effectNum = 4, viab = 0.5, maxval = 1.1)
          overallMean  <- rowMeans(exprs(allEMBL))
          nthStrongest <- apply(exprs(allEMBL), 1, function(x) sort(x)[thresh$effectNum])
          viabTab <- allEMBL[ overallMean >= thresh$viab & nthStrongest <= thresh$effectVal, ]
          viabTab <- log(exprs(viabTab))
          viabTab <- viabTab - rowMedians(viabTab, na.rm = TRUE) # center
          viabTab <- t(apply(viabTab,1,function(x) x/mad(x,na.rm = TRUE)))
          Diagnosis <- pData(allEMBL)$Diagnosis
        }
        if (input$clusterMethod == "ifTSNE") {
          tsne <- Rtsne(t(viabTab), perplexity = input$perplexity, theta = input$theta, max_iter = input$maxIter, is_distance = FALSE)
          tsne <- tsne$Y
          rownames(tsne) <- colnames(viabTab)
          colnames(tsne) <- c("x","y")
          outTab <- as.data.frame(cbind(tsne,geneAll[rownames(tsne),],Diagnosis))
          } else if (input$clusterMethod == "ifPCA") {
            pca <- prcomp(na.omit(t(viabTab)))
            varExp <- (pca$sdev^2 / sum(pca$sdev^2))*100
            pca <- data.frame(pca$x)
            names(varExp) <- colnames(pca)
            outTab <- list(data = as.data.frame(cbind(pca,geneAll[rownames(pca),],Diagnosis)), percent = varExp)
          } else if (input$clusterMethod == "ifKmeans") {
            pca <- prcomp(na.omit(t(viabTab)))
            pca <- pca$x[,1:2]
            cluster <- as.factor(kmeans(na.omit(t(viabTab)),centers = input$clusterNum)$cluster)
            outTab <- as.data.frame(cbind(pca,cluster,geneAll[rownames(pca),],Diagnosis))
          } 
    } else if (input$dataSet == "RNAseq") {
          sds <- rowSds(assay(ddsNorm))
          ddsVar <- ddsNorm[order(sds,decreasing=TRUE)[1:5000],] #only top 5000 most variant genes
          ddsTab <- assay(ddsVar)
          colnames(ddsTab) <- colData(ddsVar)$PatID
          Diagnosis <- colData(ddsVar)$diag
          if (input$clusterMethod == "ifTSNE") {
            tsne <- Rtsne(t(ddsTab), perplexity = input$perplexity, theta = input$theta, max_iter = input$maxIter, is_distance = FALSE)
            tsne <- tsne$Y
            rownames(tsne) <- colnames(ddsTab)
            colnames(tsne) <- c("x","y")
            outTab <- as.data.frame(cbind(tsne,geneAll[rownames(tsne),],Diagnosis))
          } else if (input$clusterMethod == "ifPCA") {
            pca <- prcomp(na.omit(t(ddsTab)))
            varExp <- (pca$sdev^2 / sum(pca$sdev^2))*100
            pca <- data.frame(pca$x)
            names(varExp) <- colnames(pca)
            outTab <- list(data = as.data.frame(cbind(pca,geneAll[rownames(pca),],Diagnosis)), percent = varExp)
          } else if (input$clusterMethod == "ifKmeans") {
            pca <- prcomp(na.omit(t(ddsTab)))
            pca <- pca$x[,1:2]
            cluster <- as.factor(kmeans(na.omit(t(ddsTab)),centers = input$clusterNum)$cluster)
            outTab <- as.data.frame(cbind(pca,cluster,geneAll[rownames(pca),],factor(Diagnosis)))
          } 
        } 
    outTab
  })
  
  output$tab1 <- DT::renderDataTable({
    if (input$ifGSEA == FALSE) {
    datatable(corRNA.filtered(),selection = 'single') %>% formatRound(c('coef','p.raw','p.adj.BH', 'p.adj.IHW'))} 
    else {
      resTab <- GSAsummaryTable(resGSEA())
      rownames(resTab) <- resTab$Name
      resTab$Name <- NULL
      colnames(resTab) <- c("Gene Number","Stat","p.up","p.up.adj","p.down","p.down.adj","Number up","Number down")
      datatable(resTab,selection = 'single') %>% formatRound(c('Stat','p.up','p.up.adj', 'p.down','p.down.adj'))}
    })
  
  output$tab2 <- DT::renderDataTable({
    if (input$ifGSEA==TRUE) {datatable(gseaList(),selection = 'none') %>% formatRound(c('coef','p.raw','p.adj.BH', 'p.adj.IHW'))}})
  
  output$tTable <- DT::renderDataTable(datatable(tResult()$tTab,selection = 'single' ) %>% formatRound(c('p.raw','p.adj.BH')))
  
  output$plot1 <- renderPlot({
    if (input$ifGSEA == FALSE) {
      geneID <- input$tab1_row_last_clicked
      plotTab <- data.frame(cbind(assay(ddsSea)[geneID,],seaRNA[input$seleSea,]))
      colnames(plotTab) <- c(geneID,input$seleSea)
      p <- ggplot(plotTab, aes_string(x=geneID,y=input$seleSea)) + geom_point(size=2) + geom_smooth(method = "lm") + theme(text=element_text(size=15),legend.position="bottom")
      if (input$colorGene2 == TRUE) {
        plotTab <- cbind(plotTab, geneMain[colnames(seaRNA),])
        p <- ggplot(plotTab, aes_string(x=geneID,y=input$seleSea)) + geom_point(size=2) + geom_smooth(method = "lm") + theme(text=element_text(size=15),legend.position="bottom") + geom_point(aes_string(color = input$geneBox4))
      }
      plot(p) 
    } else {
      networkPlot(resGSEA(), class="distinct", direction = 'both', significance = input$sigLevel, label='names', adjusted=TRUE)
    }
    
  })
  
  #scatter plot for the t-Test panel
  output$tPlot <- renderPlot({
    seaVec <- tResult()$seaVec
    geneVec <- geneMain[names(seaVec),input$tTable_row_last_clicked]
    plotTab <- data.frame(measure = seaVec, gene = geneVec)
    if (input$showNA2 == FALSE) plotTab <- na.omit(plotTab)
    p <- ggplot(plotTab, aes(x=factor(gene), y=measure, fill=factor(gene))) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.2) + theme(text=element_text(size=15)) + ylab(input$seleSea1) + xlab("")
    plot(p)
  })
  
  
  output$measureBox <- renderUI(
    selectInput("measureName","Select metabolic feature",unique(corFiltered()$measure),size=5, selectize = FALSE)
  )
  output$drugBox <- renderUI(
    selectInput("drugName","Select drug",sort(unique(corFiltered()[corFiltered()$measure == input$measureName,]$name)),size = 5, selectize = FALSE)
  )
  output$concBox <- renderUI(
    selectInput("drugConc","Select concentration",corFiltered()[corFiltered()$measure == input$measureName & corFiltered()$name == input$drugName,]$conc,size = 5, selectize = FALSE)
  )
  output$geneBox2 <- renderUI({
    if (input$dataSet %in% c("Main screen","EMBL screen","RNAseq")) {
      selectInput("seleGene","Select genetic feature",c("Diagnosis",colnames(geneAll)),size=5, selectize = FALSE)
    } else {
      selectInput("seleGene","Select genetic feature",colnames(geneMain),size=5, selectize = FALSE)
    }
  })
  output$geneBox3 <- renderUI({
    if (input$dataSet1 == "Main screen") {
      if (input$patSet == "All patients") {
        selectInput("seleGene1","Select genetic feature",c("Diagnosis",colnames(geneAll)),size=5, selectize = FALSE)}
      else selectInput("seleGene1","Select genetic feature",colnames(geneAll),size=5, selectize = FALSE)
    } else {
      selectInput("seleGene1","Select genetic feature",colnames(geneMain),size=5, selectize = FALSE)
    }
  })
  
  #for scatter plot of associations between seahorse measurement and viability data
  output$scatterPlot <- renderPlot({
    corrLine <- selectTab()$tabCor[selectTab()$tabCor$drugID == drugID() & selectTab()$tabCor$measure == input$measureName,,drop=FALSE]
    plotAnno = sprintf("coef = %1.2f, p = %1.2f", corrLine$coef, corrLine$pval)
    if (!input$colorGene) {
      p <- ggplot(plotTable(), aes_string(x="viab",y="meta")) + geom_point(size=2) + geom_smooth(method = "lm") + xlab(sprintf("%s(%s)",input$drugName,input$drugConc)) + ylab(input$measureName) + theme(text=element_text(size=15),legend.position="bottom") + annotate("text", x=Inf, y = Inf, label = plotAnno, vjust=1, hjust=1, size=8)
      }
    else {
      if (input$showNA) {
        plotTabNew <- plotTable()
      } else {
        plotTabNew <- plotTable()
        plotTabNew <- plotTabNew[!is.na(plotTabNew[,input$geneBox]),] #remove obeservations without certain genetic background
      }
      if (input$sperateRegression) {
        p <- ggplot(plotTabNew, aes_string(x="viab",y="meta", color = input$geneBox)) + geom_point(size=2) + geom_smooth(method = "lm") + xlab(sprintf("%s(%s)",input$drugName,input$drugConc)) + ylab(input$measureName) + theme(text=element_text(size=15),legend.position="bottom") + annotate("text", x=Inf, y = Inf, label = plotAnno, vjust=1, hjust=1, size=8)
      } else {
        p <- ggplot(plotTabNew, aes_string(x="viab",y="meta")) + geom_point(size=2, aes_string(color = input$geneBox)) + geom_smooth(method = "lm") + xlab(sprintf("%s(%s)",input$drugName,input$drugConc)) + ylab(input$measureName) + theme(text=element_text(size=15), legend.position="bottom") + annotate("text", x=Inf, y = Inf, label = plotAnno, vjust=1, hjust=1, size=8)
      }
    }
    if (input$drugName != "lymphocyte.doubling") p <- p + xlim(0,1.2)
    print(p)
  }, height = 500, width =500)
  
  #for scatter plot of t-sne results
  output$tsnePlot <- renderPlot({
    tsneTable <- tsneResult()
    if (input$clusterMethod == "ifKmeans") {
      p <- ggplot(tsneTable, aes_string(x=colnames(tsneTable)[1],y=colnames(tsneTable)[2])) + geom_point(size=2, aes_string(color = input$seleGene, shape = "cluster")) + theme(legend.position = "bottom")
    } 
    if (input$clusterMethod == "ifPCA") {
      p <- ggplot(tsneTable$data, aes_string(x=input$xaxis,y=input$yaxis)) + geom_point(size=2, aes_string(color = input$seleGene)) + theme(legend.position = "bottom") + xlab(sprintf("%s (%2.1f%%)",input$xaxis,tsneTable$percent[input$xaxis])) + ylab(sprintf("%s (%2.1f%%)",input$yaxis,tsneTable$percent[input$yaxis]))
    }
    if (input$clusterMethod == "ifTSNE") {
    p <- ggplot(tsneTable, aes_string(x=colnames(tsneTable)[1],y=colnames(tsneTable)[2])) + geom_point(size=2, aes_string(color = input$seleGene)) + theme(legend.position = "bottom") }
    print(p)
  }, height=500, width=500)
  
  #for scatter plot of all data
  output$scatterAll <- renderPlot({
    mt <- scatterTable()
    #remove data point with NA facotr
    mt <- mt[!is.na(mt[,input$seleGene1]),]
    mt$variable <- factor(mt$variable, levels = unique(mt$variable))
    #plot
    if (input$dataSet1 == "Seahorse") {
      p <- ggplot(mt, aes(x=variable,y=value)) + geom_jitter(size=2, aes_string(color = input$seleGene1), na.rm = TRUE, alpha=0.5) + theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, hjust = 1))
    } else {
      p <- ggplot(mt, aes(x=variable,y=value)) + geom_jitter(size=2, aes_string(color = input$seleGene1), na.rm = TRUE, alpha=0.5) + theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, hjust = 1)) + ylim(0,1.2)
      }
      plot(p)
  }, height = 500, width = 900)
}