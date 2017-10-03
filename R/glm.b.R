library(ggplot2)
#' @import ggplot2
glmClass <- R6::R6Class(
  "glmClass",
  inherit = glmBase,
  private=list(
    .model=NA,
    .postHocRows=NA,
    .cleanData=function() {
      
      dep <- self$options$dep
      factors <- self$options$factors
      
      covs <- NULL
      if ('covs' %in% names(self$options))
        covs <- self$options$covs
      
      data <- self$data

      if ( ! is.null(dep))
        data[[dep]] <- jmvcore::toNumeric(data[[dep]])
      
      for (factor in factors) {
        data[[factor]] <- as.factor(data[[factor]])
      }
      for (covariate in covs)
        data[[covariate]] <- jmvcore::toNumeric(data[[covariate]])
      

      for (scaling in self$options$scaling) {
        if (scaling$type=="centered") data[[scaling$var]]<-scale(data[[scaling$var]],scale = F)  
        if (scaling$type=="standardized") data[[scaling$var]]<-scale(data[[scaling$var]])
        data[[scaling$var]] <- as.numeric(data[[scaling$var]])
        
      }

      data <- na.omit(data)
      data
    },
    .init=function() {
      dep <- self$options$dep
      factors <- self$options$factors
      modelTerms <- private$.modelTerms()
      if (length(modelTerms) == 0)
        return()
      
      data <- private$.cleanData()
      
      anovaTable      <- self$results$main
      postHocTables   <- self$results$postHoc
      contrastsTables <- self$results$contrasts
      estimatesTable <- self$results$estimates
      # main table
      
      modelTerms <- private$.modelTerms()
      
      
        if (length(modelTerms) > 0) {
        for (term in modelTerms) {
        anovaTable$addRow(rowKey=term, list(name=jmvcore::stringifyTerm(term)))
        }  
        anovaTable$addFormat(col=1, rowNo=1,                  format=Cell.BEGIN_GROUP)
        anovaTable$addFormat(col=1, rowNo=length(modelTerms), format=Cell.END_GROUP)

        } else {
        anovaTable$addRow(rowKey='.', list(name='.'))
        anovaTable$addFormat(col=1, rowKey='.', format=Cell.BEGIN_END_GROUP)
      }
      
      anovaTable$addRow(rowKey='', list(name='Residuals'))
      anovaTable$addFormat(col=1, rowKey='', format=Cell.BEGIN_END_GROUP)
      

      # contrasts
      
      for (contrast in self$options$contrasts) {
        table <- contrastsTables$addItem(contrast)
        var <- data[[contrast$var]]
        if (contrast$type=="default")
          contrast$type="deviation"
        
        levels <- base::levels(var)
        labels <- private$.contrastLabels(levels, contrast$type)
        dummies<-paste(contrast$var,1:length(labels),sep="")
        groups<-paste(1:length(levels),levels,sep = "=", collapse = ", ")

        for (i in 1:length(labels)) {
          table$addRow(rowKey=labels[[i]], values=list(contrast=labels[[i]],term=dummies[[i]],groups=groups))
        }
        
      }
      
      
      # post hoc
      private$.initPostHoc(data)
      
      # descriptives
      
      descTable <- self$results$desc
      factorNames <- self$options$factors
      
      if (length(factorNames) > 0) {
        
        data <- select(data, rev(factorNames))
        al <- as.list(data)
        names(al) <- rev(paste0('f', seq_len(length(al))))
        ll <- sapply(al, base::levels, simplify=FALSE)
        ll$stringsAsFactors <- FALSE
        grid <- do.call(base::expand.grid, ll)
        grid <- rev(grid)
        
        for (i in seq_len(ncol(grid))) {
          colName <- colnames(grid)[[i]]
          descTable$addColumn(name=colName, title=factorNames[[i]], index=i)
        }
        
        for (rowNo in seq_len(nrow(grid))) {
          row <- grid[rowNo,]
          if ( ! is.list(row))
            row <- list(f1=row)
          descTable$addRow(rowKey=row, values=row)
        }
      }
      # descriptives plots
      private$.initDescPlots(data)
    },
    .run=function() {
      suppressWarnings({
        
        dep <- self$options$dep
        factors <- self$options$factors
        covs<-self$options$covs
        modelTerms <- private$.modelTerms()
        
        if (is.null(dep) ||  length(modelTerms) == 0)
          return()
        
        base::options(contrasts = c("contr.sum","contr.poly"))
        
        data <- private$.cleanData()

        # data <- lapply(data, function(x) {
        #   if (is.factor(x))
        #     levels(x) <- toB64(levels(x))
        #   return(x)
        # })
        
        if (is.factor(data[[dep]]))
          reject('Dependent variable must be numeric')
        
        for (factorName in factors) {
          lvls <- base::levels(data[[factorName]])
          if (length(lvls) == 1)
            reject("Factor '{}' contains only a single level", factorName=factorName)
          else if (length(lvls) == 0)
            reject("Factor '{}' contains no data", factorName=factorName)
        }
        
        for (contrast in self$options$contrasts) {
          levels <- base::levels(data[[contrast$var]])
          stats::contrasts(data[[contrast$var]]) <- private$.createContrasts(levels, contrast$type)
        }
        for (contrast in factors) 
            print(contrasts(data[[contrast]]))
        
        formula <- jmvcore::constructFormula(dep, modelTerms)
        formula <- stats::as.formula(formula)
        
        model <- stats::lm(formula, data)
        private$.model <- model
        self$results$.setModel(model)

        singular <- NULL
        
        if (self$options$ss == '1') {
          
          results <- try(stats::anova(private$.model), silent=TRUE)
          
        } else if (self$options$ss == '2') {
          
          results <- try(car::Anova(private$.model, type=2, singular.ok=FALSE), silent=TRUE)
          if (isError(results)) {
            message <- extractErrorMessage(results)
            if (message == 'there are aliased coefficients in the model')
              singular <- 'Singular fit encountered; one or more predictor variables are a linear combination of other predictor variables'
            results <- try(car::Anova(private$.model, type=2, singular.ok=TRUE), silent=TRUE)
          }
          
        } else {
          
          results <- try({
            r <- car::Anova(private$.model, type=3, singular.ok=FALSE, silent=TRUE)
            r <- r[-1,]
          })
          
          if (isError(results)) {
            message <- extractErrorMessage(results)
            if (message == 'there are aliased coefficients in the model')
              singular <- 'Singular fit encountered; one or more predictor variables are a linear combination of other predictor variables'
            results <- try({
              r <- car::Anova(private$.model, type=3, singular.ok=TRUE, silent=TRUE)
              r <- r[-1,]
            })
          }
        }
        
        if (isError(results)) {
          message <- extractErrorMessage(results)
          if (message == 'residual df = 0')
            reject('Residual sum of squares and/or degrees of freedom is zero, indicating a perfect fit')
        }
        
        if (results['Residuals', 'Sum Sq'] == 0 || results['Residuals', 'Df'] == 0)
          reject('Residual sum of squares and/or degrees of freedom is zero, indicating a perfect fit')
        
        anovaTable <- self$results$main
        
        if ( ! is.null(singular))
          anovaTable$setNote('singular', singular)
        
        rowCount <- dim(results)[1]
        rowNames <- dimnames(results)[[1]]
        
        errIndex <- nrow(results)
        errSS <- results[errIndex,'Sum Sq']
        errDF <- results[errIndex,'Df']
        errMS <- errSS / errDF
        totalSS <- sum(results[['Sum Sq']], na.rm=TRUE)
        
        for (i in seq_len(rowCount)) {
          rowName <- rowNames[i]
          
          ss <- results[i,'Sum Sq']
          df <- results[i,'Df']
          ms <- ss / df
          F  <- results[i,'F value']
          p  <- results[i,'Pr(>F)']
          
          if ( is.finite(F)) {
            e <- ss / totalSS
            ep <- ss / (ss + errSS)
            w <- (ss - (df * errMS)) / (totalSS + errMS)
          } else {
            e <- ''
            ep <- ''
            w <- ''
          }
          
          if ( ! is.finite(ss))
            ss <- 0
          if ( ! is.finite(ms))
            ms <- ''
          if ( ! is.finite(F))
            F <- ''
          if ( ! is.finite(p))
            p <- ''
          
          tableRow <- list(ss=ss, df=df, ms=ms, F=F, p=p, etaSq=e, etaSqP=ep, omegaSq=w)
          if (i < rowCount) {
            anovaTable$setRow(rowNo=i, tableRow)
          }
          else {
            if (rowCount < anovaTable$rowCount) {
              blankRow <- list(ss=0, df=0, ms='', F='', p='', etaSq='', etaSqP='', omegaSq='')
              for (j in seq(i, anovaTable$rowCount-1))
                anovaTable$setRow(rowNo=j, blankRow)
            }
            anovaTable$setRow(rowKey='', tableRow) # residual
          }
        }

        estimatesTable <- self$results$estimates
        
        suppressWarnings({
        eresults<-stats::summary.lm(private$.model)[['coefficients']]
#        eresults<-stats::summary.lm(mod)[['coefficients']]
        }) #end suppresseWarnings
        if ("beta" %in% self$options$effectSize)
        {
          zdata<-data
          zdata[[dep]]<-scale(zdata[[dep]])
          for (var in covs)
            zdata[[var]]<-scale(zdata[[var]])
          beta<-coef(lm(formula,data=zdata))
          beta[1]<-0
        } else beta<-1:nrow(eresults)
        
        eresults<-cbind(eresults,beta) 


        for (i in 1:nrow(eresults)) {
          tableRow=eresults[i,]
          names(tableRow)<-c("estimate","std","t","p","beta")
          estimatesTable$addRow(rowKey=rownames(eresults)[i],tableRow)
          estimatesTable$setRow(rowKey=rownames(eresults)[i],list(name=rownames(eresults)[i]))
          }
        
        private$.populateSimple(private$.model)
        private$.prepareDescPlots(private$.model)
        private$.populateLevenes(private$.model)
        private$.populatePostHoc(data)
        private$.populateDescriptives(data)
        
      }) # suppressWarnings
    },
    .initPostHoc=function(data) {
      
      bs <- self$options$factors
      phTerms <- self$options$postHoc
      
      bsLevels <- list()
      for (i in seq_along(bs))
        bsLevels[[bs[i]]] <- levels(data[[bs[i]]])
      
      tables <- self$results$postHoc
      
      postHocRows <- list()
      
      for (ph in phTerms) {
        
        table <- tables$get(key=ph)
        
        table$setTitle(paste0('Post Hoc Comparisons - ', stringifyTerm(ph)))
        
        for (i in seq_along(ph))
          table$addColumn(name=paste0(ph[i],'1'), title=ph[i], type='text', superTitle='Comparison', combineBelow=TRUE)
        
        table$addColumn(name='sep', title='', type='text', content='-', superTitle='Comparison', format='narrow')
        
        for (i in seq_along(ph))
          table$addColumn(name=paste0(ph[i],'2'), title=ph[i], type='text', superTitle='Comparison')
        
        table$addColumn(name='md', title='Mean Difference', type='number')
        table$addColumn(name='se', title='SE', type='number')
        table$addColumn(name='df', title='df', type='number')
        table$addColumn(name='t', title='t', type='number')
        
        table$addColumn(name='pnone', title='p', type='number', format='zto,pvalue', visible="(postHocCorr:none)")
        table$addColumn(name='ptukey', title='p<sub>tukey</sub>', type='number', format='zto,pvalue', visible="(postHocCorr:tukey)")
        table$addColumn(name='pscheffe', title='p<sub>scheffe</sub>', type='number', format='zto,pvalue', visible="(postHocCorr:scheffe)")
        table$addColumn(name='pbonferroni', title='p<sub>bonferroni</sub>', type='number', format='zto,pvalue', visible="(postHocCorr:bonf)")
        table$addColumn(name='pholm', title='p<sub>holm</sub>', type='number', format='zto,pvalue', visible="(postHocCorr:holm)")
        
        combin <- expand.grid(bsLevels[rev(ph)])
        combin <- sapply(combin, as.character, simplify = 'matrix')
        if (length(ph) > 1)
          combin <- combin[,rev(1:length(combin[1,]))]
        
        comp <- list()
        iter <- 1
        for (i in 1:(length(combin[,1]) - 1)) {
          for (j in (i+1):length(combin[,1])) {
            comp[[iter]] <- list()
            comp[[iter]][[1]] <- combin[i,]
            comp[[iter]][[2]] <- combin[j,]
            
            if (j == length(combin[,1]))
              comp[[iter]][[3]] <- TRUE
            else
              comp[[iter]][[3]] <- FALSE
            
            iter <- iter + 1
          }
        }
        
        postHocRows[[composeTerm(ph)]] <- comp
        
        for (i in seq_along(comp)) {
          row <- list()
          for (c in seq_along(comp[[i]][[1]]))
            row[[paste0(names(comp[[i]][[1]][c]),'1')]] <- as.character(comp[[i]][[1]][c])
          for (c in seq_along(comp[[i]][[2]]))
            row[[paste0(names(comp[[i]][[2]][c]),'2')]] <- as.character(comp[[i]][[2]][c])
          
          table$addRow(rowKey=i, row)
          if (comp[[i]][[3]] == TRUE)
            table$addFormat(rowNo=i, col=1, Cell.END_GROUP)
        }
      }
      private$.postHocRows <- postHocRows
    },
    .populatePostHoc=function(data) {
      terms <- self$options$postHoc
      
      if (length(terms) == 0)
        return()
      
      tables <- self$results$postHoc
      
      postHocRows <- list()
      
      for (ph in terms) {
        
        table <- tables$get(key=ph)
        
        term <- jmvcore::composeTerm(ph)
        termB64 <- jmvcore::composeTerm(toB64(ph))
        
        formula <- as.formula(paste('~', term))
        
        suppressWarnings({
          
          # table$setStatus('running')
          
          referenceGrid <- lsmeans::lsmeans(private$.model, formula)
          none <- summary(pairs(referenceGrid, adjust='none'))
          tukey <- summary(pairs(referenceGrid, adjust='tukey'))
          scheffe <- summary(pairs(referenceGrid, adjust='scheffe'))
          bonferroni <- summary(pairs(referenceGrid, adjust='bonferroni'))
          holm <- summary(pairs(referenceGrid, adjust='holm'))
          
        }) # suppressWarnings
        
        resultRows <- lapply(strsplit(as.character(none$contrast), ' - '), function(x) strsplit(x, ','))
        tableRows <- private$.postHocRows[[term]]
        
        for (i in seq_along(tableRows)) {
          location <- lapply(resultRows, function(x) {
            
            c1 <- identical(x[[1]], as.character(tableRows[[i]][[1]]))
            c2 <- identical(x[[1]], as.character(tableRows[[i]][[2]]))
            c3 <- identical(x[[2]], as.character(tableRows[[i]][[1]]))
            c4 <- identical(x[[2]], as.character(tableRows[[i]][[2]]))
            if (c1 && c4)
              return(list(TRUE,FALSE))
            else if (c2 && c3)
              return(list(TRUE,TRUE))
            else
              return(list(FALSE,FALSE))
          })
          index <- which(sapply(location, function(x) return(x[[1]])))
          reverse <- location[[index]][[2]]

          row <- list()
          row[['md']] <- if(reverse) -none[index,'estimate'] else none[index,'estimate']
          row[['se']] <- none[index,'SE']
          row[['df']] <- none[index,'df']
          row[['t']] <- if(reverse) -none[index,'t.ratio'] else none[index,'t.ratio']
          
          row[['pnone']] <- none[index,'p.value']
          row[['ptukey']] <- tukey[index,'p.value']
          row[['pscheffe']] <- scheffe[index,'p.value']
          row[['pbonferroni']] <- bonferroni[index,'p.value']
          row[['pholm']] <- holm[index,'p.value']
          
          table$setRow(rowNo=i, values=row)
          private$.checkpoint()
        }
        
        table$setStatus('complete')
      }
    },
    .populateLevenes=function(model) {
      
      if ( ! self$options$homo)
        return()
      data<-model$model
      data$res<-residuals(model)
      factors <- self$options$factors
      rhs <- paste0('`', factors, '`', collapse=':')
      formula <- as.formula(paste0('`res`~', rhs))
      result <- car::leveneTest(formula, data, center="mean")
      
      table <- self$results$get('assump')$get('homo')
      
      table$setRow(rowNo=1, values=list(
        F=result[1,'F value'],
        df1=result[1,'Df'],
        df2=result[2,'Df'],
        p=result[1,'Pr(>F)']))
    },
    .populateDescriptives=function(data) {
      
      if ( ! self$options$descStats)
        return()
      
      descTable <- self$results$desc
      dep <- self$options$dep
      dependent <- data[[dep]]
      factorNames <- rev(self$options$factors)
      factors <- as.list(select(data, factorNames))
      
      means <- aggregate(dependent, by=factors, base::mean, drop=FALSE)
      sds    <- aggregate(dependent, by=factors, stats::sd, drop=FALSE)
      ns <- aggregate(dependent, by=factors, base::length, drop=FALSE)
      
      stat <- data.frame(mean=means$x, sd=sds$x, n=ns$x)
      
      for (i in seq_len(nrow(stat))) {
        values <- stat[i,]
        values[is.na(values)] <- NaN
        descTable$setRow(rowNo=i, values)
      }
      
    },
    .populateSimple=function(model) {
      
      
      variable<-self$options$simpleVariable
      moderator<-self$options$simpleModerator
      threeway<-self$options$simple3way
      data<-model$model
      simpleEffectsTables<-self$results$simpleEffects
      simpleEffectsAnovas<-self$results$simpleEffectsAnovas
      
      .fillTheFTable<-function(results,aTable) {
        ftests<-results[[2]]
        ### ftests      
        for (i in seq_len(dim(ftests)[1])) {
          r<-ftests[i,]
          row<-list(variable=r$variable,
                    term=r$level,
                    ss=r$`Sum of Sq`,
                    df=r$Df,
                    F=r$`F value`,
                    p=r$`Pr(>F)`)
          aTable$addRow(rowKey=i,row)
        }
      } #### end of .fillTheFTable
      
      .fillThePTable<-function(results,aTable) {
        params<-results[[1]]
        what<-params$level
        for (i in seq_len(dim(params)[1])) {
          r<-params[i,]
          row<-list(variable=r$variable,
                    term=r$level,
                    estimate=r$Estimate,
                    std=r$`Std. Error`,
                    t=r$`t value`,
                    p=r$`Pr(>|t|)`)
          aTable$addRow(rowKey=i,row)
          if (what!=r$level)  
            aTable$addFormat(col=1, rowNo=i,format=Cell.BEGIN_GROUP)
          what<-r$level
        }
      } ##### end of .fillThePTable
      
      if (is.null(variable) | is.null(moderator)) 
         return()

        if (is.null(threeway)) {
        results<-private$.simpleEffects(model,data,variable,moderator)
        ### ftests
        title<-paste("of",variable)
        ftable<-simpleEffectsAnovas$addItem(title)
        .fillTheFTable(results,ftable)      
        ### parameters
        ptable<-simpleEffectsTables$addItem(title)
        .fillThePTable(results,ptable)
        } else {
          data$mod2<-data[,threeway]
          if (is.factor(data$mod2)) {
            levs<-levels(data$mod2)
          } else 
               levs<-c(mean(data$mod2)+sd(data$mod2),mean(data$mod2),mean(data$mod2)-sd(data$mod2))
          for(i in seq_along(levs)) {
               data[,threeway]<-data$mod2
               if (is.factor(data$mod2))
                   contrasts(data[,threeway])<-contr.treatment(length(levs),base=i)
               else
                   data[,threeway]<-data[,threeway]-levs[i]
               results<-private$.simpleEffects(model,data,variable,moderator)
               title<-paste("computed for",threeway,"at",round(levs[i],digits = 2))
               ftable<-simpleEffectsAnovas$addItem(title)
               .fillTheFTable(results,ftable)      
               ### parameters
               ptable<-simpleEffectsTables$addItem(title)
               .fillThePTable(results,ptable)
             } 
        } # end of if (is.null(threeway)) 

    },
    .simpleEffects=function(model,data,variable,moderator){
      
      atSomeLevel<-function(moderator,level,data) {
        if (is.factor(data[,moderator])) {
          .levels<-levels(data[,moderator])
          index<-which(.levels==level,arr.ind = T)
          contrasts(data[,moderator])<-contr.treatment(length(.levels),base=index)
        } else {
          data[,moderator]<-data[,moderator]-level
        }
        form<-formula(model)
        lm(form,data=data)
      }
      
      extractEffect<-function(model,variable,moderator,level) {
        ### get the parameters
        if (is.null(model$xlevels[[variable]])) {
          varname<-variable
        }  else {
          varname<-paste(variable,model$xlevels[[variable]],sep="")
        }
        ss<-summary(model)
        ss<-as.data.frame(ss$coefficients)
        a<-as.logical(apply(sapply(varname, function(a) a==rownames(ss)),1,sum))
        ss<-ss[a,] 
        if (is.numeric(level)) level<-round(level,digits=2)
        ss$level<-paste(moderator,level,sep=" at ")
        ss$variable<-rownames(ss)
        ss
      }
      extractF<-function(model,variable,moderator,level) {
        ### get the ANOVA
        ano<-drop1(model,.~.,test = "F")
        ano<-ano[rownames(ano)==variable,c(1,2,5,6)]
        if (is.numeric(level)) level<-round(level,digits=2)
        ano$level<-paste(moderator,level,sep=" at ")
        ano$variable<-variable
        ano        
      }
      
      mod<-data[,moderator]
      if (is.factor(mod)) {
        .levels<-levels(mod)
      } else {
        .levels<-c(mean(mod)+sd(mod),mean(mod),mean(mod)-sd(mod))
      }
      params<-data.frame()
      ftests<-data.frame()
      for (i in .levels) {
        mod<-atSomeLevel(moderator,i,data)
        params<-rbind(params,extractEffect(mod,variable,moderator,i))
        ftests<-rbind(ftests,extractF(mod,variable,moderator,i))
      }
      list(params,ftests)
    },
    .contrastLabels=function(levels, type) {
      nLevels <- length(levels)
      labels <- list()
      
      if (length(levels) <= 1) {
        
        # do nothing
        
      } else if (type == 'simple') {
        
        for (i in seq_len(nLevels-1))
          labels[[i]] <- paste(levels[i+1], '-', levels[1])
        
      } else if (type == 'deviation') {
        
        all <- paste(levels, collapse=', ')
        for (i in seq_len(nLevels-1))
          labels[[i]] <- paste(levels[i], '-', all)
        
      } else if (type == 'difference') {
        
        for (i in seq_len(nLevels-1)) {
          rhs <- paste0(levels[1:i], collapse=', ')
          labels[[i]] <- paste(levels[i + 1], '-', rhs)
        }
        
      } else if (type == 'helmert') {
        
        for (i in seq_len(nLevels-1)) {
          rhs <- paste(levels[(i+1):nLevels], collapse=', ')
          labels[[i]] <- paste(levels[i], '-', rhs)
        }
        
      } else if (type == 'repeated') {
        
        for (i in seq_len(nLevels-1))
          labels[[i]] <- paste(levels[i], '-', levels[i+1])
        
      } else if (type == 'polynomial') {
        
        names <- c('linear', 'quadratic', 'cubic', 'quartic', 'quintic', 'sextic', 'septic', 'octic')
        
        for (i in seq_len(nLevels-1)) {
          if (i <= length(names)) {
            labels[[i]] <- names[i]
          } else {
            labels[[i]] <- paste('degree', i, 'polynomial')
          }
        }
      }
      
      labels
    },
    .createContrasts=function(levels, type) {
      
      nLevels <- length(levels)
      
      if (type == 'simple') {
        
        contrast <- contr.treatment(levels)
        dimnames(contrast) <- NULL

      } else if (type == 'deviation') {

        contrast <- stats::contr.sum(levels)
        dimnames(contrast) <- NULL
        
      } else if (type == 'difference') {
        
        contrast <- stats::contr.helmert(levels)
        for(i in 1:ncol(contrast))
           contrast[,i]<-contrast[,i]/(i*2)
        dimnames(contrast) <- NULL
        
      } else if (type == 'helmert') {
        
        contrast <- matrix(0, nrow=nLevels, ncol=nLevels-1)
        
        for (i in seq_len(nLevels-1)) {
          p <- (1 / (nLevels - i + 1))
          contrast[i,i] <- p * (nLevels - i)
          contrast[(i+1):nLevels,i] <- -p
        }
        
        contrast
        
      } else if (type == 'polynomial') {
        
        contrast <- stats::contr.poly(levels)
        dimnames(contrast) <- NULL
        
      } else if (type == 'repeated') {
        
        contrast <- matrix(0, nrow=nLevels, ncol=nLevels-1)
        for (i in seq_len(nLevels-1)) {
          contrast[i,  i] <- 1
          contrast[i+1,i] <- -1
        }
        
      } else {
        
        contrast <- NULL
      }
      print(contrast)
      contrast
    },
    .modelTerms=function() {
      modelTerms <- self$options$modelTerms
      if (length(modelTerms) == 0)
        modelTerms <- private$.ff()
      modelTerms
    },
    .ff=function() {
      factors <- c(self$options$factors,self$options$covs)
      if (length(factors) > 1) {
        formula <- as.formula(paste('~', paste(paste0('`', factors, '`'), collapse='*')))
        terms   <- attr(stats::terms(formula), 'term.labels')
        modelTerms <- sapply(terms, function(x) as.list(strsplit(x, ':')), USE.NAMES=FALSE)
      } else {
        modelTerms <- as.list(factors)
      }
      
      for (i in seq_along(modelTerms)) {
        term <- modelTerms[[i]]
        quoted <- grepl('^`.*`$', term)
        term[quoted] <- substring(term[quoted], 2, nchar(term[quoted])-1)
        modelTerms[[i]] <- term
      }
      
      modelTerms
    },

    .initDescPlots=function(data) {
      isAxis <- ! is.null(self$options$plotHAxis)
      isMulti <- ! is.null(self$options$plotSepPlots)
      
      self$results$get('descPlot')$setVisible( ! isMulti && isAxis)
      self$results$get('descPlots')$setVisible(isMulti)
      
      if (isMulti) {
        return()    
        sepPlotsName <- self$options$plotSepPlots
        sepPlotsVar <- data[[sepPlotsName]]
        if(is.factor(sepPlotsVar))
             sepPlotsLevels <- 1:length(levels(sepPlotsVar))
        else sepPlotsLevels <- c(1,2,3)   
        array <- self$results$descPlots
        for (level in sepPlotsLevels)
          array$addItem(level)
      }
    },
    .prepareDescPlots=function(model) {
  
      depName <- self$options$dep
      groupName <- self$options$plotHAxis
      linesName <- self$options$plotSepLines
      plotsName <- self$options$plotSepPlots
      
      ciWidth   <- self$options$ciWidth
      errorBarType <- self$options$plotError
      if (length(depName) == 0 || length(groupName) == 0)
        return()
      
      plotData<-private$.preparePlotData(model)
      
      if (self$options$plotError != 'none') {
        yAxisRange <- pretty(c(plotData$lower, plotData$upper))
      } else {
        yAxisRange <- plotData$mean
      }
      if (is.null(plotsName)) {
        
        image <- self$results$get('descPlot')
        image$setState(list(data=plotData, range=yAxisRange))
        
      } else {
        sepPlotsLevels<-levels(plotData$plots)
        array <- self$results$descPlots
        for (level in sepPlotsLevels)
          array$addItem(paste(plotsName,"=",level))

        images <- self$results$descPlots
        
        for (level in images$itemKeys) {
          image <- images$get(key=level)
          real<-gsub(paste(plotsName,"= "),"",level)
          image$setState(list(data=subset(plotData,plots==real), range=yAxisRange))
          
        }
      }
    },
    .descPlot=function(image, ggtheme, theme, ...) {
      library(ggplot2)
      if (is.null(image$state))
        return(FALSE)
      
      depName <- self$options$dep
      groupName <- self$options$plotHAxis
      linesName <- self$options$plotSepLines
      plotsName <- self$options$plotSepPlots
      
      if (self$options$plotError != 'none')
        dodge <- ggplot2::position_dodge(0.2)
      else
        dodge <- ggplot2::position_dodge(0)
      
      errorType <- ''
      if (self$options$plotError != 'none') {
        if (self$options$plotError == 'ci') {
          ciWidth <- self$options$ciWidth
          errorType <- paste0('(', ciWidth, '% CI)')
        } else {
          errorType <- '(SE)'
        }
      }
      
      if ( ! is.null(linesName)) {
        p <- ggplot2::ggplot(data=image$state$data, aes(x=group, y=mean, group=factor(lines), colour=lines)) +
          geom_line(size=.8, position=dodge) +
          labs(x=groupName, y=depName, colour=paste(linesName, errorType,sep="\n")) +
          scale_y_continuous(limits=c(min(image$state$range), max(image$state$range))) +
          ggtheme

        if (is.factor(image$state$data$group)) {
            if (self$options$plotError != 'none')
                p <- p + geom_errorbar(aes(x=group, ymin=lower, ymax=upper, width=.1, group=lines), size=.8, position=dodge)
             p <- p + geom_point(shape=21, fill='white', size=3, position=dodge)
        } else {
          if (self$options$plotError != 'none')
             p <- p + geom_ribbon(aes(x=group, ymin=lower, ymax=upper,group=lines,colour=lines,fill = lines),linetype = 0,show.legend=F, alpha=.2)          
        }
#        p <- p + labs(caption = "(Pauloo, et al. 2017)")
        print(p)
        
      } else {
        
        p <- ggplot2::ggplot(data=image$state$data) +
          labs(x=groupName, y=depName, colour=paste("", errorType)) +
          scale_colour_manual(name=paste("", errorType), values=c(colour=theme$color[1]), labels='') +
          scale_y_continuous(limits=c(min(image$state$range), max(image$state$range))) +
          ggtheme
        
          
        if (is.factor(image$state$data$group)) {
            p <- p + geom_point(aes(x=group, y=mean, colour='colour'), shape=21, fill=theme$fill[1], size=3)
            if (self$options$plotError != 'none')
               p <- p + geom_errorbar(aes(x=group, ymin=lower, ymax=upper, colour='colour', width=.1), size=.8)
        }
        else { 
            p <- p+geom_line(aes(x=group,y=mean)) 
            if (self$options$plotError != 'none')
              p <- p + geom_ribbon(aes(x=group, ymin=lower, ymax=upper),show.legend=F, alpha=.3)
            
        }
        print(p)
      }
      
      TRUE
    },


    .qqPlot=function(image, ggtheme, theme, ...) {
      library(ggplot2)
      dep <- self$options$dep
      factors <- self$options$factors
      modelTerms <- private$.modelTerms()
      model<-private$.model      
      if (is.null(model) )
              return(FALSE)
      
      data <- model$model
      residuals <- rstandard(model)
      df <- as.data.frame(qqnorm(residuals, plot.it=FALSE))
      print(ggplot2::ggplot(data=df, aes(y=y, x=x)) +
              geom_abline(slope=1, intercept=0, colour=theme$color[1]) +
              geom_point(aes(x=x,y=y), size=2, colour=theme$color[1]) +
              xlab("Theoretical Quantiles") +
              ylab("Standardized Residuals") +
              ggtheme)
      
      TRUE
    },
.preparePlotData=function(model) {
  groupName <- self$options$plotHAxis
  linesName <- self$options$plotSepLines
  plotsName <- self$options$plotSepPlots
  bars<-self$options$plotError
  
  selected<-c(groupName,linesName,plotsName)  
  vars<-all.vars(model$terms)[-1]
  
  ll<-list()
  for (v in vars) {
    if (is.factor(model$model[,v])) 
      ll[v]<-list(levels(model$model[,v]))
    else {
      if (v %in% selected) {
        if (v==groupName)
            ll[v]<-list(c(-2,-1,0,1,2)*sd(model$model[,v])+mean(model$model[,v]))
        else
            ll[v]<-list(c(-1,0,1)*sd(model$model[,v])+mean(model$model[,v]))
      }
      else 
        ll[v]<-list(0)
    }
  }
  dm<-expand.grid(ll)
  for (v in names(model$contrasts)) {
    dm[,v]<-factor(dm[,v])
  }
  mm<-predict(model,dm,interval="c",level=0.95,se.fit = T)
  dm<-as.data.frame(cbind(mm,dm))
  if (length(selected)==1) by<-list(dm[,selected])
  else by<-dm[,selected]
  lnames<-c("group","lines","plots")
  dm<-aggregate(dm[,c("fit.fit","fit.lwr","fit.upr","se.fit")],by,mean)
  if (bars=="se") {
    dm$fit.lwr<-dm$fit.fit-dm$se.fit
    dm$fit.upr<-dm$fit.fit+dm$se.fit
  }
  names(dm)<-c(lnames[1:length(selected)],c("mean","lower","upper","se"))
  if (!is.null(dm$plots) & !is.factor(model$model[,plotsName])) {
        dm$plots<-factor(dm$plots)
        levels(dm$plots)<-c("-SD","Mean","+SD")
  }

  if (!is.null(dm$lines) & !is.factor(model$model[,linesName])) {
      dm$lines<-factor(dm$lines)
      levels(dm$lines)<-c("-SD","Mean","+SD")

  }
  dm
}, # end of .preparePlotData()
    .sourcifyOption = function(option) {
      
      name <- option$name
      value <- option$value
      
      if (name == 'contrasts') {
        i <- 1
        while (i <= length(value)) {
          item <- value[[i]]
          if (item$type == 'default')
            value[[i]] <- NULL
          else
            i <- i + 1
        }
        if (length(value) == 0)
          return('')
      } else if (name == 'modelTerms') {
        if (base::identical(as.list(value), private$.ff()))
          return('')
      } else if (name == 'postHoc') {
        if (length(value) == 0)
          return('')
      }
      
      super$.sourcifyOption(option)
    })
)
