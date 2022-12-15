library(mgcv)
library(itsadug)

# get onset data from GAMMs
getDivergeData<- function(df_list, k_, stimuli_, min_sig_length, measurement){
  ## initiating new dataframe 
  divergeDf<- data.frame(Type = character(),
                         Diverge = double(),
                         Triplet = character(),
                         Syl_struct = character(),
                         Speaker = character(),
                         stringsAsFactors = FALSE
  )
  ## get divergence from GAMMs 
  row<- 1
  for (d in df_list){
    triplet<- d[1, 'Triplet']
    triplet_syl<- strsplit(as.character(d[1,'Triplet']),'')
    t<- paste(triplet_syl[[1]][-1],collapse = '')
    syl<- triplet_syl[[1]][1]
    speaker<- toString(d[1,'Speaker'])
    vPair<- eval(parse(text = paste('stimuli_$', triplet, '[,1]')))
    # triplet C6 only has V pair
    if (ncol(eval(parse(text = paste("stimuli_$", triplet))))>1){
      cPair<- eval(parse(text = paste('stimuli_$', triplet, '[,2]')))
    }
    #build GAMMs model
    formula_ <- as.formula(paste(measurement, '~ Word + s(Time, by = Word, k =', k_,') + s(Time, Rep, by= Word, bs = "fs", m = 1)', sep = ' '))
    model<- bam(formula_, data = d)
    
    model_acf<- acf_resid(model)
    
    model<- bam(formula_, data = d, rho = model_acf[2], AR.start = d$start.event)
    # get the diff data frame
    vdiff<- plot_diff(model, view = 'Time', comp = list(Word = vPair), rm.ranef = TRUE, plot = FALSE)
    vsig<- find_difference(vdiff$est, vdiff$CI, f=1, xVals=vdiff$Time)
    if (length(vsig$start) <1){
      vdiverge<- NA
    }
    else{
      vsigdurations<- vector()
      for (i in 1:length(vsig$start)){
        vsigwindow<- vsig$end[i]-vsig$start[i]
        vsigdurations<- c(vsigdurations, vsigwindow)
      }
      vsigindex<- which(vsigdurations>min_sig_length)
      if (length(vsigindex)>0){
        vdiverge<- vsig$start[vsigindex[1]]
      }
      else{
        vdiverge<- NA
      }
    }
    #asemble rows
    vrow<- list('v', vdiverge, t, syl, speaker)
    divergeDf[row,]<- vrow
    row<- row + 1
    
    if (ncol(eval(parse(text = paste("stimuli_$", triplet))))>1){
      cdiff<- plot_diff(model, view = 'Time', comp = list(Word = cPair), rm.ranef = TRUE, plot = FALSE)
      csig<- find_difference(cdiff$est, cdiff$CI, f=1, xVals=cdiff$Time)
      if (length(csig$start) <1){
        cdiverge<- NA
      }
      else{
        csigdurations<- vector()
        for (i in 1:length(csig$start)){
          csigwindow<- csig$end[i]-csig$start[i]
          csigdurations<- c(csigdurations, csigwindow)
        }
        csigindex<- which(csigdurations>min_sig_length)
        if (length(csigindex)>0){
          cdiverge<- csig$start[csigindex[1]]
        }
        else{
          cdiverge<- NA
        }
      }
      crow<- list('c', cdiverge, t, syl, speaker)
      divergeDf[row,]<- crow
      row<- row + 1
    }
  }
  return(divergeDf)
}