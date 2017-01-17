dice.read.numtaxa=function(num.partitions, num.taxa){

 numtaxacaution=0

 #convert user input num.taxa to vector, assign to object numtaxa
 if(typeof(num.taxa)!='list'){
  #protect against user misspecification in the case such that there are less num.taxa than num.partitions but more num.taxa than 1, in which case, the first num.taxa is used for all partitions
  if(length(num.taxa)<num.partitions){
   numtaxa=rep(as.integer(num.taxa[1]),num.partitions)
   if(length(num.taxa)>1){
    numtaxacaution=1
   }
  } else {
   numtaxa=as.integer(num.taxa)
   if(length(num.taxa)>num.partitions){
    numtaxacaution=1
   }
  }
 } else {
  if(length(num.taxa)==1){
   if(length(num.taxa[[1]])<num.partitions){
    numtaxa=rep(as.integer(num.taxa[[1]][1]),num.partitions)
    if(length(num.taxa[[1]])>1){
     numtaxacaution=1
    }
   } else {
    numtaxa=as.integer(num.taxa[[1]])
    if(length(numtaxa)>num.partitions){
     numtaxacaution=1
    }
   }
  } else {
   numtaxa=NULL
   #iterating through the length of num.taxa protects against weird misspecifications, such as less list elements than num.partitions but multiple items in each list element so that the total length of numeric elements = num.partitions (or >=)
   for(z in 1:length(num.taxa)){
    numtaxa=c(numtaxa,as.integer(num.taxa[[z]]))
   }
   if(length(numtaxa)>num.partitions){
    numtaxacaution=1
   }
   if(length(numtaxa)<num.partitions){
    numtaxa=rep(numtaxa[1],num.partitions)
    numtaxacaution=1
   }
  }
 }

 if(sum(numtaxa<1)>0){
  return(-1)
 } else {
  return(list(numtaxa=numtaxa,numtaxacaution=numtaxacaution))
 }

}


dice.read.psiprior=function(tao.psi.prior, epsilon.psi.prior, NE.psi.prior){

 psiprior=NULL
 psicaution=0

 for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
  #if this psi is intended to be inferred
  if(z=='tao2'&&typeof(tao.psi.prior)=='list'&&length(tao.psi.prior)>1||z=='epsilon2'&&typeof(epsilon.psi.prior)=='list'&&length(epsilon.psi.prior)>1||substr(z, nchar(z), nchar(z))!='2'&&!is.null(eval(parse(text=paste(sub('2','',z),'.psi.prior',sep=''))))){
   #convert user input psi to sorted vector, assign to object psiprior; must be sorted for while loop below that determines all zeta combinations per psi, which assumes psi is in numeric order and thus allows gaps/holes and duplicates in psi distribution
   if(substr(z,nchar(z),nchar(z))=='2'){
    psiprior=append(psiprior,list(sort(as.integer(eval(parse(text=paste(sub('2','',z),'.psi.prior',sep='')))[[2]]))))
   } else {
    if(typeof(eval(parse(text=paste(z,'.psi.prior',sep=''))))=='list'){
     psiprior=append(psiprior,list(sort(as.integer(eval(parse(text=paste(z,'.psi.prior',sep='')))[[1]]))))
    } else {
     psiprior=append(psiprior,list(sort(as.integer(eval(parse(text=paste(z,'.psi.prior',sep='')))))))
    }
   }
   if(sum(psiprior[[length(psiprior)]]<0)>0){
    names(psiprior)[[length(psiprior)]]='psitest'
   } else {
    names(psiprior)[[length(psiprior)]]=z
   }
  }
  if(z=='NE'&&typeof('NE.psi.prior')=='list'&&length('NE.psi.prior')>1||substr(z,nchar(z),nchar(z))=='2'&&typeof(eval(parse(text=paste(sub('2','',z),'.psi.prior',sep=''))))=='list'&&length(eval(parse(text=paste(sub('2','',z),'.psi.prior',sep=''))))>2){
   psicaution=1
  }
 }
 
  return(append(psiprior,list(psicaution=psicaution)))

}


dice.read.zetaprior=function(num.partitions, numtaxa, idiosyncratic, psiprior, tao.zeta.prior, tao2.zeta.prior, epsilon.zeta.prior, epsilon2.zeta.prior, NE.zeta.prior, min.net.zeta.per.pulse, max.net.zeta.per.pulse, min.net.zeta.total, max.net.zeta.total){

 zetaprior=NULL
 zetatestpart=0
 minzetatestpp=0
 maxzetatestpp=0
 minzetatesttot=0
 maxzetatesttot=0
 idiosynctest=0
 zetacaution=0

 for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
  #if this psi is intended to be inferred
  if(z%in%names(psiprior)&&max(psiprior[[paste(z,sep='')]])>0){
   minzetatestpptemp=0
   maxzetatestpptemp=0
   #convert user input zeta to list, assign to object zetaprior
   zetaprior=append(zetaprior,list(NULL))
   temp=eval(parse(text=paste(z,'.zeta.prior', sep='')))
   if(typeof(temp)!='list'){
    temp=list(temp)
   }
   temporary=NULL
   for(y in 1:length(temp)){
    temporary=c(temporary,temp[[y]])
   }
   if(is.null(temporary)||sum(temporary>1)>0||sum(temporary<0)>0){
    names(zetaprior)[[length(zetaprior)]]='zetatest'
   } else {
    #protect against user misspecification in the case such that there are less zeta distributions than num.partitions but more zeta distributions than 1, in which case, the first zeta distribution is used for all partitions
    if(length(temp)<num.partitions){
     for(y in 1:num.partitions){
      #convert units of zetaprior from proportion to absolute number of taxa
      zetaprior[[length(zetaprior)]]=append(zetaprior[[length(zetaprior)]], list(as.integer(temp[[1]]*sum(numtaxa[['numtaxa']]))))
      if((min(psiprior[[paste(z,sep='')]])*min(temp[[1]]))>(numtaxa[['numtaxa']][y]/sum(numtaxa[['numtaxa']]))){
       zetatestpart=1
      }
     }
     minzetatestpptemp=minzetatestpptemp+(max(temp[[1]])*num.partitions)
     maxzetatestpptemp=maxzetatestpptemp+(min(temp[[1]])*num.partitions)
     if(length(temp)>1){
      zetacaution=1
     }
    } else {
     for(y in 1:num.partitions){
      #convert units of zetaprior from proportion to absolute number of taxa
      zetaprior[[length(zetaprior)]]=append(zetaprior[[length(zetaprior)]], list(as.integer(temp[[y]]*sum(numtaxa[['numtaxa']]))))
      if((min(psiprior[[paste(z,sep='')]])*min(temp[[y]]))>(numtaxa[['numtaxa']][y]/sum(numtaxa[['numtaxa']]))){
       zetatestpart=1
      }
      minzetatestpptemp=minzetatestpptemp+max(temp[[y]])
      maxzetatestpptemp=maxzetatestpptemp+min(temp[[y]])
     }
     if(length(temp)>num.partitions){
      zetacaution=1
     }
    }
    if(zetatestpart==1){
     names(zetaprior)[[length(zetaprior)]]='zetatestpart'
    } else {
     names(zetaprior)[[length(zetaprior)]]=z
     if(!is.null(min.net.zeta.per.pulse)&&minzetatestpptemp<min.net.zeta.per.pulse){
      minzetatestpp=1
     }
     if(!is.null(max.net.zeta.per.pulse)&&maxzetatestpptemp>max.net.zeta.per.pulse){
      maxzetatestpp=1
     }
     if(!is.null(min.net.zeta.total)&&(minzetatestpptemp*max(psiprior[[paste(z,sep='')]]))<min.net.zeta.total){
      minzetatesttot=1
     }
     if(!is.null(max.net.zeta.total)&&(maxzetatestpptemp*min(psiprior[[paste(z,sep='')]]))>max.net.zeta.total){
      maxzetatesttot=1
     }
     if(idiosyncratic==F){
      if(!is.null(min.net.zeta.total)&&min.net.zeta.total!=1){
       minzetatesttot=1
      }
      if(!is.null(max.net.zeta.total)&&max.net.zeta.total!=1){
       maxzetatesttot=1
      }
      if(!is.null(min.net.zeta.per.pulse)){
       mintemp=min.net.zeta.per.pulse
      } else {
       mintemp=NULL
      }
      if(!is.null(max.net.zeta.per.pulse)){
       maxtemp=max.net.zeta.per.pulse
      } else {
       maxtemp=NULL
      }
      marker=0
      for(y in 1:num.partitions){
       tag=0
       x=0
       while(x<length(unique(psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0]))&&tag==0){
        x=x+1
        zetacombs=NULL
        for(w in 1:unique(psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0])[x]){
         zetacombs=append(zetacombs,list(zetaprior[[paste(z,sep='')]][[y]]))
        }
        zetacombs=expand.grid(zetacombs)
        w=0
        while(w<nrow(zetacombs)&&tag==0){
         w=w+1
         if(sum(zetacombs[w,])==numtaxa[['numtaxa']][y]){
          tag=1
          marker=marker+1
         }
        }
       }
      }
      if(marker!=num.partitions){
       idiosynctest=1
      }
     }
    }
   }
  }
 }

 return(append(zetaprior,list(minzetatestpp=minzetatestpp,maxzetatestpp=maxzetatestpp,minzetatesttot=minzetatesttot,maxzetatesttot=maxzetatesttot,idiosynctest=idiosynctest,zetacaution=zetacaution)))

}



dice.read.sharedprior=function(psiprior, tao.shared.prior, tao2.shared.prior, epsilon.shared.prior, epsilon2.shared.prior, NE.shared.prior){

 sharedprior=NULL
 sharedcaution=0

 for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
  #if this psi is intended to be inferred
  if(z%in%names(psiprior)&&max(psiprior[[paste(z,sep='')]])>0){
   #convert user input shared to list, assign to object sharedprior
   sharedprior=append(sharedprior,list(NULL))
   temp=eval(parse(text=paste(z,'.shared.prior',sep='')))
   if(typeof(temp)!='list'){
    temp=list(temp)
   }
   temporary=NULL
   for(y in 1:length(temp)){
    temporary=c(temporary,temp[[y]])
   }
   if(is.null(temporary)||sum(temporary<=0)>0){
    names(sharedprior)[[length(sharedprior)]]='sharedtest'
   } else {
    names(sharedprior)[[length(sharedprior)]]=z
    #protect against user misspecification in the case such that there are less shared distributions than max psi value but more shared distributions than 1, in which case, the first shared distribution is used for all partitions
    #cannot allow user option of adding another list element to shared prior object indicating a separate idiosyncratic distribution in the case of a single shared distribution used for all pulses with max psi value > 1; potential for indistinguishability in the case of max psi value = 2, where indicating 2 different shared distributions for each pulse with no idiosyncratic distribution and indiciating 1 shared distribution for all pulses with 1 idiosyncratic distribution will both result in a list of 2 distributions in the sharedprior object; thus, to utilize a separate idiosyncratic distribution, the idiosyncratic prior specification must be used (DICE will look for this first), or can add as another list element when the number of shared distributions = max psi value, such that with an additional list element to shared prior object indicating a separate idiosyncratic distribution, length of shared prior object = max psi value + 1 (or at least >=; in cases of >, additional distributions are ignored, unless additional distributions >= num.partitions, in which case, each additional distribution will be the assigned idio distribution for its respective partition) (DICE will look for this second), or otherwise DICE defaults to use the first shared prior distribution
    if(length(temp)<max(psiprior[[paste(z,sep='')]])){
     if(substr(z,1,7)=='epsilon'){
      sharedprior[[paste(z,sep='')]]=list(temp[[1]])
     } else {
      sharedprior[[paste(z,sep='')]]=list(as.integer(temp[[1]]))
     }
     if(length(temp)>1){
      sharedcaution=1
     }
    } else {
     for(y in 1:max(psiprior[[paste(z,sep='')]])){
      if(substr(z,1,7)=='epsilon'){
       sharedprior[[paste(z,sep='')]]=append(sharedprior[[paste(z,sep='')]],list(temp[[y]]))
      } else {
       sharedprior[[paste(z,sep='')]]=append(sharedprior[[paste(z,sep='')]],list(as.integer(temp[[y]])))
      }
     }
     #no sharedcaution here since sharedpriors can appropriately be larger if including idiosyncratic distributions
    }
   }
  }
 }

 return(append(sharedprior,list(sharedcaution=sharedcaution)))

}


dice.read.buff=function(psiprior, tao.buffer=0, tao2.buffer=NULL, epsilon.buffer=0, epsilon2.buffer=NULL, NE.buffer=0){

 buff=NULL

 for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
  #if this psi is intended to be inferred
  if(z%in%names(psiprior)){
   #convert user input buffer to list, assign to object buff
   if(!is.null(eval(parse(text=paste(z,'buffer',sep='.'))))&&typeof(eval(parse(text=paste(z,'buffer',sep='.'))))=='double'&&eval(parse(text=paste(z,'buffer',sep='.')))>=0||!is.null(eval(parse(text=paste(z,'buffer',sep='.'))))&&typeof(eval(parse(text=paste(z,'buffer',sep='.'))))=='integer'&&eval(parse(text=paste(z,'buffer',sep='.')))>=0||!is.null(eval(parse(text=paste(z,'buffer',sep='.'))))&&typeof(eval(parse(text=paste(z,'buffer',sep='.'))))=='closure'){
    buff=append(buff,list(as.integer(eval(parse(text=paste(z,'buffer',sep='.'))))))
    names(buff)[[length(buff)]]=z
   } else {
    if(substr(z,nchar(z),nchar(z))=='2'){
     buff=append(buff,list(buff[[substr(z,1,(nchar(z)-1))]]))
     names(buff)[[length(buff)]]=z
    } else {
     buff=append(buff,list(-1))
     names(buff)[[length(buff)]]='bufftest'
    }
   }
  }
 }

 return(buff)

}



dice.read.numchanges=function(num.partitions, num.changes){

 numchangescaution=0

 #convert user input num.changes to vector, assign to object numchanges
 if(typeof(num.changes)!='list'){
  #protect against user misspecification in the case such that there are less num.changes than num.partitions but more num.changes than 1, in which case, the first num.changes is used for all partitions
  if(length(num.changes)<num.partitions){
   numchanges=rep(num.changes[1],num.partitions)
   if(length(num.changes)>1){
    numchangescaution=1
   }
  } else {
   numchanges=num.changes
   if(length(num.changes)>num.partitions){
    numchangescaution=1
   }
  }
 } else {
  if(length(num.changes)==1){
   if(length(num.changes[[1]])<num.partitions){
    numchanges=rep(num.changes[[1]][1],num.partitions)
    if(length(num.changes[[1]])>1){
     numchangescaution=1
    }
   } else {
    numchanges=num.changes[[1]]
    if(length(numchanges)>num.partitions){
     numchangescaution=1
    }
   }
  } else {
   numchanges=NULL
   #iterating through the length of num.changes protects against weird misspecifications, such as less list elements than num.partitions but multiple items in each list element so that the total length of numeric elements = num.partitions (or >=)
   for(z in 1:length(num.changes)){
    numchanges=c(numchanges,num.changes[[z]])
   }
   if(length(numchanges)>num.partitions){
    numchangescaution=1
   }
   if(length(numchanges)<num.partitions){
    numchanges=rep(numchanges[1],num.partitions)
    numchangescaution=1
   }
  }
 }
 if(sum(numchanges>2)>0){
  numchanges[numchanges>2]=2
  numchangescaution=1
 }

 if(sum(numchanges<1)>0){
  return(-1)
 } else {
  return(list(numchanges=numchanges,numchangescaution=numchangescaution))
 }

}


dice.read.linking=function(linked.param, attached.hyper){

 linking=NULL

 #convert user inputs linked.param and attached.hyper to vector, assign to object linking
 for(z in c('linked.param','attached.hyper')){
  if(typeof(eval(parse(text=z)))!='list'){
   linking=append(linking,list(eval(parse(text=z))))
  } else {
   linking=append(linking,list(NULL))
   #iterating through the length of linked.param and attached.hyper protects against weird misspecifications, such as multiple items in each list element
   for(y in 1:length(eval(parse(text=z)))){
    linking[[length(linking)]]=c(linking[[length(linking)]],eval(parse(text=z))[[y]])
   }
  }
 }

 if(sum(linking[[1]]=='tao')>1 || sum(linking[[1]]=='tao2')>1 || sum(linking[[1]]=='epsilon')>1 || sum(linking[[1]]=='epsilon2')>1 || sum(linking[[1]]=='NE')>1 || sum(c(sum(linking[[1]]=='tao'), sum(linking[[1]]=='tao2'), sum(linking[[1]]=='epsilon'), sum(linking[[1]]=='epsilon2'), sum(linking[[1]]=='NE'))) != length(linking[[1]]) || sum(c(sum(linking[[2]]=='tao'), sum(linking[[2]]=='tao2'), sum(linking[[2]]=='epsilon'), sum(linking[[2]]=='epsilon2'), sum(linking[[2]]=='NE'))) != length(linking[[2]]) || length(linking[[1]])!=length(linking[[2]])){
  return(-1)
 } else {
  return(linking)
 }

}


dice.read.idioprior=function(num.partitions, idiosyncratic, psiprior, tao.shared.prior, tao2.shared.prior, epsilon.shared.prior, epsilon2.shared.prior, NE.shared.prior, tao.idio.prior, tao2.idio.prior, epsilon.idio.prior, epsilon2.idio.prior, NE.idio.prior, linking, numchanges, anchor.prior, change.prior){

 idioprior=NULL
 idiocaution=0

 for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
  #if hyperparameter is being inferred and taxa are allowed to act idiosyncratically (for tao2, change.prior must also be NOT specified by user, and thus for an inferred 2nd event with idiosyncratic taxa, change.prior has 1st priority, idio.prior has 2nd priority, and then shared.prior has 3rd and last priority; anchor.prior doesn't make sense here, since anchor.prior is practically a substitute for an inferred 2nd event and cannot act as a coinciding model just for the idiosyncratic taxa as it is pulse-based for the 1st event), OR hyperparamter is not being inferred but parameter is being simulated as a nuisance (for tao2, anchor.prior and change.prior must also both be NOT specified by user, and thus for a nuisance 2nd event, anchor.prior has 1st priority, change.prior has 2nd priority, idio.prior has 3rd priority, and then shared.prior has 4th and last priority) and the parameter is not linked to another hyperparameter, then idio.prior (or a substituted shared.prior if idio.prior is not specified by user) is used for all taxa that are either idiosyncratic or nuisance for this parameter
  if(z%in%names(psiprior)&&idiosyncratic==T&&z!='tao2'||z%in%names(psiprior)&&idiosyncratic==T&&is.null(change.prior)||!z%in%names(psiprior)&&substr(z,nchar(z),nchar(z))!='2'&&!z%in%linking[[1]]||!z%in%names(psiprior)&&2%in%numchanges[['numchanges']]&&z=='epsilon2'&&!z%in%linking[[1]]||!z%in%names(psiprior)&&2%in%numchanges[['numchanges']]&&is.null(anchor.prior)&&is.null(change.prior)&&!z%in%linking[[1]]){
   #convert user input idio (or shared) to list, assign to object idioprior
   idioprior=append(idioprior,list(NULL))
   if(!is.null(eval(parse(text=paste(z,'.idio.prior',sep=''))))){
    #user input idio.prior
    temp=eval(parse(text=paste(z,'.idio.prior',sep='')))
   } else {
    #user input shared.prior
    temp=eval(parse(text=paste(z,'.shared.prior',sep='')))
    #if additional idio prior distributions are specified after true shared prior distributions in shared.prior, then taking out the true shared prior distributions; must do it here because after this if else loop, it will be impossible to distinguish between distributions from idio.prior vs. shared.prior, and thus in cases of length(idioprior) > num.partitions, whether this is due to misspecification of idio.prior (in which case, if discarding the first max(psivalue) distributions, then the latter idio.prior distributions would be used, which is not the desired default against misspecification) or additional specification of idio prior distributions in shared.prior
    if(typeof(temp)=='list'&&length(temp)>max(psiprior[[paste(z,sep='')]])){
     temporary=NULL
     for(y in (max(psiprior[[paste(z,sep='')]])+1):length(temp)){
      temporary=append(temporary,temp[[y]])
     }
     temp=temporary
    }
   }
   if(typeof(temp)!='list'){
    temp=list(temp)
   }
   temporary=NULL
   for(y in 1:length(temp)){
    temporary=c(temporary,temp[[y]])
   }
   if(is.null(temporary)||sum(temporary<=0)>0){
    names(idioprior)[[length(idioprior)]]='idiotest'
   } else {
    names(idioprior)[[length(idioprior)]]=z
    #protect against user misspecification in the case such that there are less idio distributions than num.partitions but more idio distributions than 1, in which case, the first idio distribution is used for all partitions
    if(length(temp)<num.partitions){
     for(y in 1:num.partitions){
      if(substr(z,1,7)=='epsilon'){
       idioprior[[paste(z,sep='')]]=append(idioprior[[paste(z,sep='')]],list(temp[[1]]))
      } else {
       idioprior[[paste(z,sep='')]]=append(idioprior[[paste(z,sep='')]],list(as.integer(temp[[1]])))
      }
     }
     if(length(temp)>1){
      idiocaution=1
     }
    } else {
     for(y in 1:num.partitions){
      if(substr(z,1,7)=='epsilon'){
       idioprior[[paste(z,sep='')]]=append(idioprior[[paste(z,sep='')]],list(temp[[y]]))
      } else {
       idioprior[[paste(z,sep='')]]=append(idioprior[[paste(z,sep='')]],list(as.integer(temp[[y]])))
      }
     }
     if(length(temp)>num.partitions){
      idiocaution=1
     }
    }
   }
  }
  #if hyperparamter is not being inferred but parameter is being simulated as a nuisance (for tao2, anchor.prior and change.prior must also both be NOT specified by user, and thus for a nuisance 2nd event, anchor.prior has 1st priority, change.prior has 2nd priority, idio.prior has 3rd priority, and then shared.prior has 4th and last priority) and the parameter is linked to another hyperparameter, then shared.prior (or a substituted idio.prior if shared.prior is not specified by user) is used for all taxa that are nuisance for this parameter and shared for attached hyperparameter, and idio.prior (or a substituted shared.prior if idio.prior is not specified by user) is used for all taxa that are nuisance for this parameter and not shared for attached hyperparameter
  if(!z%in%names(psiprior)&&substr(z,nchar(z),nchar(z))!='2'&&z%in%linking[[1]]||!z%in%names(psiprior)&&2%in%numchanges[['numchanges']]&&z=='epsilon2'&&z%in%linking[[1]]||!z%in%names(psiprior)&&2%in%numchanges[['numchanges']]&&is.null(anchor.prior)&&is.null(change.prior)&&z%in%linking[[1]]){
   #convert user input shared and/or idio to list, assign to object idioprior
   idioprior=append(idioprior,list(NULL))
   idioprior=append(idioprior,list(NULL))
   if(!is.null(eval(parse(text=paste(z,'.shared.prior',sep=''))))&&!is.null(eval(parse(text=paste(z,'.idio.prior',sep=''))))){
    temp1=eval(parse(text=paste(z,'.shared.prior',sep='')))
    temp2=eval(parse(text=paste(z,'.idio.prior',sep='')))
   } else {
    if(!is.null(eval(parse(text=paste(z,'.shared.prior',sep=''))))){
     temp1=eval(parse(text=paste(z,'.shared.prior',sep='')))
     temp2=eval(parse(text=paste(z,'.shared.prior',sep='')))
     #if additional idio prior distributions are specified after true shared prior distributions in shared.prior, then taking out the true shared prior distributions; must do it here because after this if else loop, it will be impossible to distinguish between distributions from idio.prior vs. shared.prior
     if(typeof(temp2)=='list'&&length(temp2)>max(psiprior[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]])){
      temporary=NULL
      for(y in (max(psiprior[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]])+1):length(temp2)){
       temporary=append(temporary,temp2[[y]])
      }
      temp2=temporary
     }
    } else {
     temp1=eval(parse(text=paste(z,'.idio.prior',sep='')))
     temp2=eval(parse(text=paste(z,'.idio.prior',sep='')))
     #if additional idio prior distributions are specified after true shared prior distributions in shared.prior, then taking out the true shared prior distributions; must do it here because after this if else loop, it will be impossible to distinguish between distributions from idio.prior vs. shared.prior
     if(typeof(temp1)=='list'&&length(temp1)>num.partitions){
      temporary=NULL
      for(y in (num.partitions+1):length(temp1)){
       temporary=append(temporary,temp1[[y]])
      }
      temp1=temporary
     }
    }
   }
   if(typeof(temp1)!='list'){
    temp1=list(temp1)
   }
   if(typeof(temp2)!='list'){
    temp2=list(temp2)
   }
   temporary=NULL
   temporarily=append(temp1,temp2)
   for(y in 1:length(temporarily)){
    temporary=c(temporary,temporarily[[y]])
   }
   if(is.null(temporary)||sum(temporary<=0)>0){
    names(idioprior)[[length(idioprior)]]='idiotest'
   } else {
    names(idioprior)[[(length(idioprior)-1)]]=paste(z,'.shared',sep='')
    names(idioprior)[[length(idioprior)]]=z
    if(max(psiprior[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]])>0){
     #protect against user misspecification in the case such that there are less shared distributions than max(psiprior) of attached hyperparameter but more shared distributions than 1, in which case, the first shared distribution is used for all partitions
     if(length(temp1)<max(psiprior[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]])){
      for(y in 1:max(psiprior[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]])){
       if(substr(z,1,7)=='epsilon'){
        idioprior[[paste(z,'.shared',sep='')]]=append(idioprior[[paste(z,'.shared',sep='')]], list(temp1[[1]]))
       } else {
        idioprior[[paste(z,'.shared',sep='')]]=append(idioprior[[paste(z,'.shared',sep='')]], list(as.integer(temp1[[1]])))
       }
      }
      if(length(temp1)>1){
       idiocaution=1
      }
     } else {
      for(y in 1:max(psiprior[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]])){
       if(substr(z,1,7)=='epsilon'){
        idioprior[[paste(z,'.shared',sep='')]]=append(idioprior[[paste(z,'.shared',sep='')]], list(temp1[[y]]))
       } else {
        idioprior[[paste(z,'.shared',sep='')]]=append(idioprior[[paste(z,'.shared',sep='')]], list(as.integer(temp1[[y]])))
       }
      }
      if(length(temp1)>max(psiprior[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]])){
       idiocaution=1
      }
     }
    }
    #protect against user misspecification in the case such that there are less idio distributions than num.partitions but more idio distributions than 1, in which case, the first idio distribution is used for all partitions
    if(length(temp2)<num.partitions){
     for(y in 1:num.partitions){
      if(substr(z,1,7)=='epsilon'){
       idioprior[[paste(z,sep='')]]=append(idioprior[[paste(z,sep='')]],list(temp2[[1]]))
      } else {
       idioprior[[paste(z,sep='')]]=append(idioprior[[paste(z,sep='')]],list(as.integer(temp2[[1]])))
      }
     }
     if(length(temp2)>1){
      idiocaution=1
     }
    } else {
     for(y in 1:num.partitions){
      if(substr(z,1,7)=='epsilon'){
       idioprior[[paste(z,sep='')]]=append(idioprior[[paste(z,sep='')]],list(temp2[[y]]))
      } else {
       idioprior[[paste(z,sep='')]]=append(idioprior[[paste(z,sep='')]],list(as.integer(temp2[[y]])))
      }
     }
     if(length(temp2)>num.partitions){
      idiocaution=1
     }
    }
   }
  }
 }

 return(append(idioprior,list(idiocaution=idiocaution)))

}


dice.read.fixing=function(linking, fix.linked.param.across.pulses, fix.linked.param.per.pulse){

 fixing=NULL
 fixingcaution=0

 if(!is.null(linking)){
  #convert user inputs fix.linked.param.across.pulses and fix.linked.param.per.pulse to vector, assign to object fixing
  for(z in c('fix.linked.param.across.pulses', 'fix.linked.param.per.pulse')){
   if(typeof(eval(parse(text=z)))!='list'){
    #protect against user misspecification in the case such that there are less fix.linked.param.across.pulses or fix.linked.param.per.pulse than length(linking) but more than 1, in which case, the first is used for all partitions
    if(length(eval(parse(text=z)))<length(linking[[1]])){
     fixing=append(fixing,list(rep(eval(parse(text=z))[1]),length(linking[[1]])))
     if(length(eval(parse(text=z)))>1){
      fixingcaution=1
     }
    } else {
     fixing=append(fixing,list(eval(parse(text=z))))
     if(length(eval(parse(text=z)))>length(linking[[1]])){
      fixingcaution=1
     }
    }
   } else {
    if(length(eval(parse(text=z)))==1){
     if(length(eval(parse(text=z))[[1]])<length(linking[[1]])){
      fixing=append(fixing,list(rep(eval(parse(text=z))[[1]][1]),length(linking[[1]])))
      if(length(eval(parse(text=z))[[1]])>1){
       fixingcaution=1
      }
     } else {
      fixing=append(fixing,list(eval(parse(text=z))[[1]]))
      if(length(eval(parse(text=z))[[1]])>length(linking[[1]])){
       fixingcaution=1
      }
     }
    } else {
     fixing=append(fixing,list(NULL))
     #iterating through the length protects against weird misspecifications, such as less list elements than length(linking) but multiple items in each list element so that the total length of numeric elements = length(linking) (or >=)
     for(y in 1:length(eval(parse(text=z)))){
      fixing[[length(fixing)]]=c(fixing[[length(fixing)]],eval(parse(text=z))[[y]])
     }
     if(length(fixing[[length(fixing)]])>length(linking[[1]])){
      fixingcaution=1
     }
     if(length(fixing[[length(fixing)]])<length(linking[[1]])){
      fixing[[length(fixing)]]=rep(fixing[[length(fixing)]][1],length(linking[[1]]))
      fixingcaution=1
     }
    }
   }
  }
 }

 return(append(fixing,list(fixingcaution=fixingcaution)))

}


dice.read.flipvector=function(num.partitions, flip){

 flipvectorcaution=0

 #convert user input flip to vector, assign to object flipvector
 if(typeof(flip)!='list'){
  #protect against user misspecification in the case such that there are less flip than num.partitions but more flip than 1, in which case, the first flip is used for all partitions
  if(length(flip)<num.partitions){
   flipvector=rep(flip[1],num.partitions)
   if(length(flipvector)>1){
    flipvectorcaution=1
   }
  } else {
   flipvector=flip
   if(length(flipvector)>num.partitions){
    flipvectorcaution=1
   }
  }
 } else {
  if(length(flip)==1){
   if(length(flip[[1]])<num.partitions){
    flipvector=rep(flip[[1]][1],num.partitions)
    if(length(flip[[1]])>1){
     flipvectorcaution=1
    }
   } else {
    flipvector=flip[[1]]
    if(length(flipvector)>num.partitions){
     flipvectorcaution=1
    }
   }
  } else {
   flipvector=NULL
   #iterating through the length of flip protects against weird misspecifications, such as less list elements than num.partitions but multiple items in each list element so that the total length of numeric elements = num.partitions (or >=)
   for(z in 1:length(flip)){
    flipvector=c(flipvector,flip[[z]])
   }
   if(length(flipvector)>num.partitions){
    flipvectorcaution=1
   }
   if(length(flipvector)<num.partitions){
    flipvector=rep(flipvector[1],num.partitions)
    flipvectorcaution=1
   }
  }
 }

 return(list(flipvector=flipvector,flipvectorcaution=flipvectorcaution))

}


#if a 2 event model is used for at least one partition, tao2 is not being inferred, and anchor.prior is specified by user, then the 'anchor' model is used, where ALL taxa will experience a 2 event model, and the timing of the older event is anchored by a set time delta; this anchor is the same across all taxa within the same pulse; each pulse thus draws from an anchor prior; the value of this anchor per pulse will be inferred; each pulse may have a different distribution defined; idiosyncratic taxa may draw from yet another distribution, though all idiosyncratic taxa will have independent draws from each other from the same distribution (in this idiosyncratic case, it acts identically as the change.prior)
dice.read.anchorprior=function(psipriortao, anchor.prior){

 anchorprior=NULL
 anchorcaution=0

 #convert user input anchor to list, assign to object anchorprior
 temp=anchor.prior
 if(typeof(temp)!='list'){
  temp=list(temp)
 }
 temporary=NULL
 for(y in 1:length(temp)){
  temporary=c(temporary,temp[[y]])
 }
 if(sum(temporary<1)>0){
  anchorprior=list(anchorprior)
  names(anchorprior)='anchortest'
 } else {
  #protect against user misspecification in the case such that there are less anchor distributions than max(psivalue) but more anchor distributions than 1, in which case, the first anchor distribution is used for all pulses
  if(max(psipriortao)>0){
   if(length(temp)<max(psipriortao)){
    for(y in 1:(max(psipriortao))){
     anchorprior=append(anchorprior,list(as.integer(temp[[1]])))
    }
    if(length(temp)>1){
     anchorcaution=1
    }
   } else {
    for(y in 1:(max(psipriortao))){
     anchorprior=append(anchorprior,list(as.integer(temp[[y]])))
    }
    if(length(temp)>(max(psipriortao)+1)){
     anchorcaution=1
    }
   }
  }
  #if an extra distribution is specified for anchor prior (i.e. length(anchorpriortemp)>=max(psivalue)+1), then the next distribution following the max(psivalue) distribution will be used for idiosyncratic taxa; otherwise, the first distribution will be used for idiosyncratic taxa
  if(length(temp)>=(max(psipriortao)+1)){
   anchorprior=append(anchorprior,list(as.integer(temp[[(max(psipriortao)+1)]])))
  } else {
   anchorprior=append(anchorprior,list(as.integer(temp[[1]])))
   anchorcaution=1
  }
 }

 return(list(anchorprior=anchorprior,anchorcaution=anchorcaution))

}


#if a 2 event model is used for at least one partition, tao2 is not being inferred, anchor.prior is NOT specified by user, and change.prior is specified by user (thus for a nuisance 2nd event, anchor.prior has 1st priority, then change.prior has 2nd priority), OR tao2 IS being inferred, taxa are allowed to act idiosyncratically, and change.prior is specified by user (thus for an inferred 2nd event with idiosyncratic taxa, change.prior has 1st priority; anchor.prior doesn't make sense here, since anchor.prior is practically a substitute for an inferred 2nd event and cannot act as a coinciding model just for the idiosyncratic taxa as it is pulse-based for the 1st event), then the "change" model is used, where all taxa that are in a partition experiencing a 2 event model and are either nuisance or idiosyncratic in the 2nd event will independently draw from a prior a value that indicates the time delta between the 1st and 2nd event
dice.read.changeprior=function(num.partitions, change.prior){

 changeprior=NULL
 changecaution=0

 #convert user input change to list, assign to object changeprior
 temp=change.prior
 if(typeof(temp)!='list'){
  temp=list(temp)
 }
 temporary=NULL
 for(y in 1:length(temp)){
  temporary=c(temporary,temp[[y]])
 }
 if(sum(temporary<1)>0){
  changeprior=list(changeprior)
  names(changeprior)='changetest'
 } else {
  #protect against user misspecification in the case such that there are less change distributions than num.partitions but more change distributions than 1, in which case, the first change distribution is used for all partitions
  if(length(temp)<num.partitions){
   for(y in 1:num.partitions){
    changeprior=append(changeprior,list(as.integer(temp[[1]])))
   }
   if(length(temp)>1){
    changecaution=1
   }
  } else {
   for(y in 1:num.partitions){
    changeprior=append(changeprior,list(as.integer(temp[[y]])))
   }
   if(length(temp)>num.partitions){
    changecaution=1
   }
  }
 }

 return(list(changeprior=changeprior,changecaution=changecaution))

}


dice.read.exp.grow=function(num.partitions, exponential.growth.rate.prior, exponential.growth.rate.prior2){

 exp.grow=NULL
 exp.growcaution=0

 exp.dists=NULL
 if(!is.null(exponential.growth.rate.prior)){
  exp.dists=c(exp.dists,1)
 }
 if(!is.null(exponential.growth.rate.prior2)){
  exp.dists=c(exp.dists,2)
 }
 exp.names=c('exponential.growth.rate.prior','exponential.growth.rate.prior2')

 #convert user input(s) exponential.growth.rate(2) to list, assign to object exp.grow
 for(z in exp.dists){
  exp.grow=append(exp.grow,list(NULL))
  names(exp.grow)[[length(exp.grow)]]=exp.names[z]
  temp=eval(parse(text=exp.names[z]))
  if(typeof(temp)!='list'){
   temp=list(temp)
  }
  #protect against user misspecification in the case such that there are less exponential.growth.rate distributions than num.partitions but more exponential.growth.rate distributions than 1, in which case, the first exponential.growth.rate distribution is used for all partitions
  if(length(temp)<num.partitions){
   for(y in 1:num.partitions){
    exp.grow[[exp.names[z]]]=append(exp.grow[[exp.names[z]]],list(temp[[1]]))
   }
   if(length(temp)>1){
    exp.growcaution=1
   }
  } else {
   for(y in 1:num.partitions){
    exp.grow[[exp.names[z]]]=append(exp.grow[[exp.names[z]]],list(temp[[y]]))
   }
   if(length(temp)>num.partitions){
    exp.growcaution=1
   }
  }
 }

 return(append(exp.grow,list(exp.growcaution=exp.growcaution)))

}



dice.read.fold=function(folded, num.partitions){

 foldcaution=0

 #convert user input folded to vector, assign to object fold
 if(typeof(folded)!='list'){
  #protect against user misspecification in the case such that there are less folded than num.partitions but more folded than 1, in which case, the first folded is used for all partitions
  if(length(folded)<num.partitions){
   fold=rep(folded[1],num.partitions)
   if(length(folded)>1){
    foldcaution=1
   }
  } else {
   fold=folded
   if(length(folded)>num.partitions){
    foldcaution=1
   }
  }
 } else {
  if(length(folded)==1){
   if(length(folded[[1]])<num.partitions){
    fold=rep(folded[[1]][1],num.partitions)
    if(length(folded[[1]])>1){
     foldcaution=1
    }
   } else {
    fold=folded[[1]]
    if(length(folded[[1]])>num.partitions){
     foldcaution=1
    }
   }
  } else {
   fold=NULL
   #iterating through the length of folded protects against weird misspecifications, such as less list elements than num.partitions but multiple items in each list element so that the total length of numeric elements = num.partitions (or >=)
   for(z in 1:length(folded)){
    fold=c(fold,folded[[z]])
   }
   if(length(fold)>num.partitions){
    foldcaution=1
   }
   if(length(fold)<num.partitions){
    fold=rep(fold[1],num.partitions)
    foldcaution=1
   }
  }
 }

 return(list(fold=fold,foldcaution=foldcaution))

}


dice.read.numhaps=function(num.partitions, num.haploid.samples){

 numhapscaution=0

 #convert user input num.haploid.samples to vector, assign to object numhaps
 if(typeof(num.haploid.samples)!='list'){
  #protect against user misspecification in the case such that there are less numsites than num.partitions but more numsites than 1, in which case, the first numsites is used for all partitions
  if(length(num.haploid.samples)<num.partitions){
   numhaps=rep(num.haploid.samples[1],num.partitions)
   if(length(num.haploid.samples)>1){
    numhapscaution=1
   }
  } else {
   numhaps=num.haploid.samples
   if(length(num.haploid.samples)>num.partitions){
    numhapscaution=1
   }
  }
 } else {
  if(length(num.haploid.samples)==1){
   if(length(num.haploid.samples[[1]])<num.partitions){
    numhaps=rep(num.haploid.samples[[1]][1],num.partitions)
    if(length(num.haploid.samples[[1]])>1){
     numhapscaution=1
    }
   } else {
    numhaps=num.haploid.samples[[1]]
    if(length(numhaps)>num.partitions){
     numhapscaution=1
    }
   }
  } else {
   numhaps=NULL
   #iterating through the length of num.haploid.samples protects against weird misspecifications, such as less list elements than num.partitions but multiple items in each list element so that the total length of numeric elements = num.partitions (or >=)
   for(z in 1:length(num.haploid.samples)){
    numhaps=c(numhaps,num.haploid.samples[[z]])
   }
   if(length(numhaps)>num.partitions){
    numhapscaution=1
   }
   if(length(numhaps)<num.partitions){
    numhapscaution=1
    numhaps=rep(numhaps[1],num.partitions)
   }
  }
 }

 if(sum(numhaps<1)>0){
  return(-1)
 } else {
  return(list(numhaps=numhaps,numhapscaution=numhapscaution))
 }

}


dice.read.samtimes=function(num.partitions, sampling.times){

 samtimescaution=0

 #convert user input sampling.times to vector, assign to object samtimes
 if(typeof(sampling.times)!='list'){
  #protect against user misspecification in the case such that there are less numsites than num.partitions but more numsites than 1, in which case, the first numsites is used for all partitions
  if(length(sampling.times)<num.partitions){
   samtimes=rep(sampling.times[1],num.partitions)
   if(length(sampling.times)>1){
    samtimescaution=1
   }
  } else {
   samtimes=sampling.times
   if(length(sampling.times)>num.partitions){
    samtimescaution=1
   }
  }
 } else {
  if(length(sampling.times)==1){
   if(length(sampling.times[[1]])<num.partitions){
    samtimes=rep(sampling.times[[1]][1],num.partitions)
    if(length(sampling.times[[1]])>1){
     samtimescaution=1
    }
   } else {
    samtimes=sampling.times[[1]]
    if(length(samtimes)>num.partitions){
     samtimescaution=1
    }
   }
  } else {
   samtimes=NULL
   #iterating through the length of num.haploid.samples protects against weird misspecifications, such as less list elements than num.partitions but multiple items in each list element so that the total length of numeric elements = num.partitions (or >=)
   for(z in 1:length(sampling.times)){
    samtimes=c(samtimes,sampling.times[[z]])
   }
   if(length(samtimes)>num.partitions){
    samtimescaution=1
   }
   if(length(samtimes)<num.partitions){
    samtimescaution=1
    samtimes=rep(samtimes[1],num.partitions)
   }
  }
 }

 if(sum(samtimes<0)>0){
  return(-1)
 } else {
  return(list(samtimes=samtimes,samtimescaution=samtimescaution))
 }

}


dice.read.numsites=function(num.partitions, num.ind.sites, num.SNPs, length.seq){

 numsitescaution=0

 #convert user input num.ind.sites, num.SNPs, or length.seq to vector, assign to object numsites; since numsites is being defined first here before reformatting, this is why the caution assessed before manipulating numsites into the proper format
 numsites=NULL
 if(!is.null(num.ind.sites)){
  numsites=num.ind.sites
 } else {
  if(!is.null(num.SNPs)){
   numsites=num.SNPs
  } else {
   if(!is.null(length.seq)){
    numsites=length.seq
   }
  }
 }
 if(typeof(numsites)!='list'){
  #protect against user misspecification in the case such that there are less numsites than num.partitions but more numsites than 1, in which case, the first numsites is used for all partitions
  if(length(numsites)<num.partitions){
   if(length(numsites)>1){
    numsitescaution=1
   }
   numsites=rep(numsites[1],num.partitions)
  } else {
   if(length(numsites)>num.partitions){
    numsitescaution=1
   }
  }
 } else {
  if(length(numsites)==1){
   if(length(numsites[[1]])<num.partitions){
    if(length(numsites)[[1]]>1){
     numsitescaution=1
    }
    numsites=rep(numsites[[1]][1],num.partitions)
   } else {
    if(length(numsites[[1]])>num.partitions){
     numsitescaution=1
    }
    numsites=numsites[[1]]
   }
  } else {
   temp=NULL
   #iterating through the length of numsites protects against weird misspecifications, such as less list elements than num.partitions but multiple items in each list element so that the total length of numeric elements = num.partitions (or >=)
   for(z in 1:length(numsites)){
    temp=c(temp,numsites[[z]])
   }
   numsites=temp
   if(length(numsites)>num.partitions){
    numsitescaution=1
   }
   if(length(numsites)<num.partitions){
    numsitescaution=1
    numsites=rep(numsites[1],num.partitions)
   }
  }
 }

 if(sum(numsites<1)>0||is.null(numsites)){
  return(-1)
 } else {
  return(list(numsites=numsites,numsitescaution=numsitescaution))
 }

}


dice.read.mutrate=function(num.partitions, mut.rate){

 mutrate=NULL
 mutratecaution=0

 #convert user input mut.rate to vector, assign to object mutrate
 temp=mut.rate
 if(typeof(temp)!='list'){
  temp=list(temp)
 }
 temporary=NULL
 for(y in 1:length(temp)){
  temporary=c(temporary,temp[[y]])
 }
 if(sum(temporary<=0)>0){
  return(-1)
 } else {
  #protect against user misspecification in the case such that there are less mut.rate than num.partitions but more mut.rate than 1, in which case, the first mut.rate is used for all partitions
  if(length(temp)<num.partitions){
   for(z in 1:num.partitions){
    mutrate=append(mutrate, list(temp[[1]]))
   }
   if(length(temp)>1){
    mutratecaution=1
   }
  } else {
   for(z in 1:num.partitions){
    mutrate=append(mutrate, list(temp[[z]]))
   }
   if(length(temp)>num.partitions){
    mutratecaution=1
   }
  }
  return(append(mutrate,list(mutratecaution=mutratecaution)))
 }

}


dice.read.gentimes=function(num.partitions, gen.times){

 gentimescaution=0

 #convert user input gen.times to vector, assign to object gentimes
 if(typeof(gen.times)!='list'){
  #user allowed to specify as many generation times as total number of taxa, in which case, ordering of generation times does not matter within partitions (since parameters values are randomly assigned and aSFS is order-independent) but does between partitions
  if(length(gen.times)>=sum(numtaxa[['numtaxa']])){
   gentimes=gen.times
   if(length(gen.times)>sum(numtaxa[['numtaxa']])){
    gentimescaution=1
   }
  } else {
   #user also allowed to specify as many generation times as number of partitions, in which case, all taxa within a partition will have the same generation time
   if(length(gen.times)>=num.partitions){
    gentimes=NULL
    for(a in 1:num.partitions){
     gentimes=c(gentimes,rep(gen.times[a],numtaxa[['numtaxa']][a]))
    }
    if(length(gen.times)>num.partitions){
     gentimescaution=1
    }
   #protect against user misspecification in the case such that there are less gen.times than num.partitions or num.taxa but more gen.times than 1, in which case, the first gen.times is used for all partitions
   } else {
    gentimes=rep(gen.times[1],sum(numtaxa[['numtaxa']]))
    if(length(gen.times)>1){
     gentimescaution=1
    }
   }
  }
 } else {
  if(length(gen.times)>=sum(numtaxa[['numtaxa']])){
   gentimes=NULL
   for(a in 1:length(gen.times)){
    gentimes=c(gentimes,gen.times[[a]][1])
    if(length(gen.times[[a]]>1)){
     gentimescaution=1
    }
   }
   if(length(gen.times)>sum(numtaxa[['numtaxa']])){
    gentimescaution=1
   }
  } else {
   if(length(gen.times)>=num.partitions){
    gentimes=NULL
    for(a in 1:num.partitions){
     if(length(gen.times[[a]])>=numtaxa[['numtaxa']][a]){
      gentimes=c(gentimes,gen.times[[a]][1:numtaxa[['numtaxa']][a]])
      if(length(gen.times[[a]]>numtaxa[['numtaxa']][a])){
       gentimescaution=1
      }
     } else {
      gentimes=c(gentimes,rep(gen.times[[a]][1],numtaxa[['numtaxa']][a]))
      if(length(gen.times[[a]]>1)){
       gentimescaution=1
      }
     }
    }
    if(length(gen.times)>num.partitions){
     gentimescaution=1
    }
   } else {
    if(length(gen.times[[1]])>=sum(numtaxa[['numtaxa']])){
     gentimes=gen.times[[1]]
     if(length(gentimes>sum(numtaxa[['numtaxa']]))){
      gentimescaution=1
     }
    } else {
     if(length(gen.times[[1]])>=num.partitions){
      gentimes=NULL
      for(a in 1:num.partitions){
       gentimes=c(gentimes,rep(gen.times[[1]][a],numtaxa[['numtaxa']][a]))
      }
      if(length(gen.times[[1]])>num.partitions){
       gentimescaution=1
      }
     } else {
      gentimes=rep(gen.times[[1]][1],sum(numtaxa[['numtaxa']]))
      if(length(gen.times[[1]])>1){
       gentimescaution=1
      }
     }
     if(length(gen.times)>1){
      gentimescaution=1
     }
    }
   }
  }
 }

 if(sum(gentimes<=0)>0){
  return(-1)
 } else {
  return(list(gentimes=gentimes,gentimescaution=gentimescaution))
 }

}



dice.read.inputdir=function(input.directory){

 #convert user input input.directory to vector, assign to object inputdir
 if(typeof(input.directory)!='list'){
  inputdir=input.directory
 } else {
  inputdir=NULL
  #iterating through the length of input.directory protects against weird misspecifications, such as multiple items in each list element
  for(z in 1:length(input.directory)){
   inputdir=c(temp,input.directory[[z]])
  }
 }

 for(z in 1:length(inputdir)){
  if(substr(inputdir[z],nchar(inputdir[z]),nchar(inputdir[z]))!='/'){
   inputdir[z]=paste(inputdir[z],'/',sep='')
  }
 }

 return(inputdir)

}


dice.read.inputs=function(input.base, input.files, numtaxa){

 inputscaution=0

 #convert user input input.base or input.files to vector, assign to object inputs
 if(!is.null(input.base)){
  if(typeof(input.base)!='list'){
   temp=input.base
  } else {
   temp=input.base[[1]]
  }
  inputs=NULL
  for(z in 1:sum(numtaxa[['numtaxa']])){
   inputs=c(inputs,paste(temp,z,sep=''))
  }
 } else {
  if(typeof(input.files)!='list'){
   inputs=input.files
  } else {
   inputs=NULL
   #iterating through the length of input.files protects against weird misspecifications, such as less list elements than numtaxa but multiple items in each list element so that the total length of numeric elements = numtaxa (or >=)
   for(z in 1:length(input.files)){
    inputs=c(inputs,input.files[[z]])
   }
  }
  #protect against user misspecification in the case such that there are less input.files than numtaxa
  if(length(inputs)<sum(numtaxa[['numtaxa']])){
   inputs=-1
  }
  if(length(inputs)>sum(numtaxa[['numtaxa']])){
   inputscaution=1
   inputs=inputs[1:sum(numtaxa[['numtaxa']])]
  }
 }

 if(-1%in%inputs){
  return(-1)
 } else {
  return(list(inputs=inputs,inputscaution=inputscaution))
 }

}


dice.read.afc=function(remove.afclasses, num.partitions){

 afc=NULL
 afccaution=0

 #convert user input remove.afclasses to vector, assign to object afc
 temp=remove.afclasses
 if(typeof(temp)!='list'){
  temp=list(temp)
 }
 temporary=NULL
 for(y in 1:length(temp)){
  temporary=c(temporary,temp[[y]])
 }
 if(sum(temporary<0)>0){
  return(-1)
 } else {
  #protect against user misspecification in the case such that there are less remove.afclasses than num.partitions but more remove.afclasses than 1, in which case, the first remove.afclasses is used for all partitions
  if(length(temp)<num.partitions){
   for(z in 1:num.partitions){
    afc=append(afc, list(as.integer(temp[[1]])))
   }
   if(length(temp)>1){
    afccaution=1
   }
  } else {
   for(z in 1:num.partitions){
    afc=append(afc, list(as.integer(temp[[z]])))
   }
   if(length(temp)>num.partitions){
    afccaution=1
   }
  }
  return(append(afc,afccaution=afccaution))
 }

}




build.dirichlet.prior=function(num.partitions=1, num.taxa, idiosyncratic=T, tao.psi.prior=NULL, epsilon.psi.prior=NULL, NE.psi.prior=NULL, tao.zeta.prior=NULL, tao2.zeta.prior=NULL, epsilon.zeta.prior=NULL, epsilon2.zeta.prior=NULL, NE.zeta.prior=NULL, min.net.zeta.per.pulse=NULL, max.net.zeta.per.pulse=NULL, min.net.zeta.total=NULL, max.net.zeta.total=NULL, fsc2template.path, fsc2path, input.directory, input.base, input.files, output.directory, folded, remove.afclasses, append.sims, keep.taxa.draws, output.hyper.draws, output.taxa.draws, keep.fsc2.files, num.sims, messages.sims, num.haploid.samples, sampling.times, num.ind.sites, num.SNPs, length.seq, mut.rate, gen.times, idiosyncratic.buffer, dirichlet.process, tao.shared.prior, tao2.shared.prior, epsilon.shared.prior, epsilon2.shared.prior, NE.shared.prior, tao.buffer, tao2.buffer, epsilon.buffer, epsilon2.buffer, NE.buffer, net.zeta.per.pulse, net.zeta.total, mean.tao.shared, disp.index.tao.shared, mean.tao, disp.index.tao, mean.tao2.shared, disp.index.tao2.shared, mean.tao2, disp.index.tao2, mean.epsilon.shared, disp.index.epsilon.shared, mean.epsilon, disp.index.epsilon, mean.epsilon2.shared, disp.index.epsilon2.shared, mean.epsilon2, disp.index.epsilon2, mean.NE.shared, disp.index.NE.shared, mean.NE, disp.index.NE, tao.idio.prior, tao2.idio.prior, epsilon.idio.prior, epsilon2.idio.prior, NE.idio.prior, linked.param, attached.hyper, fix.linked.param.across.pulses, fix.linked.param.per.pulse, num.changes, flip, anchor.prior, change.prior, exponential.growth.rate.prior, exponential.growth.rate.prior2, dirichlet.object, roll.object, play.object){

 output=NULL

 numtaxa=dice.read.numtaxa(num.partitions=num.partitions, num.taxa=num.taxa)

 if(-1%in%numtaxa){
  print("You have included at least 1 non-positive number in the 'num.taxa' object; DICE cannot continue")
 }

 psiprior=dice.read.psiprior(tao.psi.prior=tao.psi.prior, epsilon.psi.prior=epsilon.psi.prior, NE.psi.prior=NE.psi.prior)

 if(is.null(psiprior)){
  print("You have not specified any 'psi' priors; DICE cannot continue")
 }

 if('psitest'%in%names(psiprior)){
  print("At least one of the 'psi' priors that you have specified contains at least 1 value that is negative; DICE cannot continue")
 }

 zetaprior=dice.read.zetaprior(num.partitions=num.partitions, numtaxa=numtaxa, idiosyncratic=idiosyncratic, psiprior=psiprior, tao.zeta.prior=tao.zeta.prior, tao2.zeta.prior=tao2.zeta.prior, epsilon.zeta.prior=epsilon.zeta.prior, epsilon2.zeta.prior=epsilon2.zeta.prior, NE.zeta.prior=NE.zeta.prior, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total)

 if('zetatest'%in%names(zetaprior)){
  print("For at least one of the 'psi' priors that you have specified, you have not specified a 'zeta' prior and/or have specified an improper 'zeta' prior containing at least 1 value that is above 1.0 and/or negative; DICE cannot continue")
 }

 if('zetatestpart'%in%names(zetaprior)){
  print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the number of taxa that you have specified per partition; DICE cannot continue")
 }

 if(zetaprior[['minzetatestpp']]==1){
  print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the minimum zeta value per pulse that you have specified; DICE cannot continue")
 }

 if(zetaprior[['maxzetatestpp']]==1){
  print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the maximum zeta value per pulse that you have specified; DICE cannot continue")
 }

 if(zetaprior[['minzetatesttot']]==1){
  print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the minimum zeta value across all pulses (i.e. in total among the dataset) that you have specified; if you have specified 'idiosyncratic==F' (i.e. all taxa must be in a pulse), then the minimum zeta value across all pulses must either be NULL or equal to 1; DICE cannot continue")
 }

 if(zetaprior[['maxzetatesttot']]==1){
  print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the maximum zeta value across all pulses (i.e. in total among the dataset) that you have specified; if you have specified 'idiosyncratic==F' (i.e. all taxa must be in a pulse), then the maximum zeta value across all pulses must either be NULL or equal to 1; DICE cannot continue")
 }

 if(zetaprior[['idiosynctest']]==1){
  print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on you specifying 'idiosyncratic==F' (i.e. all taxa must be in a pulse), the respective 'psi' prior that you have specified, and the number of taxa that you have specified; DICE cannot continue")
 }

 if(!-1%in%numtaxa&&!is.null(psiprior)&&!'psitest'%in%names(psiprior)&&!'zetatest'%in%names(zetaprior)&&!'zetatestpart'%in%names(zetaprior)&&zetaprior[['minzetatestpp']]==0&&zetaprior[['maxzetatestpp']]==0&&zetaprior[['minzetatesttot']]==0&&zetaprior[['maxzetatesttot']]==0&&zetaprior[['idiosynctest']]==0){

  if(numtaxa[['numtaxacaution']]==1){
   print("Caution: You have specified an inappropriate number of vector elements for 'num.taxa'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'num.taxa', with x = 'num.partitions', or if not applicable, with the first element of 'num.taxa' used for all partitions")
  }

  if(psiprior[['psicaution']]==1){
   print("Caution: For at least one of the 'psi' priors that you have specified, you have included more than the allowable number of distributions (1 for Ne and 2 for tao and epsilon); DICE will continue with only the first distributions up to the allowable amount")
  }

  if(zetaprior[['zetacaution']]==1){
   print("Caution: You have specified an inappropriate number of distributions for at least one of your 'zeta' priors; the number of distributions should be = 'num.partitions' or 1; DICE will continue with only the first x distributions of 'zeta', with x = 'num.partitions', or if not applicable, with the first distribution of 'zeta' used for all partitions")
  }

  for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
   #if this psi is intended to be inferred
   if(z%in%names(psiprior)){
    #protect against user misspecification in the case such that there are more zeta distributions than appropriate (either 1 or num.partitions), as well as undoing duplication of the first element in the case of length = 1; whereas in other objects, the first element is duplicated to length of num.partitions in the case of length < num.partitions, or the latter elements will simply just be ignored in the case of length > num.partitions, here, since the length of the vector will undergo this somewhat computationally expensive combinatorics process (and it is purposely the length of the vector becuase of in the case that there is only one distribution for all partitions; whereas in other objects this one element would just get duplicated to the number of partitions, and thus the object up to the length of num.partitions would always be in consideration with latter elements simply ignored, here it is kept as just a single element before the combinatorics process so that the identical process would not have to be duplicated up to the length of num.partitions and thus creating a waste of resources), it is more efficient to restrict the vector now
    zetapriortemp=eval(parse(text=paste(z, '.zeta.prior', sep='')))
    if(typeof(zetapriortemp)!='list'){
     zetapriortemp=list(zetapriortemp)
    }
    if(length(zetapriortemp)<num.partitions){
     zetaprior[[paste(z,sep='')]]=list(zetaprior[[paste(z,sep='')]][[1]])
    }
    if(length(zetapriortemp)>num.partitions){
     for(y in (num.partitions+1):length(zetapriortemp)){
      zetaprior[[paste(z,sep='')]][[y]]=NULL
     }
    }
    combos=NULL
    #z will not be in zetaprior if psiprior==0 only
    if(paste(z,sep='')%in%names(zetaprior)){
     #construct all combinations of zeta prior draws, including for distinct zeta prior distributions per each partition, per each possible psi prior draw
     for(y in 1:length(zetaprior[[paste(z,sep='')]])){
      combos=append(combos,list(NULL))
      combinations=matrix(zetaprior[[paste(z,sep='')]][[y]],ncol=1)
      counter=1
      for(x in psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0]){
       while(counter<x){
        counter=counter+1
        combination=combinations
        combinations=NULL
        for(w in 1:nrow(combination)){
         for(v in zetaprior[[paste(z,sep='')]][[y]]){
          combinations=rbind(combinations,c(combination[w,],v))
         }
        }
       }
       combos[[y]]=append(combos[[y]],list(combinations))
      }
     }
     #allow a single zeta distribution to be used by all partitions
     if(num.partitions>length(combos)){
      for(y in 2:num.partitions){
       combos=append(combos,list(combos[[1]]))
      }
     }
    }
    #add psi=0 if part of psi prior distribution
    if(0 %in% psiprior[[paste(z,sep='')]]){
     dppprior=matrix(rep(as.integer(0),num.partitions),nrow=1)
     dpppriordraw=rep(as.integer(1),sum(psiprior[[paste(z,sep='')]]==0))
    } else {
     dppprior=NULL
     dpppriordraw=NULL
    }
    #filter out combinations of zeta prior draws that do not meet user specifications, assign valid combinations to object dppprior and the according number of draws as determined by the combinatorics to object dpppriordraw, which then becomes the dirichlet-process prior distribution
    #protect against user specification of a psi prior of only 0(s)
    if(length(psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0])>0){
     #could also catalog via a three number system (partitions number, psi number, row number); this might be slower, as it's 3x the memory, though it would reduce the memory from duplicating the combos list object in cases when there's only one zeta distribution used by multiple partitions; however, this duplication process should be quite quick, and it taking more memory isn't that big of an issue since it's not being manipulated (i.e. adding onto, re-arranging, etc.) after duplication and used just as a library, and doing the catalog via single counter number should be quick since it just uses for loops and counting to convert the single counter number to a three number catalog index, and then directly accesses matrices via that, rather than storing all the three number indices and manipulating these three number index data matrices
     counter=0
     #combination of zeta prior draws across partitions must be of equivalent psi
     for(y in 1:length(psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0])){
      dpppriorperpsi=NULL
      switch=0
      for(x in 1:num.partitions){
       dpppriorperpsi=append(dpppriorperpsi,list(NULL))
       for(w in 1:nrow(combos[[x]][[y]])){
        counter=counter+1
        #combination of zeta prior draws must be of appropriate number of taxa per partition, which is determined by the idiosyncratic logical value; by verifying every combination is of appropriate number of taxa per partition, automatically any combination of draws among partitions will be of appropriate number of total taxa
        if(idiosyncratic==T&&sum(combos[[x]][[y]][w,])<=numtaxa[['numtaxa']][x]&&switch==0||idiosyncratic==F&&sum(combos[[x]][[y]][w,])==numtaxa[['numtaxa']][x]&&switch==0){
         dpppriorperpsi[[x]]=c(dpppriorperpsi[[x]],as.integer(counter))
        }
       }
       if(is.null(dpppriorperpsi[[x]])&&switch==0){
        switch=1
       }
      }
      if(switch==0){
       dpppriorperpsi=as.matrix(expand.grid(dpppriorperpsi))
       factorial.library=NULL
       for(x in 1:nrow(dpppriorperpsi)){
        flag=0
        #combination of zeta prior draws across partitions must be of appropriate min and max net zeta per pulse and in total
        if(!is.null(min.net.zeta.per.pulse)||!is.null(max.net.zeta.per.pulse)||!is.null(min.net.zeta.total)||!is.null(max.net.zeta.total)){
         temp=NULL
         trigger=0
         clock=0
         w=0
         while(w<length(psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0])&&trigger<num.partitions){
          w=w+1
          v=0
          while(v<num.partitions&&trigger<num.partitions){
           v=v+1
           u=0
           while(u<nrow(combos[[v]][[w]])&&trigger<num.partitions){
            u=u+1
            clock=clock+1
            if(clock%in%dpppriorperpsi[x,]){
             trigger=trigger+1
             temp=append(temp,list(combos[[v]][[w]][u,]))
            }
           }
          }
         }
         #min and max net zeta per pulse
         temporary=0
         for(w in 1:psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0][y]){
          tentative=0
          for(v in 1:num.partitions){
           tentative=tentative+temp[[v]][w]
          }
          #convert user input min/max.net.zeta.per.pulse to a single numeric value
          for(v in c('min','max')){
           if(!is.null(eval(parse(text=paste(v,'.net.zeta.per.pulse',sep=''))))){
            interim=eval(parse(text=paste(v,'.net.zeta.per.pulse',sep='')))
            if(typeof(interim)!='list'){
             if(length(interim)==1){
              interim=interim*sum(numtaxa[['numtaxa']])
             } else {
              interim=interim[w]*sum(numtaxa[['numtaxa']])
             }
            } else {
             if(length(interim)==1){
              if(length(interim[[1]])==1){
               interim=interim[[1]]*sum(numtaxa[['numtaxa']])
              } else {
               interim=interim[[1]][w]*sum(numtaxa[['numtaxa']])
              }
             } else {
              interim=interim[[w]][1]*sum(numtaxa[['numtaxa']])
             }
            }
            if(v=='min'&&tentative<interim||v=='max'&&tentative>interim){
             flag=1
            }
           }
          }
          temporary=temporary+tentative
         }
         #min and max net zeta in total across all pulses
         if(!is.null(min.net.zeta.total)&&temporary<(min.net.zeta.total*sum(numtaxa[['numtaxa']]))){
          flag=1
         }
         if(!is.null(max.net.zeta.total)&&temporary>(max.net.zeta.total*sum(numtaxa[['numtaxa']]))){
          flag=1
         }
         #remove draws that are already in factorial.library from object temp to match object temp below in else statement
         if(flag==0){
          temporary=NULL
          for(w in c(1:num.partitions)[!dpppriorperpsi[x,]%in%factorial.library[,1]]){
           temporary=append(temporary,list(temp[[w]]))
          }
          temp=temporary
         }
        } else {
         #assign zeta values for draws that are not in factorial.library to object temp
         card=sum(!dpppriorperpsi[x,]%in%factorial.library[,1])
         temp=NULL
         trigger=0
         clock=0
         w=0
         while(w<length(psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0])&&trigger<card){
          w=w+1
          v=0
          while(v<num.partitions&&trigger<card){
           v=v+1
           u=0
           while(u<nrow(combos[[v]][[w]])&&trigger<card){
            u=u+1
            clock=clock+1
            if(clock%in%dpppriorperpsi[x,][!dpppriorperpsi[x,]%in%factorial.library[,1]]){
             trigger=trigger+1
             temp=append(temp,list(combos[[v]][[w]][u,]))
            }
           }
          }
         }
        }
        if(flag==0){
         #assign valid combinations to object dppprior
         dppprior=rbind(dppprior,dpppriorperpsi[x,])
         #determine the according number of draws via combinatorics; utilize and build factorial.library for efficiency
         temporary=NULL
         v=0
         for(w in 1:num.partitions){
          if(dpppriorperpsi[x,w]%in%factorial.library[,1]){
           temporary=c(temporary,factorial.library[factorial.library[,1]==dpppriorperpsi[x,w],2])
          } else {
           v=v+1
           temporary=c(temporary,factorial(numtaxa[['numtaxa']][w]))
           for(u in 1:length(temp[[v]])){
            temporary[w]=temporary[w]/(factorial(temp[[v]][u]))
           }
           temporary[w]=temporary[w]/(factorial(numtaxa[['numtaxa']][w]-sum(temp[[v]])))
           factorial.library=rbind(factorial.library,c(dpppriorperpsi[x,w],temporary[w]))
          }
         }
         temporary=prod(temporary)
         #assign according number of draws to object dpppriordraw, which then becomes the dirichlet-process prior distribution
         dpppriordraw=c(dpppriordraw,rep(as.integer(nrow(dppprior)),temporary))
        }
       }
      }
     }
    }
    output=append(output,list(dpppriordraw,dppprior,combos))
    #dirichlet-process prior distribution; each draw specifies a row in the indices matrix
    names(output)[[length(output)-2]]=paste('dpp', z, sep='.')
    #indices matrix that lists the appropriate indexed zeta value(s) per partition; each row is a different draw, and each column is an index number corresponding to a set of zeta values in the zeta matrix for that respective partition
    names(output)[[length(output)-1]]=paste('indices', z, sep='.')
    #zeta matrix that lists the zeta value(s) for each valid draw, organized by partition (list elements), then psi (list elements within list elements), then combinations of zeta values (rows of ncol = value of respective psi), but indexed by psi, then partition, then combinations of zeta values (i.e. counting up by combinations of zeta values (smallest unit), then by partition, then by psi (largest unit))
    names(output)[[length(output)]]=paste('zeta', z, sep='.')
   }
  }
 return(output)
 }
 
}



roll.dice=function(num.sims=1000000, num.partitions=1, num.taxa, idiosyncratic=T, dirichlet.process=F, tao.psi.prior=NULL, epsilon.psi.prior=NULL, NE.psi.prior=NULL, tao.zeta.prior=NULL, tao2.zeta.prior=NULL, epsilon.zeta.prior=NULL, epsilon2.zeta.prior=NULL, NE.zeta.prior=NULL, min.net.zeta.per.pulse=NULL, max.net.zeta.per.pulse=NULL, min.net.zeta.total=NULL, max.net.zeta.total=NULL, tao.shared.prior=NULL, tao2.shared.prior=NULL, epsilon.shared.prior=NULL, epsilon2.shared.prior=NULL, NE.shared.prior=NULL, tao.buffer=0, tao2.buffer=NULL, epsilon.buffer=0, epsilon2.buffer=NULL, NE.buffer=0, dirichlet.object=NULL, fsc2template.path, fsc2path, input.directory, input.base, input.files, output.directory, folded, remove.afclasses, append.sims, keep.taxa.draws, output.hyper.draws, output.taxa.draws, keep.fsc2.files, messages.sims, num.haploid.samples, sampling.times, num.ind.sites, num.SNPs, length.seq, mut.rate, gen.times, idiosyncratic.buffer, net.zeta.per.pulse, net.zeta.total, mean.tao.shared, disp.index.tao.shared, mean.tao, disp.index.tao, mean.tao2.shared, disp.index.tao2.shared, mean.tao2, disp.index.tao2, mean.epsilon.shared, disp.index.epsilon.shared, mean.epsilon, disp.index.epsilon, mean.epsilon2.shared, disp.index.epsilon2.shared, mean.epsilon2, disp.index.epsilon2, mean.NE.shared, disp.index.NE.shared, mean.NE, disp.index.NE, tao.idio.prior, tao2.idio.prior, epsilon.idio.prior, epsilon2.idio.prior, NE.idio.prior, linked.param, attached.hyper, fix.linked.param.across.pulses, fix.linked.param.per.pulse, num.changes, flip, anchor.prior, change.prior, exponential.growth.rate.prior, exponential.growth.rate.prior2, roll.object, play.object){

 crush=0
 messages=0
 if(dirichlet.process==T&&is.null(dirichlet.object)){
  dirichlet.object=build.dirichlet.prior(num.partitions=num.partitions, num.taxa=num.taxa, idiosyncratic=idiosyncratic, tao.psi.prior=tao.psi.prior, epsilon.psi.prior=epsilon.psi.prior, NE.psi.prior=NE.psi.prior, tao.zeta.prior=tao.zeta.prior, tao2.zeta.prior=tao2.zeta.prior, epsilon.zeta.prior=epsilon.zeta.prior, epsilon2.zeta.prior=epsilon2.zeta.prior, NE.zeta.prior=NE.zeta.prior, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total)
  if(is.null(dirichlet.object)){
   crush=1
  }
 messages=1
 }

 if(crush==0){

  numtaxa=dice.read.numtaxa(num.partitions=num.partitions, num.taxa=num.taxa)

  if(-1%in%numtaxa){
   print("You have included at least 1 non-positive number in the 'num.taxa' object; DICE cannot continue")
  }

  psiprior=dice.read.psiprior(tao.psi.prior=tao.psi.prior, epsilon.psi.prior=epsilon.psi.prior, NE.psi.prior=NE.psi.prior)

  if(is.null(psiprior)){
   print("You have not specified any 'psi' priors; DICE cannot continue")
  }

  if('psitest'%in%names(psiprior)){
   print("At least one of the 'psi' priors that you have specified contains at least 1 value that is negative; DICE cannot continue")
  }

  zetaprior=dice.read.zetaprior(num.partitions=num.partitions, numtaxa=numtaxa, idiosyncratic=idiosyncratic, psiprior=psiprior, tao.zeta.prior=tao.zeta.prior, tao2.zeta.prior=tao2.zeta.prior, epsilon.zeta.prior=epsilon.zeta.prior, epsilon2.zeta.prior=epsilon2.zeta.prior, NE.zeta.prior=NE.zeta.prior, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total)

  if('zetatest'%in%names(zetaprior)){
   print("For at least one of the 'psi' priors that you have specified, you have not specified a 'zeta' prior and/or have specified an improper 'zeta' prior containing at least 1 value that is above 1.0 and/or negative; DICE cannot continue")
  }

  if('zetatestpart'%in%names(zetaprior)){
   print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the number of taxa that you have specified per partition; DICE cannot continue")
  }

  if(zetaprior[['minzetatestpp']]==1){
   print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the minimum zeta value per pulse that you have specified; DICE cannot continue")
  }

  if(zetaprior[['maxzetatestpp']]==1){
   print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the maximum zeta value per pulse that you have specified; DICE cannot continue")
  }

  if(zetaprior[['minzetatesttot']]==1){
   print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the minimum zeta value across all pulses (i.e. in total among the dataset) that you have specified; if you have specified 'idiosyncratic==F' (i.e. all taxa must be in a pulse), then the minimum zeta value across all pulses must either be NULL or equal to 1; DICE cannot continue")
  }

  if(zetaprior[['maxzetatesttot']]==1){
   print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the maximum zeta value across all pulses (i.e. in total among the dataset) that you have specified; if you have specified 'idiosyncratic==F' (i.e. all taxa must be in a pulse), then the maximum zeta value across all pulses must either be NULL or equal to 1; DICE cannot continue")
  }

  if(zetaprior[['idiosynctest']]==1){
   print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on you specifying 'idiosyncratic==F' (i.e. all taxa must be in a pulse), the respective 'psi' prior that you have specified, and the number of taxa that you have specified; DICE cannot continue")
  }

  sharedprior=dice.read.sharedprior(psiprior=psiprior, tao.shared.prior=tao.shared.prior, tao2.shared.prior=tao2.shared.prior, epsilon.shared.prior=epsilon.shared.prior, epsilon2.shared.prior=epsilon2.shared.prior, NE.shared.prior=NE.shared.prior)

  if('sharedtest'%in%names(sharedprior)){
   print("For at least one of the 'psi' priors that you have specified, you have not specified a 'shared' prior and/or have specified an improper 'shared' prior containing at least 1 value that is non-positive; DICE cannot continue")
  }

  buff=dice.read.buff(psiprior=psiprior, tao.buffer=tao.buffer, tao2.buffer=tao2.buffer, epsilon.buffer=epsilon.buffer, epsilon2.buffer=epsilon2.buffer, NE.buffer=NE.buffer)

  sign=0
  if('bufftest'%in%names(buff)){
   print("For at least one of the 'psi' priors that you have specified, you have not specified a 'buffer' and/or have specified an improper 'buffer' of the wrong type (i.e. not a numeric value or function) and/or containing at least 1 value that is negative; DICE cannot continue")
  } else {
   for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
    if(z%in%names(psiprior)){
     if(max(psiprior[[paste(z,sep='')]])>1&&typeof(buff[[paste(z,sep='')]])=='double'||max(psiprior[[paste(z,sep='')]])>1&&typeof(buff[[paste(z,sep='')]])=='integer'){
      if(length(sharedprior[[paste(z,sep='')]])==1){
       if(((max(psiprior[[paste(z,sep='')]])*buff[[paste(z,sep='')]])+1)>=length(sharedprior[[paste(z,sep='')]][[1]])){
        sign=1
       }
      } else {
       temp=sharedprior[[paste(z,sep='')]][[1]][0:buff[[paste(z,sep='')]]]
       for(y in 2:length(sharedprior[[paste(z,sep='')]])){
        if(length(sharedprior[[paste(z,sep='')]][[y]][!sharedprior[[paste(z,sep='')]][[y]]%in%temp])==0){
         sign=1
        } else {
         temp=c(temp,sharedprior[[paste(z,sep='')]][[y]][!sharedprior[[paste(z,sep='')]][[y]]%in%temp][0:buff[[paste(z,sep='')]]])
        }
       }
      }
      if(sign==1){
       print("For at least one of the 'psi' priors that you have specified, you have specified a 'shared' prior distribution that is too small given your upper 'psi' boundary and 'buffer' size; DICE cannot continue")
      }
     }
    }
   }
  }

  if(!-1%in%numtaxa&&!is.null(psiprior)&&!'psitest'%in%names(psiprior)&&!'zetatest'%in%names(zetaprior)&&!'zetatestpart'%in%names(zetaprior)&&zetaprior[['minzetatestpp']]==0&&zetaprior[['maxzetatestpp']]==0&&zetaprior[['minzetatesttot']]==0&&zetaprior[['maxzetatesttot']]==0&&zetaprior[['idiosynctest']]==0&&!'sharedtest'%in%names(sharedprior)&&!'bufftest'%in%names(buff)&&sign!=1){

   if(messages==0){

    if(numtaxa[['numtaxacaution']]==1){
     print("Caution: You have specified an inappropriate number of vector elements for 'num.taxa'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'num.taxa', with x = 'num.partitions', or if not applicable, with the first element of 'num.taxa' used for all partitions")
    }

    if(psiprior[['psicaution']]==1){
     print("Caution: For at least one of the 'psi' priors that you have specified, you have included more than the allowable number of distributions (1 for Ne and 2 for tao and epsilon); DICE will continue with only the first distributions up to the allowable amount")
    }

    if(zetaprior[['zetacaution']]==1){
     print("Caution: You have specified an inappropriate number of distributions for at least one of your 'zeta' priors; the number of distributions should be = 'num.partitions' or 1; DICE will continue with only the first x distributions of 'zeta', with x = 'num.partitions', or if not applicable, with the first distribution of 'zeta' used for all partitions")
    }

   }

   if(sharedprior[['sharedcaution']]==1){
    print("Caution: You have specified an inappropriate number of distributions for at least one of your 'shared' priors; the number of distributions should be = maximum value in the corresponding '.psi.prior' (may include more than the max 'psi' value if you also included 'idiosyncratic' distributions) or 1; DICE will continue with only the first 'shared' distribution for each pulse")
   }

   draws.psi=NULL
   draws.zeta=NULL
   draws.pulse.values=NULL

   for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
    #if this psi is intended to be inferred
    if(z%in%names(psiprior)){
     #protect against user specification of a psi prior of only 0(s)
     if(length(psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0])>0){
      for(y in 1:num.partitions){
       draws.zeta=append(draws.zeta,list(matrix(rep(as.integer(0),(num.sims*max(psiprior[[paste(z,sep='')]]))),ncol=max(psiprior[[paste(z,sep='')]]))))
       names(draws.zeta)[[length(draws.zeta)]]=paste(z,y,sep='.')
      }
      draws.pulse.values=append(draws.pulse.values,list(draws.zeta[[length(draws.zeta)]]))
      names(draws.pulse.values)[[length(draws.pulse.values)]]=z
      if(dirichlet.process==T){
       if(length(dirichlet.object[[paste('dpp',z,sep='.')]])==1){
        draws.coded=rep(dirichlet.object[[paste('dpp',z,sep='.')]],num.sims)
       } else {
        draws.coded=sample(dirichlet.object[[paste('dpp',z,sep='.')]],num.sims,replace=T)
       }
       draws.psi=append(draws.psi,list(matrix(rep(as.integer(0),num.sims),ncol=1)))
       for(y in unique(draws.coded)){
        if(sum(dirichlet.object[[paste('indices',z,sep='.')]][y,])!=0){
         for(x in 1:num.partitions){
          click=0
          counter=0
          w=0
          while(w<length(psiprior[[paste(z,sep='')]][psiprior[[paste(z,sep='')]]>0])&&click==0){
           w=w+1
           v=0
           while(v<num.partitions&&click==0){
            v=v+1
            u=0
            while(u<nrow(dirichlet.object[[paste('zeta',z,sep='.')]][[v]][[w]])&&click==0){
             u=u+1
             counter=counter+1
             if(counter==dirichlet.object[[paste('indices',z,sep='.')]][y,x]){
              click=1
              temp=dirichlet.object[[paste('zeta',z,sep='.')]][[v]][[w]][u,]
              if(x==1){
               draws.psi[[length(draws.psi)]][draws.coded==y,]=length(temp)
              }
              for(t in c(1:num.sims)[draws.coded==y]){
               draws.zeta[[((length(draws.zeta)-num.partitions)+x)]][t,1:length(temp)]=temp
              }
             }
            }
           }
          }
         }
        }
       }
      } else {
       if(length(psiprior[[paste(z,sep='')]])==1){
        draws.psi=append(draws.psi,list(matrix(rep(psiprior[[paste(z,sep='')]],num.sims),ncol=1)))
       } else {
        draws.psi=append(draws.psi, list(matrix(sample(psiprior[[paste(z,sep='')]], num.sims, replace=T),ncol=1)))
       }
       for(y in 1:num.sims){
        if(draws.psi[[length(draws.psi)]][y,]!=0){
         #check that zeta prior draws meet user specifications, assign valid combinations to object draws.zeta
         operator=0
         while(operator==0){
          twitch=0
          temp=NULL
          x=0
          while(x<num.partitions&&twitch==0){
           x=x+1
           if(length(zetaprior[[paste(z,sep='')]][[x]])==1){
            temp=append(temp,list(rep(zetaprior[[paste(z,sep='')]][[x]],draws.psi[[length(draws.psi)]][y,])))
           } else {
            temp=append(temp,list(sample(zetaprior[[paste(z,sep='')]][[x]],draws.psi[[length(draws.psi)]][y,],replace=T)))
           }
           #combination of zeta prior draws must be of appropriate number of taxa per partition, which is determined by the idiosyncratic logical value; by verifying combination is of appropriate number of taxa per partition, automatically combination is of appropriate number of total taxa across partitions
           if(idiosyncratic==T&&sum(temp[[length(temp)]])>numtaxa[['numtaxa']][x]||idiosyncratic==F&&sum(temp[[length(temp)]])!=numtaxa[['numtaxa']][x]){
            twitch=1
           }
          }
          if(twitch==0){
           #combination of zeta prior draws across partitions must be of appropriate min and max net zeta per pulse and in total
           if(!is.null(min.net.zeta.per.pulse)||!is.null(max.net.zeta.per.pulse)||!is.null(min.net.zeta.total)||!is.null(max.net.zeta.total)){
            #min and max net zeta per pulse
            temporary=0
            for(x in 1:draws.psi[[length(draws.psi)]][y,]){
             tentative=0
             for(w in 1:num.partitions){
              tentative=tentative+temp[[w]][x]
             }
             #convert user input min/max.net.zeta.per.pulse to a single numeric value
             for(w in c('min','max')){
              if(!is.null(eval(parse(text=paste(w,'.net.zeta.per.pulse',sep=''))))){
               interim=eval(parse(text=paste(w,'.net.zeta.per.pulse',sep='')))
               if(typeof(interim)!='list'){
                if(length(interim)==1){
                 interim=interim*sum(numtaxa[['numtaxa']])
                } else {
                 interim=interim[x]*sum(numtaxa[['numtaxa']])
                }
               } else {
                if(length(interim)==1){
                 if(length(interim[[1]])==1){
                  interim=interim[[1]]*sum(numtaxa[['numtaxa']])
                 } else {
                  interim=interim[[1]][x]*sum(numtaxa[['numtaxa']])
                 }
                } else {
                 interim=interim[[x]][1]*sum(numtaxa[['numtaxa']])
                }
               }
               if(w=='min'&&tentative<interim||w=='max'&&tentative>interim){
                twitch=1
               }
              }
             }
            temporary=temporary+tentative
            }
            #min and max net zeta in total across all pulses
            if(!is.null(min.net.zeta.total)&&temporary<(min.net.zeta.total*sum(numtaxa[['numtaxa']]))){
             twitch=1
            }
            if(!is.null(max.net.zeta.total)&&temporary>(max.net.zeta.total*sum(numtaxa[['numtaxa']]))){
             twitch=1
            }
           }
          }
          if(twitch==0){
           operator=1
          }
         }
         #assign valid combinations to object draws.zeta
         for(x in 1:num.partitions){
          draws.zeta[[((length(draws.zeta)-num.partitions)+x)]][y,1:draws.psi[[length(draws.psi)]][y,]]=temp[[x]]
         }
        }
       }
      }
      names(draws.psi)[[length(draws.psi)]]=z
      if(length(sharedprior[[paste(z,sep='')]])==1){
       for(y in 1:num.sims){
        if(draws.psi[[length(draws.psi)]][y,]!=0){
         temp1=NULL
         temp2=NULL
         for(x in 1:draws.psi[[length(draws.psi)]][y,]){
          if(length(sharedprior[[paste(z,sep='')]][[1]][!sharedprior[[paste(z,sep='')]][[1]]%in%temp2])==1){
           temp1=c(temp1,sharedprior[[paste(z,sep='')]][[1]][!sharedprior[[paste(z,sep='')]][[1]]%in%temp2])
          } else {
           temp1=c(temp1,sample(sharedprior[[paste(z,sep='')]][[1]][!sharedprior[[paste(z,sep='')]][[1]]%in%temp2],1,replace=T))
          }
          if(typeof(buff[[paste(z,sep='')]])=='closure'){
           buff.temp=buff[[paste(z,sep='')]](temp1[x])
           temp2=unique(c(temp2,buff.temp))
          } else {
           temp2=unique(c(temp2,(temp1[x]-buff[[paste(z,sep='')]]):(temp1[x]+buff[[paste(z,sep='')]])))
          }
         }
         draws.pulse.values[[length(draws.pulse.values)]][y,1:draws.psi[[length(draws.psi)]][y,]]=sort(temp1)
        }
       }
      } else {
       #not buffering out former draws iteratively as done before since the draws are being done explicitly in the order of pulses, and blocking out buffered zones in this order may place a drawing bias; thus, taking the approach of taking random draws until the draws satisfy the requirements, which is less efficient but less biased
       for(y in 1:num.sims){
        if(draws.psi[[length(draws.psi)]][y,]!=0){
         option=0
         while(option==0){
          temp=NULL
          for(x in 1:draws.psi[[length(draws.psi)]][y,]){
           if(length(sharedprior[[paste(z,sep='')]][[x]])==1){
            temp=c(temp,sharedprior[[paste(z,sep='')]][[x]])
           } else {
            temp=c(temp,sample(sharedprior[[paste(z,sep='')]][[x]],1,replace=T))
           }
          }
          #confirming draws are in chronological/ascending order
          if(sum(sort(temp)==temp)!=length(temp)){
           option=1
          }
          #confirming draws are outside of buffer zone
          if(length(temp)>1&&option==0){
           if(typeof(buff[[paste(z,sep='')]])=='closure'){
            for(x in 1:(length(temp))){
             if(x>1){
              if(temp[x]%in%buff[[paste(z,sep='')]](temp1[(x-1)])){
               option=1
              }
             }
             if(x<length(temp)){
              if(temp[x]%in%buff[[paste(z,sep='')]](temp1[(x+1)])){
               option=1
              }
             }
            }
           } else {
            for(x in 1:(length(temp)-1)){
             if((temp[(x+1)]-temp[x])<buff[[paste(z,sep='')]]){
              option=1
             }
            }
           }
          }
          if(option==0){
           option=1
          } else {
           option=0
          }
         }
         draws.pulse.values[[length(draws.pulse.values)]][y,1:draws.psi[[length(draws.psi)]][y,]]=temp
        }
       }
      }
     } else {
      for(y in 1:num.partitions){
       draws.zeta=append(draws.zeta,list(matrix(rep(as.integer(0),num.sims),ncol=1)))
       names(draws.zeta)[[length(draws.zeta)]]=paste(z,y,sep='.')
      }
      draws.pulse.values=append(draws.pulse.values,list(draws.zeta[[length(draws.zeta)]]))
      names(draws.pulse.values)[[length(draws.pulse.values)]]=z
      draws.psi=append(draws.psi,list(matrix(rep(as.integer(0),num.sims),ncol=1)))
      names(draws.psi)[[length(draws.psi)]]=z
     }
    }
   }

   #these list elements consist of the hyperparameters to be inferred
   output=list(draws.psi,draws.zeta,draws.pulse.values)
   names(output)=c('draws.psi','draws.zeta','draws.pulse.values')
   return(output)

  }

 }

}



play.dice=function(num.sims=1000000, num.partitions=1, num.taxa, idiosyncratic=T, idiosyncratic.buffer=T, dirichlet.process=F, tao.psi.prior=NULL, epsilon.psi.prior=NULL, NE.psi.prior=NULL, tao.zeta.prior=NULL, tao2.zeta.prior=NULL, epsilon.zeta.prior=NULL, epsilon2.zeta.prior=NULL, NE.zeta.prior=NULL, min.net.zeta.per.pulse=NULL, max.net.zeta.per.pulse=NULL, min.net.zeta.total=NULL, max.net.zeta.total=NULL, tao.shared.prior=NULL, tao2.shared.prior=NULL, epsilon.shared.prior=NULL, epsilon2.shared.prior=NULL, NE.shared.prior=NULL, tao.buffer=0, tao2.buffer=NULL, epsilon.buffer=0, epsilon2.buffer=NULL, NE.buffer=0, net.zeta.per.pulse=F, net.zeta.total=F, mean.tao.shared=F, disp.index.tao.shared=F, mean.tao=F, disp.index.tao=F, mean.tao2.shared=F, disp.index.tao2.shared=F, mean.tao2=F, disp.index.tao2=F, mean.epsilon.shared=F, disp.index.epsilon.shared=F, mean.epsilon=F, disp.index.epsilon=F, mean.epsilon2.shared=F, disp.index.epsilon2.shared=F, mean.epsilon2=F, disp.index.epsilon2=F, mean.NE.shared=F, disp.index.NE.shared=F, mean.NE=F, disp.index.NE=F, tao.idio.prior=NULL, tao2.idio.prior=NULL, epsilon.idio.prior=NULL, epsilon2.idio.prior=NULL, NE.idio.prior=NULL, linked.param=NULL, attached.hyper=NULL, fix.linked.param.across.pulses=F, fix.linked.param.per.pulse=F, num.changes=1, flip=F, anchor.prior=NULL, change.prior=NULL, exponential.growth.rate.prior=NULL, exponential.growth.rate.prior2=NULL, dirichlet.object=NULL, roll.object=NULL, fsc2template.path, fsc2path, input.directory, input.base, input.files, output.directory, folded, remove.afclasses, append.sims, keep.taxa.draws, output.hyper.draws, output.taxa.draws, keep.fsc2.files, messages.sims, num.haploid.samples, sampling.times, num.ind.sites, num.SNPs, length.seq, mut.rate, gen.times, play.object, skipper=NULL){

 messages=0
 if(is.null(roll.object)){
  roll.object=roll.dice(num.sims=num.sims, num.partitions=num.partitions, num.taxa=num.taxa, idiosyncratic=idiosyncratic, dirichlet.process=dirichlet.process, tao.psi.prior=tao.psi.prior, epsilon.psi.prior=epsilon.psi.prior, NE.psi.prior=NE.psi.prior, tao.zeta.prior=tao.zeta.prior, tao2.zeta.prior=tao2.zeta.prior, epsilon.zeta.prior=epsilon.zeta.prior, epsilon2.zeta.prior=epsilon2.zeta.prior, NE.zeta.prior=NE.zeta.prior, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total, tao.shared.prior=tao.shared.prior, tao2.shared.prior=tao2.shared.prior, epsilon.shared.prior=epsilon.shared.prior, epsilon2.shared.prior=epsilon2.shared.prior, NE.shared.prior=NE.shared.prior, tao.buffer=tao.buffer, tao2.buffer=tao2.buffer, epsilon.buffer=epsilon.buffer, epsilon2.buffer=epsilon2.buffer, NE.buffer=NE.buffer, dirichlet.object=dirichlet.object)
  messages=1
 }

 if(!is.null(roll.object)){

  numtaxa=dice.read.numtaxa(num.partitions=num.partitions, num.taxa=num.taxa)

  if(-1%in%numtaxa){
   print("You have included at least 1 non-positive number in the 'num.taxa' object; DICE cannot continue")
  }

  psiprior=dice.read.psiprior(tao.psi.prior=tao.psi.prior, epsilon.psi.prior=epsilon.psi.prior, NE.psi.prior=NE.psi.prior)

  if(is.null(psiprior)){
   print("You have not specified any 'psi' priors; DICE cannot continue")
  }

  if('psitest'%in%names(psiprior)){
   print("At least one of the 'psi' priors that you have specified contains at least 1 value that is negative; DICE cannot continue")
  }

  #nuisance draws are never buffered, and thus although an idiosyncratic prior is needed for nuisance draws, a buffer is not needed 
  if(idiosyncratic==T){
   buff=dice.read.buff(psiprior=psiprior, tao.buffer=tao.buffer, tao2.buffer=tao2.buffer, epsilon.buffer=epsilon.buffer, epsilon2.buffer=epsilon2.buffer, NE.buffer=NE.buffer)
  }

  numchanges=dice.read.numchanges(num.partitions=num.partitions, num.changes=num.changes)

  if(-1%in%numchanges){
   print("You have included at least 1 non-positive number in the 'num.changes' object; DICE cannot continue")
  }

  halt=0

  #if idiosyncratic=T, then an idiosyncratic prior is needed; if any of the first event hyperparameters are not being inferred, or if any of the second event hyperparameters are not being inferred and a 2 event model is used for at least one partition, then an idiosyncratic prior is needed for these nuisance parameter draws (for tao2, anchor.prior and change.prior must also both be NOT specified by user)
  if(idiosyncratic==T||!'tao'%in%names(roll.object$draws.psi)||!'epsilon'%in%names(roll.object$draws.psi)||!'NE'%in%names(roll.object$draws.psi)||2%in%numchanges[['numchanges']]&&!'tao2'%in%names(roll.object$draws.psi)&&is.null(anchor.prior)&&is.null(change.prior)||2%in%numchanges[['numchanges']]&&!'epsilon2'%in%names(roll.object$draws.psi)){

   if(!is.null(linked.param)||!is.null(attached.hyper)){
    linking=dice.read.linking(linked.param=linked.param, attached.hyper=attached.hyper)
   } else {
    linking=NULL
   }

   if(-1%in%linking){
    halt=1
    print("You have misspecified the 'linked.param' and/or 'attached.hyper' objects; 'linked.param' cannot have more than one of any parameter and each of these objects may contain only the allowable parameter names (i.e. 'tao', 'tao2', 'epsilon', 'epsilon2', 'NE') and must be of the same length; DICE cannot continue")
   }

   if(halt==0){

    idioprior=dice.read.idioprior(num.partitions=num.partitions, idiosyncratic=idiosyncratic, psiprior=psiprior, tao.shared.prior=tao.shared.prior, tao2.shared.prior=tao2.shared.prior, epsilon.shared.prior=epsilon.shared.prior, epsilon2.shared.prior=epsilon2.shared.prior, NE.shared.prior=NE.shared.prior, tao.idio.prior=tao.idio.prior, tao2.idio.prior=tao2.idio.prior, epsilon.idio.prior=epsilon.idio.prior, epsilon2.idio.prior=epsilon2.idio.prior, NE.idio.prior=NE.idio.prior, linking=linking, numchanges=numchanges, anchor.prior=anchor.prior, change.prior=change.prior)

    if('idiotest'%in%names(idioprior)){
     print("For at least one of the taxon-specific parameters, you have not specified an 'idiosyncratic' prior and/or have specified an improper 'idiosyncratic' prior containing at least 1 value that is non-positive; DICE cannot continue")
     halt=1
    }

   }

   if(halt==0){

    fixing=dice.read.fixing(linking=linking, fix.linked.param.across.pulses=fix.linked.param.across.pulses, fix.linked.param.per.pulse=fix.linked.param.per.pulse)

    if(idiosyncratic==T){
     if('bufftest'%in%names(buff)){
      halt=1
      print("For at least one of the 'psi' priors that you have specified, you have not specified a 'buffer' and/or have specified an improper 'buffer' of the wrong type (i.e. not a numeric value or function) and/or containing at least 1 value that is negative; DICE cannot continue")
     } else {
      if(idiosyncratic.buffer==T){
       zetaprior=dice.read.zetaprior(num.partitions=num.partitions, numtaxa=numtaxa, idiosyncratic=idiosyncratic, psiprior=psiprior, tao.zeta.prior=tao.zeta.prior, tao2.zeta.prior=tao2.zeta.prior, epsilon.zeta.prior=epsilon.zeta.prior, epsilon2.zeta.prior=epsilon2.zeta.prior, NE.zeta.prior=NE.zeta.prior, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total)
       if('zetatest'%in%names(zetaprior)){
        halt=1
        print("For at least one of the 'psi' priors that you have specified, you have not specified a 'zeta' prior and/or have specified an improper 'zeta' prior containing at least 1 value that is above 1.0 and/or negative; DICE cannot continue")
       }
       if('zetatestpart'%in%names(zetaprior)){
        halt=1
        print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the number of taxa that you have specified per partition; DICE cannot continue")
       }
       if(zetaprior[['minzetatestpp']]==1){
        halt=1
        print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the minimum zeta value per pulse that you have specified; DICE cannot continue")
       }
       if(zetaprior[['maxzetatestpp']]==1){
        halt=1
        print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the maximum zeta value per pulse that you have specified; DICE cannot continue")
       }
       if(zetaprior[['minzetatesttot']]==1){
        halt=1
        print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the minimum zeta value across all pulses (i.e. in total among the dataset) that you have specified; if you have specified 'idiosyncratic==F' (i.e. all taxa must be in a pulse), then the minimum zeta value across all pulses must either be NULL or equal to 1; DICE cannot continue")
       }
       if(zetaprior[['maxzetatesttot']]==1){
        halt=1
        print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on the maximum zeta value across all pulses (i.e. in total among the dataset) that you have specified; if you have specified 'idiosyncratic==F' (i.e. all taxa must be in a pulse), then the maximum zeta value across all pulses must either be NULL or equal to 1; DICE cannot continue")
       }
       if(zetaprior[['idiosynctest']]==1){
        halt=1
        print("At least one of the 'zeta' priors that you have specified does not allow for any valid draws based on you specifying 'idiosyncratic==F' (i.e. all taxa must be in a pulse), the respective 'psi' prior that you have specified, and the number of taxa that you have specified; DICE cannot continue")
       }
      }
      for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
       if(z%in%names(psiprior)){
        if(typeof(buff[[paste(z,sep='')]])=='double'||typeof(buff[[paste(z,sep='')]])=='integer'){
         #this is a bit conservative as it assumes all the sharedprior draws + their buffer would overlap into the idio distribution, but it also doesn't consider, in the case of inferring psi2/tao2, the chunk of distribution taken out due to the first event tao draw (though this isn't a huge concern, since it is coded that if there is no idio distribution left in psi2/tao2, the lowest available tao2 will be used)
         if(idiosyncratic.buffer==T){
          temp=NULL
          for(y in 1:num.partitions){
           if(z%in%names(zetaprior)){
            minzeta=min(zetaprior[[paste(z,sep='')]][[y]])
           } else {
            minzeta=0
           }
           if(minzeta==0){
            zetacorrect=0
           } else {
            zetacorrect=1
           }
           #this was intended to consider the maximum number of idiosyncratic taxa there would be; however, this statement originally was problematic with zeta distributions that include 0 and min(psiprior)>0, since this statement assumed each pulse would have at least 1 taxon, thus the zetacorrect now
           if(((((min(psiprior[[paste(z,sep='')]])*zetacorrect)+(numtaxa[['numtaxa']][y]-(min(psiprior[[paste(z,sep='')]])*minzeta)))*buff[[paste(z,sep='')]])+1)>=length(idioprior[[paste(z,sep='')]][[y]][!idioprior[[paste(z,sep='')]][[y]]%in%temp])){
            halt=1
            print("For at least one of the 'psi' priors that you have specified, you have specified an 'idiosyncratic' prior distribution that is too small given your number of taxa, 'buffer' size, and 'psi'/'zeta'/'idiosyncratic.buffer' specifications; DICE cannot continue")
           } else {
            temp=c(temp,idioprior[[paste(z,sep='')]][[y]][!idioprior[[paste(z,sep='')]][[y]]%in%temp][0:((min(psiprior[[paste(z,sep='')]])+(numtaxa[['numtaxa']][y]-(min(psiprior[[paste(z,sep='')]])*minzeta)))*buff[[paste(z,sep='')]])])
           }
          }
         } else {
          for(y in 1:num.partitions){
           if(((max(psiprior[[paste(z,sep='')]])*buff[[paste(z,sep='')]])+1)>=length(idioprior[[paste(z,sep='')]][[y]])){
            halt=1
            print("For at least one of the 'psi' priors that you have specified, you have specified an 'idiosyncratic' prior distribution that is too small given your number of taxa, 'buffer' size, and 'psi'/'zeta'/'idiosyncratic.buffer' specifications; DICE cannot continue")
           }
          }
         }
        }
       }
      }
     }
    }

   }

  }

  if(2%in%numchanges[['numchanges']]&&!'tao2'%in%names(psiprior)&&is.null(anchor.prior)){
   flipvector=dice.read.flipvector(num.partitions, flip)
  }

  if(!'tao2'%in%names(psiprior)&&2%in%numchanges[['numchanges']]&&!is.null(anchor.prior)){

   anchorprior=dice.read.anchorprior(psipriortao=psiprior$tao, anchor.prior=anchor.prior)

   if('anchortest'%in%names(anchorprior)){
    print("You have specified an improper 'anchor' prior containing at least 1 value that is non-positive; DICE cannot continue")
    halt=1
   }

  }

  if(!'tao2'%in%names(psiprior)&&2%in%numchanges[['numchanges']]&&is.null(anchor.prior)&&!is.null(change.prior)||'tao2'%in%names(roll.object$draws.psi)&&idiosyncratic==T&&!is.null(change.prior)){

   changeprior=dice.read.changeprior(num.partitions=num.partitions, change.prior=change.prior)

   if('changetest'%in%names(changeprior)){
    print("You have specified an improper 'change' prior containing at least 1 value that is non-positive; DICE cannot continue")
    halt=1
   }

  }

  #it would be difficult to buffer around the beginning and end exponential growth times because this rate draw would need to be done in the roll.dice function when drawing hyperparameters since buffering occurs during those hyperparameter draws; instead, buffering occurs around the beginning time
  if(!is.null(exponential.growth.rate.prior)||!is.null(exponential.growth.rate.prior2)){

   exp.grow=dice.read.exp.grow(num.partitions=num.partitions, exponential.growth.rate.prior=exponential.growth.rate.prior, exponential.growth.rate.prior2=exponential.growth.rate.prior2)

  }

  if(!-1%in%numtaxa&&!is.null(psiprior)&&!'psitest'%in%names(psiprior)&&!-1%in%numchanges&&halt!=1){

   if(is.null(skipper)){

    if(messages==0){

     if(numtaxa[['numtaxacaution']]==1){
      print("Caution: You have specified an inappropriate number of vector elements for 'num.taxa'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'num.taxa', with x = 'num.partitions', or if not applicable, with the first element of 'num.taxa' used for all partitions")
     }

     if(psiprior[['psicaution']]==1){
      print("Caution: For at least one of the 'psi' priors that you have specified, you have included more than the allowable number of distributions (1 for Ne and 2 for tao and epsilon); DICE will continue with only the first distributions up to the allowable amount")
     }

    }

    if(numchanges[['numchangescaution']]==1){
     print("Caution: You have specified an inappropriate number of vector elements for 'num.changes', and/or specified an inappropriate number within 'num.changes'; this should be of length = 'num.partitions' or 1, and include values of only 1 or 2 (i.e. DICE currently allows only up to two demographic events per taxon); DICE will continue with only the first x elements of 'num.changes', with x = 'num.partitions', or if not applicable, with the first element of 'num.changes' used for all partitions, as well as any number >2 converted to 2")
    }

    if(idiosyncratic==T||!'tao'%in%names(roll.object$draws.psi)||!'epsilon'%in%names(roll.object$draws.psi)||!'NE'%in%names(roll.object$draws.psi)||2%in%numchanges[['numchanges']]&&!'tao2'%in%names(roll.object$draws.psi)&&is.null(anchor.prior)&&is.null(change.prior)||2%in%numchanges[['numchanges']]&&!'epsilon2'%in%names(roll.object$draws.psi)){

     if(idioprior[['idiocaution']]==1){
      print("Caution: You have specified an inappropriate number of distributions for at least one of your 'idiosyncratic' priors, contained in the corresponding '.idio.prior' and/or '.shared.prior' objects; the number of distributions should be = 'num.partitions' or 1, or in the case of a parameter linked to another hyperparameter, the number of distributions should be = maximum value in the '.psi.prior' of the attached hyperparameter or 1; DICE will continue with only the first x 'idiosyncratic' distributions, with x = 'num.partitions' or maximum value in the '.psi.prior' of the attached hyperparameter, or if not applicable, with the first 'idiosyncratic' distribution used for all partitions or pulses in the attached hyperparameter")
     }

     if(fixing[['fixingcaution']]==1){
      print("Caution: You have specified an inappropriate number of vector elements for 'fix.linked.param.across.pulses' and/or 'fix.linked.param.per.pulse'; these should be of the same length as 'linked.param' and 'attached.hyper' or of length = 1; DICE will continue with only the first x elements, with x = the length of 'linked.param' and 'attached.hyper', or if not applicable, with the first element used for all linked parameters")
     }

     if(idiosyncratic==T&&idiosyncratic.buffer==T){
      if(zetaprior[['zetacaution']]==1){
       print("Caution: You have specified an inappropriate number of distributions for at least one of your 'zeta' priors; the number of distributions should be = 'num.partitions' or 1; DICE will continue with only the first x distributions of 'zeta', with x = 'num.partitions', or if not applicable, with the first distribution of 'zeta' used for all partitions")
      }
     }

    }

    if(2%in%numchanges[['numchanges']]&&!'tao2'%in%names(psiprior)&&is.null(anchor.prior)){
     if(flipvector[['flipvectorcaution']]==1){
      print("Caution: You have specified an inappropriate number of vector elements for 'flip'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'flip', with x = 'num.partitions', or if not applicable, with the first element of 'flip' used for all partitions")
     }
    }

    if(!'tao2'%in%names(psiprior)&&2%in%numchanges[['numchanges']]&&!is.null(anchor.prior)){
     if(anchorprior[['anchorcaution']]==1){
      print("Caution: You have specified an inappropriate number of distributions for your 'anchor' prior; the number of distributions should be = maximum value in the corresponding '.psi.prior' (may include an additional distribution for idiosyncratic taxa) or 1; DICE will continue with only the first x 'anchor' distributions, with x = maximum value in the corresponding '.psi.prior' + 1, or if not applicable, with the first 'anchor' distribution used for all pulses and idiosyncratic taxa")
     }
    }

    if(!'tao2'%in%names(psiprior)&&2%in%numchanges[['numchanges']]&&is.null(anchor.prior)&&!is.null(change.prior)||'tao2'%in%names(roll.object$draws.psi)&&idiosyncratic==T&&!is.null(change.prior)){
     if(changeprior[['changecaution']]==1){
      print("Caution: You have specified an inappropriate number of distributions for your 'change' prior; the number of distributions should be = 'num.partitions' or 1; DICE will continue with only the first x 'change' distributions, with x = 'num.partitions', or if not applicable, with the first 'change' distribution used for all partitions")
     }
    }

    if(!is.null(exponential.growth.rate.prior)||!is.null(exponential.growth.rate.prior2)){
     if(exp.grow[['exp.growcaution']]==1){
      print("Caution: For at least one of the 'exponential.growth.rate' priors that you have specified, you have included an inappropriate number of distributions; the number of distributions should be = 'num.partitions' or 1; DICE will continue with only the first x distributions, with x = 'num.partitions', or if not applicable, with the first distribution used for all partitions")
     }
    }

   }

   sim.specs=NULL

   if(!is.null(exponential.growth.rate.prior)||!is.null(exponential.growth.rate.prior2)){
    for(z in c('exponential.growth.rate.prior','exponential.growth.rate.prior2')){
     if(z%in%names(exp.grow)){
      sim.specs=append(sim.specs,list(matrix(rep(0,(sum(numtaxa[['numtaxa']])*num.sims)),ncol=sum(numtaxa[['numtaxa']]))))
      names(sim.specs)[[length(sim.specs)]]=z
      for(y in 1:num.sims){
       for(x in 1:num.partitions){
        if(z=='exponential.growth.rate.prior'||numchanges[['numchanges']][x]==2){
         if(length(exp.grow[[paste(z,sep='')]][[x]])==1){
          sim.specs[[length(sim.specs)]][y,(sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x])]=rep(exp.grow[[paste(z,sep='')]][[x]],numtaxa[['numtaxa']][x])
         } else {
          sim.specs[[length(sim.specs)]][y,(sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x])]=sample(exp.grow[[paste(z,sep='')]][[x]],numtaxa[['numtaxa']][x],replace=T)
         }
        }
       }
      }
     }
    }
   }

   if(dirichlet.process==T){
    starter=-1
   } else {
    starter=NULL
   }

   #assign shared pulse values to individual taxa; can't assign idiosyncratic draws here since, for a 2-event model with psi1/zeta1 and psi2/zeta2 being inferred, the idiosyncratic taxa in the first event need to be unassigned when assigning shared pulse values to individual taxa for the second event (so as not to interfere timing-wise)
   for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
    #if this parameter will be simulated
    if(substr(z,nchar(z),nchar(z))!='2'||substr(z,nchar(z),nchar(z))=='2'&&2%in%numchanges[['numchanges']]){
     sim.specs=append(sim.specs,list(matrix(rep(as.integer(0),(sum(numtaxa[['numtaxa']])*num.sims)),ncol=sum(numtaxa[['numtaxa']]))))
     names(sim.specs)[[length(sim.specs)]]=z
    }
    #if this psi is intended to be inferred
    if(z%in%names(psiprior)){
     #if mean of shared pulse values of this parameter is to be inferred
     if(eval(parse(text=paste('mean.',z,'.shared',sep='')))==T){
      roll.object=append(roll.object,list(matrix(rep(0,num.sims),ncol=1)))
      names(roll.object)[[length(roll.object)]]=paste('mean.',z,'.shared',sep='')
     }
     #if dispersion index of shared pulse values of this parameter is to be inferred
     if(eval(parse(text=paste('disp.index.',z,'.shared',sep='')))==T){
      roll.object=append(roll.object,list(matrix(rep(0,num.sims),ncol=1)))
      names(roll.object)[[length(roll.object)]]=paste('disp.index.',z,'.shared',sep='')
     }
     #if net.zeta.per.pulse is to be inferred
     if(net.zeta.per.pulse==T){
      roll.object=append(roll.object,list(matrix(rep(as.integer(0),(num.sims*ncol(roll.object[['draws.zeta']][[paste(z,'.1',sep='')]]))),ncol=ncol(roll.object[['draws.zeta']][[paste(z,'.1',sep='')]]))))
      names(roll.object)[[length(roll.object)]]=paste('net.zeta.per.pulse.',z,sep='')
      for(y in 1:num.sims){
       for(x in 1:ncol(roll.object[['draws.zeta']][[paste(z,'.1',sep='')]])){
        temp=0
        for(w in 1:num.partitions){
         temp=temp+roll.object[['draws.zeta']][[paste(z,w,sep='.')]][y,x]
        }
        roll.object[[length(roll.object)]][y,x]=temp
       }
      }
     }
     #if net.zeta.total is to be inferred
     if(net.zeta.total==T){
      roll.object=append(roll.object,list(matrix(rep(as.integer(0),num.sims,ncol=1))))
      names(roll.object)[[length(roll.object)]]=paste('net.zeta.total.',z,sep='')
      for(y in 1:num.sims){
       temp=0
       for(x in 1:ncol(roll.object[['draws.zeta']][[paste(z,'.1',sep='')]])){
        for(w in 1:num.partitions){
         temp=temp+roll.object[['draws.zeta']][[paste(z,w,sep='.')]][y,x]
        }
       }
       roll.object[[length(roll.object)]][y,]=temp
      }
     }
     for(y in 1:num.sims){
      if(roll.object[['draws.psi']][[paste(z,sep='')]][y,]>0){
       #randomly select taxa assignment
       if(z!='tao2'){
        for(x in 1:num.partitions){
         for(w in 1:roll.object[['draws.psi']][[paste(z,sep='')]][y,]){
          if(roll.object[['draws.zeta']][[paste(z,x,sep='.')]][y,w]>0){
           if(sum(sim.specs[[length(sim.specs)]][y,(sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x])]==0)==1){
            sim.specs[[length(sim.specs)]][y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))[sim.specs[[length(sim.specs)]][y,(sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x])]==0]]=roll.object[['draws.pulse.values']][[paste(z,sep='')]][y,w]
           } else {
            sim.specs[[length(sim.specs)]][y,combn(c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))[sim.specs[[length(sim.specs)]][y,(sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x])]==0],roll.object[['draws.zeta']][[paste(z,x,sep='.')]][y,w])[,sample(1:ncol(combn(c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))[sim.specs[[length(sim.specs)]][y,(sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x])]==0],roll.object[['draws.zeta']][[paste(z,x,sep='.')]][y,w])),1,replace=T)]]=roll.object[['draws.pulse.values']][[paste(z,sep='')]][y,w]
           }
          }
         }
        }
       } else {
        detonator=0
        button=0
        #as soon as a partition shows to have no valid tao2 assignments that wouldn't interfere with tao1 draws, draw is invalidated via detonator and loop repeats with a new draw via roll.dice function (hence the button activated after the first pass, thus allowing new draws in every loop after the first that failed); since this is built to re-draw tao2 if incompatible with tao1, there is, so to speak, an unaffected tao1, while tao2 distributions/taxa assignments are influenced by tao1 draws
        while(detonator==0){
         #re-draw if last draw did not allow valid tao1 draws
         if(button==1){
          if(-1%in%starter){
           starter=list(dirichlet.object$dpp.tao.prior2,as.matrix(dirichlet.object$zeta.tao.prior2),as.matrix(dirichlet.object$indices.tao.prior2))
           names(starter)=c('dpp.tao.prior','zeta.tao.prior','indices.tao.prior')
          }
          #treat user specified tao2-related distributions as tao1 distributions for this particular use of the roll.dice function since the function is built so that if tao2-related hyperparameters are to be inferred, then so must too tao1-related hyperparameters, which is not of interest in this case, so just treat tao2 as tao1 since this doesn't affect the drawing procedure; NULL everything else
          new.draw=roll.dice(num.sims=1, num.partitions=num.partitions, num.taxa=num.taxa, idiosyncratic=idiosyncratic, dirichlet.process=dirichlet.process, tao.psi.prior=tao.psi.prior[[2]], epsilon.psi.prior=NULL, NE.psi.prior=NULL, tao.zeta.prior=tao2.zeta.prior, tao2.zeta.prior=NULL, epsilon.zeta.prior=NULL, epsilon2.zeta.prior=NULL, NE.zeta.prior=NULL, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total, tao.shared.prior=tao2.shared.prior, tao2.shared.prior=NULL, epsilon.shared.prior=NULL, epsilon2.shared.prior=NULL, NE.shared.prior=NULL, tao.buffer=buff$tao2, tao2.buffer=NULL, epsilon.buffer=NULL, epsilon2.buffer=NULL, NE.buffer=NULL, dirichlet.object=starter)
          roll.object$draws.psi$tao2[y,]=new.draw$draws.psi$tao
          for(x in 1:num.partitions){
           roll.object[['draws.zeta']][[paste('tao2.',x,sep='')]][y,]=new.draw[['draws.zeta']][[paste('tao.',x,sep='')]]
          }
          roll.object$draws.pulse.values$tao2[y,]=new.draw$draws.pulse.values$tao
         }
         button=1
         #must replicate if statement because of additional loops with redraws that might draw psi=0
         if(roll.object[['draws.psi']][[paste(z,sep='')]][y,]>0){
          #if the first event will need to do a future idiosyncratic draw
          if(idiosyncratic==T){
           buffer=NULL
           if(roll.object$draws.psi$tao[y,]>0){
            for(x in 1:roll.object$draws.psi$tao[y,]){
             if(typeof(buff$tao)=='closure'){
              buff.temp=buff$tao(roll.object$draws.pulse.values$tao[y,x])
              buffer=unique(c(buffer,buff.temp))
             } else {
              buffer=unique(c(buffer,(roll.object$draws.pulse.values$tao[y,x]-buff$tao):(roll.object$draws.pulse.values$tao[y,x]+buff$tao)))
             }
            }
           }
          }
          tao1.draws=sim.specs[['tao']][y,]
          tao2.psi=roll.object[['draws.psi']][['tao2']][y,]
          tao2.pv=roll.object[['draws.pulse.values']][['tao2']][y,]
          indi.combos=NULL
          x=0
          while(x<num.partitions&&detonator==0){
           x=x+1
           ind.count.pulse=1:numtaxa[['numtaxa']][x]
           ind.count.total=ind.count.pulse+sum(numtaxa[['numtaxa']][0:(x-1)])
           tao2.zeta=roll.object[['draws.zeta']][[paste(z,x,sep='.')]][y,]
           tao2.zeta.real.id=c(1:tao2.psi)[tao2.zeta>0]
           indi.combos=append(indi.combos,list(NULL))
           if(length(tao2.zeta.real.id)>0){
            #assign the minimum tao2 synchronous draw allowed per taxon based on tao1 assignment (i.e. tao1 must be less than tao2; if tao1 is 0 i.e. idiosyncratic, tao2 will certainly be greater than tao1 i.e. 0 and tao2 will serve as the upper bound for downstream tao1 idiosyncratic draw, thus tao1 idiosyncratic prior (which is modified by the buffer(s) around the tao1 shared pulse draws) must be sufficient to accommodate tao2)
            ind.allowed=rep((tao2.psi+1),numtaxa[['numtaxa']][x])
            for(w in tao2.psi:1){
             ind.allowed[ind.count.pulse[ind.allowed==(w+1)][tao1.draws[ind.count.total[ind.allowed==(w+1)]]<tao2.pv[w]]]=w
            }
            if(sum(tao1.draws[ind.count.total]==0)>0){
             tick=0
             w=tao2.psi
             while(w>0&&tick==0){
              if(length(idioprior$tao[[x]][!idioprior$tao[[x]]%in%buffer][idioprior$tao[[x]][!idioprior$tao[[x]]%in%buffer]<tao2.pv[w]])==0){
               ind.allowed[ind.count.pulse[tao1.draws[ind.count.total]==0][ind.allowed[ind.count.pulse[tao1.draws[ind.count.total]==0]]<=w]]=w+1
               tick=1
              }
              w=w-1
             }
            }
            w=0
            while(w<tao2.psi&&detonator==0){
             w=w+1
             if(sum(ind.allowed<=w)<sum(tao2.zeta[(1:w)])){
              detonator=1
             }
            }
            if(detonator==0){
             combo.inds=NULL
             selected.combos=NULL
             tot.zeta=NULL
             for(w in tao2.zeta.real.id){
              allowed.ind=c(1:length(ind.allowed))[ind.allowed<=w]
              if(length(allowed.ind[!(allowed.ind+sum(numtaxa[['numtaxa']][0:(x-1)]))%in%tot.zeta])==1){
               combo.inds=append(combo.inds,list(matrix((allowed.ind[!(allowed.ind+sum(numtaxa[['numtaxa']][0:(x-1)]))%in%tot.zeta]+sum(numtaxa[['numtaxa']][0:(x-1)])),ncol=1)))
              } else {
               combo.inds=append(combo.inds,list(combn((allowed.ind[!(allowed.ind+sum(numtaxa[['numtaxa']][0:(x-1)]))%in%tot.zeta]+sum(numtaxa[['numtaxa']][0:(x-1)])),tao2.zeta[w])))
              }
              selected.combos=c(selected.combos,sample(1:ncol(combo.inds[[length(combo.inds)]]),1,replace=T))
              tot.zeta=c(tot.zeta,combo.inds[[length(combo.inds)]][,selected[length(selected)]])
             }
             if(idiosyncratic==T&&idiosyncratic.buffer==T){
              circle=0
              while(circle==0){
               tot.zeta=NULL
               lever=0
               #at one point, this was coded with an object of unavailable tao1 idio distribution, which was added to at each loop from the lower bound of the leftover tao1 idio distribution of length = buff$tao * # of tao1 idiosyncratic taxa with this tao2 pulse; this determined if there is still available tao1 idioprior space given the added idiosyncratic buffering of previously looped taxa; this same function is conserved here, since this goes in order of pulses (and thus smallest to largest available distribution) and the available amount of idio distribution in earlier pulses is completely subsumed by latter pulses; if the while loop continues without breaking, then earlier pulses guaranteed to have a large enough distribution, and now for the current pulse in the loop, if there is enough extra distribution exclusive from the earlier pulses for the idiosyncratic taxa with this tao2 pulse value
               w=0
               while(w<length(tao2.zeta.real.id)&&lever==0){
                w=w+1
                tot.zeta=c(tot.zeta,combo.inds[[w]][,selected.combos[w]])
                #if the idiosyncratic prior of tao1 (which is modified by the buffer(s) around the shared pulse draws, and needs to be of a certain size to accommodate the number of idiosyncratic taxa if idiosyncratic.buffer=T) is sufficient to accommodate this combination of taxa assignment for tao2; loop through all tao2 pulses with >= 2 tao1 idiosyncratic taxa; 1 idiosyncratic taxon should be no concern for idiosyncratic buffer; if tao1 idio distribution can accommodate # of tao1 idiosyncratic taxa with this tao2 pulse value and more recent pulse values; since looping is done in order of pulses, looping goes from smallest available idiosyncratic distribution to largest since the order of pulses is from most recent to most ancient (which imposes the upper bounds; thus the most recent tao2 will impose the greatest restriction) and the lower bounds is shared across all pulses
                #this approach is a bit more conservative since it's based on the length of the tao1 idio.prior rather than the range, since if the tao1 idio.prior is not wholly continuous, then it could possibly fit more tao1 idio draws than this statement would indicate, meaning this tao2 taxa assignment could be valid (i.e. fail if statement and not flagged) yet still pass this if statement; however, given how exactly the distribution is discontinued, using the range only may be too liberal, since for example, a distribution with a very large gulf in between two extreme ends, both ends of which have very small limited available prior space, would seem to be a very large distribution due to the range but may only have very limited draws
                if(((sum(tao1.draws[tot.zeta]==0)*buff$tao)+1)>=length(idioprior$tao[[x]][!idioprior$tao[[x]]%in%buffer][idioprior$tao[[x]][!idioprior$tao[[x]]%in%buffer]<tao2.pv[tao2.zeta.real.id[w]]])&&sum(tao1.draws[tot.zeta]==0)>1){
                 lever=1
                }
               }
               if(lever==1){
                w=length(selected.combos)
                while(lever==1){
                 if(w>0){
                  combo.inds[[w]]=combo.inds[[w]][,-selected.combos[w]]
                  if(length(combo.inds[[w]])>0){
                   lever=-1
                  } else {
                   w=w-1
                  }
                 } else {
                  lever=-1
                 }
                }
                if(w==0){
                 circle=1
                 detonator=1
                } else {
                 selected.combos[w]=sample(1:ncol(combo.inds[[w]]),1,replace=T)
                 if(w<length(selected.combos)){
                  tot.zeta=NULL
                  for(v in 1:w){
                   tot.zeta=c(tot.zeta,combo.inds[[v]][,selected[v]])
                  }
                  for(v in (w+1):length(selected.combos)){
                   if(length(allowed.ind[!(allowed.ind+sum(numtaxa[['numtaxa']][0:(x-1)]))%in%tot.zeta])==1){
                    combo.inds[[v]]=matrix((allowed.ind[!(allowed.ind+sum(numtaxa[['numtaxa']][0:(x-1)]))%in%tot.zeta]+sum(numtaxa[['numtaxa']][0:(x-1)])),ncol=1)
                   } else {
                    combo.inds[[v]]=combn((allowed.ind[!(allowed.ind+sum(numtaxa[['numtaxa']][0:(x-1)]))%in%tot.zeta]+sum(numtaxa[['numtaxa']][0:(x-1)])),tao2.zeta[tao2.zeta.real.id[w]])
                   }
                   selected.combos[v]=sample(1:ncol(combo.inds[[v]]),1,replace=T)
                   tot.zeta=c(tot.zeta,combo.inds[[v]][,selected[v]])
                  }
                 }
                }
               } else {
                circle=1
               }
              }
             }
            }
            if(detonator==0){
             tock=0
             for(w in 1:tao2.psi){
              if(w%in%tao2.zeta.real.id){
               tock=tock+1
               indi.combos[[length(indi.combos)]]=append(indi.combos[[length(indi.combos)]], list(combo.inds[[tock]][,selected.combos[tock]]))
              } else {
               indi.combos[[length(indi.combos)]]=append(indi.combos[[length(indi.combos)]], list(NULL))
              }
             }
            }
           }
          }
          if(detonator==0){
           for(x in 1:tao2.psi){
            for(w in 1:num.partitions){
             if(!is.null(indi.combos[[w]][[x]])){
              sim.specs[['tao2']][y,indi.combos[[w]][[x]]]=tao2.pv[x]
             }
            }
           }
          }
         }
         if(detonator==0){
          detonator=1
         } else {
          detonator=0
         }
        }
       }
       if(roll.object[['draws.psi']][[paste(z,sep='')]][y,]>0){
        #mean of shared pulse values of this parameter
        if(eval(parse(text=paste('mean.',z,'.shared',sep='')))==T){
         roll.object[[paste('mean.',z,'.shared',sep='')]][y,]=mean(sim.specs[[length(sim.specs)]][y,][sim.specs[[length(sim.specs)]][y,]!=0])
        }
        #dispersion index of shared pulse values of this parameter
        if(eval(parse(text=paste('disp.index.',z,'.shared',sep='')))==T){
         roll.object[[paste('disp.index.',z,'.shared',sep='')]][y,]=var(sim.specs[[length(sim.specs)]][y,][sim.specs[[length(sim.specs)]][y,]!=0])/mean(sim.specs[[length(sim.specs)]][y,][sim.specs[[length(sim.specs)]][y,]!=0])
        }
       }
      }
     }
    }
   }
   #assign idiosyncratic/nuisance values to individual taxa
   for(z in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
    #if this parameter exists in this model
    if(substr(z,nchar(z),nchar(z))!='2'||2%in%numchanges[['numchanges']]){
     #if mean of this parameter is to be inferred
     if(eval(parse(text=paste('mean.',z,sep='')))==T){
      roll.object=append(roll.object,list(matrix(rep(0,num.sims),ncol=1)))
      names(roll.object)[[length(roll.object)]]=paste('mean.',z,sep='')
     }
     #if dispersion index of this parameter is to be inferred
     if(eval(parse(text=paste('disp.index.',z,sep='')))==T){
      roll.object=append(roll.object,list(matrix(rep(0,num.sims),ncol=1)))
      names(roll.object)[[length(roll.object)]]=paste('disp.index.',z,sep='')
     }
     #if this parameter has idiosyncratic taxa
     if(z%in%names(roll.object$draws.psi)){
      #if statement remains separate from previous if statement, because the attached else statment only applies to the previous clause
      if(idiosyncratic==T){
       for(y in 1:num.sims){
        #if this simulation has idiosyncratic taxa
        if(sum(sim.specs[[paste(z,sep='')]][y,1:sum(numtaxa[['numtaxa']])]==0)>0){
         buffer=NULL
         if(roll.object[['draws.psi']][[paste(z,sep='')]][y,]>0){
          for(x in 1:roll.object[['draws.psi']][[paste(z,sep='')]][y,]){
           if(typeof(buff$tao)=='closure'){
            buff.temp=buff[[paste(z,sep='')]](roll.object$draws.pulse.values[[paste(z,sep='')]][y,x])
            buffer=unique(c(buffer,buff.temp))
           } else {
            buffer=unique(c(buffer,(roll.object[['draws.pulse.values']][[paste(z,sep='')]][y,x]-buff[[paste(z,sep='')]]):(roll.object[['draws.pulse.values']][[paste(z,sep='')]][y,x]+buff[[paste(z,sep='')]])))
           }
          }
         }
         check=0
         while(check==0){
          idio.draws=NULL
          for(x in 1:num.partitions){
           #idiosyncratic taxa
           inds=c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))[sim.specs[[paste(z,sep='')]][y,(sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x])]==0]
           #if this partition has idiosyncratic taxa
           if(length(inds)>0){
            #change.prior outranks idio.prior for tao2
            if(z!='tao2'||is.null(change.prior)){
             idio.dist=idioprior[[paste(z,sep='')]][[x]][!idioprior[[paste(z,sep='')]][[x]]%in%buffer]
            }
            #if tao1 is bounded by tao2 shared pulse values
            if(z=='tao'&&'tao2'%in%names(roll.object$draws.psi)){
             for(w in unique(sim.specs$tao2[y,inds])[unique(sim.specs$tao2[y,inds])!=0]){
              if(length(idio.dist[idio.dist<w])==1){
               idio.draws=rbind(idio.draws,cbind(c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))[sim.specs$tao2[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]==w],rep(idio.dist[idio.dist<w],sum(sim.specs$tao2[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]==w))))
              } else {
               idio.draws=rbind(idio.draws,cbind(c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))[sim.specs$tao2[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]==w],sample(idio.dist[idio.dist<w],sum(sim.specs$tao2[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]==w),replace=T)))
              }
             }
            }
            if(z=='tao2'){
             idio.lengths=NULL
             if(is.null(change.prior)){
              for(w in inds){
               idio.lengths=c(idio.lengths,sum(idio.dist>sim.specs$tao[y,w]))
              }
             } else {
              idio.dists=NULL
              for(w in inds){
               idio.dists=append(idio.dists,list((changeprior[['changeprior']][[x]]+sim.specs$tao[y,w])[!(changeprior[['changeprior']][[x]]+sim.specs$tao[y,w])%in%buffer]))
               idio.lengths=c(idio.lengths,sum(idio.dists[[length(idio.dists)]]>sim.specs$tao[y,w]))
              }
             }
             idio.buffer=NULL
             #smallest to largest distributions
             for(w in inds[order(idio.lengths)]){
              if(!is.null(change.prior)){
               idio.dist=idio.dists[[c(1:length(inds))[inds==w]]]
              }
              if(length(idio.dist[idio.dist>sim.specs$tao[y,w]][!idio.dist[idio.dist>sim.specs$tao[y,w]]%in%idio.buffer])>0){
               #tao2 is bounded by tao1
               if(length(idio.dist[idio.dist>sim.specs$tao[y,w]][!idio.dist[idio.dist>sim.specs$tao[y,w]]%in%idio.buffer])==1){
                idio.draws=rbind(idio.draws,cbind(w,idio.dist[idio.dist>sim.specs$tao[y,w]][!idio.dist[idio.dist>sim.specs$tao[y,w]]%in%idio.buffer]))
               } else {
                idio.draws=rbind(idio.draws,cbind(w,sample(idio.dist[idio.dist>sim.specs$tao[y,w]][!idio.dist[idio.dist>sim.specs$tao[y,w]]%in%idio.buffer],1,replace=T)))
               }
               #for efficiency sake, idiosyncratic.buffer for tao2 operates iteratively from smallest to largest distribution i.e. idiosyncratic taxa later in the loop are restricted in their distribution based on idiosyncratic draws earlier in the loop rather than all idiosyncratic drawing randomly and independently from distributions and checking afterward if fitting idiosyncratic.buffer, continuing to re-draw in a while loop until idiosyncratic.buffer condition is met
               if(idiosyncratic.buffer==T){
                if(typeof(buff$tao)=='closure'){
                 buff.temp=buff$tao2(idio.draws[nrow(idio.draws),2])
                 idio.buffer=unique(c(idio.buffer,buff.temp))
                } else {
                 idio.buffer=unique(c(idio.buffer,(idio.draws[nrow(idio.draws),2]-buff$tao2):(idio.draws[nrow(idio.draws),2]+buff$tao2)))
                }
               }
              #account for possibility of drawing from beyond bounds
              } else {
               v=max(idio.dist)
               #lowest value > max(idio distribution) not in a buffer; protect against tao2 distribution having a smaller max than tao1
               while(v%in%idio.buffer||v%in%buffer||v<=sim.specs$tao[y,w]){
                v=v+1
               }
               idio.draws=rbind(idio.draws,cbind(w,v))
               if(idiosyncratic.buffer==T){
                if(typeof(buff$tao)=='closure'){
                 buff.temp=buff$tao2(v)
                 idio.buffer=unique(c(idio.buffer,buff.temp))
                } else {
                 idio.buffer=unique(c(idio.buffer,(v-buff$tao2):(v+buff$tao2)))
                }
               }
              }
             }
            }
            if(z!='tao2'){
             if(length(idio.dist)==1){
              idio.draws=rbind(idio.draws,cbind(inds[!inds%in%idio.draws[,1]],rep(idio.dist,length(inds[!inds%in%idio.draws[,1]]))))
             } else {
              idio.draws=rbind(idio.draws,cbind(inds[!inds%in%idio.draws[,1]],sample(idio.dist,length(inds[!inds%in%idio.draws[,1]]),replace=T)))
             }
            }
           }
          }
          #if idiosyncratic draws satisfy idiosyncratic.buffer; tao2 has already been checked, and cannot go through this loop in case draws were made beyond bounds
          if(idiosyncratic.buffer==T&&z!='tao2'){
           #if only 1 idiosyncratic taxon in total across all taxa in this simulation, should be no concern for idiosyncratic buffer even if T
           if(nrow(idio.draws)>1){
            x=0
            while(x<nrow(idio.draws)&&check==0){
             x=x+1
             if(typeof(buff[[paste(z,sep='')]])=='closure'){
              if(x>1){
               if(sort(idio.draws[,2])[x]%in%buff[[paste(z,sep='')]](sort(idio.draws[,2])[x-1])){
                check=1
               }
              }
              if(x<nrow(idio.draws)){
               if(sort(idio.draws[,2])[x]%in%buff[[paste(z,sep='')]](sort(idio.draws[,2])[x+1])){
                check=1
               }
              }
             } else {
              if(x<nrow(idio.draws)){
               if((sort(idio.draws[,2])[x+1]-sort(idio.draws[,2])[x])<buff[[paste(z,sep='')]]){
                check=1
               }
              }
             }
            }
           }
          }
          if(check==0){
           check=1
           for(x in 1:nrow(idio.draws)){
            sim.specs[[paste(z,sep='')]][y,idio.draws[x,1]]=idio.draws[x,2]
           }
          } else {
           check=0
          }
         }
        }
       }
      }
     #if this parameter has nuisance taxa
     } else {
      if(z!='tao2'){
       if(substr(z,nchar(z),nchar(z))!='2'||2%in%numchanges[['numchanges']]){
        for(y in 1:num.sims){
         if(z%in%linking[[1]]){
          if(fixing[[1]][c(1:length(linking[[1]])[linking[[1]]==z])]==T){
           if(length(idioprior[[paste(z,'.shared',sep='')]][[1]])==1){
            linked.nuisance=idioprior[[paste(z,'.shared',sep='')]][[1]]
           } else {
            linked.nuisance=sample(idioprior[[paste(z,'.shared',sep='')]][[1]],1,replace=T)
           }
          } else {
           if(fixing[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]==T){
            linked.nuisance=NULL
            for(x in 1:roll.object$draws.psi[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,]){
             if(length(idioprior[[paste(z,'.shared',sep='')]][[x]])==1){
              linked.nuisance=c(linked.nuisance,idioprior[[paste(z,'.shared',sep='')]][[x]])
             } else {
              linked.nuisance=c(linked.nuisance,sample(idioprior[[paste(z,'.shared',sep='')]][[x]],1,replace=T))
             }
            }
           }
          }
         }
         for(x in 1:num.partitions){
          if(substr(z,nchar(z),nchar(z))!='2'||numchanges[['numchanges']][x]==2){
           for(w in (sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x])){
            if(z%in%linking[[1]]){
             if(sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w]%in%roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,]){
              if(fixing[[1]][c(1:length(linking[[1]])[linking[[1]]==z])]==T){
               sim.specs[[paste(z,sep='')]][y,w]=linked.nuisance
              } else {
               if(fixing[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]==T){
                sim.specs[[paste(z,sep='')]][y,w]=linked.nuisance[c(1:roll.object$draws.psi[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,])(roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,]==sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w])]
               } else {
                if(length(idioprior[[paste(z,'.shared',sep='')]][[which(sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w]==roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,])]])==1){
                 sim.specs[[paste(z,sep='')]][y,w]=idioprior[[paste(z,'.shared',sep='')]][[which(sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w]==roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,])]]
                } else {
                 sim.specs[[paste(z,sep='')]][y,w]=sample(idioprior[[paste(z,'.shared',sep='')]][[which(sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w]==roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,])]],1,replace=T)
                }
               }
              }
             }
            }
            if(sim.specs[[paste(z,sep='')]][y,w]==0){
             if(length(idioprior[[paste(z,sep='')]][[x]])==1){
              sim.specs[[paste(z,sep='')]][y,w]=idioprior[[paste(z,sep='')]][[x]]
             } else {
              sim.specs[[paste(z,sep='')]][y,w]=sample(idioprior[[paste(z,sep='')]][[x]],1,replace=T)
             }
            }
           }
          }
         }
        }
       }
      } else {
       if(2%in%numchanges[['numchanges']]){
        #cannot flip anchor because the anchor makes the timing of both events a function of each other, thus flipping would have no effect; should use change if user wants to flip some partitions and not others
        if(!is.null(anchor.prior)){
         roll.object=append(roll.object,list(matrix(rep(0,(num.sims*max(psiprior$tao))),ncol=max(psiprior$tao))))
         names(roll.object)[[length(roll.object)]]='draws.anchors'
         for(y in 1:num.sims){
          for(x in 1:max(psiprior$tao)){
           if(length(anchorprior[['anchorprior']][[x]])==1){
            roll.object$draws.anchors[y,x]=anchorprior[['anchorprior']][[x]]
           } else {
            roll.object$draws.anchors[y,x]=sample(anchorprior[['anchorprior']][[x]],1,replace=T)
           }
           sim.specs$tao2[y,c(1:sum(numtaxa[['numtaxa']]))[sim.specs$tao[y,]==roll.object$draws.pulse.values$tao[y,x]]]=as.integer(roll.object$draws.pulse.values$tao[y,x]+roll.object$draws.anchors[y,x])
          }
          if(length(anchorprior[['anchorprior']][[(max(psiprior$tao)+1)]])==1){
           sim.specs$tao2[y,c(1:sum(numtaxa[['numtaxa']]))[sim.specs$tao2[y,]==0]]=as.integer(sim.specs$tao[y,c(1:sum(numtaxa[['numtaxa']]))[sim.specs$tao2[y,]==0]]+rep(anchorprior[['anchorprior']][[(max(psiprior$tao)+1)]],sum(sim.specs$tao2[y,]==0)))
          } else {
           sim.specs$tao2[y,c(1:sum(numtaxa[['numtaxa']]))[sim.specs$tao2[y,]==0]]=as.integer(sim.specs$tao[y,c(1:sum(numtaxa[['numtaxa']]))[sim.specs$tao2[y,]==0]]+sample(anchorprior[['anchorprior']][[(max(psiprior$tao)+1)]],sum(sim.specs$tao2[y,]==0),replace=T))
          }
         }
        } else {
         if(!is.null(change.prior)){
          for(y in 1:num.sims){
           for(x in 1:num.partitions){
            if(numchanges[['numchanges']][x]==2){
             if(flipvector[['flipvector']][x]==F){
              if(length(changeprior[['changeprior']][[x]])==1){
               sim.specs$tao2[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]=as.integer(sim.specs$tao[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]+rep(changeprior[['changeprior']][[x]],numtaxa[['numtaxa']][x]))
              } else {
               sim.specs$tao2[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]=as.integer(sim.specs$tao[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]+sample(changeprior[['changeprior']][[x]],numtaxa[['numtaxa']][x],replace=T))
              }
             } else {
              if(length(changeprior[['changeprior']][[x]])==1){
               sim.specs$tao2[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]=as.integer(sim.specs$tao[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]-rep(changeprior[['changeprior']][[x]],numtaxa[['numtaxa']][x]))
              } else {
               sim.specs$tao2[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]=as.integer(sim.specs$tao[y,c((sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x]))]-sample(changeprior[['changeprior']][[x]],numtaxa[['numtaxa']][x],replace=T))
              }
             }
            }
           }
          }
         } else {
          for(y in 1:num.sims){
           if(z%in%linking[[1]]){
            if(fixing[[1]][c(1:length(linking[[1]])[linking[[1]]==z])]==T){
             if(length(idioprior[[paste(z,'.shared',sep='')]][[1]])==1){
              linked.nuisance=idioprior[[paste(z,'.shared',sep='')]][[1]]
             } else {
              linked.nuisance=sample(idioprior[[paste(z,'.shared',sep='')]][[1]],1,replace=T)
             }
            } else {
             if(fixing[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]==T){
              linked.nuisance=NULL
              for(x in 1:roll.object$draws.psi[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,]){
               if(length(idioprior[[paste(z,'.shared',sep='')]][[x]])==1){
                linked.nuisance=c(linked.nuisance,idioprior[[paste(z,'.shared',sep='')]][[x]])
               } else {
                linked.nuisance=c(linked.nuisance,sample(idioprior[[paste(z,'.shared',sep='')]][[x]],1,replace=T))
               }
              }
             }
            }
           }
           for(x in 1:num.partitions){
            if(numchanges[['numchanges']][x]==2){
             for(w in (sum(numtaxa[['numtaxa']][0:(x-1)])+1):sum(numtaxa[['numtaxa']][0:x])){
              if(z%in%linking[[1]]){
               if(sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w]%in%roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,]){
                if(fixing[[1]][c(1:length(linking[[1]])[linking[[1]]==z])]==T){
                 sim.specs[[paste(z,sep='')]][y,w]=linked.nuisance
                } else {
                 if(fixing[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]==T){
                  sim.specs[[paste(z,sep='')]][y,w]=linked.nuisance[c(1:roll.object$draws.psi[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,])(roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,]==sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w])]
                 } else {
                  if(length(idioprior[[paste(z,'.shared',sep='')]][[which(sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w]==roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,])]])==1){
                   sim.specs[[paste(z,sep='')]][y,w]=idioprior[[paste(z,'.shared',sep='')]][[which(sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w]==roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,])]]
                  } else {
                   sim.specs[[paste(z,sep='')]][y,w]=sample(idioprior[[paste(z,'.shared',sep='')]][[which(sim.specs[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,w]==roll.object$draws.pulse.values[[linking[[2]][c(1:length(linking[[1]])[linking[[1]]==z])]]][y,])]],1,replace=T)
                  }
                 }
                }
               }
              }
              if(sim.specs[[paste(z,sep='')]][y,w]!=0){
               if(flipvector[['flipvector']][x]==F&&sim.specs$tao[y,w]>=sim.specs$tao2[y,w]){
                sim.specs$tao2[y,w]=sim.specs$tao[y,w]+1
               }
               if(flipvector[['flipvector']][x]==T&&sim.specs$tao[y,w]<=sim.specs$tao2[y,w]){
                sim.specs$tao2[y,w]=sim.specs$tao[y,w]-1
               }
              } else {
               if(flipvector[['flipvector']][x]==F){
                if(sim.specs$tao[y,w]>=max(idioprior$tao2[[x]])){
                 sim.specs$tao2[y,w]=sim.specs$tao[y,w]+1
                } else {
                 if(length(idioprior$tao2[[x]][idioprior$tao2[[x]]>sim.specs$tao[y,w]])==1){
                  sim.specs$tao2[y,w]=idioprior$tao2[[x]][idioprior$tao2[[x]]>sim.specs$tao[y,w]]
                 } else {
                  sim.specs$tao2[y,w]=sample(idioprior$tao2[[x]][idioprior$tao2[[x]]>sim.specs$tao[y,w]],1,replace=T)
                 }
                }
               } else {
                if(sim.specs$tao[y,w]<=min(idioprior$tao2[[x]])){
                 sim.specs$tao2[y,w]=sim.specs$tao[y,w]-1
                } else {
                 if(length(idioprior$tao2[[x]][idioprior$tao2[[x]]<sim.specs$tao[y,w]])==1){
                  sim.specs$tao2[y,w]=idioprior$tao2[[x]][idioprior$tao2[[x]]<sim.specs$tao[y,w]]
                 } else {
                  sim.specs$tao2[y,w]=sample(idioprior$tao2[[x]][idioprior$tao2[[x]]<sim.specs$tao[y,w]],1,replace=T)
                 }
                }
               }
              }
             }
            }
           }
          }
         }
        }
       }
      }
     }
     #mean of this parameter
     if(eval(parse(text=paste('mean.',z,sep='')))==T){
      for(y in 1:num.sims){
       roll.object[[paste('mean.',z,sep='')]][y,]=mean(sim.specs[[paste(z,sep='')]][y,])
      }
     }
     #dispersion index of this parameter
     if(eval(parse(text=paste('disp.index.',z,sep='')))==T){
      for(y in 1:num.sims){
       roll.object[[paste('disp.index.',z,sep='')]][y,]=var(sim.specs[[paste(z,sep='')]][y,])/mean(sim.specs[[paste(z,sep='')]][y,])
       if(is.na(roll.object[[paste('disp.index.',z,sep='')]][y,])){
        roll.object[[paste('disp.index.',z,sep='')]][y,]=0
       }
      }
     }
    }
   }

   output=list(roll.object,sim.specs)
   for(z in names(output[['roll.object']][['draws.zeta']])){
    output[['roll.object']][['draws.zeta']][[paste(z,sep='')]]=output[['roll.object']][['draws.zeta']][[paste(z,sep='')]]/sum(numtaxa[['numtaxa']])
   }
   names(output)=c('roll.object','sim.specs')
   return(output)

  }

 }

}



dice.sims=function(fsc2template.path='fsc2template', fsc2path='fsc25211', output.directory='.', folded=T, append.sims=F, keep.taxa.draws=F, output.hyper.draws=T, output.taxa.draws=F, keep.fsc2.files=F, num.sims=1000000, messages.sims=10000, num.partitions=1, num.taxa, num.haploid.samples, sampling.times=NULL, num.ind.sites=NULL, num.SNPs=NULL, length.seq=NULL, mut.rate, gen.times=NULL, idiosyncratic=T, idiosyncratic.buffer=T, dirichlet.process=F, tao.psi.prior=NULL, epsilon.psi.prior=NULL, NE.psi.prior=NULL, tao.zeta.prior=NULL, tao2.zeta.prior=NULL, epsilon.zeta.prior=NULL, epsilon2.zeta.prior=NULL, NE.zeta.prior=NULL, min.net.zeta.per.pulse=NULL, max.net.zeta.per.pulse=NULL, min.net.zeta.total=NULL, max.net.zeta.total=NULL, tao.shared.prior=NULL, tao2.shared.prior=NULL, epsilon.shared.prior=NULL, epsilon2.shared.prior=NULL, NE.shared.prior=NULL, tao.buffer=0, tao2.buffer=NULL, epsilon.buffer=0, epsilon2.buffer=NULL, NE.buffer=0, net.zeta.per.pulse=F, net.zeta.total=F, mean.tao.shared=F, disp.index.tao.shared=F, mean.tao=F, disp.index.tao=F, mean.tao2.shared=F, disp.index.tao2.shared=F, mean.tao2=F, disp.index.tao2=F, mean.epsilon.shared=F, disp.index.epsilon.shared=F, mean.epsilon=F, disp.index.epsilon=F, mean.epsilon2.shared=F, disp.index.epsilon2.shared=F, mean.epsilon2=F, disp.index.epsilon2=F, mean.NE.shared=F, disp.index.NE.shared=F, mean.NE=F, disp.index.NE=F, tao.idio.prior=NULL, tao2.idio.prior=NULL, epsilon.idio.prior=NULL, epsilon2.idio.prior=NULL, NE.idio.prior=NULL, linked.param=NULL, attached.hyper=NULL, fix.linked.param.across.pulses=F, fix.linked.param.per.pulse=F, num.changes=1, flip=F, anchor.prior=NULL, change.prior=NULL, exponential.growth.rate.prior=NULL, exponential.growth.rate.prior2=NULL, dirichlet.object=NULL, roll.object=NULL, play.object=NULL, input.directory, input.base, input.files, remove.afclasses){

 messages1=0
 if(is.null(roll.object)&&is.null(play.object)){
  roll.object=roll.dice(num.sims=num.sims, num.partitions=num.partitions, num.taxa=num.taxa, idiosyncratic=idiosyncratic, dirichlet.process=dirichlet.process, tao.psi.prior=tao.psi.prior, epsilon.psi.prior=epsilon.psi.prior, NE.psi.prior=NE.psi.prior, tao.zeta.prior=tao.zeta.prior, tao2.zeta.prior=tao2.zeta.prior, epsilon.zeta.prior=epsilon.zeta.prior, epsilon2.zeta.prior=epsilon2.zeta.prior, NE.zeta.prior=NE.zeta.prior, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total, tao.shared.prior=tao.shared.prior, tao2.shared.prior=tao2.shared.prior, epsilon.shared.prior=epsilon.shared.prior, epsilon2.shared.prior=epsilon2.shared.prior, NE.shared.prior=NE.shared.prior, tao.buffer=tao.buffer, tao2.buffer=tao2.buffer, epsilon.buffer=epsilon.buffer, epsilon2.buffer=epsilon2.buffer, NE.buffer=NE.buffer, dirichlet.object=dirichlet.object)
  messages1=1
 }

 temp=NULL
 messages2=0
 if(!is.null(roll.object)&&is.null(play.object)){
  if(keep.taxa.draws==T){
   play.object=play.dice(num.sims=num.sims, num.partitions=num.partitions, num.taxa=num.taxa, idiosyncratic=idiosyncratic, idiosyncratic.buffer=idiosyncratic.buffer, dirichlet.process=dirichlet.process, tao.psi.prior=tao.psi.prior, epsilon.psi.prior=epsilon.psi.prior, NE.psi.prior=NE.psi.prior, tao.zeta.prior=tao.zeta.prior, tao2.zeta.prior=tao2.zeta.prior, epsilon.zeta.prior=epsilon.zeta.prior, epsilon2.zeta.prior=epsilon2.zeta.prior, NE.zeta.prior=NE.zeta.prior, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total, tao.shared.prior=tao.shared.prior, tao2.shared.prior=tao2.shared.prior, epsilon.shared.prior=epsilon.shared.prior, epsilon2.shared.prior=epsilon2.shared.prior, NE.shared.prior=NE.shared.prior, tao.buffer=tao.buffer, tao2.buffer=tao2.buffer, epsilon.buffer=epsilon.buffer, epsilon2.buffer=epsilon2.buffer, NE.buffer=NE.buffer, net.zeta.per.pulse=net.zeta.per.pulse, net.zeta.total=net.zeta.total, mean.tao.shared=mean.tao.shared, disp.index.tao.shared=disp.index.tao.shared, mean.tao=mean.tao, disp.index.tao=disp.index.tao, mean.tao2.shared=mean.tao2.shared, disp.index.tao2.shared=disp.index.tao2.shared, mean.tao2=mean.tao2, disp.index.tao2=disp.index.tao2, mean.epsilon.shared=mean.epsilon.shared, disp.index.epsilon.shared=disp.index.epsilon.shared, mean.epsilon=mean.epsilon, disp.index.epsilon=disp.index.epsilon, mean.epsilon2.shared=mean.epsilon2.shared, disp.index.epsilon2.shared=disp.index.epsilon2.shared, mean.epsilon2=mean.epsilon2, disp.index.epsilon2=disp.index.epsilon2, mean.NE.shared=mean.NE.shared, disp.index.NE.shared=disp.index.NE.shared, mean.NE=mean.NE, disp.index.NE=disp.index.NE, tao.idio.prior=tao.idio.prior, tao2.idio.prior=tao2.idio.prior, epsilon.idio.prior=epsilon.idio.prior, epsilon2.idio.prior=epsilon2.idio.prior, NE.idio.prior=NE.idio.prior, num.changes=num.changes, flip=flip, anchor.prior=anchor.prior, change.prior=change.prior, exponential.growth.rate.prior=exponential.growth.rate.prior, exponential.growth.rate.prior2=exponential.growth.rate.prior2, dirichlet.object=dirichlet.object, roll.object=roll.object)
  } else {
   temp=play.dice(num.sims=1, num.partitions=num.partitions, num.taxa=num.taxa, idiosyncratic=idiosyncratic, idiosyncratic.buffer=idiosyncratic.buffer, dirichlet.process=dirichlet.process, tao.psi.prior=tao.psi.prior, epsilon.psi.prior=epsilon.psi.prior, NE.psi.prior=NE.psi.prior, tao.zeta.prior=tao.zeta.prior, tao2.zeta.prior=tao2.zeta.prior, epsilon.zeta.prior=epsilon.zeta.prior, epsilon2.zeta.prior=epsilon2.zeta.prior, NE.zeta.prior=NE.zeta.prior, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total, tao.shared.prior=tao.shared.prior, tao2.shared.prior=tao2.shared.prior, epsilon.shared.prior=epsilon.shared.prior, epsilon2.shared.prior=epsilon2.shared.prior, NE.shared.prior=NE.shared.prior, tao.buffer=tao.buffer, tao2.buffer=tao2.buffer, epsilon.buffer=epsilon.buffer, epsilon2.buffer=epsilon2.buffer, NE.buffer=NE.buffer, net.zeta.per.pulse=net.zeta.per.pulse, net.zeta.total=net.zeta.total, mean.tao.shared=mean.tao.shared, disp.index.tao.shared=disp.index.tao.shared, mean.tao=mean.tao, disp.index.tao=disp.index.tao, mean.tao2.shared=mean.tao2.shared, disp.index.tao2.shared=disp.index.tao2.shared, mean.tao2=mean.tao2, disp.index.tao2=disp.index.tao2, mean.epsilon.shared=mean.epsilon.shared, disp.index.epsilon.shared=disp.index.epsilon.shared, mean.epsilon=mean.epsilon, disp.index.epsilon=disp.index.epsilon, mean.epsilon2.shared=mean.epsilon2.shared, disp.index.epsilon2.shared=disp.index.epsilon2.shared, mean.epsilon2=mean.epsilon2, disp.index.epsilon2=disp.index.epsilon2, mean.NE.shared=mean.NE.shared, disp.index.NE.shared=disp.index.NE.shared, mean.NE=mean.NE, disp.index.NE=disp.index.NE, tao.idio.prior=tao.idio.prior, tao2.idio.prior=tao2.idio.prior, epsilon.idio.prior=epsilon.idio.prior, epsilon2.idio.prior=epsilon2.idio.prior, NE.idio.prior=NE.idio.prior, num.changes=num.changes, flip=flip, anchor.prior=anchor.prior, change.prior=change.prior, exponential.growth.rate.prior=exponential.growth.rate.prior, exponential.growth.rate.prior2=exponential.growth.rate.prior2, dirichlet.object=dirichlet.object, roll.object=roll.object)
  }
  messages2=1
 }

 if(!is.null(play.object)||!is.null(temp)){

  if(substr(output.directory,nchar(output.directory),nchar(output.directory))=='/'){
   output=paste(output.directory,'dice.sims',sep='')
   output.temp=paste(output.directory,'dice.simulations',sep='')
  } else {
   output=paste(output.directory,'/dice.sims',sep='')
   output.temp=paste(output.directory,'/dice.simulations',sep='')
  }

  if(is.null(play.object)){
   for(a in names(temp$roll.object)){
    if(substr(a,1,4)!='draw'){
     roll.object=append(roll.object,list(matrix(rep(0,(num.sims*length(temp$roll.object[[paste(a,sep='')]]))),ncol=length(temp$roll.object[[paste(a,sep='')]]))))
     names(roll.object)[[length(roll.object)]]=a
    }
   }
  }

  numtaxa=dice.read.numtaxa(num.partitions=num.partitions, num.taxa=num.taxa)

  if(-1%in%numtaxa){
   print("You have included at least 1 non-positive number in the 'num.taxa' object; DICE cannot continue")
  }

  numchanges=dice.read.numchanges(num.partitions=num.partitions, num.changes=num.changes)

  if(-1%in%numchanges){
   print("You have included at least 1 non-positive number in the 'num.changes' object; DICE cannot continue")
  }

  fold=dice.read.fold(folded, num.partitions)

  numhaps=dice.read.numhaps(num.partitions, num.haploid.samples)

  if(-1%in%numhaps){
   print("You have included at least 1 non-positive number in the 'num.haploid.samples' object; DICE cannot continue")
  }

  trap=0

  if(!is.null(sampling.times)){

   samtimes=dice.read.samtimes(num.partitions, sampling.times)

   if(-1%in%samtimes){
    trap=1
    print("You have included at least 1 negative number in the 'sampling.times' object; DICE cannot continue")
   }

  }

  numsites=dice.read.numsites(num.partitions, num.ind.sites, num.SNPs, length.seq)

  if(-1%in%numsites){
   print("You have included at least 1 non-positive number in the 'num.ind.sites', 'num.SNPs', or 'length.seq' object, or have not specified any of these objects at all; DICE cannot continue")
  }

  if(is.null(num.ind.sites)&&is.null(num.SNPs)&&!is.null(length.seq)){

   mutrate=dice.read.mutrate(num.partitions, mut.rate)

   if(-1%in%mutrate){
    trap=1
    print("You have specified an improper 'mut.rate' prior containing at least 1 value that is non-positive; DICE cannot continue")
   }

  }

  if(!is.null(gen.times)){

   gentimes=dice.read.gentimes(num.partitions, gen.times)

   if(-1%in%gentimes){
    trap=1
    print("You have included at least 1 non-positive number in the 'gen.times' object; DICE cannot continue")
   }

  }

  if(append.sims==F){
   a=0
   while(a<sum(numtaxa[['numtaxa']])&&trap==0){
    a=a+1
    if(file.exists(paste(output.temp,a,sep=''))){
     trap=1
     print("You have at least 1 file in your output.directory with the same filename as one of the simulation output files (to append to the existing files, specify the option 'append.sims' to TRUE); DICE cannot continue")
    }
   }
  }

  if(!-1%in%numtaxa&&!-1%in%numchanges&&!-1%in%numhaps&&!-1%in%numsites&&trap==0){

   if(messages1==0){

    if(numtaxa[['numtaxacaution']]==1){
     print("Caution: You have specified an inappropriate number of vector elements for 'num.taxa'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'num.taxa', with x = 'num.partitions', or if not applicable, with the first element of 'num.taxa' used for all partitions")
    }

   }

   if(messages2==0){

    if(numchanges[['numchangescaution']]==1){
     print("Caution: You have specified an inappropriate number of vector elements for 'num.changes', and/or specified an inappropriate number within 'num.changes'; this should be of length = 'num.partitions' or 1, and include values of only 1 or 2 (i.e. DICE currently allows only up to two demographic events per taxon); DICE will continue with only the first x elements of 'num.changes', with x = 'num.partitions', or if not applicable, with the first element of 'num.changes' used for all partitions, as well as any number >2 converted to 2")
    }

   }

   if(fold[['foldcaution']]==1){
    print("Caution: You have specified an inappropriate number of vector elements for 'folded'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'folded', with x = 'num.partitions', or if not applicable, with the first element of 'folded' used for all partitions")
   }

   if(numhaps[['numhapscaution']]==1){
    print("Caution: You have specified an inappropriate number of vector elements for 'num.haploid.samples'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'num.haploid.samples', with x = 'num.partitions', or if not applicable, with the first element of 'num.haploid.samples' used for all partitions")
   }

   if(!is.null(sampling.times)){
    if(samtimes[['samtimescaution']]==1){
     print("Caution: You have specified an inappropriate number of vector elements for 'sampling.times'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'sampling.times', with x = 'num.partitions', or if not applicable, with the first element of 'sampling.times' used for all partitions")
    }
   }

   if(numsites[['numsitescaution']]==1){
    print("Caution: You have specified an inappropriate number of vector elements for 'num.ind.sites', 'num.SNPs', or 'length.seq'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'num.ind.sites', 'num.SNPs', or 'length.seq', with x = 'num.partitions', or if not applicable, with the first element of 'num.ind.sites', 'num.SNPs', or 'length.seq' used for all partitions")
   }

   if(is.null(num.ind.sites)&&is.null(num.SNPs)&&!is.null(length.seq)){
    if(mutrate[['mutratecaution']]==1){
     print("Caution: You have specified an inappropriate number of distributions for your 'mut.rate' prior; the number of distributions should be = 'num.partitions' or 1; DICE will continue with only the first x 'mut.rate' distributions, with x = 'num.partitions', or if not applicable, with the first 'mut.rate' distribution used for all partitions")
    }
   }

   if(!is.null(gen.times)){
    if(gentimes[['gentimescaution']]==1){
     print("Caution: You have specified an inappropriate number of vector elements for 'gen.times'; this should be of length = total number of taxa, 'num.partitions', or 1; DICE will continue with only the first x elements of 'gen.times', with x = total number of taxa, 'num.partitions' (with each element representing an entire partition), or 1 (with the first element representing all taxa)")
    }
   }

   if(output.hyper.draws==T&&!is.null(play.object)){
    for(a in names(play.object$roll.object)){
     if(substr(a,1,4)=='draw'){
      for(b in names(play.object$roll.object[[paste(a,sep='')]])){
       if(substr(a,7,10)=='zeta'){
        write.table(play.object$roll.object[[paste(a,sep='')]][[paste(b,sep='')]]/sum(numtaxa[['numtaxa']]), paste(output,'hyper',a,b,'temp',sep='.'), row.names=F, col.names=F)
       } else {
        write.table(play.object$roll.object[[paste(a,sep='')]][[paste(b,sep='')]], paste(output,'hyper',a,b,'temp',sep='.'), row.names=F, col.names=F)
       }
       system(paste('cat ',output,'.hyper.',a,'.',b,'.temp >>',output,'.hyper.',a,'.',b,sep=''))
       system(paste('rm ',output,'.hyper.',a,'.',b,'.temp',sep=''))
      }
     } else {
      if(sub('zeta','',a)==a){
       write.table(play.object$roll.object[[paste(a,sep='')]], paste(output,'hyper.draws',a,'temp',sep='.'), row.names=F, col.names=F)
      } else {
       write.table(play.object$roll.object[[paste(a,sep='')]]/sum(numtaxa[['numtaxa']]), paste(output,'hyper.draws',a,'temp',sep='.'), row.names=F, col.names=F)
      }
      system(paste('cat ',output,'.hyper.draws.',a,'.temp >>',output,'.hyper.draws.',a,sep=''))
      system(paste('rm ',output,'.hyper.draws.',a,'.temp',sep=''))
     }
    }
   }

   if(output.taxa.draws==T&&!is.null(play.object)){
    for(a in names(play.object$sim.specs)){
     write.table(play.object$sim.specs[[paste(a,sep='')]], paste(output,'taxa.draws',a,'temp',sep='.'), row.names=F, col.names=F)
     system(paste('cat ',output,'.taxa.draws.',a,'.temp >>',output,'.taxa.draws.',a,sep=''))
     system(paste('rm ',output,'.taxa.draws.',a,'.temp',sep=''))
    }
   }

   messages=NULL
   if(!is.null(messages.sims)){
    if(messages.sims>0){
     for(a in 1:round(num.sims/messages.sims)){
      messages=c(messages,(a*messages.sims))
     }
     messages=messages[messages!=num.sims]
    }
   }

   for(a in 1:num.sims){
    if((a-1)%in%messages){
     print(paste('Simulation #',(a-1),' has been completed',sep=''))
    }
    if(is.null(temp)){
     play.object.temp=NULL
     for(b in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE', 'exponential.growth.rate', 'exponential.growth.rate2')){
      if(b %in% names(play.object$sim.specs)){
       play.object.temp=append(play.object.temp,list(play.object$sim.specs[[paste(b,sep='')]][a,]))
       names(play.object.temp)[length(play.object.temp)]=b
      }
     }
    } else {
     roll.object.temp=NULL
     for(b in c('draws.psi', 'draws.pulse.values', 'draws.zeta')){
      roll.object.temp=append(roll.object.temp,list(NULL))
      names(roll.object.temp)[length(roll.object.temp)]=b
      for(c in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE')){
       if(c%in%names(roll.object[['draws.psi']])){
        if(b=='draws.zeta'){
         for(d in 1:num.partitions){
          roll.object.temp[[paste(b,sep='')]]=append(roll.object.temp[[paste(b,sep='')]], list(matrix(roll.object[[paste(b,sep='')]][[paste(c,d,sep='.')]][a,],nrow=1)))
          names(roll.object.temp[[paste(b,sep='')]])[length(roll.object.temp[[paste(b,sep='')]])]=paste(c,d,sep='.')
         }
        } else {
         roll.object.temp[[paste(b,sep='')]]=append(roll.object.temp[[paste(b,sep='')]], list(matrix(roll.object[[paste(b,sep='')]][[paste(c,sep='')]][a,],nrow=1)))
         names(roll.object.temp[[paste(b,sep='')]])[length(roll.object.temp[[paste(b,sep='')]])]=c
        }
       }
      }
     }
     play.object=play.dice(num.sims=1, num.partitions=num.partitions, num.taxa=num.taxa, idiosyncratic=idiosyncratic, idiosyncratic.buffer=idiosyncratic.buffer, dirichlet.process=dirichlet.process, tao.psi.prior=tao.psi.prior, epsilon.psi.prior=epsilon.psi.prior, NE.psi.prior=NE.psi.prior, tao.zeta.prior=tao.zeta.prior, tao2.zeta.prior=tao2.zeta.prior, epsilon.zeta.prior=epsilon.zeta.prior, epsilon2.zeta.prior=epsilon2.zeta.prior, NE.zeta.prior=NE.zeta.prior, min.net.zeta.per.pulse=min.net.zeta.per.pulse, max.net.zeta.per.pulse=max.net.zeta.per.pulse, min.net.zeta.total=min.net.zeta.total, max.net.zeta.total=max.net.zeta.total, tao.shared.prior=tao.shared.prior, tao2.shared.prior=tao2.shared.prior, epsilon.shared.prior=epsilon.shared.prior, epsilon2.shared.prior=epsilon2.shared.prior, NE.shared.prior=NE.shared.prior, tao.buffer=tao.buffer, tao2.buffer=tao2.buffer, epsilon.buffer=epsilon.buffer, epsilon2.buffer=epsilon2.buffer, NE.buffer=NE.buffer, net.zeta.per.pulse=net.zeta.per.pulse, net.zeta.total=net.zeta.total, mean.tao.shared=mean.tao.shared, disp.index.tao.shared=disp.index.tao.shared, mean.tao=mean.tao, disp.index.tao=disp.index.tao, mean.tao2.shared=mean.tao2.shared, disp.index.tao2.shared=disp.index.tao2.shared, mean.tao2=mean.tao2, disp.index.tao2=disp.index.tao2, mean.epsilon.shared=mean.epsilon.shared, disp.index.epsilon.shared=disp.index.epsilon.shared, mean.epsilon=mean.epsilon, disp.index.epsilon=disp.index.epsilon, mean.epsilon2.shared=mean.epsilon2.shared, disp.index.epsilon2.shared=disp.index.epsilon2.shared, mean.epsilon2=mean.epsilon2, disp.index.epsilon2=disp.index.epsilon2, mean.NE.shared=mean.NE.shared, disp.index.NE.shared=disp.index.NE.shared, mean.NE=mean.NE, disp.index.NE=disp.index.NE, tao.idio.prior=tao.idio.prior, tao2.idio.prior=tao2.idio.prior, epsilon.idio.prior=epsilon.idio.prior, epsilon2.idio.prior=epsilon2.idio.prior, NE.idio.prior=NE.idio.prior, num.changes=num.changes, flip=flip, anchor.prior=anchor.prior, change.prior=change.prior, exponential.growth.rate.prior=exponential.growth.rate.prior, exponential.growth.rate.prior2=exponential.growth.rate.prior2, dirichlet.object=dirichlet.object, roll.object=roll.object.temp, skipper=1)
     play.object.temp=NULL
     for(b in c('tao', 'tao2', 'epsilon', 'epsilon2', 'NE', 'exponential.growth.rate', 'exponential.growth.rate2')){
      if(b %in% names(play.object$sim.specs)){
       play.object.temp=append(play.object.temp,list(play.object$sim.specs[[paste(b,sep='')]]))
       names(play.object.temp)[length(play.object.temp)]=b
      }
     }
     for(b in names(play.object$roll.object)){
      if(substr(b,1,4)!='draw'){
       roll.object[[paste(b,sep='')]][a,]=play.object$roll.object[[paste(b,sep='')]][1,]
      }
     }
     if(output.hyper.draws==T){
      for(b in names(play.object$roll.object)){
       if(substr(b,1,4)=='draw'){
        for(c in names(play.object$roll.object[[paste(b,sep='')]])){
         if(substr(b,7,10)=='zeta'){
          write.table(play.object$roll.object[[paste(b,sep='')]][[paste(c,sep='')]]/sum(numtaxa[['numtaxa']]), paste(output,'hyper',b,c,'temp',sep='.'), row.names=F, col.names=F)
         } else {
          write.table(play.object$roll.object[[paste(b,sep='')]][[paste(c,sep='')]], paste(output,'hyper',b,c,'temp',sep='.'), row.names=F, col.names=F)
         }
         system(paste('cat ',output,'.hyper.',b,'.',c,'.temp >>',output,'.hyper.',b,'.',c,sep=''))
         system(paste('rm ',output,'.hyper.',b,'.',c,'.temp',sep=''))
        }
       } else {
        if(sub('zeta','',b)==b){
         write.table(play.object$roll.object[[paste(b,sep='')]], paste(output,'hyper.draws',b,'temp',sep='.'), row.names=F, col.names=F)
        } else {
         write.table(play.object$roll.object[[paste(b,sep='')]]/sum(numtaxa[['numtaxa']]), paste(output,'hyper.draws',b,'temp',sep='.'), row.names=F, col.names=F)
        }
        system(paste('cat ',output,'.hyper.draws.',b,'.temp >>',output,'.hyper.draws.',b,sep=''))
        system(paste('rm ',output,'.hyper.draws.',b,'.temp',sep=''))
       }
      }
     }
     if(output.taxa.draws==T){
      for(b in names(play.object$sim.specs)){
       write.table(play.object$sim.specs[[paste(b,sep='')]], paste(output,'taxa.draws',b,'temp',sep='.'), row.names=F, col.names=F)
       system(paste('cat ',output,'.taxa.draws.',b,'.temp >>',output,'.taxa.draws.',b,sep=''))
       system(paste('rm ',output,'.taxa.draws.',b,'.temp',sep=''))
      }
     }
    }
    for(b in 1:num.partitions){
     if(is.null(num.ind.sites)&&is.null(num.SNPs)&&!is.null(length.seq)){
      asub1=0
      asub2=0
      for(d in 1:(numhaps[['numhaps']][b]-1)){
       asub1=asub1+(1/d)
       asub2=asub2+(1/(d^2))
      }
      bsub1=(numhaps[['numhaps']][b]+1)/(3*(numhaps[['numhaps']][b]-1))
      bsub2=(2*((numhaps[['numhaps']][b]^2)+numhaps[['numhaps']][b]+3))/((9*numhaps[['numhaps']][b])*(numhaps[['numhaps']][b]-1))
      csub1=bsub1-(1/asub1)
      csub2=bsub2-((numhaps[['numhaps']][b]+2)/(asub1*numhaps[['numhaps']][b]))+(asub2/(asub1^2))
      esub1=csub1/asub1
      esub2=csub2/((asub1^2)+asub2)
     }
     for(c in (sum(numtaxa[['numtaxa']][0:(b-1)])+1):sum(numtaxa[['numtaxa']][0:b])){
      system(paste('cp ',fsc2template.path,' ',output,c,'.par',sep=''))
      if(!is.null(sampling.times)){
       system(paste('cat ',output,c,'.par|awk "NR<6" >',output,c,'.par.top',sep=''))
       system(paste('cat ',output,c,'.par|awk "NR==6" >',output,c,'.par.mid0',sep=''))
       system(paste('echo "',samtimes[['samtimes']][b],'" >',output,c,'.par.mid1',sep=''))
       system(paste('paste -d " " ',output,c,'.par.mid0 ',output,c,'.par.mid1 >',output,c,'.par.mid',sep=''))
       system(paste('cat ',output,c,'.par|awk "NR>6" >',output,c,'.par.bot',sep=''))
       system(paste('cat ',output,c,'.par.top ',output,c,'.par.mid ',output,c,'.par.bot >',output,c,'.par',sep=''))
      }
      if(is.null(num.ind.sites)&&!is.null(num.SNPs)){
       system(paste('cat ',output,c,'.par|awk "NR<15" >',output,c,'.par.top',sep=''))
       system(paste('echo "',numsites[['numsites']][b],' 0" >',output,c,'.par.mid',sep=''))
       system(paste('cat ',output,c,'.par|awk "NR>15"|sed -e "s/FREQ/SNP/g" -e "s/2.5e-8 OUTEXP/0/g" >',output,c,'.par.bot',sep=''))
       system(paste('cat ',output,c,'.par.top ',output,c,'.par.mid ',output,c,'.par.bot >',output,c,'.par',sep=''))
      }
      if(is.null(num.ind.sites)&&is.null(num.SNPs)&&!is.null(length.seq)){
       if(length(mutrate[[b]])==1){
        system(paste('cat ',output,c,'.par|sed -e "s/FREQ 1/SNP ',numsites[['numsites']][b],'/g" -e "s/2.5e-8 OUTEXP/',mutrate[[b]],'/g" >',output,c,'.par.temp',sep=''))
       } else {
        system(paste('cat ',output,c,'.par|sed -e "s/FREQ 1/SNP ',numsites[['numsites']][b],'/g" -e "s/2.5e-8 OUTEXP/',sample(mutrate[[b]],1,replace=T),'/g" >',output,c,'.par.temp',sep=''))
       }
       system(paste('mv ',output,c,'.par.temp ',output,c,'.par',sep=''))
      }
      for(d in c('tao','epsilon','NE','num.haploid.samples')){
       if(d=='tao'&&!is.null(gen.times)){
        play.object.temp$tao[c]=play.object.temp$tao[c]/gentimes[['gentimes']][c]
        if(numchanges[['numchanges']][b]==2){
         play.object.temp$tao2[c]=play.object.temp$tao2[c]/gentimes[['gentimes']][c]
        }
       }
       if(d=='NE'){
        find='NPOP'
       } else {
        find=d
       }
       if(d=='epsilon'&&'exponential.growth.rate'%in%names(play.object.temp)){
        replace=1
        system(paste('cat ',output,c,'.par|awk "NR<12" >',output,c,'.par.top',sep=''))
        system(paste('cat ',output,c,'.par|awk "NR==12"|cut -d " " -f 1 >',output,c,'.par.mid0',sep=''))
        hist=read.table(paste(output,c,'.par.mid0',sep=''))+1
        system(paste('echo ',hist,' historical events >',output,c,'.par.mid1',sep=''))
        system(paste('echo "',(round(play.object.temp$tao[c]-(log(play.object.temp$epsilon[c])/play.object.temp$exponential.growth.rate[c]))),' 0 0 0 1 ',play.object.temp$exponential.growth.rate[c],' 0" >',output,c,'.par.mid2',sep=''))
        system(paste('cat ',output,c,'.par|awk "NR>12" >',output,c,'.par.bot',sep=''))
        system(paste('cat ',output,c,'.par.top ',output,c,'.par.mid1 ',output,c,'.par.mid2 ',output,c,'.par.bot >',output,c,'.par',sep=''))
       } else {
        if(d=='num.haploid.samples'){
         replace=numhaps[['numhaps']][b]
        } else {
         replace=play.object.temp[[paste(d,sep='')]][c]
        }
       }
       system(paste('cat ',output,c,'.par|sed -e "s/',find,'/',replace,'/g" >',output,c,'.temp; mv ',output,c,'.temp ',output,c,'.par',sep=''))
       if(d=='tao'&&numchanges[['numchanges']][b]==2){
        if('exponential.growth.rate2'%in%names(play.object.temp)){
         replace=1
         system(paste('cat ',output,c,'.par|awk "NR<14"|sed -e "s/1 historical event/3 historical events/g" >',output,c,'.par.top',sep=''))
         system(paste('echo "',round((play.object.temp$tao2[c]-(log(play.object.temp$epsilon2[c])/play.object.temp$exponential.growth.rate2[c]))),' 0 0 0 1 ',play.object.temp$exponential.growth.rate2[c],' 0" >',output,c,'.par.mid0',sep=''))
         system(paste('cat ',output,c,'.par.top ',output,c,'.par.mid0 >',output,c,'.par.mid; mv ',output,c,'.par.mid ',output,c,'.par.top',sep=''))
        } else {
         replace=play.object.temp$epsilon2[c]
         system(paste('cat ',output,c,'.par|awk "NR<14"|sed -e "s/1 historical event/2 historical events/g" >',output,c,'.par.top',sep=''))
        }
        system(paste('echo "',play.object.temp$tao2[c],' 0 0 0 ',replace,' 0 0" >',output,c,'.par.mid',sep=''))
        system(paste('cat ',output,c,'.par|awk "NR>13" >',output,c,'.par.bot',sep=''))
        system(paste('cat ',output,c,'.par.top ',output,c,'.par.mid ',output,c,'.par.bot >',output,c,'.par',sep=''))
       }
      }
      if(fold[['fold']][b]==T){
       freqtype1='m'
       freqtype2='M'
      } else {
       freqtype1='d'
       freqtype2='D'
      }
      work.dr=getwd()
      setwd(output.directory)
      if(!is.null(num.ind.sites)){
       if(substr(fsc2path,1,1)=='/'){
        system(paste(fsc2path,' -i dice.sims',c,'.par -n ',numsites[['numsites']][b],' -',freqtype1,' -c 0 >dice.log 2>&1',sep=''))
       } else {
        system(paste('./',fsc2path,' -i dice.sims',c,'.par -n ',numsites[['numsites']][b],' -',freqtype1,' -c 0 >dice.log 2>&1',sep=''))
       }
       system(paste('cat dice.sims',c,'/dice.sims',c,'_',freqtype2,'AFpop0.txt|awk "NR==2" >>dice.simulations',c,sep=''))
      } else {
       if(!is.null(num.SNPs)){
        if(substr(fsc2path,1,1)=='/'){
         system(paste(fsc2path,' -i dice.sims',c,'.par -n 1 -',freqtype1,' -c 0 >dice.log 2>&1',sep=''))
        } else {
         system(paste('./',fsc2path,' -i dice.sims',c,'.par -n 1 -',freqtype1,' -c 0 >dice.log 2>&1',sep=''))
        }
        system(paste('cat dice.sims',c,'/dice.sims',c,'_',freqtype2,'AFpop0.obs|awk "NR==3" >dice.sims',c,'.par.top',sep=''))
        convert=read.table(paste('dice.sims',c,'.par.top',sep=''))
        convert=convert/sum(convert)
        write.table(convert,paste('dice.sims',c,'.sfs',sep=''),row.names=F,col.names=F)
        system(paste('cat dice.sims',c,'.sfs >>dice.simulations',c,sep=''))
       } else {
        if(!is.null(length.seq)){
         if(substr(fsc2path,1,1)=='/'){
          system(paste(fsc2path,' -i dice.sims',c,'.par -n 1 -c 0 >dice.log 2>&1',sep=''))
         } else {
          system(paste('./',fsc2path,' -i dice.sims',c,'.par -n 1 -c 0 >dice.log 2>&1',sep=''))
         }
         system(paste('cat dice.sims',c,'/dice.sims',c,'_1_1.arp|awk "NR>23"|awk "NR<',(numhaps[['numhaps']][b]+1),'"|cut -f 3|tr -d " " >dice.sims',c,'.hap',sep=''))
         haplotypes=c(read.table(paste('dice.sims',c,'.hap',sep=''))[[1]])
         sumstats=length(unique(haplotypes))
         hapdiv=0
         for(d in unique(haplotypes)){
          hapdiv=hapdiv+((sum(haplotypes==d)/length(haplotypes))^2)
         }
         sumstats=c(sumstats,((length(haplotypes)/(length(haplotypes)-1))*(1-hapdiv)))
         system(paste("cat dice.sims",c,".hap|sed -e 's/\\(.\\)/\\1 /g' >dice.sims",c,".seq",sep=""))
         sequences=as.matrix(read.table(paste('dice.sims',c,'.seq',sep='')))
         pair.combos=combn(1:length(unique(haplotypes)),2)
         nucdiv=0
         for(d in 1:ncol(pair.combos)){
          nucdiv=nucdiv+((sum(haplotypes==unique(haplotypes)[pair.combos[1,d]])/length(haplotypes))*(sum(haplotypes==unique(haplotypes)[pair.combos[2,d]])/length(haplotypes))*(sum(sequences[which(haplotypes==unique(haplotypes)[pair.combos[1,d]])[1],]!=sequences[which(haplotypes==unique(haplotypes)[pair.combos[2,d]])[1],])/numsites[['numsites']][b]))
         }
         sumstats=c(sumstats,(2*nucdiv))
         pair.combos=combn(1:length(haplotypes),2)
         pwd=0
         for(d in 1:ncol(pair.combos)){
          pwd=pwd+sum(sequences[pair.combos[1,d],]!=sequences[pair.combos[2,d],])
         }
         pwd=pwd/ncol(pair.combos)
         ess=0
         for(d in 1:ncol(sequences)){
          if((max(sequences[,d])-min(sequences[,d]))>0){
           ess=ess+1
          }
         }
         sumstats=c(sumstats,((pwd-(ess/asub1))/(((esub1*ess)+(esub2*ess*(ess-1)))^.5)))
         write.table(matrix(sumstats,nrow=1), sep='\t', paste('dice.sims',c,'.ss',sep=''), row.names=F, col.names=F)
         system(paste('cat dice.sims',c,'.ss >>dice.simulations',c,sep=''))
        }
       }
      }
      if(keep.fsc2.files==T){
       system(paste('cat seed.txt|sed -e "s/ :/:/g"|tr ":" "\n" >>dice.sims',c,'.seed',sep=''))
       if(!is.null(num.ind.sites)||!is.null(num.SNPs)){
        system(paste('paste dice.sims',c,'/dice.sims',c,'.lhoods dice.sims',c,'.seed >>dice.sims',c,'.fsc2',sep=''))
       } else {
        if(!is.null(length.seq)){
         system(paste('cat dice.sims',c,'.par|grep "SNP"|cut -d " " -f 4 >>dice.sims',c,'.fsc2',sep=''))
        }
       }
      }
      setwd(work.dr)
     }
    }
   }
   print(paste('Simulation #',num.sims,' has been completed',sep=''))
   if(substr(output.directory,nchar(output.directory),nchar(output.directory))=='/'){
    system(paste('rm ',output.directory,'seed.txt ',output.directory,'MRCAs.txt ',output.directory,'dice.log 2>/dev/null',sep=''))
   } else {
    system(paste('rm ',output.directory,'/seed.txt ',output.directory,'/MRCAs.txt ',output.directory,'/dice.log 2>/dev/null',sep=''))
   }
   for(a in 1:sum(numtaxa[['numtaxa']])){
    system(paste('rm ',output,a,'.*',sep=''))
    system(paste('rm -r ',output,a,sep=''))
   }

   if(keep.taxa.draws==T){
    return(play.object)
   } else {
    if(!is.null(temp)){
     for(z in names(roll.object[['draws.zeta']])){
      roll.object[['draws.zeta']][[paste(z,sep='')]]=roll.object[['draws.zeta']][[paste(z,sep='')]]/sum(numtaxa[['numtaxa']])
     }
     return(list(roll.object=roll.object))
    } else {
     return(list(roll.object=play.object$roll.object))
    }
   }

  }

 }

}



if(!'bigmemory'%in%rownames(installed.packages())){
 install.packages('bigmemory')
}
library(bigmemory)

dice.aSFS=function(input.directory=NULL, input.base=NULL, input.files=NULL, output.directory='.', folded=T, remove.afclasses=NULL, num.sims=1000000, num.partitions=1, num.taxa, num.haploid.samples, fsc2template.path, fsc2path, append.sims, keep.taxa.draws, output.hyper.draws, output.taxa.draws, keep.fsc2.files, messages.sims, sampling.times, num.ind.sites, num.SNPs, length.seq, mut.rate, gen.times, idiosyncratic, idiosyncratic.buffer, dirichlet.process, tao.psi.prior, epsilon.psi.prior, NE.psi.prior, tao.zeta.prior, tao2.zeta.prior, epsilon.zeta.prior, epsilon2.zeta.prior, NE.zeta.prior, min.net.zeta.per.pulse, max.net.zeta.per.pulse, min.net.zeta.total, max.net.zeta.total, tao.shared.prior, tao2.shared.prior, epsilon.shared.prior, epsilon2.shared.prior, NE.shared.prior, tao.buffer, tao2.buffer, epsilon.buffer, epsilon2.buffer, NE.buffer, net.zeta.per.pulse, net.zeta.total, mean.tao.shared, disp.index.tao.shared, mean.tao, disp.index.tao, mean.tao2.shared, disp.index.tao2.shared, mean.tao2, disp.index.tao2, mean.epsilon.shared, disp.index.epsilon.shared, mean.epsilon, disp.index.epsilon, mean.epsilon2.shared, disp.index.epsilon2.shared, mean.epsilon2, disp.index.epsilon2, mean.NE.shared, disp.index.NE.shared, mean.NE, disp.index.NE, tao.idio.prior, tao2.idio.prior, epsilon.idio.prior, epsilon2.idio.prior, NE.idio.prior, linked.param, attached.hyper, fix.linked.param.across.pulses, fix.linked.param.per.pulse, num.changes, flip, anchor.prior, change.prior, exponential.growth.rate.prior, exponential.growth.rate.prior2, dirichlet.object, roll.object, play.object){

 numtaxa=dice.read.numtaxa(num.partitions=num.partitions, num.taxa=num.taxa)

 if(-1%in%numtaxa){
  print("You have included at least 1 non-positive number in the 'num.taxa' object; DICE cannot continue")
 }

 fold=dice.read.fold(folded, num.partitions)

 numhaps=dice.read.numhaps(num.partitions, num.haploid.samples)

 if(-1%in%numhaps){
  print("You have included at least 1 non-positive number in the 'num.haploid.samples' object; DICE cannot continue")
 }

 onoff=0
 if(!is.null(input.directory)){
  inputdir=dice.read.inputdir(input.directory)
  if(!is.null(input.base)||!is.null(input.files)){
   inputs=dice.read.inputs(input.base, input.files, numtaxa)
   if(-1%in%inputs){
    onoff=1
    print("You have specified less vector elements for 'input.files' than the total number of taxa that you specified in 'num.taxa'; DICE cannot continue")
   }
  } else {
   inputs=NULL
   for(a in 1:sum(numtaxa[['numtaxa']])){
    inputs=c(inputs,paste('dice.simulations',a,sep=''))
   }
   inputs=list(inputs=inputs,inputscaution=0)
  }
 }

 if(!is.null(remove.afclasses)){
  afc=dice.read.afc(remove.afclasses, num.partitions)
 } else {
  afc=NULL
  for(a in 1:num.partitions){
   afc=append(afc,list(0))
  }
  afc=append(afc,list(afccaution=0))
 }

 if(-1%in%afc){
  print("You have included at least 1 negative number in the 'remove.afclasses' object; DICE cannot continue")
 }

 if(!-1%in%numtaxa&&!-1%in%numhaps&&onoff==0&&!-1%in%afc){

  if(numtaxa[['numtaxacaution']]==1){
     print("Caution: You have specified an inappropriate number of vector elements for 'num.taxa'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'num.taxa', with x = 'num.partitions', or if not applicable, with the first element of 'num.taxa' used for all partitions")
  }

  if(fold[['foldcaution']]==1){
   print("Caution: You have specified an inappropriate number of vector elements for 'folded'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'folded', with x = 'num.partitions', or if not applicable, with the first element of 'folded' used for all partitions")
  }

  if(numhaps[['numhapscaution']]==1){
   print("Caution: You have specified an inappropriate number of vector elements for 'num.haploid.samples'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'num.haploid.samples', with x = 'num.partitions', or if not applicable, with the first element of 'num.haploid.samples' used for all partitions")
  }

  if(!is.null(input.directory)){
   if(inputs[['inputscaution']]==1){
    print("You have specified more vector elements for 'input.files' than the total number of taxa that you specified in 'num.taxa'; DICE will continue with only the first x elements of 'input.files', with x = the total number of taxa that you have specified in 'num.taxa'")
   }
  }

  if(afc[['afccaution']]==1){
     print("Caution: You have specified an inappropriate number of lists for 'remove.afclasses'; this should be of length = 'num.partitions' or 1 (a non-list object will be treated as a list of length = 1); DICE will continue with only the first x lists of 'remove.afclasses', with x = 'num.partitions', or if not applicable, with the first list of 'remove.afclasses' used for all partitions")
  }

  SFSs=NULL
  if(!is.null(input.directory)){
   for(a in inputs$inputs){
    SFSs=append(SFSs,list(NULL))
    for(b in inputdir){
     SFSs[[length(SFSs)]]=rbind(SFSs[[length(SFSs)]],read.big.matrix(paste(b,a,sep=''), sep='\t', type='double')[])
    }
   }
  } else {
   for(a in 1:sum(numtaxa[['numtaxa']])){
    if(substr(output.directory,nchar(output.directory),nchar(output.directory))=='/'){
     SFSs=append(SFSs,list(read.big.matrix(paste(output.directory,'dice.simulations',a,sep=''), sep='\t', type='double')))
    } else {
     SFSs=append(SFSs,list(read.big.matrix(paste(output.directory,'/dice.simulations',a,sep=''), sep='\t', type='double')))
    }
   }
  }
  for(a in 1:num.partitions){
   for(b in (sum(numtaxa[['numtaxa']][0:(a-1)])+1):sum(numtaxa[['numtaxa']][0:a])){
    if(fold[['fold']][a]==T){
     SFSs[[b]]=matrix(SFSs[[b]][,1:(round(numhaps[['numhaps']][a]/2)+1)],ncol=length(1:(round(numhaps[['numhaps']][a]/2)+1)))
    } else {
     SFSs[[b]]=matrix(SFSs[[b]][,1:numhaps[['numhaps']][a]],ncol=length(1:numhaps[['numhaps']][a]))
    }
    SFSs[[b]]=matrix(SFSs[[b]][,-(unique(c(0,afc[[a]]))+1)],ncol=(ncol(SFSs[[b]])-length(unique(c(0,afc[[a]])))))
   }
  }
  bins=NULL
  for(a in 1:num.partitions){
   bins=c(bins,(numtaxa[['numtaxa']][a]*ncol(SFSs[[sum(numtaxa[['numtaxa']][1:a])]])))
  }
  aSFSs=matrix(rep(0,(sum(bins)*num.sims)),ncol=sum(bins))
  for(a in 1:num.sims){
   for(b in 1:num.partitions){
    SFS=NULL
    for(c in (sum(numtaxa[['numtaxa']][0:(b-1)])+1):sum(numtaxa[['numtaxa']][0:b])){
     SFS=rbind(SFS,SFSs[[c]][a,]/sum(SFSs[[c]][a,]))
    }
    aSFS=NULL
    for(c in 1:ncol(SFS)){
     aSFS=c(aSFS, sort(SFS[,c], decreasing=T))
    }
    aSFSs[a,(sum(bins[0:(b-1)])+1):sum(bins[0:b])]=aSFS
   }
  }

  return(aSFSs)

 }

}



if(!'fBasics'%in%rownames(installed.packages())){
 install.packages('fBasics')
}
library(fBasics)

dice.sumstats=function(input.directory=NULL, input.base=NULL, input.files=NULL, output.directory='.', num.sims=1000000, num.partitions=1, num.taxa, fsc2template.path, fsc2path, folded, remove.afclasses, append.sims, keep.taxa.draws, output.hyper.draws, output.taxa.draws, keep.fsc2.files, messages.sims, num.haploid.samples, sampling.times, num.ind.sites, num.SNPs, length.seq, mut.rate, gen.times, idiosyncratic, idiosyncratic.buffer, dirichlet.process, tao.psi.prior, epsilon.psi.prior, NE.psi.prior, tao.zeta.prior, tao2.zeta.prior, epsilon.zeta.prior, epsilon2.zeta.prior, NE.zeta.prior, min.net.zeta.per.pulse, max.net.zeta.per.pulse, min.net.zeta.total, max.net.zeta.total, tao.shared.prior, tao2.shared.prior, epsilon.shared.prior, epsilon2.shared.prior, NE.shared.prior, tao.buffer, tao2.buffer, epsilon.buffer, epsilon2.buffer, NE.buffer, net.zeta.per.pulse, net.zeta.total, mean.tao.shared, disp.index.tao.shared, mean.tao, disp.index.tao, mean.tao2.shared, disp.index.tao2.shared, mean.tao2, disp.index.tao2, mean.epsilon.shared, disp.index.epsilon.shared, mean.epsilon, disp.index.epsilon, mean.epsilon2.shared, disp.index.epsilon2.shared, mean.epsilon2, disp.index.epsilon2, mean.NE.shared, disp.index.NE.shared, mean.NE, disp.index.NE, tao.idio.prior, tao2.idio.prior, epsilon.idio.prior, epsilon2.idio.prior, NE.idio.prior, linked.param, attached.hyper, fix.linked.param.across.pulses, fix.linked.param.per.pulse, num.changes, flip, anchor.prior, change.prior, exponential.growth.rate.prior, exponential.growth.rate.prior2, dirichlet.object, roll.object, play.object){

 numtaxa=dice.read.numtaxa(num.partitions=num.partitions, num.taxa=num.taxa)

 if(-1%in%numtaxa){
  print("You have included at least 1 non-positive number in the 'num.taxa' object; DICE cannot continue")
 }

 onoff=0
 if(!is.null(input.directory)){
  inputdir=dice.read.inputdir(input.directory)
  if(!is.null(input.base)||!is.null(input.files)){
   inputs=dice.read.inputs(input.base, input.files, numtaxa)
   if(-1%in%inputs){
    onoff=1
    print("You have specified less vector elements for 'input.files' than the total number of taxa that you specified in 'num.taxa'; DICE cannot continue")
   }
  } else {
   inputs=NULL
   for(a in 1:sum(numtaxa[['numtaxa']])){
    inputs=c(inputs,paste('dice.simulations',a,sep=''))
   }
   inputs=list(inputs=inputs,inputscaution=0)
  }
 }

 if(!-1%in%numtaxa&&onoff==0){

  if(numtaxa[['numtaxacaution']]==1){
     print("Caution: You have specified an inappropriate number of vector elements for 'num.taxa'; this should be of length = 'num.partitions' or 1; DICE will continue with only the first x elements of 'num.taxa', with x = 'num.partitions', or if not applicable, with the first element of 'num.taxa' used for all partitions")
  }

  if(!is.null(input.directory)){
   if(inputs[['inputscaution']]==1){
    print("You have specified more vector elements for 'input.files' than the total number of taxa that you specified in 'num.taxa'; DICE will continue with only the first x elements of 'input.files', with x = the total number of taxa that you have specified in 'num.taxa'")
   }
  }

  sumstats=NULL
  if(!is.null(input.directory)){
   for(a in inputs$inputs){
    sumstats=append(sumstats,list(NULL))
    for(b in inputdir){
     sumstats[[length(sumstats)]]=rbind(sumstats[[length(sumstats)]],read.big.matrix(paste(b,a,sep=''), sep='\t', type='double'))
    }
   }
  } else {
   for(a in 1:sum(numtaxa[['numtaxa']])){
    if(substr(output.directory,nchar(output.directory),nchar(output.directory))=='/'){
     sumstats=append(sumstats,list(read.big.matrix(paste(output.directory,'dice.simulations',a,sep=''), sep='\t', type='double')))
    } else {
     sumstats=append(sumstats,list(read.big.matrix(paste(output.directory,'/dice.simulations',a,sep=''), sep='\t', type='double')))
    }
   }
  }
  sumstatsvectors=matrix(rep(0,(16*num.partitions*num.sims)),ncol=(16*num.partitions))
  for(a in 1:num.sims){
   for(b in 1:num.partitions){
    for(c in 1:4){
     sumstat=NULL
     for(d in (sum(numtaxa[['numtaxa']][0:(b-1)])+1):sum(numtaxa[['numtaxa']][0:b])){
      sumstat=c(sumstat,sumstats[[d]][a,c])
     }
     sumstatsvectors[a,(((b-1)*16)+((c-1)*4)+1)]=mean(sumstat)
     sumstatsvectors[a,(((b-1)*16)+((c-1)*4)+2)]=var(sumstat)
     sumstatsvectors[a,(((b-1)*16)+((c-1)*4)+3)]=skewness(sumstat, method=c('moment'))
     sumstatsvectors[a,(((b-1)*16)+((c-1)*4)+4)]=kurtosis(sumstat, method=c('moment'))
    }
   }
  }

  return(sumstatsvectors)

 }

}
