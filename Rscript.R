library(dplyr)
library(rlist)
library(openxlsx)
library(stringr)
readT1map = function(path){
  con <- file(path, open = "rt", raw = TRUE)
  textraw <- readLines(con, skipNul = TRUE)
  close(con)
  for(i in 1:length(textraw)){
    textraw[i] = toString(textraw[[i]])
  }
  textraw = textraw[-which(str_detect(textraw,'-------------'))]
  textraw = textraw[-which(textraw == "")]
  return(textraw)
}

Info_table = function(textdata){
  sortinformation = sort(c(which(textdata == 'Native T1'),which(textdata == 'CA T1'),which (textdata =='ECV'),which(textdata=='\xbb\003'),length(textdata)))
  if(exists(c('maininfo','T1Info','ECV_LambdaInfo'))){rm(maininfo,T1Info,ECV_LambdaInfo)}
  for(i in sortinformation){
    if(i == sortinformation[1]){
      key = rep(NA,i-1)
      value = rep(NA, i-1)
      for(j in 1:(i-1)){
        key[j] = strsplit(textdata[j],'\t',useBytes = T)[[1]][1]
        value[j] = strsplit(textdata[j],'\t',useBytes = T)[[1]][2]
      }
      maininfo = as.data.frame(value[-1],key[-1])
      colnames(maininfo) = "Main Info"
    }
    else{
      if(textdata[sortinformation[match(i,sortinformation)-1]] %in% c('Native T1','CA T1')){
        textinfo = textdata[sortinformation[match(i,sortinformation)-1]:(i-1)]
        tableNum = which(str_detect(textinfo,'Global Myo'))
        key = rep(NA,tableNum)
        value = rep(NA, tableNum)
        for(j in 1:(tableNum)){
          key[j] = strsplit(textinfo[j],'\t',useBytes = T)[[1]][1]
          value[j] = strsplit(textinfo[j],'\t',useBytes = T)[[1]][2]
        }
        if(!exists('T1Info')){
          T1Info = as.data.frame(value[-1],key[-1])
          colnames(T1Info)[match(i,sortinformation)-1] = key[1]
        }
        else{
          newcolumn = rep(NA,nrow(T1Info))
          for(k in 1:(length(key))){
            if(key[k] %in% rownames(T1Info)){
              newcolumn[k-1] = value[k]
            }
          }
          
          if(length(newcolumn) != nrow(T1Info)){
            newcolumn = newcolumn[which(newcolumn!='NA')]
          }
          T1Info = cbind(T1Info,newcolumn)
          colnames(T1Info)[match(i,sortinformation)-1] = key[1]
        }
      }
      else{
        textinfo = textdata[sortinformation[match(i,sortinformation)-1]:(i-1)]
        tableNum = which(str_detect(textinfo,'Global Myo'))
        key = rep(NA,tableNum)
        value = rep(NA, tableNum)
        for(j in 1:(tableNum)){
          key[j] = strsplit(textinfo[j],'\t',useBytes = T)[[1]][1]
          value[j] = strsplit(textinfo[j],'\t',useBytes = T)[[1]][2]
        }
        if(!exists('ECV_LambdaInfo')){
          ECV_LambdaInfo = as.data.frame(value[-1],key[-1])
          colnames(ECV_LambdaInfo)[match(i,sortinformation)-3] = key[1]
        }
        else{
          newcolumn = rep(NA,nrow(ECV_LambdaInfo))
          for(k in 2:(length(key))){
            if(key[k] %in% rownames(ECV_LambdaInfo)){
              newcolumn[k-1] = value[k]
            }
          }
          newcolumn[length(key)] = value[length(key)]
          ECV_LambdaInfo = cbind(ECV_LambdaInfo,newcolumn)
          if(key[1]=='\xbb\003'){
            key[1] = 'Lambda'
          }
          colnames(ECV_LambdaInfo)[match(i,sortinformation)-3] = key[1]
        }
      }
    }
  }
  rtn = list()
  if(exists('maininfo')){
    rtn=list.append(rtn,maininfo)
  }
  if(exists('T1Info')){
    rtn=list.append(rtn,T1Info)
  }
  if(exists('ECV_LambdaInfo')){
    rtn=list.append(rtn,ECV_LambdaInfo)
  }
  return(rtn)
}

Segment_Table = function(textdata){
  sortinformation = sort(c(which(textdata == 'Native T1'),which(textdata == 'CA T1'),which (textdata =='ECV'),which(textdata =='\xbb\003'),length(textdata)))
  if(exists(c('SegmentTable'))){rm(SegmentTable)}
  
  for(i in 1:(length(sortinformation)-1)){
    name = textdata[sortinformation[i]]
    testdata = textdata[sortinformation[i]:(sortinformation[i+1]-1)]
    segmentnum = which(str_detect(testdata,'Segment	Mean'))
    cname = paste(testdata[segmentnum-2],testdata[segmentnum],sep = " ")
    key = rep(NA, 16)
    for(j in (segmentnum+1):(segmentnum+16)){
      key[j-segmentnum] = strsplit(testdata[j],'\t')[[1]][2]
    }
    if(!exists('SegmentTable')){
      SegmentTable = as.data.frame(key)
      colnames(SegmentTable)[i] = cname
    }
    else{
      SegmentTable = cbind(SegmentTable,key)
      colnames(SegmentTable)[i] = cname
    }
  }
  if(ncol(SegmentTable)>2){
    colnames(SegmentTable)[4] = 'Regional Lambda (AHA Segmentation) Segment\tMean Lambda (%)'
  }
  
  return(SegmentTable)
}

Regional_Table = function(textdata){
  sortinformation = sort(c(which(textdata == 'Native T1'),which(textdata == 'CA T1'),which (textdata =='ECV'),which(textdata =='\xbb\003'),length(textdata)))
  rtn = list()
  for(i in 1:(length(sortinformation)-1)){
    name = textdata[sortinformation[i]]
    testdata = textdata[sortinformation[i]:(sortinformation[i+1])]
    checkname = paste('Regional', name, 'Slice',sep = ' ')
    dataIdx = which(str_detect(testdata, checkname))
    if(name == '\xbb\003'){
      name = 'Lambda'
    }
    for(k in 1:length(dataIdx)){
      for(j in (dataIdx[k]+2):(dataIdx[k]+9)){
        if(k == 1 && j == dataIdx[k]+2){
          tbl = data.frame(paste('Slice_',k,'Segment',strsplit(testdata[j],'\t')[[1]][1]),
                           strsplit(testdata[j],'\t')[[1]][2],
                           strsplit(testdata[j],'\t')[[1]][3],
                           strsplit(testdata[j],'\t')[[1]][4],
                           strsplit(testdata[j],'\t')[[1]][5],
                           strsplit(testdata[j],'\t')[[1]][6],
                           strsplit(testdata[j],'\t')[[1]][7])
        }
        else{
          tbl = rbind(tbl,data.frame(paste('Slice_',k,'Segment',strsplit(testdata[j],'\t')[[1]][1]),
                                     strsplit(testdata[j],'\t')[[1]][2],
                                     strsplit(testdata[j],'\t')[[1]][3],
                                     strsplit(testdata[j],'\t')[[1]][4],
                                     strsplit(testdata[j],'\t')[[1]][5],
                                     strsplit(testdata[j],'\t')[[1]][6],
                                     strsplit(testdata[j],'\t')[[1]][7]))
        }
        
      }
    }
    colnames(tbl) = strsplit(testdata[dataIdx[1]+1],'\t',useBytes = T)[[1]]
    colnames(tbl)[7] = 'Area(mmsq2)'
    colnames(tbl)[1]  = paste(name,'_',colnames(tbl)[1])
    rtn = list.append(rtn,tbl)
    
  }
  return(rtn)
}

patientReport = function(path){
  textdata = readT1map(path)
  infoTable = Info_table(textdata)
  segmentTable = Segment_Table(textdata)
  regionalTable = Regional_Table(textdata)
  report = createWorkbook()
  addWorksheet(report,'Main Info')
  addWorksheet(report, 'Segment record')
  addWorksheet(report,'T1 Info')
  addWorksheet(report, 'Native T1 Region Segment')
  addWorksheet(report, 'CA T1 Region Segment')
  if(length(infoTable) == 3){
    addWorksheet(report,'ECV Labmda Info')
  }
  if(length(regionalTable) == 4){
    addWorksheet(report, 'ECV Region Segment')
    addWorksheet(report, 'LAMBDA Region Segment')
  }
  
  writeData(report, sheet = 'Main Info',x = data.frame(infoTable[[1]]),rowNames = T)
  writeData(report, sheet = 'Segment record',x = as.data.frame(segmentTable),rowNames = T)
  writeData(report, sheet = 'T1 Info',x = data.frame(infoTable[[2]]),rowNames = T)
  writeData(report, sheet = 'Native T1 Region Segment',x = as.data.frame(regionalTable[[1]]),rowNames = T)
  writeData(report, sheet = 'CA T1 Region Segment',x = as.data.frame(regionalTable[[2]]),rowNames = T)
  if(length(infoTable) == 3){
    writeData(report, sheet = 'ECV Labmda Info',x = as.data.frame(infoTable[[3]]),rowNames = T)
  }
  if(length(regionalTable) == 4){
    writeData(report, sheet = 'ECV Region Segment',x = as.data.frame(regionalTable[[3]]),rowNames = T)
    writeData(report, sheet = 'LAMBDA Region Segment',x = as.data.frame(regionalTable[[4]]),rowNames = T)
  }
  PatientID = as.data.frame(infoTable[[1]])[3,]
  saveWorkbook(report,paste('Patient',PatientID,'report.xlsx'))
}
