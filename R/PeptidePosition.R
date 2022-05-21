
#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export


PeptidePosition = function(ProteinList,PeptideList,database){
  Peptide_query<-ProteinList
  
  
  mouseUp <- database
  
  keys <- Peptide_query$Uniprot
  columns <- c("PROTEIN-NAMES","SEQUENCE")
  kt <- "UNIPROTKB"
  res <- UniProt.ws::select(mouseUp, keys, columns, kt)
  Peptide_ProteinSeq=res
  
  
  data_protein_seq<-Peptide_ProteinSeq# change to use itself information 
  # get the peptide location
  data_peptide=PeptideList
  for(n in 1: nrow(data_peptide)){# begin 0
    
    data_peptide[n,"Protein.Id"]
    data_protein_seq_test<-subset(data_protein_seq,UNIPROTKB==data_peptide[n,"Protein.Id"])
    if(!is.na(data_protein_seq_test[1,"SEQUENCE"])&&nchar(data_protein_seq_test[1,"SEQUENCE"])>10){# begin1
      a=data_protein_seq_test[1,"SEQUENCE"]
      data_peptide[n,"PeptideLocation"]=indexOf(a,data_peptide[n,"PeptideSeq"])
      data_peptide[n,"PeptideLen"]=nchar(data_peptide[n,"PeptideSeq"])
      data_peptide[n,"ProteinLen"]=nchar(data_protein_seq_test[1,"SEQUENCE"])
      data_peptide[n,"RTPosition"]=indexOf(a,data_peptide[n,"PeptideSeq"])/nchar(data_protein_seq_test[1,"SEQUENCE"])
      ##print(data_peptide[n,1])
      ##print(n)
      
    }# final1
  }
  return(data_peptide)
  
}
indexOf = function(str,str2){
  cd=nchar(str);
  cd2=nchar(str2);
  if(cd==0||cd2==0){
    return(0);
  }
  for(i in 1:cd){
    t=substr(str,i,i);
    for(j in 1:cd2){
      if(t==substr(str2,j,j)&&j==1){
        if(cd2==1){
          return(i);
        }else{
          c=TRUE;
          for(k in 1:(cd2-1)){
            if(substr(str,i+k,i+k)!=substr(str2,j+k,j+k)){
              c=FALSE;
              break;
            }
          }
          if(c==TRUE){
            return(i);
          }
        }
      }else{
        break;
      }
    }
  }
  return(0);
}