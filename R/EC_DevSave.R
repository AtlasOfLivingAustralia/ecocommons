EC_DevSave <- function(fileroot, ext=".pdf") {
  if (ext==".eps") {dev.copy2eps(file=paste(fileroot,ext,sep="."))
  } else {
      dev.copy2pdf(file=paste(fileroot,"pdf",sep="."))}
}
