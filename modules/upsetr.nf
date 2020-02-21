process upsetr {
      publishDir "${params.output}/upsetr", mode: 'copy', pattern: "upsetr.svg"
      label 'upsetr'
      errorStrategy{task.exitStatus=1 ?'ignore':'terminate'}

    input:
      file(files)

    output:
      file("upsetr.svg")

    script:
      """
      #!/usr/bin/env Rscript
      library(UpSetR)

      files_for_upset <- list.files(path="./", full.names=T, pattern=".txt")
      sets <- lapply(files_for_upset, readLines)
      sets_names <- sapply(gsub("_RENAMED_subset", "", files_for_upset), function(x){
        splits <- strsplit(strsplit(x,"//")[[1]][2],".txt")[[1]]
        return(splits)
      })

      names(sets) <- sets_names
      sets <- sets[order(sapply(sets,length),decreasing=T)]
      svg(filename="upsetr.svg", 
          width=${params.width}, 
          height=${params.height}, 
          pointsize=12)
      upset(fromList(sets), sets = names(sets),
          mainbar.y.label = "No. of common genes", sets.x.label = "No. of identified genes", 
          nsets = 20, nintersects = 40,
          order.by = "freq", sets.bar.color = "#56B4E9", keep.order = F, 
          text.scale = 1.4, point.size = 2.6, line.size = 0.8, set_size.show = TRUE)
      dev.off()
      """
}


process upsetr_subset {
      publishDir "${params.output}/upsetr", mode: 'copy', pattern: "upsetr_subset.svg"
      label 'upsetr'
      errorStrategy{task.exitStatus=1 ?'ignore':'terminate'}

    input:
      file(files)

    output:
      file("upsetr_subset.svg")

    script:
        """
      #!/usr/bin/env Rscript
      library(UpSetR)

      files_for_upset <- list.files(path="./", full.names=T, pattern=".txt")
      sets <- lapply(files_for_upset, readLines)
      sets_names <- sapply(gsub("_RENAMED_subset", "", files_for_upset), function(x){
        splits <- strsplit(strsplit(x,"//")[[1]][2],".txt")[[1]]
        return(splits)
      })

      names(sets) <- sets_names
      sets <- sets[order(sapply(sets,length),decreasing=T)]
      svg(filename="upsetr_subset.svg", 
          width=${params.width}, 
          height=${params.height}, 
          pointsize=12)
      upset(fromList(sets), sets = c(${params.sets}),
          mainbar.y.label = "No. of common genes", sets.x.label = "No. of identified genes", 
          nsets = 20, nintersects = 40,
          order.by = "freq", sets.bar.color = "#56B4E9", keep.order = T, 
          text.scale = 1.4, point.size = 2.6, line.size = 0.8, set_size.show = TRUE)
      dev.off()
        """
}