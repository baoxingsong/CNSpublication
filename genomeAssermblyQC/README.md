run 
``` 
perl ~/software/mypipeline/song_tabasco_pipeline.pl Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa 0.6
```

you will have a result file named `tabasco_summary`

run R code:
```
data =read.table("tabasco_summary")                                                                                                                                                                       
mean(data$V2)
```

It will report the number of copies for each TABASCO orthogroup.




song_tabasco_pipeline.r is used to draw a plot for the song_tabasco_pipeline.pl output(s)
and could generate a BUSCO similar plot
