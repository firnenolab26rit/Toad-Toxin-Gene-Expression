# Generating a Better Quality Transcriptome By Reducing Redundant Transcripts



## Methods: 

### DeNovo Run Stats

**Busco Results***
```
# tertapoda_odb12 dataset 
        ***** Results: *****

        C:81.4%[S:19.4%,D:62.0%],F:4.1%,M:14.5%,n:5623
        4578    Complete BUSCOs (C)
        1089    Complete and single-copy BUSCOs (S)
        3489    Complete and duplicated BUSCOs (D)
        228     Fragmented BUSCOs (F)
        817     Missing BUSCOs (M)
        5623    Total BUSCO groups searched
```

**N_Sequences**
```622,491 ```

### 1.) EvidentialGene tr2aacds.pl 
Ran tr2acds.pl pipeline and retained the ```/okeyset*.okay.mrna``` transcripts. 

**Busco Results***
```
# tertapoda_0db12 dataset
        ***** Results: *****

        C:81.3%[S:44.1%,D:37.1%],F:4.0%,M:14.7%,n:5623
        4570    Complete BUSCOs (C)
        2482    Complete and single-copy BUSCOs (S)
        2088    Complete and duplicated BUSCOs (D)
        224     Fragmented BUSCOs (F)
        829     Missing BUSCOs (M)
        5623    Total BUSCO groups searched
```

**N_Sequences**
```159,857 ```
* This is the version of the transcriptome we used.

### 2.) CD-Hit + TD2 


**Busco Results***
```
# BUSCO was run in mode: euk_tran

        ***** Results: *****

        C:77.2%[S:53.2%,D:23.9%],F:5.4%,M:17.4%,n:5623
        4339    Complete BUSCOs (C)
        2994    Complete and single-copy BUSCOs (S)
        1345    Complete and duplicated BUSCOs (D)
        306     Fragmented BUSCOs (F)
        978     Missing BUSCOs (M)
        5623    Total BUSCO groups searched
```

**N_Sequences**
```72,553 ```
