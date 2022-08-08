#!/bin/sh
plotsr --genomes genomes_lex.txt -o ~/plotsr_lex.png\
	--sr an1_c24syri.out\
	--sr c24_colsyri.out\
	--sr col_cvisyri.out\
	--sr cvi_erisyri.out\
	--sr eri_kyosyri.out\
	--sr kyo_lersyri.out\
	--sr ler_shasyri.out

plotsr --genomes genomes_lex.txt -o ~/plotsr_pansyri.png\
	--sr kyo_cvisyri.out\
	--sr cvi_an1syri.out\
	--sr an1_colsyri.out\
	--sr col_lersyri.out\
	--sr ler_c24syri.out\
	--sr c24_erisyri.out\
	--sr eri_shasyri.out
