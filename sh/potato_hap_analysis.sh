# Potato hap assemblies copied from  copied from /dss/dsslegfs01/pn29fi/pn29fi-dss-0012/02_reference/Haplotype_resolved_potato_assemblies_v2_20240126/

# remove scaffolds and only keep chromosomes
chars=({A..J} O)
for g in ${chars[@]}; do
  {
  for i in 1 2 3 4; do
    f=${g}_hap${i}_genome.fasta
    echo $f
    b=$(echo $f | sed 's/\.fasta//g')
    scafs=$(grep 'PGA' $f |  sed 's/>//g')
    hometools getchr -v --chrs ${scafs[@]} -o ${b}_filtered.fasta $f
  done
  } &
done
