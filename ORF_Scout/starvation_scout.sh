#!/usr/bin/env bash
# Comparación de marcadores de escasez de carbono por identidad proteica y cobertura.
set -euo pipefail
export LC_ALL=C
base=$(dirname -- "$(realpath -- "$0")")
catalog="$base/catalogo/marcadores.tsv"
if [[ ${1:-} == --help ]]; then
  echo 'Uso: bash starvation_scout.sh [carpeta_salida opción muestra1.fna muestra2.fna ...]'
  echo 'Desde ORF: bash starvation_scout.sh --orf carpeta_muestras carpeta_salida opción [otra_carpeta ...]'
  echo '1 todo, 2 marcadores, 3 sistemas, 4 variación del ADN, 0 salir.'
  exit
fi
batch=""
roots=(); labels=()
explicit_orf=0
if [[ ${1:-} == --orf ]]; then
  (($# >= 4)) || { echo 'Uso: --orf carpeta_muestras carpeta_salida opción [otra_carpeta ...]' >&2; exit 2; }
  explicit_orf=1
  batch=$(realpath -- "$2"); parent=$3; choice=$4; shift 4; roots=("$batch" "$@")
elif (($#)); then
  (($# >= 4)) || { echo 'Se requieren salida, opción y al menos dos .fna.' >&2; exit 2; }
  parent=$1; choice=$2; shift 2; samples=("$@")
else
  samples=()
  read -r -p 'Carpeta de una muestra o carpeta con todas las muestras (Enter para ingresar ensamblajes): ' batch
  if [[ -n $batch ]]; then
    batch=$(realpath -- "$batch")
  else
  while true; do
    read -r -p 'Archivo .fna (Enter al terminar): ' file
    [[ -n $file ]] || break
    samples+=("$file")
  done
  fi
  read -r -p 'Carpeta de salida [resultados]: ' parent
  printf '%s\n' 'Starvation-Scout' '1. Ejecutar todo' '2. Detectar marcadores' '3. Comparar sistemas' '4. Comparar variantes de ADN' '0. Salir'
  read -r -p 'Opción [1]: ' choice
fi
parent=${parent:-resultados}; choice=${choice:-1}
case "$choice" in 0) exit;; 1|2|3|4) ;; *) echo 'Opción inválida.' >&2; exit 2;; esac
if [[ -n $batch ]]; then
  samples=(); dirs=(); declare -A used=()
  add_sample(){
    local dir=$1 label=$2 original=$3
    [[ -n $label && $label != */* && $label != .* && $label != *"|"* && ! ${used[$label]+x} && -f "$dir/orf.fna" && -s "$dir/orf.tsv" ]] || { echo "Muestra repetida o incompleta: $label" >&2; exit 2; }
    used[$label]=1; labels+=("$label"); samples+=("$original"); dirs+=("$dir")
  }
  (( ${#roots[@]} )) || roots=("$batch")
  for root in "${roots[@]}"; do
    root=$(realpath -- "$root")
    if [[ -f "$root/origen.txt" ]]; then
      add_sample "$root/01_orf" "$(basename -- "$root")" "$(<"$root/origen.txt")"
    elif [[ -f "$root/muestras.tsv" ]]; then
      while IFS=$'\t' read -r id original relative; do
        add_sample "$root/$relative" "$id" "$original"
      done < <(tail -n +2 "$root/muestras.tsv")
    else
      found=0
      for source in "$root"/*/origen.txt; do
        [[ -f $source ]] || continue
        dir=${source%/*}; add_sample "$dir/01_orf" "$(basename -- "$dir")" "$(<"$source")"
        found=1
      done
      ((found)) || { echo "No hay muestras terminadas en $root." >&2; exit 2; }
    fi
  done
  (( ${#samples[@]} )) || { echo 'Sin muestras.' >&2; exit 2; }
else
(( ${#samples[@]} >= 2 )) || { echo 'Ingresa al menos dos ensamblajes.' >&2; exit 2; }
declare -A seen=() named=()
for i in "${!samples[@]}"; do
  path=$(realpath -- "${samples[i]}")
  [[ -f $path && -r $path && ! ${seen[$path]+x} ]] || { echo "Archivo repetido o ilegible: $path" >&2; exit 2; }
  label=$(basename -- "$path"); label=${label%.*}
  [[ -n $label && $label != .* && $label != *"|"* && $label != *$'\t'* && $label != *$'\n'* && ! ${named[$label]+x} ]] || { echo "Nombre inválido o repetido: $label" >&2; exit 2; }
  labels+=("$label"); named[$label]=1
  samples[i]=$path; seen[$path]=1
done
fi
if (( explicit_orf )); then
  [[ ! -e $parent ]] || { echo "La carpeta de resultados ya existe: $parent" >&2; exit 2; }
  mkdir -p -- "$(dirname -- "$parent")"
  mkdir -- "$parent"
  run=$(realpath -- "$parent")
else
  mkdir -p -- "$parent"
  run=$(mktemp -d -- "$(realpath -- "$parent")/starvation_scout.XXXXXXXX")
fi
complete=0
trap 'if (( ! complete )); then rm -rf -- "$run"; fi' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
cd -- "$run"; SECONDS=0
mkdir -p 01_marcadores/secuencias tmp
printf 'id\tarchivo\n' > muestras.tsv
printf 'muestra\tgen\taccesion\tcontig\tinicio\tfin\thebra\torf\tinicio_elegido\tsistemas\tidentidad_aa_pct\tcobertura_ref_aa_pct\tidentidad_adn_ref_pct\tcobertura_orf_aa_pct\taa_comparados\tasignacion\n' > 01_marcadores/hallazgos.tsv
printf 'muestra\torf\tcontig\treferencia\testado\tidentidad_aa_pct\tcobertura_orf_pct\tcobertura_ref_pct\taa_comparados\n' > 01_marcadores/revision.tsv
touch tmp/alelos.tsv
printf '[1/3] Entradas validadas. Muestras: %d. Marcadores en catálogo: %d.\n' "${#samples[@]}" "$(awk 'END{print NR-1}' "$catalog")"
for i in "${!samples[@]}"; do
  sample=$(printf 'M%03d' "$((i+1))")
  sample=${labels[i]:-$sample}
  printf '%s\t%s\n' "$sample" "${samples[i]}" >> muestras.tsv
  mkdir "tmp/$sample"
  printf '[2/3] Buscando marcadores en %s (%s)\n' "$sample" "$(basename -- "${samples[i]}")"
  if [[ -n $batch ]]; then
    files=("${dirs[i]}")
  else
  bash "$base/orf_scout.sh" "${samples[i]}" "tmp/$sample" 2 1 > "tmp/$sample/ejecucion.log"
  files=(tmp/"$sample"/*/01_orf)
  fi
  awk -F '\t' -v sample="$sample" -v refs="$catalog" -v coords="${files[0]}/orf.tsv" '
  BEGIN{
    OFS="\t"; bases="TCAG"; amino="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    for(i=1;i<=4;i++) for(j=1;j<=4;j++) for(k=1;k<=4;k++)
      code[substr(bases,i,1) substr(bases,j,1) substr(bases,k,1)]=substr(amino,++n,1)
    getline < refs
    while((getline < refs)>0){systemList[$1]=$2; acc[$1]=$3; protein[$1]=$4; referenceDNA[$1]=$6;index_ref($1)}
    close(refs)
    getline < coords
    while((getline < coords)>0){location[$5]=$1; left[$5]=$2;right[$5]=$3;strand[$5]=$4; allowedStarts[$5]=$7}
    close(coords)
  }
  # Indexar palabras de cinco aminoácidos para proponer desplazamientos.
  function index_ref(g,   j,word,key) {
    for(j=1;j+4<=length(protein[g]);j++){
      word=substr(protein[g],j,5); freq[word]++; key=word SUBSEP freq[word]
      seedgene[key]=g; seedpos[key]=j
    }
  }
  function evaluate(g,delta,    ref,nr,nq,j,q0,r0,len,bestDist,p,first,c,dist,h,pid,qcov,rcov,tier,score){
    ref=protein[g]; nr=length(ref); nq=length(peptide)
    # Inicio válido más próximo al comienzo proyectado de la referencia.
    bestDist=nq+nr+1
    for(j=1;j<=nStarts;j++){
      p=starts[j]; dist=p-(1+delta); if(dist<0) dist=-dist
      if(dist<bestDist){first=p;bestDist=dist}
    }
    q0=(first>1+delta?first:1+delta); r0=q0-delta
    len=nq-q0+1; if(nr-r0+1<len) len=nr-r0+1
    if(len<30) return
    qcov=100*len/(nq-first+1); rcov=100*len/nr
    if(qcov<30 || rcov<30) return
    h=0
    for(j=0;j<len;j++){
      c=(q0+j==first?"M":substr(peptide,q0+j,1))
      if(c==substr(ref,r0+j,1)) h++
    }
    pid=100*h/len
    if(pid<70) return
    tier=(pid>=90 && qcov>=80 && rcov>=80 ? 2 : 1)
    score=2*h-len
    if(tier>level[g] || (tier==level[g] && score>points[g])){
      level[g]=tier;points[g]=score;ident[g]=pid;covQ[g]=qcov;covR[g]=rcov
      startAA[g]=first;span[g]=len;queryBegin[g]=q0;refBegin[g]=r0
    }
  }
  function process(    j,k,p,word,key,g,nq,delta,i,winner,strong,nt,a,b,kind,sets,total,path,dnaMatches,qnt,rnt,basesCompared,status){
    peptide=""; nStarts=0; delete visited;delete level;delete points
    for(p=1;p<length(dna)-2;p+=3){
      codon=substr(dna,p,3); peptide=peptide code[codon]
      if(p==1) starts[++nStarts]=1
    }
    if(allowedStarts[id]!="."){
      total=split(allowedStarts[id],anchors,",")
      for(j=1;j<=total;j++){
        delta=(strand[id]=="+" ? anchors[j]-left[id] : right[id]-anchors[j])
        starts[++nStarts]=1+delta/3
      }
    }
    for(j=1;j+4<=length(peptide);j+=3){
      word=substr(peptide,j,5)
      if(!(word in freq)) continue
      for(k=1;k<=freq[word];k++){
        key=word SUBSEP k;g=seedgene[key];delta=j-seedpos[key]
        if(!((g SUBSEP delta) in visited)){visited[g,delta]=1;evaluate(g,delta)}
      }
    }
    strong=0
    for(g in level) if(level[g]==2) strong++
    for(g in level){
      if(level[g]<1) continue
      if(level[g]==1 || strong>1){
        status=(level[g]==2?"ambiguo":"revisar")
        print sample,id,location[id],g,status,sprintf("%.2f",ident[g]),sprintf("%.2f",covQ[g]),sprintf("%.2f",covR[g]),span[g] >> "01_marcadores/revision.tsv"
        continue
      }
      p=3*(startAA[g]-1);nt=substr(dna,p+1)
      a=left[id];b=right[id];if(strand[id]=="+") a+=p;else b-=p
      kind=(p?"alternativo":"principal")
      qnt=3*(queryBegin[g]-1)+1;rnt=3*(refBegin[g]-1)+1;basesCompared=3*span[g]
      if(queryBegin[g]+span[g]-1==length(peptide) && refBegin[g]+span[g]-1==length(protein[g])) basesCompared+=3
      dnaMatches=0
      for(j=0;j<basesCompared;j++) if(substr(dna,qnt+j,1)==substr(referenceDNA[g],rnt+j,1)) dnaMatches++
      status=(ident[g]==100 && covQ[g]==100 && covR[g]==100?"identica":"probable")
      if(covQ[g]<100 || covR[g]<100) status="probable_cobertura_parcial"
      print sample,g,acc[g],location[id],a,b,strand[id],id,kind,systemList[g],sprintf("%.2f",ident[g]),sprintf("%.2f",covR[g]),sprintf("%.2f",100*dnaMatches/basesCompared),sprintf("%.2f",covQ[g]),span[g],status >> "01_marcadores/hallazgos.tsv"
      print sample,g,nt >> "tmp/alelos.tsv"
      total=split(systemList[g],sets,",")
      for(i=1;i<=total;i++){
        path="01_marcadores/secuencias/" sets[i] ".fna"
        print ">" sample "|" g "|" id "|" location[id] ":" a "-" b "(" strand[id] ")" >> path
        for(j=1;j<=length(nt);j+=80) print substr(nt,j,80) >> path
      }
    }
  }
  /^>/{if(id!="")process();split(substr($0,2),header," ");id=header[1];dna="";next}
  {dna=dna $0}
  END{if(id!="")process()}' "${files[0]}/orf.fna"
  rm -rf -- "tmp/$sample"
done
awk -F '\t' -v choice="$choice" -v catfile="$catalog" '
BEGIN{
 OFS="\t"; matrix="01_marcadores/matriz.tsv"; systems="sistemas.tsv"; variants="variacion.tsv"
 getline < catfile
 while((getline < catfile)>0){
   genes[++ng]=$1; n=split($2,parts,",")
   for(i=1;i<=n;i++){s=parts[i]; if(!(s in size)) order[++ns]=s; size[s]++; member[s,$1]=1}
 }
 close(catfile); getline < "muestras.tsv"
 while((getline < "muestras.tsv")>0) samples[++nm]=$1
 close("muestras.tsv")
}
{ key=$1 SUBSEP $2; copies[key]++
  if(!(key in dna)) dna[key]=$3; else if(dna[key]!=$3) mixed[key]=1
}
END{
 printf "muestra" > matrix
 for(g=1;g<=ng;g++) printf "\t%s",genes[g] > matrix
 print "" > matrix
 if(choice==1 || choice==3) print "muestra","sistema","detectados","esperados","genes_detectados","componentes_pct" > systems
 for(i=1;i<=nm;i++){
   m=samples[i]; found=0;printf "%s",m > matrix
   for(g=1;g<=ng;g++){v=copies[m,genes[g]]+0; printf "\t%d",v > matrix;found+=(v>0)}
   print "" > matrix; print m,found "/" ng " marcadores",sprintf("%.2f%%",100*found/ng)
   if(choice==1 || choice==3) for(z=1;z<=ns;z++){
     s=order[z];count=0;names=""
     for(g=1;g<=ng;g++) if(member[s,genes[g]] && copies[m,genes[g]]){
       count++;names=names (names?",":"") genes[g]
     }
     print m,s,count,size[s],(names?names:"."),sprintf("%.2f",100*count/size[s]) > systems
   }
 }
 if(choice==1 || choice==4){
   print "muestra_A","muestra_B","gen","resultado" > variants
   for(i=1;i<=nm;i++) for(j=i+1;j<=nm;j++) for(g=1;g<=ng;g++){
     a=samples[i] SUBSEP genes[g]; b=samples[j] SUBSEP genes[g]
     if(!copies[a] || !copies[b]) state="no_evaluable"
     else if(mixed[a] || mixed[b]) state="multiples_alelos"
     else state=(dna[a]==dna[b]?"identico":"diferente")
     print samples[i],samples[j],genes[g],state > variants
   }
 }
}' tmp/alelos.tsv | tee resumen.tsv
printf '[3/3] Matrices, tablas por sistema y secuencias FASTA generadas.\n'
if [[ -f sistemas.tsv ]]; then mkdir 02_sistemas; mv sistemas.tsv 02_sistemas/; fi
if [[ -f variacion.tsv ]]; then mkdir 03_variacion; mv variacion.tsv 03_variacion/; fi
mkdir por_muestra
while IFS=$'\t' read -r name original; do
  dir="por_muestra/$name"; mkdir -p -- "$dir/secuencias"
  for table in 01_marcadores/hallazgos.tsv 01_marcadores/revision.tsv 02_sistemas/sistemas.tsv; do
    [[ -f $table ]] || continue
    awk -F '\t' -v name="$name" 'NR==1 || $1==name' "$table" > "$dir/${table##*/}"
  done
done < <(tail -n +2 muestras.tsv)
for fasta in 01_marcadores/secuencias/*.fna; do
  [[ -f $fasta ]] || continue
  awk -v file="${fasta##*/}" '/^>/{if(out)close(out); split(substr($0,2),h,"[|]");out="por_muestra/" h[1] "/secuencias/" file} {print >> out}' "$fasta"
done
rm -rf -- tmp
printf 'tiempo_segundos\t%s\n' "$SECONDS" | tee -a resumen.tsv
printf 'Resultados: %s\n' "$run"
printf '  %s\n' */ resumen.tsv muestras.tsv
complete=1
