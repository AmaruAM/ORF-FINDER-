#!/usr/bin/env bash
# Coincidencias de proteínas de respuesta a limitación nutricional en E. coli.
set -euo pipefail
export LC_ALL=C
base=$(dirname -- "$(realpath -- "$0")")
catalog="$base/catalogo/marcadores.tsv"
if [[ ${1:-} == --help ]]; then
  echo 'Uso: bash starvation_scout.sh [carpeta_salida opción muestra.fna ...]'
  echo 'Desde ORF: bash starvation_scout.sh --orf carpeta_muestras carpeta_salida opción'
  echo '1 coincidencias por sistema, 2 coincidencias y comparación de variantes, 0 salir.'
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
  (($# >= 3)) || { echo 'Se requieren salida, opción y al menos un .fna.' >&2; exit 2; }
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
  printf '%s\n' 'Starvation-Scout' '1. Buscar coincidencias por sistema' '2. Buscar coincidencias y comparar variantes' '0. Salir'
  read -r -p 'Opción [1]: ' choice
fi
parent=${parent:-resultados}; choice=${choice:-1}
case "$choice" in 0) exit;; 1|2) ;; *) echo 'Opción inválida.' >&2; exit 2;; esac
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
(( ${#samples[@]} >= 1 )) || { echo 'Ingresa al menos un ensamblaje.' >&2; exit 2; }
if [[ $choice == 2 && ${#samples[@]} -lt 2 ]]; then
  echo 'Para comparar variantes se requieren al menos dos muestras.' >&2
  exit 2
fi
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
mkdir -p tmp
printf 'id\tarchivo\n' > muestras.tsv
printf 'muestra\tsistema\tproteina\taccesion_referencia\torf\tcontig\tinicio\tfin\thebra\tinicio_usado\tidentidad_aa_pct\tcobertura_referencia_pct\tcobertura_orf_pct\tidentidad_adn_pct\taa_comparados\testado\n' > coincidencias.tsv
touch tmp/alelos.tsv
printf '[1/3] Entradas validadas. Muestras: %d. Proteínas de referencia: %d.\n' "${#samples[@]}" "$(awk 'END{print NR-1}' "$catalog")"
for i in "${!samples[@]}"; do
  sample=$(printf 'M%03d' "$((i+1))")
  sample=${labels[i]:-$sample}
  printf '%s\t%s\n' "$sample" "${samples[i]}" >> muestras.tsv
  mkdir "tmp/$sample"
  printf '[2/3] Comparando proteínas con el catálogo para %s (%s)\n' "$sample" "$(basename -- "${samples[i]}")"
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
  function process(    j,k,p,word,key,g,nq,delta,i,winner,strong,nt,a,b,kind,sets,total,dnaMatches,qnt,rnt,basesCompared,status){
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
        total=split(systemList[g],sets,",")
        for(i=1;i<=total;i++) print sample,sets[i],g,acc[g],id,location[id],left[id],right[id],strand[id],"no_asignado",sprintf("%.2f",ident[g]),sprintf("%.2f",covR[g]),sprintf("%.2f",covQ[g]),"NA",span[g],status >> "coincidencias.tsv"
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
      total=split(systemList[g],sets,",")
      for(i=1;i<=total;i++) print sample,sets[i],g,acc[g],id,location[id],a,b,strand[id],kind,sprintf("%.2f",ident[g]),sprintf("%.2f",covR[g]),sprintf("%.2f",covQ[g]),sprintf("%.2f",100*dnaMatches/basesCompared),span[g],status >> "coincidencias.tsv"
      print sample,g,nt >> "tmp/alelos.tsv"
    }
  }
  /^>/{if(id!="")process();split(substr($0,2),header," ");id=header[1];dna="";next}
  {dna=dna $0}
  END{if(id!="")process()}' "${files[0]}/orf.fna"
  rm -rf -- "tmp/$sample"
done
if [[ $choice == 2 && ${#samples[@]} -lt 2 ]]; then
  echo 'Para comparar variantes se requieren al menos dos muestras.' >&2
  exit 2
fi
sorted=$(mktemp)
{ head -n 1 coincidencias.tsv; tail -n +2 coincidencias.tsv | sort -t $'\t' -k2,2 -k1,1 -k3,3; } > "$sorted"
mv -- "$sorted" tmp/coincidencias_todas.tsv
awk -F '\t' -v total="$(awk 'END{print NR-1}' "$catalog")" '
BEGIN{OFS="\t"; print "muestra","proteinas_asignadas","proteinas_en_catalogo","asignadas_pct","candidatos_revision"}
NR>1 {key=$1 SUBSEP $3; if($16=="revisar" || $16=="ambiguo") review[key]=1; else assigned[key]=1; samples[$1]=1}
END{
  while((getline line < "muestras.tsv")>0){
    if(line !~ /^id\t/){split(line,f,"\t");sample=f[1];a=0;r=0
      for(k in assigned)if(index(k,sample SUBSEP)==1)a++
      for(k in review)if(index(k,sample SUBSEP)==1)r++
      print sample,a,total,sprintf("%.2f",100*a/total),r
    }
  }
}' tmp/coincidencias_todas.tsv | tee resumen.tsv
mkdir -p por_cepa
while IFS=$'\t' read -r sample original; do
  [[ $sample == id ]] && continue
  mkdir -p -- "por_cepa/$sample"
  awk -F '\t' -v sample="$sample" 'NR==1 || $1==sample' tmp/coincidencias_todas.tsv > "por_cepa/$sample/coincidencias.tsv"
done < muestras.tsv
if [[ $choice == 2 ]]; then
  awk -F '\t' -v catfile="$catalog" '
  BEGIN{
    OFS="\t"; variants="variantes.tsv"
    getline < catfile
    while((getline < catfile)>0) genes[++ng]=$1
    close(catfile)
    getline < "muestras.tsv"
    while((getline < "muestras.tsv")>0){split($0,f,"\t");samples[++nm]=f[1]}
    close("muestras.tsv")
  }
  {key=$1 SUBSEP $2;copies[key]++;if(!(key in dna))dna[key]=$3;else if(dna[key]!=$3)mixed[key]=1}
  END{
    print "muestra_A","muestra_B","proteina","estado" > variants
    for(i=1;i<=nm;i++)for(j=i+1;j<=nm;j++)for(g=1;g<=ng;g++){
      a=samples[i] SUBSEP genes[g];b=samples[j] SUBSEP genes[g]
      if(!copies[a]||!copies[b])state="no_evaluable"
      else if(mixed[a]||mixed[b])state="multiples_alelos"
      else state=(dna[a]==dna[b]?"identico":"diferente")
      print samples[i],samples[j],genes[g],state > variants
    }
  }' tmp/alelos.tsv
fi
printf '[3/3] Tablas de coincidencias separadas por cepa y ordenadas por sistema.\n'
rm -rf -- tmp
printf 'Tiempo total: %s s\n' "$SECONDS"
printf 'Resultados: %s\n' "$run"
printf '  por_cepa/CEPA/coincidencias.tsv\n  resumen.tsv\n  muestras.tsv\n'
[[ $choice != 2 ]] || printf '  variantes.tsv\n'
complete=1
