#!/usr/bin/env bash
# ORF-Scout: prototipo para ensamblajes bacterianos con tabla de traducción 11.
set -euo pipefail
export LC_ALL=C
if [[ ${1:-} == --help ]]; then
  cat <<'HELP'
Uso: bash orf_scout.sh [ensamblaje.fna [carpeta_padre [opción [modo]]]]
Por muestra: bash orf_scout.sh --muestras carpeta_salida opción modo archivo.fna [...]
Opciones: 1 todo, 2 ORF, 3 posibles operones, 4 GC, 0 salir.
Modo de inicio: 1 habitual (ATG GTG TTG), 2 ampliado (+ CTG ATT ATC ATA).
Coordenadas inclusivas, desde 1. El FASTA contiene ADN orientado 5′ a 3′.
Un registro representa el inicio válido más lejano a su STOP en ese marco.
Los inicios alternativos listados producen también ORF de al menos 150 pb.
Los contigs se tratan como lineales; solo se consideran ORF completos.
Una base ambigua interrumpe el ORF. El GC global excluye bases ambiguas.
HELP
  exit
fi
self=$(realpath -- "$0")
if (($# == 0)); then
  printf '%s\n' 'ORF-Scout' '1. Un ensamblaje' '2. Varios ensamblajes'
  read -r -p 'Modo [1]: ' input_mode
  case ${input_mode:-1} in
    1) ;;
    2)
      inputs=()
      while true; do
        read -r -p 'Ruta .fna (Enter para terminar): ' path
        [[ -n $path ]] || break
        inputs+=("$path")
      done
      read -r -p 'Carpeta de resultados [resultados_orf]: ' dest
      read -r -p '1 todo, 2 ORF, 3 operones, 4 GC [1]: ' opt
      read -r -p '1 inicios habituales, 2 ampliados [1]: ' ini
      set -- --muestras "${dest:-resultados_orf}" "${opt:-1}" "${ini:-1}" "${inputs[@]}"
      ;;
    *) echo 'Modo inválido.' >&2; exit 2;;
  esac
fi
if [[ ${1:-} == --muestras || ${1:-} == --lote ]]; then
  shift
  (($# >= 4)) || { echo 'Uso: --muestras salida opción modo archivo.fna [...]' >&2; exit 2; }
  parent=$1; choice=$2; mode=$3; shift 3
  jobs=${ORF_JOBS:-$(nproc)}
  [[ $jobs =~ ^[1-9][0-9]*$ && $choice =~ ^[1-4]$ && $mode =~ ^[12]$ ]] || { echo 'Revisa ORF_JOBS, opción y modo.' >&2; exit 2; }
  mkdir -p -- "$parent"
  declare -A names=()
  paths=()
  for path in "$@"; do
    abs=$(realpath -- "$path")
    name=$(basename -- "$abs"); name=${name%.*}
    [[ -f $abs && -r $abs && -n $name && $name != .* && $name != *"|"* && $abs != *$'\t'* && $abs != *$'\n'* ]] || { echo "Entrada inválida: $path" >&2; exit 2; }
    [[ ! ${names[$name]+x} && ! -e "$parent/$name" ]] || { echo "Nombre repetido o salida existente: $name. Usa nombres distintos." >&2; exit 2; }
    paths+=("$abs"); names[$name]=1
  done
  printf '%s\0' "${paths[@]}" | xargs -0 -n 1 -P "$jobs" bash -c 'bash "$1" "$5" "$2" "$3" "$4"' _ "$self" "$parent" "$choice" "$mode"
  printf 'Resultados por genoma: %s\n' "$(realpath -- "$parent")"
  exit
fi
if (($# > 4)); then echo 'Demasiados argumentos. Usa --help.' >&2; exit 2; fi
input=${1:-}
if [[ -z $input ]]; then read -r -p 'Ruta del ensamblaje .fna: ' input; fi
[[ -f $input && -r $input ]] || { echo 'Archivo inexistente o ilegible.' >&2; exit 1; }
input=$(realpath -- "$input")
parent=${2:-}
if [[ -z $parent ]]; then read -r -p 'Carpeta de resultados [resultados]: ' parent; fi
parent=${parent:-resultados}
choice=${3:-}
if [[ -z $choice ]]; then
  printf '%s\n' 'ORF-Scout' '1. Ejecutar todo' '2. Detectar ORF' '3. Posibles operones' '4. Composición GC' '0. Salir'
  read -r -p 'Opción [1]: ' choice
fi
choice=${choice:-1}
case $choice in 0) exit;; 1|2|3|4) ;; *) echo 'Opción inválida.' >&2; exit 2;; esac
mode=${4:-}
if [[ -z $mode ]]; then read -r -p 'Inicios: 1 habituales, 2 ampliados [1]: ' mode; fi
mode=${mode:-1}
[[ $mode == 1 || $mode == 2 ]] || { echo 'Modo inválido.' >&2; exit 2; }
mkdir -p -- "$parent"
name=$(basename -- "$input"); name=${name%.*}
[[ -n $name && $name != .* && $name != *"|"* && $input != *$'\t'* && $input != *$'\n'* ]] || { echo 'Nombre de muestra inválido.' >&2; exit 2; }
run="$(realpath -- "$parent")/$name"
mkdir -- "$run" || { echo "Elige otra carpeta de salida para $name." >&2; exit 2; }
complete=0
trap 'if (( ! complete )); then rm -rf -- "$run"; fi' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
SECONDS=0
cd -- "$run"
mkdir 01_orf
awk -v mode="$mode" '
function fail(msg) {bad=1; print msg > "/dev/stderr"; exit 1}
function revcomp(s,    i,t) {
  for(i=length(s);i>0;i--) t=t comp[substr(s,i,1)]
  return t
}
function scan(strand,frame,    p,step,c,first,nalt,alts,j,a,b,nt,gc,len,anchor,id) {
  step=(strand=="+" ? 3 : -3)
  for(p=(step>0 ? 1+frame : length(seq)-2-frame); p>=1 && p+2<=length(seq); p+=step) {
    c=substr(seq,p,3)
    if(c ~ /[^ACGT]/) {first=0; nalt=0; continue}
    if(step<0) c=reversecodon[c]
    if(c in start) {
      if(!first) first=p
      else positions[++nalt]=p
    }
    if(c!="TAA" && c!="TAG" && c!="TGA") continue
    if(first && (p-first)*(step>0?1:-1)+3>=150) {
      a=(step>0?first:p); b=(step>0?p:first)+2
      nt=substr(seq,a,b-a+1)
      if(step<0) nt=revcomp(nt)
      len=length(nt); alts=""
      for(j=1;j<=nalt;j++) {
        if((p-positions[j])*(step>0?1:-1)+3<150) break
        anchor=positions[j]+(step<0?2:0)
        alts=alts (alts?",":"") anchor
      }
      id=sprintf("ORF_%07d",++count)
      print ">" id " " contig ":" a "-" b "(" strand ")" > fasta
      for(j=1;j<=len;j+=80) print substr(nt,j,80) > fasta
      gc=gsub(/[GC]/,"",nt)
      print contig,a,b,strand,id,len,(alts?alts:"."),sprintf("%.8f",100*gc/len) > rows
    }
    first=0; nalt=0
  }
}
function finish(    f,copy,i) {
  if(contig=="") return
  seq=""
  for(i=1;i<=chunks;i++) {seq=seq chunk[i]; delete chunk[i]}
  seq=seq buffer; buffer=""; chunks=0
  if(!length(seq)) fail("Contig vacío: " contig)
  copy=seq; global_gc+=gsub(/[GC]/,"",copy)
  copy=seq; canonical+=gsub(/[ACGT]/,"",copy)
  for(f=0;f<3;f++) {scan("+",f); scan("-",f)}
}
BEGIN {
  OFS="\t"; rows="01_orf/orf.tsv"; fasta="01_orf/orf.fna"
  printf "" > fasta
  print "contig","inicio","fin","hebra","id","longitud_pb","inicios_alt_5p","gc_pct" > rows
  comp["A"]="T";comp["C"]="G";comp["G"]="C";comp["T"]="A"
  for(i in comp) for(j in comp) for(k in comp) reversecodon[i j k]=comp[k] comp[j] comp[i]
  split("ATG GTG TTG" (mode==2?" CTG ATT ATC ATA":""),list," ")
  for(i in list) start[list[i]]=1
}
{ sub(/\r$/,"") }
/^>/ {
  finish(); header=substr($0,2); sub(/^[[:space:]]+/,"",header)
  split(header,fields,/[[:space:]]+/); contig=fields[1]
  if(!length(contig) || seen[contig]++) fail("Identificador FASTA vacío o repetido: " contig)
  ncontigs++; seq=""; next
}
{
  gsub(/[[:space:]]/,""); if(!length($0)) next
  if(!length(contig)) fail("Secuencia sin encabezado FASTA.")
  line=toupper($0)
  if(line ~ /[^ACGTRYSWKMBDHVN]/) fail("Símbolo de ADN inválido en " contig)
  buffer=buffer line
  if(length(buffer)>=32768) {chunk[++chunks]=buffer; buffer=""}
}
END {
  if(bad) exit 1
  finish()
  if(!canonical) fail("Ensamblaje sin bases A/C/G/T.")
  print "contigs",ncontigs
  print "orf_candidatos",count+0
  print "gc_ensamblaje_pct",sprintf("%.8f",100*global_gc/canonical)
}' "$input" > resumen.tsv
printf 'opcion\t%s\nmodo_inicios\t%s\nmin_orf_pb\t150\n' "$choice" "$mode" >> resumen.tsv
if [[ $choice == 1 || $choice == 3 ]]; then
  mkdir 02_operones
  printf 'grupo\tcontig\thebra\torf_orden_genomico\tdistancias_pb\n' > 02_operones/grupos.tsv
  # Seleccionar vecinos compatibles, priorizando ORF largos en zonas solapadas.
  # Un índice por ventanas evita comparar cada candidato contra todo el genoma.
  tail -n +2 01_orf/orf.tsv | sort -t $'\t' -k1,1 -k6,6nr -k2,2n -k5,5 |
    awk -F '\t' '
    $1 != contig {delete windows; contig=$1}
    {
      left=int($2/150); right=int($3/150); conflict=0
      for(w=left;w<=right && !conflict;w++){
        count=split(windows[w],near," ")
        for(j=1;j<=count;j++){
          k=near[j]
          overlap=($3<finish[k]?$3:finish[k])-($2>begin[k]?$2:begin[k])+1
          if(overlap>4){conflict=1;break}
        }
      }
      if(conflict) next
      begin[NR]=$2; finish[NR]=$3
      for(w=left;w<=right;w++) windows[w]=windows[w] (windows[w]?" ":"") NR
      print
    }' | sort -t $'\t' -k1,1 -k2,2n |
    awk -F '\t' 'BEGIN{OFS="\t"}
    function flush(){if(n>=2){print ++groups,chr,strand,ids,gaps >> "02_operones/grupos.tsv"; members+=n}}
    {d=$2-end-1
     if(n && $1==chr && $4==strand && d>=-4 && d<=75 && $3>end){
       ids=ids "," $5; gaps=gaps (n>1?",":"") d; n++
     }else{flush(); n=1; ids=$5; gaps=""}
     chr=$1; strand=$4; end=$3
    }
    END{flush();print "posibles_operones",groups+0;print "orf_en_grupos",members+0}' >> resumen.tsv
fi
if [[ $choice == 1 || $choice == 4 ]]; then
  mkdir 03_gc
  awk -F '\t' 'NR>1 && $6>=300 {print $8}' 01_orf/orf.tsv | sort -n > gc_ordenado.tmp
  genome=$(awk -F '\t' '$1=="gc_ensamblaje_pct"{print $2}' resumen.tsv)
  awk -F '\t' -v genome="$genome" 'BEGIN{OFS="\t";out="03_gc/gc.tsv"
    print "orf","longitud_pb","gc_pct","delta_pp","estado" > out
    while((getline x < "gc_ordenado.tmp")>0) v[++n]=x
    if(n>=20){q1=q(.25);q3=q(.75); lo=q1-1.5*(q3-q1);hi=q3+1.5*(q3-q1)}
  }
  function q(p,  i,j){i=1+(n-1)*p;j=int(i);return v[j]+(i-j)*(v[j+1]-v[j])}
  NR>1{
    if($6<300) state="corto"
    else if(n<20) state="muestra_insuficiente"
    else if($8+0<lo){state="bajo";low++}
    else if($8+0>hi){state="alto";high++}
    else state="dentro_del_intervalo"
    print $5,$6,$8,sprintf("%.4f",$8-genome),state > out
  }
  END{print "orf_elegibles_gc",n+0
    if(n>=20){print "limite_gc_inferior",lo;print "limite_gc_superior",hi;print "gc_bajo",low+0;print "gc_alto",high+0;print "gc_atipico",low+high}
    else print "gc_atipico","no evaluable (menos de 20 ORF de 300 pb)"
  }' 01_orf/orf.tsv >> resumen.tsv
  rm gc_ordenado.tmp
fi
printf 'tiempo_segundos\t%s\n' "$SECONDS" >> resumen.tsv
cat resumen.tsv
printf '\nResultados: %s\n' "$run"
printf '  %s\n' */ resumen.tsv
printf '%s\n' "$input" > origen.txt
complete=1
