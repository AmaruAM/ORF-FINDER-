# ORF-Scout

ORF-Scout procesa ensamblajes bacterianos en FASTA de ADN. Para esta prueba, usa la tabla de traducción bacteriana 11 y ORF completos de al menos 150 pb. La identificación de ORF es una predicción por codones de inicio y terminación, no una anotación funcional definitiva.

## Inicio interactivo y corridas por grupo

Ejecuta `./orf_scout.sh` desde `ORF_Scout/` o `bash orf_scout.sh`. El menú solicita un nombre de grupo, acepta archivos `.fna` individuales o carpetas que los contienen, y crea una corrida nueva en `corridas/`. Copia las entradas para conservar la procedencia. Un nombre de grupo puede repetirse porque cada ejecución recibe una carpeta única.

Al concluir ORF-Scout, presiona Enter para terminar la primera etapa. La anotación se inicia solo si respondes `s`. Para ejecutarla más tarde, copia el comando que muestra ORF-Scout. Cada mensaje `[n/N]` indica la etapa que está avanzando.

```text
corridas/FECHA_GRUPO_ID/
├── genomas/                 copias de los ensamblajes de entrada
├── muestras.tsv             nombre de muestra y ruta original
├── ORF/
│   ├── MG1655/
│   │   ├── 01_orf/
│   │   │   ├── orf.fna      secuencias de ADN de ORF candidatos
│   │   │   ├── orf.faa      traducción de los ORF a proteínas
│   │   │   └── orf.tsv      coordenadas, hebra, longitud e inicio
│   │   ├── 02_operones/      grupos de ORF vecinos, si se calcularon
│   │   ├── 03_gc/            composición GC, si se calculó
│   │   └── resumen.tsv       parámetros y conteos de la muestra
│   └── otra_cepa/            misma estructura para cada genoma
└── Starvation/               aparece al ejecutar el segundo análisis
```

## Archivos de ORF

- `orf.fna` contiene ADN orientado de cada ORF candidato, incluido el codón STOP.
- `orf.faa` contiene la traducción proteica, sin el STOP terminal. Un inicio bacteriano alternativo se traduce como metionina inicial.
- `orf.tsv` contiene una fila por ORF. `contig`, `inicio`, `fin` y `hebra` ubican la secuencia en el ensamblaje. `longitud_pb`, `inicios_alt_5p` y `gc_pct` describen el candidato.
- `02_operones/grupos.tsv` recoge grupos de ORF consecutivos, en la misma hebra y separados entre −4 y 75 pb. Son candidatos basados en vecindad, no operones confirmados por promotor o expresión.
- `03_gc/gc.tsv` registra el GC de cada ORF, la diferencia frente al ensamblaje en puntos porcentuales y el estado estadístico. Un GC atípico no demuestra transferencia horizontal.
- `resumen.tsv` presenta contigs, GC del ensamblaje, conteo de ORF, grupos de vecindad, datos GC, parámetros y tiempo.

## Parámetros del detector

En el modo habitual se usan ATG, GTG y TTG. El modo ampliado añade CTG, ATT, ATC y ATA. Se buscan STOP TAA, TAG y TGA en los seis marcos. Los contigs se tratan como lineales y las bases ambiguas interrumpen la búsqueda. El detector guarda el inicio válido más lejano al STOP de cada marco y registra otros inicios elegibles. Los ORF solapados pueden aparecer como candidatos independientes.

Para reproducir el análisis de una muestra sin preguntas:

```bash
bash orf_scout.sh ejemplos/MG1655.fna resultados 1 1
```

El tercer argumento ejecuta todo (1), solo ORF (2), posibles operones (3) o GC (4). El cuarto selecciona inicios habituales (1) o ampliados (2). Se requiere Bash 4+, awk y utilidades GNU habituales.

## Referencias

1. NCBI ORFfinder. https://www.ncbi.nlm.nih.gov/orffinder/
2. NCBI Genetic Codes. https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
