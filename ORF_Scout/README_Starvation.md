# Starvation-Scout

Starvation-Scout compara ORF predichos con un catálogo local de **79 marcadores candidatos** organizados en nueve módulos de limitación nutricional en *Escherichia coli*. El catálogo incluye referencias proteicas de MG1655. Está pensado para comparar cepas de *E. coli*, pero no garantiza detectar genes ausentes, muy divergentes o no incluidos en el catálogo.

Los módulos son captación de carbono, reservas de carbono, respuesta estricta, supervivencia general/fase estacionaria, captación de nitrógeno, captación de fosfato, captación de azufre, adquisición de hierro y regulación asociada a limitación de aminoácidos. La respuesta real incluye muchos genes adicionales y depende del nutriente, la cepa y las condiciones. Este catálogo prioriza reguladores, componentes de captación y algunos efectores de supervivencia. No equivale a un perfil completo de expresión génica.

## Ejecución independiente

Desde la carpeta del proyecto, usa el comando que entrega ORF-Scout, o ejecuta:

```bash
bash starvation_scout.sh --orf corridas/MI_GRUPO/ORF corridas/MI_GRUPO/Starvation 1
```

La carpeta de salida debe ser nueva. La opción 1 genera todos los análisis. Las opciones 2, 3 y 4 generan respectivamente marcadores, cobertura por sistema o comparación de variantes de ADN. El segundo script funciona sin internet y no llama a anotadores externos.

## Cómo leer las salidas

```text
Starvation/
├── muestras.tsv
├── resumen.tsv
├── 01_marcadores/
│   ├── hallazgos.tsv
│   ├── revision.tsv
│   ├── matriz.tsv
│   └── secuencias/
├── 02_sistemas/
│   └── sistemas.tsv
├── 03_variacion/
│   └── variacion.tsv
└── por_muestra/
    └── NOMBRE_MUESTRA/
        ├── hallazgos.tsv
        ├── revision.tsv
        ├── sistemas.tsv
        └── secuencias/
```

| Archivo | Qué contiene | Cómo interpretarlo |
|---|---|---|
| `muestras.tsv` | Identificador de muestra y archivo de origen | Permite rastrear qué ensamblaje se analizó |
| `resumen.tsv` | Conteo de marcadores distintos por muestra y duración | Un conteo bajo puede reflejar ausencia, divergencia o límites del catálogo |
| `hallazgos.tsv` | Asignaciones que cumplen el perfil fuerte de identidad y cobertura | Cada fila es una coincidencia de secuencia con una referencia, no evidencia de expresión |
| `revision.tsv` | Coincidencias candidatas que no cumplen el perfil fuerte o tienen competencia | Requieren revisar secuencia, longitud, contexto y referencia |
| `matriz.tsv` | Filas de muestras por columnas de genes, con cantidad de copias | `0` indica que el método no asignó una copia, no prueba ausencia biológica |
| `sistemas.tsv` | Genes detectados y genes esperados por módulo | El porcentaje es cobertura del catálogo del módulo, no nivel de actividad |
| `variacion.tsv` | Comparación por pares de la secuencia de ADN de cada marcador detectado | “No evaluable” suele indicar que falta el marcador en una de las muestras |
| `secuencias/*.fna` | Secuencias de ADN de los ORF asignados, separadas por módulo | `.fna` es ADN en FASTA; se conserva el ORF detectado y sus coordenadas en la tabla |

### Columnas principales de `hallazgos.tsv`

`muestra`, `gen` y `accesion` identifican la muestra y la proteína de referencia. `contig`, `inicio`, `fin` y `hebra` ubican el ORF. `orf` enlaza la fila con los archivos `orf.fna` y `orf.faa` de ORF-Scout. `inicio_elegido` indica si se usó el inicio principal o uno alternativo. `sistemas` lista los módulos asociados al gen. `identidad_aa_pct` es la proporción de aminoácidos iguales en la región comparada. `cobertura_ref_aa_pct` y `cobertura_orf_aa_pct` indican cuánto de cada proteína fue comparado. `identidad_adn_ref_pct` compara los codones alineados y puede diferir de la identidad proteica por la degeneración del código genético. `asignacion` distingue coincidencia idéntica, probable o probable con cobertura parcial.

El perfil fuerte actual requiere al menos 90 % de identidad de aminoácidos y 80 % de cobertura de la referencia y del ORF. Se conserva un nivel exploratorio desde 70 % de identidad y 30 % de ambas coberturas en `revision.tsv`. Son umbrales iniciales del proyecto y deben evaluarse con cepas de varios grupos filogenéticos. Un cero significa que la secuencia no superó el criterio con la referencia local.

## Catálogo, alcance y fuentes

Cada fila de `catalogo/marcadores.tsv` contiene el símbolo del gen, uno o más módulos, acceso RefSeq de proteína, secuencia de aminoácidos, evidencia bibliográfica y CDS de referencia de MG1655. `catalogo/referencias.tsv` describe la base y el alcance de cada fuente. El catálogo usa una secuencia de referencia por marcador, así que las cepas que tengan proteínas muy divergentes o genes accesorios ausentes pueden producir coincidencias parciales o no asignadas.

Las respuestas de carbono, nitrógeno y fosfato incluyen sistemas de búsqueda específicos; si la búsqueda no resuelve la limitación, se induce una respuesta de transición a fase estacionaria y supervivencia. `rpoS` se trata como respuesta general de *E. coli*. El módulo de hierro representa privación de un micronutriente y se informa por separado de la respuesta estricta general. Las respuestas a nutrientes se interconectan, así que un regulador puede aparecer en más de un módulo.

Fuentes de apoyo:

1. Peterson CN, Mandel MJ, Silhavy TJ. *Escherichia coli* starvation diets: essential nutrients weigh in distinctly. J Bacteriol. 2005. PMID: 16267278. https://pmc.ncbi.nlm.nih.gov/articles/PMC1280323/
2. Shimizu K. Regulation Systems of Bacteria such as *Escherichia coli* in Response to Nutrient Limitation and Environmental Stresses. Metabolites. 2014. PMID: 24958385. https://pmc.ncbi.nlm.nih.gov/articles/PMC4018673/
3. Irving SE, Choudhury NR, Corrigan RM. The stringent response and physiological roles of (pp)pGpp in bacteria. Nat Rev Microbiol. 2021. doi:10.1038/s41579-020-00470-y.
4. Moreau PL. Regulation of phosphate starvation-specific responses in *Escherichia coli*. Microbiology. 2023. doi:10.1099/mic.0.001312. https://pmc.ncbi.nlm.nih.gov/articles/PMC10191379/
5. RpoS and the bacterial general stress response. Microbiol Mol Biol Rev. 2024. doi:10.1128/mmbr.00151-22. https://pmc.ncbi.nlm.nih.gov/articles/PMC10966952/

La procedencia individual de cada marcador aparece en la columna `evidencia`. La secuencia de referencia procede de CDS de MG1655 RefSeq GCF_000005845.2, ensamblaje GCA/RefSeq de NCBI. La presencia genómica no demuestra transcripción, traducción, actividad ni supervivencia bajo una condición concreta.

## Prueba exploratoria entre cepas

Se procesaron MG1655, W3110, REL606 y CFT073 con los mismos parámetros, sin conexión de red durante la anotación. Starvation-Scout asignó 79/79 marcadores a MG1655, 79/79 a W3110, 78/79 a REL606 y 77/79 a CFT073. En REL606 y CFT073 quedaron candidatos parciales de `cstA` o `rpoS` en `revision.tsv`; el resultado requiere inspección y no se interpreta automáticamente como pérdida del gen. Esta prueba muestra ejecución en cuatro ensamblajes, pero no valida todas las cepas ni todos los alelos de *E. coli*.
