# ORF-Scout
Ejecuta `bash orf_scout.sh`. Escribe un nombre de corrida y agrega rutas de archivos .fna o carpetas que los contengan. Pulsa Enter para iniciar el análisis. Cada corrida recibe su propio directorio en corridas/, con copias de los genomas, una tabla de procedencia y resultados por muestra. Incluso si repites el nombre, ambas corridas se conservan. Al terminar puedes iniciar la anotación del grupo con Enter. Acepta FASTA de ADN sin comprimir, con encabezados únicos, múltiples contigs y secuencias en varias líneas. Se requieren Bash 4+, awk y utilidades GNU habituales. Las rutas relativas parten de la carpeta donde abres la terminal.

Para reproducir la prueba sin preguntas:

```bash
bash orf_scout.sh ejemplos/MG1655.fna resultados 1 1
```

El primer número selecciona todo (1), ORF (2), vecindad (3) o GC (4). El segundo selecciona ATG/GTG/TTG (1) o añade CTG/ATT/ATC/ATA (2). La tabla 11 presupone TAA/TAG/TGA como STOP. Comprueba el código genético de tu organismo. Los extremos se tratan como lineales y las bases ambiguas cortan la búsqueda de ORF.

1. **ORF**. Por contig y hebra se recorren tres marcos en pasos de tres. Cada registro usa el inicio elegible más lejano al STOP en ese marco y guarda las coordenadas 5′ de los inicios alternativos elegibles. Mínimo 150 pb incluyendo STOP. En `01_orf/`, el FASTA contiene ADN orientado y el TSV da posiciones genómicas inclusivas desde 1.
2. **Posibles operones**. Solo para este análisis se ordenan los candidatos por longitud descendente dentro del contig. Se seleccionan sucesivamente los compatibles con los ya aceptados (solapamiento máximo de 4 pb); en caso de empate se prioriza la menor coordenada. En orden genómico se agrupan vecinos consecutivos de la misma hebra con distancia de −4 a 75 pb. Ver `02_operones/grupos.tsv`. En la hebra negativa se lee la lista de derecha a izquierda para obtener el orden de transcripción. La prioridad por longitud puede desplazar genes cortos reales; los grupos son hipótesis de vecindad.
3. **GC**. En `03_gc/gc.tsv` se reportan el GC de cada candidato y su diferencia respecto al ensamblaje en puntos porcentuales. El cálculo global usa A+C+G+T. Para ORF de ≥300 pb se usan límites Q1−1,5×RIC y Q3+1,5×RIC entre todos los elegibles, con cuantiles interpolados y al menos 20 observaciones. Se informa por separado si faltan datos. La etiqueta indica composición atípica.

Cada corrida crea su propia carpeta y un resumen. En MG1655, modo habitual: 24 056 candidatos, 1 021 grupos, 373 candidatos con GC atípico. El grupo 494 corresponde a los ORF de mglC, mglA y mglB en orden de coordenadas. El script 1 asigna identificadores de ORF; los nombres de genes proceden del contraste posterior.

## Fuentes
1. National Center for Biotechnology Information. ORFfinder. https://www.ncbi.nlm.nih.gov/orffinder/
2. National Center for Biotechnology Information. The genetic codes. https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
3. Moreno-Hagelsieb G, Collado-Vides J. A powerful non-homology method for the prediction of operons in prokaryotes. Bioinformatics. 2002;18 Suppl 1:S329–36.
4. García-Vallvé S, Romeu A, Palau J. Horizontal gene transfer in bacterial and archaeal complete genomes. Genome Res. 2000;10:1719–25. Los umbrales de vecindad, RIC y mínimo de observaciones son decisiones del proyecto.
5. Guía de clase, sesiones 6, 10, 13 y 14. https://github.com/1ABL0014/Principios_Programacion_Bioinformatica

## Una carpeta por genoma
En modo de argumentos, ejecuta:
```bash
bash orf_scout.sh --muestras resultados_orf 1 1 ejemplos/*.fna
bash starvation_scout.sh --orf resultados_orf anotaciones 1
```
Cada archivo genera resultados_orf/NOMBRE/ con sus análisis, resumen y origen.txt. Puedes añadir nuevos genomas a la misma carpeta usando nombres distintos. Una muestra existente se protege contra sobrescritura. La extracción acepta cualquier número de entradas sujeto a recursos del sistema. Usa automáticamente los procesadores disponibles con nproc; ORF_JOBS permite ajustar el paralelismo según CPU y memoria.
