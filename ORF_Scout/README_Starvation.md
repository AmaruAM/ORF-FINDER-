# Starvation-Scout

Starvation-Scout busca coincidencias entre ORF predichos y 79 proteínas de referencia asociadas con limitación nutricional en *Escherichia coli*. Las referencias provienen de MG1655 y están organizadas en nueve sistemas: captación de carbono, reservas de carbono, respuesta estricta, supervivencia general/fase estacionaria, nitrógeno, fosfato, azufre, hierro y limitación de aminoácidos.

El catálogo es un panel compacto de reguladores, componentes de captación y efectores de supervivencia. No representa todos los genes de respuesta ni mide expresión. Una coincidencia sugiere potencial genético y no demuestra que el gen esté activo.

## Ejecución

Desde ORF-Scout, presiona Enter al terminar para dejar la anotación para después. El comando que aparece inicia el análisis. También puedes ejecutarlo así:

```bash
bash starvation_scout.sh --orf corridas/MI_GRUPO/ORF corridas/MI_GRUPO/Starvation 1
```

La opción `1` genera la tabla de coincidencias y el resumen. La opción `2` requiere al menos dos muestras y añade la comparación de variantes de ADN. La carpeta de salida debe ser nueva. El programa trabaja localmente sin conexión ni anotadores externos.

## Archivos de salida

```text
Starvation/
├── muestras.tsv
├── por_cepa/
│   ├── MG1655/coincidencias.tsv
│   └── W3110/coincidencias.tsv
├── resumen.tsv
└── variantes.tsv       solo con la opción 2
```

- `por_cepa/CEPA/coincidencias.tsv` contiene una sola cepa por tabla y ordena sus coincidencias por sistema. Cada fila describe una proteína de referencia encontrada o candidata. Si una proteína pertenece a más de un sistema, aparece una fila por sistema. El identificador `orf` permite recuperar su secuencia en `ORF/CEPA/01_orf/orf.fna` o `orf.faa`.
- `resumen.tsv` reúne en una tabla breve el número de proteínas asignadas y de candidatos que requieren revisión para cada cepa.
- `muestras.tsv` conserva el nombre y la procedencia de cada ensamblaje.
- `variantes.tsv` compara la secuencia de ADN de las proteínas asignadas entre pares de cepas. Se crea solo con la opción 2.

### Columnas de `coincidencias.tsv`

`muestra`, `sistema`, `proteina` y `accesion_referencia` identifican la cepa, el módulo funcional y la proteína comparada. Como cada archivo está dentro de la carpeta de su cepa, todas sus filas pertenecen a esa cepa. `orf`, `contig`, `inicio`, `fin` y `hebra` ubican el ORF. `inicio_usado` indica si se usó el inicio principal, uno alternativo o si la asignación quedó pendiente. `identidad_aa_pct` informa identidad de aminoácidos. `cobertura_referencia_pct` y `cobertura_orf_pct` indican qué proporción de cada proteína se comparó. `identidad_adn_pct` compara los codones alineados. `aa_comparados` muestra el tamaño de la región evaluada. `estado` puede ser `identica`, `probable`, `probable_cobertura_parcial`, `revisar` o `ambiguo`.

El perfil de asignación fuerte exige ≥90 % de identidad de aminoácidos y ≥80 % de cobertura de la referencia y del ORF. Coincidencias desde 70 % de identidad y 30 % de ambas coberturas se muestran como `revisar`. Son umbrales exploratorios para secuencias cercanas a la referencia. Las inserciones o deleciones pueden dificultar la detección. Un cero en el resumen significa que no hubo una asignación fuerte. Conviene revisar también las filas `revisar` y `ambiguo`.

## Catálogo y fuentes

`catalogo/marcadores.tsv` contiene el símbolo del gen, sistema, acceso RefSeq, secuencia proteica, evidencia y CDS de referencia de MG1655. `catalogo/referencias.tsv` describe las fuentes que sustentan los sistemas.

1. Peterson CN, Mandel MJ, Silhavy TJ. *Escherichia coli* starvation diets: essential nutrients weigh in distinctly. J Bacteriol. 2005. PMID: 16267278. https://pmc.ncbi.nlm.nih.gov/articles/PMC1280323/
2. Shimizu K. Regulation Systems of Bacteria such as *Escherichia coli* in Response to Nutrient Limitation and Environmental Stresses. Metabolites. 2014. PMID: 24958385. https://pmc.ncbi.nlm.nih.gov/articles/PMC4018673/
3. Irving SE, Choudhury NR, Corrigan RM. The stringent response and physiological roles of (pp)pGpp in bacteria. Nat Rev Microbiol. 2021. doi:10.1038/s41579-020-00470-y.
4. Moreau PL. Regulation of phosphate starvation-specific responses in *Escherichia coli*. Microbiology. 2023. doi:10.1099/mic.0.001312. https://pmc.ncbi.nlm.nih.gov/articles/PMC10191379/
5. RpoS and the bacterial general stress response. Microbiol Mol Biol Rev. 2024. doi:10.1128/mmbr.00151-22. https://pmc.ncbi.nlm.nih.gov/articles/PMC10966952/

La prueba local previa incluyó MG1655, W3110, REL606 y CFT073. Las coincidencias parciales de `cstA` o `rpoS` se conservaron para revisión, en lugar de contarlas como ausencia confirmada.
