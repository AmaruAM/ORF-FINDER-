# ORF-Scout

Predicción de ORF bacterianos y búsqueda local de proteínas asociadas con respuestas a limitación nutricional en *Escherichia coli*.

```bash
wget -O ORF_Scout.zip https://github.com/AmaruAM/ORF-FINDER-/raw/main/ORF_Scout.zip
unzip ORF_Scout.zip
cd ORF_Scout
./orf_scout.sh
```

ORF-Scout crea una carpeta por corrida y por muestra. La anotación se inicia solo si se responde `s`. Starvation-Scout guarda una tabla de coincidencias por cepa, agrupada por sistema, y un resumen conjunto. La comparación de variantes es opcional. Consulta `ORF_Scout/README_ORF.md` y `ORF_Scout/README_Starvation.md` para conocer parámetros, interpretación de las tablas y fuentes bibliográficas.
