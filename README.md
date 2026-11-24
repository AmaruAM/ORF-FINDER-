# ORF-Scout

Predicción de ORF bacterianos y cribado local de marcadores de respuesta a limitación de nutrientes en *Escherichia coli*.

```bash
wget -O ORF_Scout.zip https://github.com/AmaruAM/ORF-FINDER-/raw/main/ORF_Scout.zip
unzip ORF_Scout.zip
cd ORF_Scout
./orf_scout.sh
```

ORF-Scout crea una carpeta por corrida y por muestra. Al finalizar, Starvation-Scout se inicia solo si se responde `s`. El catálogo local contiene 79 marcadores candidatos organizados en nueve módulos. Consulta `ORF_Scout/README_ORF.md` y `ORF_Scout/README_Starvation.md` para conocer parámetros, interpretación de las tablas y fuentes bibliográficas.
