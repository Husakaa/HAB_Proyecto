# Análisis Integrado: Propagación en Redes Biológicas y Enriquecimiento Funcional

**Autores:** Hugo Salas Calderón, Aissa Omar El Hammouti Chachoui y Yussef Barakat Nieto  
**Institución:** Universidad de Málaga - Ingeniería de la Salud  
**Asignatura:** Bioinformática Aplicada  
**Fecha:** Noviembre 2025

---

## Introducción

La identificación de genes asociados a enfermedades complejas es uno de los desafíos fundamentales en bioinformática moderna. Los estudios de expresión diferencial proporcionan listas de genes candidatos, pero estos análisis ignoran las interacciones moleculares que definen los procesos biológicos subyacentes. La topología de las redes de interacción proteína-proteína contiene información crucial sobre las relaciones funcionales entre genes, permitiendo identificar candidatos adicionales que no muestran cambios de expresión significativos pero que participan en los mismos módulos funcionales.

Este proyecto implementa un pipeline computacional completo que integra dos aproximaciones complementarias para el análisis de redes biológicas: propagación de información mediante algoritmos topológicos y estadísticos, seguido de caracterización funcional mediante análisis de enriquecimiento. El objetivo es proporcionar una herramienta automatizada, reproducible y robusta para la identificación y caracterización de módulos de enfermedad a partir de datos experimentales.

---

## Descripción del Proyecto

### Objetivo General

Desarrollar un sistema automatizado de análisis bioinformático que, a partir de una lista inicial de genes semilla (típicamente genes diferencialmente expresados), identifique genes candidatos adicionales mediante propagación en redes de interacción proteína-proteína y caracterice su función biológica mediante análisis de enriquecimiento funcional.

### Flujo de Trabajo

El pipeline implementado sigue una arquitectura modular que comprende las siguientes etapas:

**1. Validación y Normalización de Identificadores**

La primera etapa del análisis consiste en la validación de los genes de entrada. Los identificadores de genes pueden provenir de diferentes sistemas (símbolos HGNC, Entrez IDs, Ensembl, UniProt) y frecuentemente contienen errores tipográficos o están desactualizados. El sistema consulta automáticamente la base de datos MyGene.info para validar cada identificador, convertirlo a símbolos oficiales HGNC y obtener su correspondiente Entrez ID. Esta normalización es crítica para garantizar la compatibilidad entre las diferentes bases de datos utilizadas en etapas posteriores.

El sistema implementa además una búsqueda inteligente de variantes mitocondriales, intentando automáticamente agregar el prefijo MT- cuando un gen no se encuentra en la base de datos principal. Esto resuelve inconsistencias comunes en la nomenclatura de genes mitocondriales.

**2. Carga y Procesamiento de Redes PPI**

El sistema acepta redes de interacción proteína-proteína en múltiples formatos. Mediante análisis de la estructura del archivo, detecta automáticamente si se trata de un formato STRING (TSV con encabezados), DIAMOnD (CSV sin encabezados) o GUILD (formato espaciado con pesos). Esta detección automática elimina la necesidad de configuración manual y reduce errores de formato.

Las redes se representan internamente como grafos no dirigidos utilizando la biblioteca NetworkX, permitiendo operaciones eficientes sobre estructuras de gran tamaño. El sistema valida que los genes semilla estén presentes en la red antes de proceder con los análisis de propagación.

**3. Algoritmos de Propagación en Redes**

Se implementan dos algoritmos complementarios que explotan diferentes propiedades de la topología de red:

**GUILD (Random Walk with Restart)**

GUILD implementa un algoritmo de difusión global basado en caminatas aleatorias con reinicio. El proceso simula un "caminante" que se mueve por la red desde los genes semilla, con una probabilidad constante de "reiniciar" y volver a cualquier gen semilla. Este mecanismo permite propagar información a través de toda la red, asignando puntuaciones más altas a genes bien conectados con el módulo semilla.

La implementación utiliza el algoritmo PageRank de NetworkX, que resuelve el sistema iterativamente hasta convergencia. El parámetro de reinicio controla el equilibrio entre exploración local (valores altos favorecen vecindad inmediata) y global (valores bajos permiten alcanzar nodos distantes). Matemáticamente, el algoritmo converge a la distribución estacionaria de una cadena de Markov personalizada.

**DIAMOnD (Disease Module Detection)**

DIAMOnD adopta una estrategia de expansión iterativa basada en significancia estadística. En cada iteración, el algoritmo evalúa todos los nodos adyacentes al módulo actual y calcula la probabilidad de observar su conectividad bajo un modelo hipergeométrico nulo. El nodo con el p-valor más significativo se añade al módulo, y el proceso se repite hasta alcanzar el tamaño deseado.

La implementación optimizada utiliza una cola de prioridad (min-heap) para mantener los candidatos ordenados por p-valor, con actualizaciones incrementales que evitan recalcular valores innecesariamente. Esto reduce la complejidad computacional de O(K×N²) a O(K×N×log N), donde K es el número de genes a añadir y N el tamaño de la red.

La principal ventaja de DIAMOnD es su interpretación estadística clara: cada gen añadido tiene una conexión con el módulo que no puede explicarse por azar, considerando el grado del nodo y el tamaño del módulo.

**4. Análisis de Enriquecimiento Funcional**

Una vez identificados los genes candidatos mediante propagación, el sistema caracteriza su función biológica consultando tres bases de datos ortogonales:

- Gene Ontology Biological Process (GO:BP): Procesos biológicos genéricos como metabolismo, señalización celular, apoptosis
- KEGG (Kyoto Encyclopedia of Genes and Genomes): Vías metabólicas y rutas de señalización específicas
- Reactome: Vías de reacción biológica con mayor granularidad que KEGG

El análisis se realiza mediante g:Profiler, que ejecuta tests de sobrerepresentación hipergeométricos para identificar términos enriquecidos significativamente. El sistema aplica corrección de múltiples pruebas mediante False Discovery Rate (FDR) para controlar errores tipo I. Solo se reportan términos con p-valor ajustado inferior al umbral especificado (por defecto 0.05).

La ventaja de g:Profiler sobre otras herramientas es su integración directa con Ensembl, garantizando anotaciones actualizadas automáticamente con cada release de la base de datos.

**5. Visualización y Exportación**

El sistema genera visualizaciones de alta calidad (300 DPI) adecuadas para publicación científica:

- Comparación de resultados de propagación entre GUILD y DIAMOnD
- Ranking de términos enriquecidos por significancia estadística
- Distribución de términos por base de datos de origen

Todos los resultados se exportan en formatos estructurados (CSV, TSV) con metadata completa en JSON, garantizando la reproducibilidad del análisis. Se genera además un resumen ejecutivo en texto plano con las principales conclusiones.

---

## Métodos Implementados

### Validación de Identificadores de Genes

**Herramienta:** MyGene.info API  
**Proceso:** Consulta automática de identificadores contra base de datos integrada (NCBI, Ensembl, UniProt), normalización a símbolos HGNC oficiales y conversión a Entrez IDs para compatibilidad con herramientas downstream.

### Algoritmo GUILD

**Referencia:** Guney et al., Nature Communications 7:10331 (2016)  
**Fundamento:** Random Walk with Restart sobre grafo de interacciones. Implementación mediante PageRank con vector de personalización sobre genes semilla.  
**Parámetros configurables:** Probabilidad de reinicio (r), máximo de iteraciones, criterio de convergencia.

### Algoritmo DIAMOnD

**Referencia:** Ghiassian et al., PLOS Computational Biology 11(4):e1004120 (2015)  
**Fundamento:** Expansión iterativa mediante test hipergeométrico. En cada paso se añade el nodo con conectividad más significativa al módulo actual.  
**Optimización:** Cola de prioridad con actualizaciones incrementales para eficiencia en redes grandes.

### Análisis de Enriquecimiento Funcional

**Herramienta:** g:Profiler (Raudvere et al., Nucleic Acids Research 2019)  
**Bases de datos consultadas:** GO:BP, KEGG, Reactome  
**Método estadístico:** Test de sobrerepresentación hipergeométrico con corrección FDR  
**Umbral de significancia:** p-valor ajustado < 0.05 (configurable)

---

## Estructura del Proyecto
```
│
├── data/                          # Datos de entrada
│   ├── genes_seed.txt             # Genes semilla (usuario)
│   └── string_network_filtered_hugo-400.tsv  # Red PPI
│
├── results/                       # Resultados generados automáticamente
│   ├── metadata.json              # Parámetros de ejecución
│   ├── 01_validacion_genes.csv    # Genes validados con Entrez IDs
│   ├── 02_guild_propagation.tsv   # Resultados GUILD
│   ├── 03_diamond_propagation.tsv # Resultados DIAMOnD
│   ├── 04_enrichment_complete.csv # Enriquecimiento completo
│   ├── 05_enrichment_GO_BP.csv    # Procesos biológicos
│   ├── 06_enrichment_KEGG.csv     # Vías metabólicas
│   ├── 07_enrichment_REAC.csv     # Vías Reactome
│   ├── 08_resumen_ejecucion.txt   # Resumen ejecutivo
│   └── *.png                      # Visualizaciones
│
├── scripts/
│   └── network_functional_analysis.py  # Script principal
│
├── setup.sh                       # Script de instalación (Linux/macOS)
├── setup.bat                      # Script de instalación (Windows CMD)
├── requirements.txt               # Dependencias Python
└── README.md                      # Este archivo
```

---

## Instalación

### Requisitos Previos

- Python 3.8 o superior
- Gestor de paquetes pip

### Configuración del Entorno Virtual

#### En Linux/macOS:

Ejecute el script `setup.sh`, o manualmente ejecute los siguientes comandos:
```bash
# Crear entorno virtual
python3 -m venv .venv

# Activar entorno virtual
source .venv/bin/activate

# Instalar dependencias
pip install -r requirements.txt
```

#### En Windows:

Ejecute el script `setup.bat`, o manualmente ejecute los siguientes comandos:
```batch
# Crear entorno virtual
python -m venv .venv

# Activar entorno virtual
.venv\Scripts\activate

# Instalar dependencias
pip install -r requirements.txt
```

### Desactivación del Entorno Virtual

Cuando haya terminado de trabajar en el proyecto:
```bash
deactivate
```
## Uso del Sistema

### Sintaxis General
```bash
python scripts/network_functional_analysis.py -s <archivo_genes> -n <archivo_red> [opciones]
```

### Argumentos Obligatorios

- `-s, --seed-file`: Archivo de texto con genes semilla (un gen por línea o separados por comas)
- `-g, --genes`: Alternativamente, lista de genes directamente desde línea de comandos
- `-n, --network`: Archivo de red de interacción proteína-proteína

### Opciones Principales

- `-a, --algorithm`: Algoritmo a ejecutar: `guild`, `diamond`, o `both` (por defecto: both)
- `--top-genes N`: Número de genes candidatos a usar para enriquecimiento (por defecto: 100)
- `--guild-restart R`: Probabilidad de reinicio en GUILD (por defecto: 0.5)
- `--diamond-k K`: Número de genes a añadir en DIAMOnD (por defecto: 200)
- `--enrichment-threshold P`: Umbral de p-valor para enriquecimiento (por defecto: 0.05)
- `--organism ORG`: Organismo objetivo (por defecto: hsapiens)
- `-o, --output DIR`: Directorio de salida (por defecto: results)
- `--no-plots`: Omitir generación de visualizaciones

### Formatos de Entrada Soportados

**Archivo de genes:**
```
ENO1
PGK1
HK2
```

o
```
ENO1, PGK1, HK2
```

**Red PPI:**
- STRING: TSV con encabezados protein1, protein2
- DIAMOnD: CSV sin encabezados (nodo1,nodo2)
- GUILD: Formato espaciado (nodo1 peso nodo2)

### Ejemplos de Uso

**Análisis completo con ambos algoritmos:**
```bash
python scripts/network_functional_analysis.py \
    -s data/genes_seed.txt \
    -n data/string_network_filtered_hugo-400.tsv \
    -o results/analisis_completo
```

**Solo GUILD con parámetros personalizados:**
```bash
python scripts/network_functional_analysis.py \
    -s data/genes_seed.txt \
    -n data/string_network_filtered_hugo-400.tsv \
    -a guild \
    --guild-restart 0.7 \
    --top-genes 50
```

**Genes directos desde línea de comandos:**
```bash
python scripts/network_functional_analysis.py \
    -g ENO1 PGK1 HK2 GAPDH \
    -n data/string_network_filtered_hugo-400.tsv
```

**Análisis rápido sin visualizaciones:**
```bash
python scripts/network_functional_analysis.py \
    -s data/genes_seed.txt \
    -n data/string_network_filtered_hugo-400.tsv \
    --no-plots
```

---

## Interpretación de Resultados

### Archivo 01_validacion_genes.csv

Contiene el mapeo entre identificadores de entrada y símbolos validados, junto con los Entrez IDs correspondientes. Los genes marcados como "no_encontrado" no pudieron validarse en MyGene.info.

### Archivos 02_guild_propagation.tsv y 03_diamond_propagation.tsv

Rankings de genes ordenados por relevancia según cada algoritmo. GUILD proporciona scores de propagación (valores más altos indican mayor relevancia), mientras que DIAMOnD proporciona p-valores (valores más bajos indican mayor significancia estadística).

### Archivos 04-07: Enriquecimiento funcional

Términos significativamente enriquecidos en la lista de genes candidatos. Cada fila representa un proceso biológico o vía metabólica, con su p-valor ajustado, tamaño del término, y genes de intersección.

**Columnas principales:**
- source: Base de datos de origen (GO:BP, KEGG, REAC)
- name: Descripción del término
- p_value: p-valor ajustado (FDR)
- intersection_size: Número de genes en común
- term_size: Tamaño total del término
- query_size: Número de genes en la consulta

### Archivo 08_resumen_ejecucion.txt

Resumen ejecutivo en texto plano con los principales hallazgos del análisis, incluyendo estadísticas de propagación y los procesos biológicos más enriquecidos.

---

## Justificación de Decisiones Técnicas

### Elección de Algoritmos

**GUILD vs DIAMOnD:**

Ambos algoritmos ofrecen aproximaciones complementarias al problema de identificación de módulos de enfermedad. GUILD explota la estructura global de la red mediante difusión, siendo capaz de identificar genes alejados topológicamente pero funcionalmente relevantes. DIAMOnD adopta una estrategia más conservadora, identificando módulos localmente cohesivos con interpretación estadística rigurosa.

La implementación de ambos permite validación cruzada: genes que aparecen en ambos rankings tienen mayor probabilidad de ser verdaderos positivos.

### Elección de g:Profiler

Entre las herramientas disponibles para enriquecimiento funcional (DAVID, Enrichr, g:Profiler), se seleccionó g:Profiler por:

1. Integración directa con Ensembl, garantizando anotaciones actualizadas
2. API programática estable y bien documentada
3. Soporte nativo para múltiples organismos
4. Corrección de múltiples pruebas mediante FDR
5. Mantenimiento activo y desarrollo continuo

### Automatización Completa

El diseño del sistema prioriza la automatización para reducir errores humanos:

- Detección automática de formatos de entrada
- Conversión automática entre sistemas de identificadores
- Validación de datos en cada etapa
- Generación automática de metadata para reproducibilidad

---

## Limitaciones y Trabajo Futuro

**Limitaciones actuales:**

1. La propagación asume que la red PPI es completa y libre de errores, cuando en realidad contiene falsos positivos y negativos
2. El análisis no considera pesos diferenciales en las interacciones (ej. fuerza de la interacción, contexto celular)
3. La selección de genes candidatos para enriquecimiento es arbitraria (top N)

**Mejoras propuestas:**

1. Integración de redes específicas de tejido o tipo celular
2. Incorporación de información de expresión diferencial como pesos en la propagación
3. Implementación de análisis de estabilidad mediante bootstrap
4. Extensión a otras especies modelo (ratón, rata, pez cebra)

---

## Referencias

1. Guney, E., Menche, J., Vidal, M., & Barabasi, A. L. (2016). Network-based in silico drug efficacy screening. Nature Communications, 7, 10331.

2. Ghiassian, S. D., Menche, J., & Barabasi, A. L. (2015). A DIseAse MOdule Detection (DIAMOnD) Algorithm Derived from a Systematic Analysis of Connectivity Patterns of Disease Proteins in the Human Interactome. PLOS Computational Biology, 11(4), e1004120.

3. Raudvere, U., Kolberg, L., Kuzmin, I., et al. (2019). g:Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update). Nucleic Acids Research, 47(W1), W191-W198.

4. Xin, J., Mark, A., Afrasiabi, C., et al. (2016). High-performance web services for querying gene and variant annotation. Genome Biology, 17, 91.

5. Szklarczyk, D., Gable, A. L., Nastou, K. C., et al. (2023). The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest. Nucleic Acids Research, 51(D1), D638-D646.

---

## Licencia

Este proyecto se distribuye bajo licencia MIT. Ver archivo LICENSE para detalles completos.

---

**Universidad de Málaga - Ingeniería de la Salud**  
**Noviembre 2025**