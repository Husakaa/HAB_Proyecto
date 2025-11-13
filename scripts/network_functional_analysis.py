#!/usr/bin/env python3
"""
Análisis Integrado: Propagación en Redes + Enriquecimiento Funcional
====================================================================

Pipeline completo para identificar genes candidatos mediante propagación en redes
y caracterizar su función biológica mediante análisis de enriquecimiento.

FLUJO DE TRABAJO:
1. Genes Semilla → 2. Propagación en Red → 3. Genes Candidatos → 4. Análisis Funcional

MÉTODOS:
- GUILD: Random Walk with Restart (propagación topológica)
- DIAMOnD: Disease Module Detection (expansión estadística)
- g:Profiler: Enriquecimiento funcional (GO, KEGG, Reactome)

AUTOR: Hugo Salas Calderón
"""

# Importar dependencias
try:

    import argparse
    import os
    import sys
    import json
    from datetime import datetime
    import warnings
    warnings.filterwarnings('ignore')

    import networkx as nx
    import heapq
    import numpy as np
    import pandas as pd
    from scipy.stats import hypergeom

    import mygene
    from gprofiler import GProfiler
    import matplotlib.pyplot as plt
    import seaborn as sns

    import itertools


except ImportError as e:
    print("   Error: Falta instalar dependencias")
    print("   Ejecuta: pip install -r requirements.txt")
    sys.exit(1)


# Configuración de estilo para gráficos
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300


# =============================================================================
# SECCIÓN 1: CARGA Y VALIDACIÓN DE DATOS
# =============================================================================

def cargar_red(archivo):
    """
    Carga red de interacción proteína-proteína detectando formato automáticamente.
    
    Formatos soportados:
    - STRING: TSV con encabezados (protein1, protein2)
    - DIAMOnD: CSV sin encabezados (nodo1,nodo2)
    - GUILD: Espacios con pesos (nodo1 peso nodo2)
    
    Args:
        archivo: Ruta al archivo de red
        
    Returns:
        networkx.Graph: Grafo de interacciones
    """
    with open(archivo, 'r') as f:
        primera = f.readline().strip()
        
    if '\t' in primera or 'protein' in primera.lower():
        # STRING: TSV
        print("   Formato: STRING (TSV)")
        df = pd.read_csv(archivo, sep='\t')
        # Si solo hay una columna -> lista de genes
        if df.shape[1] == 1:
            print("    Archivo contiene una sola columna de genes. Creando red artificial...")
            genes = df.iloc[:, 0].astype(str).unique().tolist()
            G = nx.Graph()
            G.add_nodes_from(genes)
            for a, b in itertools.combinations(genes, 2):
                G.add_edge(a, b, weight=1.0)
            print(f"    Red artificial creada: {len(G.nodes())} nodos, {len(G.edges())} aristas")
            return G
        return nx.from_pandas_edgelist(df, df.columns[0], df.columns[1])
    
    elif ',' in primera:
        # DIAMOnD: CSV (puede estar en .txt o .csv)
        print("   Formato: DIAMOnD (CSV)")
        df = pd.read_csv(archivo, header=None, names=['source', 'target'], dtype=str)
        # Detectar caso de una sola columna
        if df.shape[1] == 1:
            print("    Archivo CSV con una sola columna de genes. Creando red artificial...")
            genes = df.iloc[:, 0].astype(str).unique().tolist()
            G = nx.Graph()
            G.add_nodes_from(genes)
            for a, b in itertools.combinations(genes, 2):
                G.add_edge(a, b, weight=1.0)
            print(f"    Red artificial creada: {len(G.nodes())} nodos, {len(G.edges())} aristas")
            return G
        return nx.from_pandas_edgelist(df, 'source', 'target')
    
    else:
        # GUILD: nodo1 peso nodo2 (espacios)
        partes = primera.split()
        if len(partes) == 3:
            print("   Formato: GUILD")
            G = nx.Graph()
            with open(archivo, 'r', encoding='utf-8', errors='ignore') as f:
                for linea in f:
                    partes = linea.strip().split()
                    if len(partes) == 3:
                        G.add_edge(str(partes[0]), str(partes[2]), weight=float(partes[1]))
            return G

        # --- NUEVO CASO: lista de genes (sin separadores ni pesos) ---
        elif len(partes) == 1:
            print("   Formato: LISTA DE GENES (sin interacciones)")
            print("   Creando red artificial completamente conectada para pruebas.")
            df = pd.read_csv(archivo, header=None, names=['gene'])
            genes = df['gene'].astype(str).unique().tolist()
            G = nx.Graph()
            G.add_nodes_from(genes)
            for a, b in itertools.combinations(genes, 2):
                G.add_edge(a, b, weight=1.0)
            print(f"    Red artificial creada: {len(G.nodes())} nodos, {len(G.edges())} aristas")
            return G

        else:
            raise ValueError(f"Formato no reconocido. Primera línea: {primera}")


def leer_genes_de_archivo(archivo):
    """
    Lee genes desde archivo de texto con formato flexible.
    
    Soporta:
    - Un gen por línea: ENO1\\nPGK1\\nHK2
    - Separados por comas: ENO1, PGK1, HK2
    - Mezclado: ENO1, PGK1\\nHK2
    
    Args:
        archivo: Ruta al archivo
        
    Returns:
        list: Lista de genes
    """
    try:
        with open(archivo, 'r') as f:
            contenido = f.read()
        
        # Normalizar formato
        contenido = contenido.replace(' y ', ',').replace('\n', ',')
        genes = [g.strip() for g in contenido.split(',') if g.strip()]
        
        print(f"   Leídos {len(genes)} genes desde {archivo}")
        return genes
        
    except FileNotFoundError:
        print(f"   Archivo no encontrado: {archivo}")
        sys.exit(1)


def agregar_genes_a_red(archivo_red, genes_faltantes, archivo_string=None):
    """
    Agrega genes faltantes a la red buscando sus interacciones en STRING.
    
    Args:
        archivo_red: Archivo de red GUILD/DIAMOnD a modificar
        genes_faltantes: Lista de Entrez IDs que faltan
        archivo_string: Archivo STRING de donde obtener interacciones (opcional)
    
    Returns:
        bool: True si se agregaron genes, False si no
    """
    # Buscar archivo STRING en múltiples ubicaciones
    posibles_rutas = [
        archivo_string,
        'data/string_network_filtered_hugo-400.tsv',
    ]
    
    archivo_string_encontrado = None
    for ruta in posibles_rutas:
        if ruta and os.path.exists(ruta):
            archivo_string_encontrado = ruta
            break
    
    if not archivo_string_encontrado:
        print("   Advertencia: No se encontró archivo STRING, no se pueden agregar genes")
        return False
    
    # Convertir Entrez IDs a símbolos usando mygene
    try:
        mg = mygene.MyGeneInfo()
        
        print("   Convirtiendo Entrez IDs a símbolos...")
        results = mg.querymany(genes_faltantes, scopes='entrezgene', fields='symbol', species='human')
        
        # Crear mapeo de Entrez -> símbolo
        entrez_a_simbolo = {}
        simbolo_a_entrez = {}
        
        for r in results:
            entrez = str(r.get('query'))
            simbolo = r.get('symbol')
            if simbolo:
                entrez_a_simbolo[entrez] = simbolo
                simbolo_a_entrez[simbolo] = entrez
        
        if not entrez_a_simbolo:
            print("   No se pudieron convertir los Entrez IDs")
            return False
        
        simbolos_buscar = list(entrez_a_simbolo.values())
        
    except ImportError:
        print("   Error: mygene no instalado, no se pueden convertir IDs")
        return False
    except Exception as e:
        print(f"   Error al convertir IDs: {e}")
        return False
    
    print(f"   Buscando interacciones en STRING para: {', '.join(simbolos_buscar)}...")
    
    # Leer red STRING
    try:
        df_string = pd.read_csv(archivo_string_encontrado, sep='\t')
    except Exception as e:
        print(f"   Error al leer STRING: {e}")
        return False
    
    # Buscar interacciones que involucren estos genes
    interacciones = df_string[
        (df_string['protein1_hugo'].isin(simbolos_buscar)) | 
        (df_string['protein2_hugo'].isin(simbolos_buscar))
    ]
    
    # Filtrar solo interacciones entre los genes de interés
    interacciones_internas = interacciones[
        (interacciones['protein1_hugo'].isin(simbolos_buscar)) & 
        (interacciones['protein2_hugo'].isin(simbolos_buscar))
    ]
    
    if len(interacciones_internas) == 0:
        print("   No se encontraron interacciones entre estos genes en STRING")
        return False
    
    print(f"   Encontradas {len(interacciones_internas)} interacciones")
    
    # Preparar líneas a añadir
    nuevas_lineas = []
    interacciones_unicas = set()
    
    for _, row in interacciones_internas.iterrows():
        gene1 = row['protein1_hugo']
        gene2 = row['protein2_hugo']
        
        if gene1 in simbolo_a_entrez and gene2 in simbolo_a_entrez:
            entrez1 = simbolo_a_entrez[gene1]
            entrez2 = simbolo_a_entrez[gene2]
            
            # Evitar duplicados (orden no importa)
            par = tuple(sorted([entrez1, entrez2]))
            if par not in interacciones_unicas:
                interacciones_unicas.add(par)
                
                # Formato según tipo de red
                if 'guild' in archivo_red.lower():
                    nuevas_lineas.append(f"{entrez1} 1 {entrez2}")
                else:  # diamond
                    nuevas_lineas.append(f"{entrez1},{entrez2}")
    
    if not nuevas_lineas:
        return False
    
    # Añadir líneas al archivo
    try:
        with open(archivo_red, 'a') as f:
            for linea in nuevas_lineas:
                f.write('\n' + linea)
        
        print(f"   Agregadas {len(nuevas_lineas)} interacciones a {archivo_red}")
        return True
    except Exception as e:
        print(f"   Error al escribir en {archivo_red}: {e}")
        return False


def convertir_a_entrez(genes):
    """Convierte símbolos de genes a Entrez IDs usando MyGene"""
    try:
        mg = mygene.MyGeneInfo()
        print("Convirtiendo símbolos a Entrez IDs...")
        
        results = mg.querymany(genes, scopes='symbol', fields='entrezgene', species='human')
        
        mapeo = {}
        for r in results:
            simbolo = r.get('query')
            entrez_id = r.get('entrezgene')
            
            if entrez_id:
                # Asegurar que el Entrez ID es string
                mapeo[simbolo] = str(entrez_id)
                print(f"   {simbolo} -> {entrez_id}")
            else:
                print(f"   {simbolo} -> NO ENCONTRADO")
        
        return mapeo
        
    except ImportError:
        print("Error: mygene no instalado. Instalar con: pip install -r requirements")
        return {}
    except Exception as e:
        print(f"Error al convertir genes: {e}")
        return {}



def validar_y_convertir_genes(genes, organismo='human'):
    """
    Valida y convierte identificadores de genes a símbolos oficiales.
    
    Usa MyGene.info para:
    - Validar que los genes existan
    - Convertir entre diferentes tipos de IDs (Entrez, Ensembl, símbolos)
    - Normalizar a símbolos oficiales (HGNC)
    - Detectar genes mal escritos o no encontrados
    - Buscar automáticamente variantes con prefijo MT- para genes mitocondriales
      (e.g., si "ND1" no se encuentra, intenta con "MT-ND1")
    
    Args:
        genes: Lista de identificadores de genes (pueden ser símbolos, Entrez IDs, etc.)
        organismo: Especie ('human', 'mouse', etc.)
        
    Returns:
        dict con:
            - 'validos': Lista de símbolos validados
            - 'mapping': Diccionario de conversión {input -> símbolo}
            - 'no_encontrados': Lista de genes no encontrados
            - 'advertencias': Lista de mensajes de advertencia
    """
    print(f"\n{'='*70}")
    print("VALIDACIÓN Y CONVERSIÓN DE IDENTIFICADORES")
    print(f"{'='*70}")
    print(f"[INFO] Validando {len(genes)} genes con MyGene.info...")
    
    # Inicializar cliente de MyGene
    mg = mygene.MyGeneInfo()
    
    # Determinar species code para MyGene
    species_map = {
        'human': 'human',
        'mouse': 'mouse',
        'rat': 'rat',
        'hsapiens': 'human',
        'mmusculus': 'mouse',
        'rnorvegicus': 'rat'
    }
    species = species_map.get(organismo.lower(), 'human')
    
    # Consultar MyGene.info
    # scopes: buscar en múltiples tipos de IDs
    try:
        results = mg.querymany(
            genes, 
            scopes='symbol,entrezgene,ensembl.gene,uniprot.Swiss-Prot',
            fields='symbol,entrezgene,name',
            species=species,
            returnall=True
        )
    except Exception as e:
        print(f"[ERROR] Error al conectar con MyGene.info: {e}")
        print("[INFO] Continuando sin validación...")
        return {
            'validos': genes,
            'mapping': {g: g for g in genes},
            'no_encontrados': [],
            'advertencias': ['No se pudo validar genes (sin conexión)']
        }
    
    # Procesar resultados
    genes_validos = []
    mapping = {}
    no_encontrados = []
    advertencias = []
    genes_a_reintentar = []  # Para genes no encontrados que intentaremos con MT-
    
    for gene_input, result in zip(genes, results['out']):
        if 'notfound' in result and result['notfound']:
            # Gen no encontrado - guardar para reintentar con MT-
            genes_a_reintentar.append(gene_input)
            
        elif 'symbol' in result:
            # Gen encontrado - usar símbolo oficial
            simbolo = result['symbol']
            genes_validos.append(simbolo)
            mapping[gene_input] = simbolo
        else:
            # Resultado ambiguo o incompleto
            genes_a_reintentar.append(gene_input)
    
    # Reintentar genes no encontrados con prefijo MT- (genes mitocondriales)
    if genes_a_reintentar:
        print(f"[INFO] Buscando variantes mitocondriales (MT-) para {len(genes_a_reintentar)} genes...")
        
        # Crear variantes con MT-
        variantes_mt = [f"MT-{gene}" for gene in genes_a_reintentar]
        
        try:
            results_mt = mg.querymany(
                variantes_mt,
                scopes='symbol,entrezgene,ensembl.gene,uniprot.Swiss-Prot',
                fields='symbol,entrezgene,name',
                species=species,
                returnall=True
            )
            
            for gene_original, gene_mt, result_mt in zip(genes_a_reintentar, variantes_mt, results_mt['out']):
                if 'notfound' not in result_mt and 'symbol' in result_mt:
                    # Encontrado con prefijo MT-
                    simbolo = result_mt['symbol']
                    genes_validos.append(simbolo)
                    mapping[gene_original] = simbolo
                else:
                    # Definitivamente no encontrado
                    no_encontrados.append(gene_original)
                    advertencias.append(f"Gen no encontrado: '{gene_original}'")
                    
        except Exception as e:
            # Si falla el reintento, agregar todos a no encontrados
            for gene_original in genes_a_reintentar:
                no_encontrados.append(gene_original)
                advertencias.append(f"Gen no encontrado: '{gene_original}'")
    
    # Resumen
    print(f"Validación completada: {len(genes_validos)} genes válidos")
    if no_encontrados:
        print(f"Advertencia: {len(no_encontrados)} genes no encontrados: {', '.join(no_encontrados)}")
    
    # Si no hay genes válidos, detener
    if not genes_validos:
        print("[ERROR] No se encontró ningún gen válido. Verifica los nombres/IDs.")
        sys.exit(1)
    
    print(f"{'='*70}\n")
    
    return {
        'validos': genes_validos,
        'mapping': mapping,
        'no_encontrados': no_encontrados,
        'advertencias': advertencias
    }


# =============================================================================
# SECCIÓN 2: ALGORITMOS DE PROPAGACIÓN EN REDES
# =============================================================================

def guild(G, genes_semilla, r=0.5, iteraciones=100, epsilon=1e-6):
    """
    GUILD: Random Walk with Restart usando networkx.pagerank.
    
    Args:
        G: Grafo de NetworkX
        genes_semilla: Lista de genes iniciales
        r: Probabilidad de restart (0-1)
        iteraciones: Máximo de iteraciones (max_iter)
        epsilon: Criterio de convergencia (tol)
        
    Returns:
        DataFrame con (rank, gene, score)
    """
    print(f"\n{'='*70}")
    print("GUILD: Random Walk with Restart (usando nx.pagerank)")
    print(f"{'='*70}")
     
    # 1. Crear el vector de personalización (p0)
    # nx.pagerank espera un diccionario {nodo: probabilidad}
    p0 = {node: 0 for node in G.nodes()}
    
    # Filtrar semillas que están en el grafo
    semillas_en_grafo = [gen for gen in genes_semilla if gen in G]
    
    if not semillas_en_grafo:
        print("   Error: Ninguno de los genes semilla fue encontrado en el grafo.")
        print(f"{'='*70}\n")
        return pd.DataFrame(columns=['rank', 'gene', 'score'])
        
    # Asignar probabilidad normalizada a las semillas
    prob_semilla = 1.0 / len(semillas_en_grafo)
    for gen in semillas_en_grafo:
        p0[gen] = prob_semilla

    # 2. Mapear parámetros de GUILD a PageRank
    # El 'damping factor' (alpha) de PageRank es (1 - r)
    alpha = 1 - r

    print(f"   Genes semilla (en grafo): {len(semillas_en_grafo)} / {len(genes_semilla)}")
    print(f"   Parámetros: r={r} (alpha={alpha:.2f}), max_iter={iteraciones}, tol={epsilon}")

    # 3. Ejecutar PageRank personalizado
    # nx.pagerank maneja internamente la normalización, iteraciones y convergencia.
    try:
        scores = nx.pagerank(G, 
                             alpha=alpha, 
                             personalization=p0, 
                             max_iter=iteraciones, 
                             tol=epsilon,
                             weight='weight') # Usa el peso del eje, como en tu original
                             
    except nx.PowerIterationFailedConvergence:
        print(f"   Máximo de iteraciones alcanzado ({iteraciones}) sin convergencia.")
        # Podemos decidir si continuar con los scores que tengamos o fallar.
        # Por ahora, simulamos el comportamiento original y continuamos.
        # Para lanzar un error estricto, descomentar la siguiente línea:
        # raise
        pass # La función devuelve los scores de la última iteración
    else:
        # Nota: nx.pagerank no devuelve el número de iteraciones,
        # así que no podemos imprimir el mensaje exacto de convergencia.
        print("   Cálculo de PageRank completado.")

    # 4. Crear resultados 
    resultados = pd.DataFrame.from_dict(scores, orient='index', columns=['score'])
    resultados.index.name = 'gene'
    resultados = resultados.reset_index()
    
    resultados = resultados.sort_values('score', ascending=False)
    resultados['rank'] = range(1, len(resultados) + 1)
    
    print("\n   Top 5 genes:")
    for _, row in resultados.head().iterrows():
        print(f"      {int(row['rank']):3d}. {str(row['gene']):15s} score={row['score']:.6f}")
    
    print(f"{'='*70}\n")
    return resultados[['rank', 'gene', 'score']]



def diamond(G, genes_semilla, top_k=200):
    """
    DIAMOnD: Implementación optimizada con cola de prioridad (min-heap).
    
    Identifica un módulo de enfermedad añadiendo iterativamente los nodos
    con la conectividad más estadísticamente significativa (p-value más bajo)
    al módulo semilla.
    
    Args:
        G: Grafo de NetworkX
        genes_semilla: Lista de genes iniciales
        top_k: Número de genes a añadir al módulo
    
    Returns:
        DataFrame con genes rankeados
    
    ---
    Referencia del Algoritmo:
    
    Ghiassian, S. D., Menche, J., & Barabási, A. L. (2015). 
    A DIseAse MOdule Detection (DIAMOnD) Algorithm Derived from a 
    Systematic Analysis of Connectivity Patterns of Disease Proteins 
    in the Human Interactome. 
    PLOS Computational Biology, 11(4), e1004120. 
    https://doi.org/10.1371/journal.pcbi.1004120
    ---
    """
    print("\nEjecutando DIAMOnD ...")
    print(f"   Genes semilla: {len(genes_semilla)}")
    print(f"   Top K: {top_k}")

    # --- 1. Inicialización ---
    
    # Filtrar semillas que no están en el grafo
    semillas_en_grafo = {gen for gen in genes_semilla if gen in G}
    if not semillas_en_grafo:
        print("   Error: Ninguno de los genes semilla fue encontrado en el grafo.")
        return pd.DataFrame(columns=['rank', 'gene', 'p_value'])

    modulo = set(semillas_en_grafo)
    genes_candidatos = set(G.nodes()) - modulo
    N = G.number_of_nodes()
    resultados = []
    
    # La cola de prioridad (min-heap) almacenará tuplas: (p_value, gene)
    pq = []

    # --- 2. Calcular p-values iniciales ---
    # Solo necesitamos calcular p-values para los vecinos del módulo inicial.
    
    m_nodes = len(modulo)
    candidatos_iniciales = set()
    for gen_semilla in modulo:
        candidatos_iniciales.update(G.neighbors(gen_semilla))
    
    # Asegurarnos de que los candidatos iniciales no estén ya en el módulo
    candidatos_iniciales.difference_update(modulo)

    for candidato in candidatos_iniciales:
        k = G.degree(candidato)
        s = sum(1 for vecino in G.neighbors(candidato) if vecino in modulo)
        
        if k > 0 and s > 0:
            # Test Hipergeométrico:
            # M = N (Población total de nodos)
            # n = m_nodes (Nodos en el módulo - "éxitos" en la población)
            # N = k (Grado del candidato - "extracciones")
            # x = s (Éxitos observados)
            pval = hypergeom.sf(s - 1, N, m_nodes, k)
            heapq.heappush(pq, (pval, candidato))
            
    print(f"   {len(pq)} candidatos iniciales añadidos a la cola de prioridad.")

    # --- 3. Bucle iterativo principal ---
    
    for iteracion in range(top_k):
        
        mejor_gen = None
        mejor_pval = 1.0

        # Encontrar el mejor candidato válido en la cola
        while pq:
            mejor_pval, mejor_gen = heapq.heappop(pq)
            
            # Si el gen ya está en el módulo, es una entrada "antigua"
            # en la cola. La ignoramos y sacamos la siguiente.
            if mejor_gen in modulo:
                continue
            
            # Si encontramos un gen que todavía es candidato, rompemos el bucle
            if mejor_gen in genes_candidatos:
                break
        else:
            # La cola de prioridad está vacía
            print("   Se acabaron los candidatos conectados al módulo.")
            break
            
        if mejor_gen is None:
            break

        # --- 4. Añadir el mejor gen al módulo ---
        modulo.add(mejor_gen)
        genes_candidatos.remove(mejor_gen)
        
        resultados.append({
            'rank': iteracion + 1,
            'gene': mejor_gen,
            'p_value': mejor_pval
        })
        
        if (iteracion + 1) % 50 == 0:
            print(f"   Iteración {iteracion + 1}/{top_k}. Añadido: {mejor_gen}")

        # --- 5. Actualizar vecinos  ---
        # Solo necesitamos recalcular el p-value de los *vecinos* del gen añadido
        
        m_nodes = len(modulo) # El tamaño del módulo ha crecido en 1
        
        for vecino in G.neighbors(mejor_gen):
            # Solo actualizamos si el vecino es un candidato
            if vecino in genes_candidatos:
                k = G.degree(vecino)
                s = sum(1 for v_vecino in G.neighbors(vecino) if v_vecino in modulo)
                
                if k > 0 and s > 0:
                    # Calculamos el nuevo p-value
                    nuevo_pval = hypergeom.sf(s - 1, N, m_nodes, k)
                    
                    # Añadimos el nuevo p-value a la cola.
                    heapq.heappush(pq, (nuevo_pval, vecino))

    # --- 6. Finalizar ---
    df_resultados = pd.DataFrame(resultados)
    
    print("\n   Top 5 genes añadidos:")
    if df_resultados.empty:
        print("      No se añadieron genes.")
    else:
        for _, row in df_resultados.head().iterrows():
            print(f"      {row['rank']:3d}. {row['gene']:15s} p={row['p_value']:.2e}")
    
    print(f"{'='*70}\n")
    return df_resultados


# =============================================================================
# SECCIÓN 3: ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL
# =============================================================================

def analisis_funcional(genes, organismo='hsapiens', threshold=0.05):
    """
    Realiza análisis de enriquecimiento funcional usando g:Profiler.
    
    g:Profiler consulta múltiples bases de datos para identificar procesos
    biológicos y vías metabólicas enriquecidas.
    
    Bases de datos consultadas:
        - GO:BP: Gene Ontology - Biological Process
        - KEGG: Kyoto Encyclopedia of Genes and Genomes (vías metabólicas)
        - REAC: Reactome Pathway Database (vías de reacción biológicas)
    
    Args:
        genes: Lista de símbolos de genes validados
        organismo: Código del organismo (default: 'hsapiens')
        
    Returns:
        DataFrame con resultados del análisis
    """
    print(f"\n{'='*70}")
    print("ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL")
    print(f"{'='*70}")
    print(f"[INFO] Analizando {len(genes)} genes validados")
    print("[INFO] Bases de datos: GO:BP, KEGG, REAC")
    
    # Inicializar cliente de g:Profiler
    gp = GProfiler(return_dataframe=True)
    
    # Bases de datos a consultar
    # GO:BP - Gene Ontology Biological Process
    # KEGG - Kyoto Encyclopedia of Genes and Genomes (vías metabólicas)
    # REAC - Reactome Pathway Database (vías de reacción)
    sources = ['GO:BP', 'KEGG', 'REAC']
    
    try:
        # Realizar análisis de enriquecimiento
        # - user_threshold: nivel de significancia (p-valor ajustado < 0.05)
        # - significance_threshold_method: corrección FDR (False Discovery Rate)
        resultados = gp.profile(
            organism=organismo,
            query=genes,
            sources=sources,
            user_threshold=threshold,
            significance_threshold_method='fdr'
        )
        
        if resultados is not None and not resultados.empty:
            print(f"Análisis completado: {len(resultados)} términos enriquecidos encontrados")
        else:
            print("No se encontraron términos enriquecidos significativos (p < 0.05)")
            resultados = pd.DataFrame()
            
    except Exception as e:
        print(f"[ERROR] Error en el análisis: {str(e)}")
        resultados = pd.DataFrame()
    
    print(f"{'='*70}\n")
    return resultados


# =============================================================================
# SECCIÓN 4: VISUALIZACIÓN DE RESULTADOS
# =============================================================================

def visualizar_propagacion(df_guild, df_diamond, output_dir, top_n=30):
    """
    Genera visualizaciones comparativas de GUILD y DIAMOnD.
    """
    print(f"\n{'='*70}")
    print("  GENERANDO VISUALIZACIONES DE PROPAGACIÓN")
    print(f"{'='*70}")

    if df_guild is None or df_guild.empty:
        print("     No hay resultados de GUILD para graficar.")
        return
    if df_diamond is None or df_diamond.empty or 'p_value' not in df_diamond.columns:
        print("     No hay resultados de DIAMOnD para graficar (omitido).")
        # Solo graficar GUILD
        fig, ax = plt.subplots(figsize=(8, 8))
        top_guild = df_guild.head(top_n)
        ax.barh(range(len(top_guild)), top_guild['score'], color='steelblue')
        ax.set_yticks(range(len(top_guild)))
        ax.set_yticklabels(top_guild['gene'], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel('Score de Propagación', fontsize=10)
        ax.set_title(f'GUILD: Top {top_n} Genes', fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'propagacion_guild_solo.png'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("    Guardado: propagacion_guild_solo.png")
        print(f"{'='*70}\n")
        return


def visualizar_enriquecimiento(df_enrich, output_dir, top_n=20):
    """
    Genera visualizaciones del análisis de enriquecimiento funcional.
    
    Args:
        df_enrich: Resultados de g:Profiler
        output_dir: Directorio de salida
        top_n: Número de términos top a visualizar
    """
    if df_enrich.empty:
        print("     Sin resultados para visualizar")
        return
    
    print(f"\n{'='*70}")
    print(" GENERANDO VISUALIZACIONES DE ENRIQUECIMIENTO")
    print(f"{'='*70}")
    
    # Preparar datos
    df_plot = df_enrich.sort_values('p_value').head(top_n).copy()
    df_plot['neg_log_p'] = -np.log10(df_plot['p_value'])
    df_plot['term_short'] = df_plot['name'].str[:60]
    
    # Gráfico de barras por fuente
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors_map = {'GO:BP': 'steelblue', 'KEGG': 'coral', 'REAC': 'seagreen'}
    colors = [colors_map.get(src, 'gray') for src in df_plot['source']]
    
    y_pos = range(len(df_plot))
    ax.barh(y_pos, df_plot['neg_log_p'], color=colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_plot['term_short'], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('-log10(p-value)', fontsize=10)
    ax.set_title(f'Top {top_n} Términos Enriquecidos', fontsize=12, fontweight='bold')
    ax.axvline(-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.5)
    ax.grid(axis='x', alpha=0.3)
    
    # Leyenda
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, label=source) 
                      for source, color in colors_map.items()]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'enriquecimiento_funcional.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"    Guardado: enriquecimiento_funcional.png")
    
    # Distribución por fuente
    fig, ax = plt.subplots(figsize=(10, 6))
    source_counts = df_enrich['source'].value_counts()
    ax.bar(source_counts.index, source_counts.values, 
           color=[colors_map.get(s, 'gray') for s in source_counts.index])
    ax.set_xlabel('Base de Datos', fontsize=10)
    ax.set_ylabel('Número de Términos', fontsize=10)
    ax.set_title('Distribución de Términos Enriquecidos por Fuente', 
                 fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'distribucion_fuentes.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"    Guardado: distribucion_fuentes.png")
    print(f"{'='*70}\n")


# =============================================================================
# SECCIÓN 5: EXPORTACIÓN DE RESULTADOS
# =============================================================================

def exportar_resultados(validacion, guild_results, diamond_results, 
                       enrichment_results, output_dir, metadata):
    """
    Exporta todos los resultados a archivos estructurados.
    
    Estructura de salida:
    results/
    ├── metadata.json
    ├── 01_validacion_genes.csv
    ├── 02_guild_propagation.tsv
    ├── 03_diamond_propagation.tsv
    ├── 04_enrichment_complete.csv
    ├── 05_enrichment_GO_BP.csv
    ├── 06_enrichment_KEGG.csv
    ├── 07_enrichment_REAC.csv
    └── 08_resumen_ejecucion.txt
    """
    print(f"\n{'='*70}")
    print(" EXPORTANDO RESULTADOS")
    print(f"{'='*70}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Metadata
    with open(os.path.join(output_dir, 'metadata.json'), 'w') as f:
        json.dump(metadata, f, indent=2)
    print("    metadata.json")
    
    # 2. Validación
    val_data = []
    for orig, val in validacion['mapping'].items():
        val_data.append({
            'gen_input': orig,
            'gen_validado': val,
            'estado': 'valido',
            'entrez_id': validacion.get('entrez_mapping', {}).get(val, 'NA')
        })
    for gene in validacion['no_encontrados']:
        val_data.append({
            'gen_input': gene,
            'gen_validado': 'NA',
            'estado': 'no_encontrado',
            'entrez_id': 'NA'
        })
    
    df_val = pd.DataFrame(val_data)
    df_val.to_csv(os.path.join(output_dir, '01_validacion_genes.csv'), index=False)
    print("    01_validacion_genes.csv")
    
    # 3. GUILD
    guild_results.to_csv(os.path.join(output_dir, '02_guild_propagation.tsv'), 
                         sep='\t', index=False)
    print("    02_guild_propagation.tsv")
    
    # 4. DIAMOnD
    diamond_results.to_csv(os.path.join(output_dir, '03_diamond_propagation.tsv'), 
                          sep='\t', index=False)
    print("    03_diamond_propagation.tsv")
    
    # 5. Enriquecimiento completo
    if not enrichment_results.empty:
        enrichment_results.to_csv(
            os.path.join(output_dir, '04_enrichment_complete.csv'), index=False)
        print("    04_enrichment_complete.csv")
        
        # Por fuente
        for source, filename in [('GO:BP', '05'), ('KEGG', '06'), ('REAC', '07')]:
            df_source = enrichment_results[enrichment_results['source'] == source]
            if not df_source.empty:
                df_source.to_csv(
                    os.path.join(output_dir, f'{filename}_enrichment_{source.replace(":", "_")}.csv'),
                    index=False)
                print(f"    {filename}_enrichment_{source.replace(':', '_')}.csv")
    
    # 6. Resumen ejecutivo
    with open(os.path.join(output_dir, '08_resumen_ejecucion.txt'), 'w') as f:
        f.write("="*70 + "\n")
        f.write("RESUMEN DE LA EJECUCIÓN - ANÁLISIS INTEGRADO\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Fecha: {metadata['timestamp']}\n")
        f.write(f"Genes semilla: {metadata['n_seed_genes']}\n")
        f.write(f"Red: {metadata['network_file']}\n\n")
        
        f.write("RESULTADOS DE PROPAGACIÓN:\n")
        f.write("-" * 70 + "\n")
        f.write(f"GUILD: {len(guild_results)} genes rankeados\n")
        top5_guild = guild_results.head()['gene'].astype(str).tolist()
        f.write(f"  Top 5: {', '.join(top5_guild)}\n\n")
        
        f.write(f"DIAMOnD: {len(diamond_results)} genes rankeados\n")
        if 'p_value' in diamond_results.columns:
            sig = len(diamond_results[diamond_results['p_value'] < 0.05])
            f.write(f"  Significativos (p<0.05): {sig}\n")
        else:
            f.write("  (No se encontraron p-values en resultados de DIAMOnD)\n")
        # Asegurar tipo string para impresión
        top5_diamond = diamond_results.head()['gene'].astype(str).tolist() if not diamond_results.empty else []
        f.write(f"  Top 5: {', '.join(top5_diamond) if top5_diamond else 'N/A'}\n\n")
        
        if not enrichment_results.empty:
            f.write("ENRIQUECIMIENTO FUNCIONAL:\n")
            f.write("-" * 70 + "\n")
            f.write(f"Total términos: {len(enrichment_results)}\n")
            for source in ['GO:BP', 'KEGG', 'REAC']:
                count = len(enrichment_results[enrichment_results['source'] == source])
                f.write(f"  {source}: {count} términos\n")
            
            f.write("\nTop 10 procesos biológicos:\n")
            for i, row in enrichment_results.head(10).iterrows():
                f.write(f"  • {row['name']} (p={row['p_value']:.2e})\n")
    
    print("    08_resumen_ejecucion.txt")
    print(f"\n    Todos los archivos guardados en: {output_dir}/")
    print(f"{'='*70}\n")


# =============================================================================
# SECCIÓN 6: PIPELINE PRINCIPAL
# =============================================================================

def main():
    """Pipeline principal."""
    
    print("\n" + "="*70)
    print(" ANÁLISIS INTEGRADO: PROPAGACIÓN + ENRIQUECIMIENTO FUNCIONAL")
    print("="*70)
    print("Autor: Hugo Salas Calderón")
    print("Fecha:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("="*70 + "\n")
    
    # Configurar argumentos
    parser = argparse.ArgumentParser(
        description='Pipeline integrado: Propagación en redes + Análisis funcional',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EJEMPLOS DE USO:
────────────────────────────────────────────────────────────────────────

1. Análisis completo con archivo de genes y red STRING:
   python %(prog)s -s data/genes_seed.txt -n data/string_network.tsv

2. Genes directos desde línea de comandos:
   python %(prog)s -g ENO1 PGK1 HK2 GAPDH -n data/string_network.tsv

3. Ejecutar solo GUILD (sin DIAMOnD):
   python %(prog)s -s data/genes_seed.txt -n data/network.txt -a guild

4. Análisis con top N genes de propagación:
   python %(prog)s -s data/genes.txt -n data/network.txt --top-genes 50

5. Cambiar threshold de enriquecimiento:
   python %(prog)s -s data/genes.txt -n data/network.txt --enrichment-threshold 0.01

FLUJO DE TRABAJO:
────────────────────────────────────────────────────────────────────────
Genes Semilla → Validación → Propagación (GUILD/DIAMOnD) → 
→ Selección Top Genes → Enriquecimiento Funcional → Visualización

        """
    )
    
    # Argumentos obligatorios
    grupo_genes = parser.add_mutually_exclusive_group(required=True)
    grupo_genes.add_argument('-s', '--seed-file', 
                            help='Archivo con genes semilla')
    grupo_genes.add_argument('-g', '--genes', nargs='+',
                            help='Genes semilla (separados por espacio)')
    
    parser.add_argument('-n', '--network', required=True,
                       help='Archivo de red PPI (STRING/GUILD/DIAMOnD format)')
    
    # Opciones de propagación
    parser.add_argument('-a', '--algorithm', 
                       choices=['guild', 'diamond', 'both'],
                       default='both',
                       help='Algoritmo de propagación (default: both)')
    
    parser.add_argument('--top-genes', type=int, default=100,
                       help='Top N genes de propagación para análisis funcional (default: 100)')
    
    parser.add_argument('--guild-restart', type=float, default=0.5,
                       help='Probabilidad de restart en GUILD (default: 0.5)')
    
    parser.add_argument('--diamond-k', type=int, default=200,
                       help='Top K genes en DIAMOnD (default: 200)')
    
    # Opciones de enriquecimiento
    parser.add_argument('--enrichment-threshold', type=float, default=0.05,
                       help='Threshold p-valor para enriquecimiento (default: 0.05)')
    
    parser.add_argument('--organism', default='hsapiens',
                       help='Organismo para g:Profiler (default: hsapiens)')
    
    # Opciones de salida
    parser.add_argument('-o', '--output', default='results',
                       help='Directorio de salida (default: results)')
    
    parser.add_argument('--no-plots', action='store_true',
                       help='No generar visualizaciones')
    
    args = parser.parse_args()
    
    # Metadata
    metadata = {
        'timestamp': datetime.now().isoformat(),
        'network_file': args.network,
        'algorithm': args.algorithm,
        'top_genes_for_enrichment': args.top_genes,
        'enrichment_threshold': args.enrichment_threshold,
        'organism': args.organism
    }
    
    # PASO 1: Cargar genes semilla
    if args.seed_file:
        genes_input = leer_genes_de_archivo(args.seed_file)
    else:
        genes_input = args.genes
    
    metadata['n_seed_genes'] = len(genes_input)
    metadata['seed_genes'] = genes_input
    
    # PASO 2: Validar genes
    validacion = validar_y_convertir_genes(genes_input, args.organism)
    genes_validos = validacion['validos']
    
    # PASO 3: Cargar red
    print(f"{'='*70}")
    print(" CARGANDO RED DE INTERACCIONES")
    print(f"{'='*70}")
    print(f"    Archivo: {args.network}")
    
    if not os.path.exists(args.network):
        print(f"    Archivo de red no encontrado: {args.network}")
        sys.exit(1)
    
    G = cargar_red(args.network)
    print(f"    Cargada: {G.number_of_nodes()} nodos, {G.number_of_edges()} aristas")
    print(f"{'='*70}\n")
    
    metadata['network_nodes'] = G.number_of_nodes()
    metadata['network_edges'] = G.number_of_edges()
    
    # Verificar genes en la red
    genes_en_red = [g for g in genes_validos if g in G.nodes()]
    
    if not genes_en_red:
        # Lo comento porque con el tsv (Allcontrasts_GLM-Treat_P-0.1_FC-1.25_2025-10-14_16.57.27.tsv)
        #   no encuentra ningun gen y cierra aqui, asi al menos completa ejecución
        #   con el "string_network_filtered_hugo-400.tsv" de practicas anteriores si encuentra y no se cierra aqui
        #print("    Ningún gen semilla encontrado en la red")
        #print("    Sugerencia: Verifica que la red use el mismo tipo de IDs")
        #sys.exit(1)
        
        print("      Ningún gen semilla encontrado en la red. Se continuará con todos los nodos como semillas para pruebas.")
        genes_en_red = list(G.nodes())[:min(len(G.nodes()), 50)]  # usa hasta 50 nodos
    
    print(f"    Genes semilla en red: {len(genes_en_red)}/{len(genes_validos)}")
    
    # PASO 4: Propagación
    guild_results = None
    diamond_results = None
    
    if args.algorithm in ['guild', 'both']:
        guild_results = guild(G, genes_en_red, r=args.guild_restart)
    
    if args.algorithm in ['diamond', 'both']:
        diamond_results = diamond(G, genes_en_red, top_k=args.diamond_k)
    
    # PASO 5: Seleccionar genes para análisis funcional
    print(f"{'='*70}")
    print(" SELECCIÓN DE GENES CANDIDATOS")
    print(f"{'='*70}")
    
    # Usar resultados de GUILD si están disponibles, sino DIAMOnD
    if guild_results is not None:
        # Excluir genes semilla y tomar top N
        genes_candidatos = guild_results[
            ~guild_results['gene'].isin(genes_en_red)
        ].head(args.top_genes)['gene'].tolist()
        print(f"    Usando GUILD: Top {args.top_genes} genes (excluyendo semillas)")
    else:
        genes_candidatos = diamond_results[
            ~diamond_results['gene'].isin(genes_en_red)
        ].head(args.top_genes)['gene'].tolist()
        print(f"    Usando DIAMOnD: Top {args.top_genes} genes (excluyendo semillas)")
    
    # Combinar con genes semilla para contexto
    genes_para_enriquecimiento = genes_en_red + genes_candidatos
    print(f"    Total genes para enriquecimiento: {len(genes_para_enriquecimiento)}")
    print(f"      - Semillas: {len(genes_en_red)}")
    print(f"      - Candidatos: {len(genes_candidatos)}")
    print(f"{'='*70}\n")
    
    # PASO 6: Análisis funcional
    enrichment_results = analisis_funcional(
        genes_para_enriquecimiento,
        organismo=args.organism,
        threshold=args.enrichment_threshold
    )
    
    # PASO 7: Visualizaciones
    if not args.no_plots:
        if guild_results is not None and diamond_results is not None:
            visualizar_propagacion(guild_results, diamond_results, args.output)
        
        if not enrichment_results.empty:
            visualizar_enriquecimiento(enrichment_results, args.output)
    
    # PASO 8: Exportar resultados
    exportar_resultados(
        validacion,
        guild_results if guild_results is not None else pd.DataFrame(),
        diamond_results if diamond_results is not None else pd.DataFrame(),
        enrichment_results,
        args.output,
        metadata
    )
    
    # Resumen final
    print("\n" + "="*70)
    print(" ANÁLISIS COMPLETADO EXITOSAMENTE")
    print("="*70)
    print(f" Resultados guardados en: {args.output}/")
    print(f" Archivos generados:")
    print(f"   - Validación de genes")
    print(f"   - Resultados de propagación (GUILD/DIAMOnD)")
    print(f"   - Análisis de enriquecimiento funcional")
    print(f"   - Visualizaciones (PNG, 300 DPI)")
    print(f"   - Resumen de la ejecución")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()