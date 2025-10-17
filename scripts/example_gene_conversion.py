import argparse
import mygene
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Conversión de símbolos de genes a UniProt IDs usando MyGene.info")
    parser.add_argument("--input", required=True, help="Ruta al archivo de entrada con una lista de genes (uno por línea)")
    parser.add_argument("--output", required=True, help="Ruta al archivo de salida para guardar los resultados")
    args = parser.parse_args()

    # Leer genes
    with open(args.input, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]

    # Inicializar cliente
    mg = mygene.MyGeneInfo()
    results = mg.querymany(genes, scopes='symbol', fields='uniprot', species='human')

    # Procesar resultados
    output_data = []
    for r in results:
        query = r.get('query', '')
        uniprot = r.get('uniprot', {})
        if isinstance(uniprot, dict):
            uniprot_id = uniprot.get('Swiss-Prot', 'N/A')
        elif isinstance(uniprot, list):
            uniprot_id = ", ".join([u.get('Swiss-Prot', '') for u in uniprot if isinstance(u, dict)])
        else:
            uniprot_id = 'N/A'
        output_data.append({"Gene Symbol": query, "UniProt ID": uniprot_id})

    # Guardar resultados
    df = pd.DataFrame(output_data)
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Conversión completada. Resultados guardados en {args.output}")

if __name__ == "__main__":
    main()
