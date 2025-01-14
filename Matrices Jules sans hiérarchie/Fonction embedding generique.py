#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:59:50 2024

@author: massambadiop
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sentence_transformers import SentenceTransformer

DF = pd.read_csv('Chemin vers votre fichier') # Si autre format, adapter

# Charger le modèle RoBERTa
model_RoBERTa = SentenceTransformer('stsb-roberta-large')
model_SBERT = SentenceTransformer('all-MiniLM-L6-v2')



def embedding_analysis(text, model, data, target_column, metric='cosine', output_results=True):
    """
    Analyze text embeddings, calculate distance matrices, and find most similar/distant pairs.

    Args:
        text (list[str]): List of sentences or texts.
        model (object): Embedding model (e.g., SentenceTransformer) with an `encode` method.
        data (pd.DataFrame): DataFrame containing the original texts.
        target_column (str): Name of the column in `data` containing the texts to analyze.
        metric (str): Metric for distance calculation ('cosine', 'euclidean', etc.).
        output_results (bool): If True, prints results to the console.

    Returns:
        dict: A dictionary containing:
            - embeddings (pd.DataFrame): DataFrame of embeddings.
            - distance_matrix (np.ndarray): Pairwise distance matrix.
            - most_similar (list): Top 5 most similar pairs with distances and indices.
            - most_distant (list): Top 5 most distant pairs with distances and indices.
            - similar_content (list): Content of most similar text pairs.
            - distant_content (list): Content of most distant text pairs.
    """
    # Input validation
    if not isinstance(text, list) or not all(isinstance(t, str) for t in text):
        raise ValueError("`text` must be a list of strings.")
    if not hasattr(model, "encode"):
        raise ValueError("`model` must have an `encode` method.")
    if not isinstance(data, pd.DataFrame):
        raise ValueError("`data` must be a pandas DataFrame.")
    if target_column not in data.columns:
        raise ValueError(f"DataFrame must contain a `{target_column}` column.")

    # Generate embeddings
    embeddings = np.array(model.encode(text, batch_size=32, show_progress_bar=True))

    # Calculate distance matrix
    distance_matrix = cdist(embeddings, embeddings, metric=metric)

    # Find the 5 most similar pairs
    np.fill_diagonal(distance_matrix, np.inf)  # Ignore self-distances for minimum search
    flat_min_indices = np.argsort(distance_matrix.flatten())[:5]
    row_min_indices, col_min_indices = np.unravel_index(flat_min_indices, distance_matrix.shape)
    most_similar = [(distance_matrix[i, j], (i, j)) for i, j in zip(row_min_indices, col_min_indices)]
    similar_content = [(data.iloc[i][target_column], data.iloc[j][target_column]) for i, j in zip(row_min_indices, col_min_indices)]

    # Find the 5 most distant pairs
    np.fill_diagonal(distance_matrix, -np.inf)  # Ignore self-distances for maximum search
    flat_max_indices = np.argsort(distance_matrix.flatten())[::-1][:5]
    row_max_indices, col_max_indices = np.unravel_index(flat_max_indices, distance_matrix.shape)
    most_distant = [(distance_matrix[i, j], (i, j)) for i, j in zip(row_max_indices, col_max_indices)]
    distant_content = [(data.iloc[i][target_column], data.iloc[j][target_column]) for i, j in zip(row_max_indices, col_max_indices)]
    
    # Reinitialise distance matrix
    distance_matrix = cdist(embeddings, embeddings, metric=metric)

    # Output results if requested
    if output_results:
        print("Distance Matrix:")
        print(distance_matrix)
        print("\n5 Most Similar Pairs:")
        print(pd.DataFrame(most_similar, columns=['Distance', 'Indices']))
        print("\nContent of Similar Pairs:")
        print(pd.DataFrame(similar_content, columns=['Text 1', 'Text 2']))
        print("\n5 Most Distant Pairs:")
        print(pd.DataFrame(most_distant, columns=['Distance', 'Indices']))
        print("\nContent of Distant Pairs:")
        print(pd.DataFrame(distant_content, columns=['Text 1', 'Text 2']))

    # Return results
    return {
        'embeddings': pd.DataFrame(embeddings),
        'distance_matrix': distance_matrix,
        'most_similar': most_similar,
        'most_distant': most_distant,
        'similar_content': similar_content,
        'distant_content': distant_content
    }


### Usage


model = model_SBERT

# Run analysis
results = embedding_analysis(DF['liste_phenotype'].tolist(), model, DF, target_column='liste_phenotype', metric='cosine', output_results=True)
result['embeddings']
# Save results to CSV (optional)
results['embeddings'].to_csv('embeddings_pheno_SBERT.csv', index=False)


