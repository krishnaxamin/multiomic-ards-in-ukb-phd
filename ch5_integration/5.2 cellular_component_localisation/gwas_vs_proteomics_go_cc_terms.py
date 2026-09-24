""" Plot GO-CC terms enriched in GWAS vs proteomics genes by counts in broader themes. """
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.use('TkAgg')

""" Read in data """
all_olink = pd.read_csv('~/data/internal/proteomics/enrichment_analysis/olink_vs_all_nom_sig.csv')
high_effect_proteins = pd.read_csv('~/data/internal/proteomics/enrichment_analysis/high_effect_cox_pathways_gocc_enrichment.csv')
high_effect_genes = pd.read_csv('~/data/internal/genomics/magma/lifestyles/high_effect_pathways_gocc_enrichment.csv')

collated = pd.concat([
    all_olink[(all_olink['fdr_sig'] == 1) & (all_olink['ontology'] == 'go_C')][['annotation_name']].assign(cat='All Olink'),
    high_effect_proteins[high_effect_proteins['fdr_sig'] == 1][['annotation_name']].assign(cat='High-effect proteomics\n pathways proteins'),
    high_effect_genes[high_effect_genes['fdr_sig'] == 1][['annotation_name']].assign(cat='High-effect GWAS\n pathways genes')
])

""" Define broader themes """
theme_map = {
    # Secretion / circulation
    "specific granule": "Secretion / circulation",
    "tertiary granule": "Secretion / circulation",
    "secretory granule": "Secretion / circulation",
    "azurophil granule": "Secretion / circulation",
    "platelet alpha granule": "Secretion / circulation",
    "platelet dense granule": "Secretion / circulation",
    "platelet alpha granule lumen": "Secretion / circulation",
    "platelet dense granule lumen": "Secretion / circulation",
    "secretory granule lumen": "Secretion / circulation",
    "specific granule lumen": "Secretion / circulation",
    "tertiary granule lumen": "Secretion / circulation",
    "azurophil granule lumen": "Secretion / circulation",
    "Golgi lumen": "Secretion / circulation",
    "endoplasmic reticulum lumen": "Secretion / circulation",
    "secretory granule membrane": "Secretion / circulation",
    "lamellar body": "Secretion / circulation",
    "blood microparticle": "Secretion / circulation",
    "chylomicron": "Secretion / circulation",
    "COPII-coated ER to Golgi transport vesicle": "Secretion / circulation",
    "ER to Golgi transport vesicle membrane": "Secretion / circulation",
    "endoplasmic reticulum-Golgi intermediate compartment": "Secretion / circulation",
    "endoplasmic reticulum-Golgi intermediate compartment membrane": "Secretion / circulation",

    # ECM / cell surface
    "external encapsulating structure": "ECM / cell surface",
    "extracellular matrix": "ECM / cell surface",
    "external side of plasma membrane": "ECM / cell surface",
    "extrinsic component of plasma membrane": "ECM / cell surface",
    "protein complex involved in cell adhesion": "ECM / cell surface",
    "basement membrane": "ECM / cell surface",
    "perisynaptic extracellular matrix": "ECM / cell surface",
    "synapse-associated extracellular matrix": "ECM / cell surface",
    "specialized extracellular matrix": "ECM / cell surface",
    "collagen trimer": "ECM / cell surface",
    "focal adhesion": "ECM / cell surface",
    "cell-substrate junction": "ECM / cell surface",
    "adherens junction": "ECM / cell surface",
    "axon": "ECM / cell surface",

    # Membrane receptors / signalling
    "receptor complex": "Membrane receptors / signalling",
    "plasma membrane raft": "Membrane receptors / signalling",
    "membrane raft": "Membrane receptors / signalling",
    "membrane microdomain": "Membrane receptors / signalling",
    "plasma membrane signaling receptor complex": "Membrane receptors / signalling",
    "basal plasma membrane": "Membrane receptors / signalling",
    "basal part of cell": "Membrane receptors / signalling",
    "apical part of cell": "Membrane receptors / signalling",
    "integrin complex": "Membrane receptors / signalling",
    "catenin complex": "Membrane receptors / signalling",
    "semaphorin receptor complex": "Membrane receptors / signalling",
    "voltage-gated potassium channel complex": "Membrane receptors / signalling",
    "potassium channel complex": "Membrane receptors / signalling",
    "acetylcholine-gated channel complex": "Membrane receptors / signalling",

    # Lipid complexes
    "protein-lipid complex": "Lipid complexes",
    "plasma lipoprotein particle": "Lipid complexes",
    "high-density lipoprotein particle": "Lipid complexes",
    "lipoprotein particle": "Lipid complexes",
    "very-low-density lipoprotein particle": "Lipid complexes",
    "triglyceride-rich plasma lipoprotein particle": "Lipid complexes",

    # Immune
    "immunological synapse": "Immune",
    "immunoglobulin complex": "Immune",
    "immunoglobulin complex, circulating": "Immune",
    "Bcl-2 family protein complex": "Immune",

    # Endocytosis / non-secretory vesicular transport
    "vesicle lumen": "Endocytosis / non-secretory-specific\nvesicular transport",
    "cytoplasmic vesicle lumen": "Endocytosis / non-secretory-specific\nvesicular transport",
    "lysosomal lumen": "Endocytosis / non-secretory-specific\nvesicular transport",
    "endocytic vesicle": "Endocytosis / non-secretory-specific\nvesicular transport",
    "endocytic vesicle membrane": "Endocytosis / non-secretory-specific\nvesicular transport",
    "endocytic vesicle lumen": "Endocytosis / non-secretory-specific\nvesicular transport",
    "coated vesicle membrane": "Endocytosis / non-secretory-specific\nvesicular transport",
    "clathrin-coated vesicle membrane": "Endocytosis / non-secretory-specific\nvesicular transport",
    "clathrin-coated endocytic vesicle membrane": "Endocytosis / non-secretory-specific\nvesicular transport",
    "clathrin-coated endocytic vesicle": "Endocytosis / non-secretory-specific\nvesicular transport",
    "transport vesicle membrane": "Endocytosis / non-secretory-specific\nvesicular transport",
    "multivesicular body": "Endocytosis / non-secretory-specific\nvesicular transport",
    "early endosome": "Endocytosis / non-secretory-specific\nvesicular transport",
    "endosome membrane": "Endocytosis / non-secretory-specific\nvesicular transport",
    "endosome lumen": "Endocytosis / non-secretory-specific\nvesicular transport",
    "vacuolar lumen": "Endocytosis / non-secretory-specific\nvesicular transport",
    "peroxisomal matrix": "Endocytosis / non-secretory-specific\nvesicular transport",
    "microbody lumen": "Endocytosis / non-secretory-specific\nvesicular transport",
    "endolysosome": "Endocytosis / non-secretory-specific\nvesicular transport",
    "endolysosome membrane": "Endocytosis / non-secretory-specific\nvesicular transport",

    # Nuclear / chromatin
    "nuclear inclusion body": "Nuclear / chromatin",
    "nuclear periphery": "Nuclear / chromatin",
    "telomere cap complex": "Nuclear / chromatin",
    "nuclear telomere cap complex": "Nuclear / chromatin",
    "protein-DNA complex": "Nuclear / chromatin",
    "nucleosome": "Nuclear / chromatin",

    # Other
    "other organism part": "Other",
    "symbiont cell surface": "Other",
    "A band": "Other",
    "serine-type peptidase complex": "Other",
    "peptidase complex": "Other",
    "proteasome regulatory particle, base subcomplex": "Other",
}

""" Map and count """
collated['theme'] = collated['annotation_name'].map(theme_map)
totals = collated.value_counts('cat').reset_index().set_axis(['cat', 'cat_count'], axis=1)
theme_totals = collated.value_counts(['cat', 'theme']).reset_index().set_axis(['cat', 'Theme', 'theme_count'], axis=1)
collated_counts = theme_totals.merge(totals)
collated_counts['proportion'] = collated_counts['theme_count'] * 100 / collated_counts['cat_count']

collated.sort_values(by=['cat', 'theme', 'annotation_name'], ascending=True).to_csv(
    '~/data/internal/integration/gwas_vs_proteomics_go_cc_terms.csv', index=False)

""" Plot """
seen = set()
seq = list(theme_map.values())
seen_add = seen.add
ordered_themes = [x for x in seq if not (x in seen or seen_add(x))]

collated_counts['cat'] = pd.Categorical(collated_counts['cat'], ['All Olink', 'High-effect proteomics\n pathways proteins', 'High-effect GWAS\n pathways genes'])

plt.figure(figsize=(8, 6))
sns.barplot(data=collated_counts, x='cat', y='proportion', hue='Theme', hue_order=ordered_themes)
plt.xlabel('')
plt.tick_params(axis='x', labelsize=8)
plt.ylabel('%')
plt.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', fontsize=8)
plt.gca().yaxis.grid(True, which='major', color='grey', alpha=0.6)
plt.gca().set_axisbelow(True)
plt.tight_layout()
plt.savefig('~/ch5_integration/5.2 cellular_component_localisation/plots/gwas_vs_proteomics_go_cc_terms.png')
plt.savefig('~/ch5_integration/5.2 cellular_component_localisation/plots/gwas_vs_proteomics_go_cc_terms.svg')
plt.close()
