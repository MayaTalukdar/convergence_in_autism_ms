#conda activate pdfkit_env

import pandas as pd
import holoviews as hv
from holoviews import opts, dim

hv.extension('matplotlib')
hv.output(fig='svg', size=500)
output_path = "/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Final_Submission/Code/Figure_1/"

nodes_df = pd.read_csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/circos_plot_cross_enrichment/Input_Files/nodes.csv") #In Input_Files_Not_Generated_By_Scripts

#all genes 
chords = pd.read_csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/circos_plot_cross_enrichment/Input_Files/adj_matrix_all.csv")
chords.index = chords.type
chords["color"] = ["red" if x == "opp" else "green" for x in chords.index]
min_value = chords['value'].min()
max_value = chords['value'].max()
chords['normalized_value'] = 5 * (chords['value'] - min_value) / (max_value - min_value)
nodes = hv.Dataset(nodes_df)
chord_plot = hv.Chord((chords, nodes)).opts(
opts.Chord(edge_color=chords.color, labels='name', edge_linewidth = chords.normalized_value))
hv.save(chord_plot, filename=output_path + 'circos_all_genes.pdf')

#all sfari genes
chords = pd.read_csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/circos_plot_cross_enrichment/Input_Files/adj_matrix_sfari.csv")
chords.index = chords.type
chords["color"] = ["red" if x == "opp" else "green" for x in chords.index]
chords['normalized_value'] = 5 * (chords['value'] - min_value) / (max_value - min_value)
nodes = hv.Dataset(nodes_df)
chord_plot = hv.Chord((chords, nodes)).opts(
opts.Chord(edge_color=chords.color, labels='name', edge_linewidth = chords.normalized_value))
hv.save(chord_plot, filename=output_path + 'circos_all_sfari_genes.pdf')

#all high confidence sfari genes
chords = pd.read_csv("/n/groups/walsh/indData/Maya/Finalized_Convergence_In_ASD/Cell_Lines/figures/circos_plot_cross_enrichment/Input_Files/adj_matrix_high_conf_sfari.csv")
chords.index = chords.type
chords["color"] = ["red" if x == "opp" else "green" for x in chords.index]
chords['normalized_value'] = 5 * (chords['value'] - min_value) / (max_value - min_value)
nodes = hv.Dataset(nodes_df)
chord_plot = hv.Chord((chords, nodes)).opts(
opts.Chord(edge_color=chords.color, labels='name', edge_linewidth = chords.normalized_value))
hv.save(chord_plot, filename=output_path + 'circos_all_high_conf_sfari_genes.pdf')
