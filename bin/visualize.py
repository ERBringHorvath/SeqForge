

def run_visualization(filtered_df, args):
    import os
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    if filtered_df.empty:
        print("\n\033[91mNo data to visualize. Exiting\033[0m")
        return

    required_cols = {'database', "query_file_name", "pident"}
    if not required_cols.issubset(set(filtered_df.columns)):
        print("\033[91mError: filtered_df must contain 'query_file_name', 'database', and 'pident' columns\033[0m")
        return

    # Pivot table: average percent identity per (query, genome)
    heatmap_data = (
        filtered_df
        .groupby(['query_file_name', 'database'])['pident']
        .mean()
        .unstack(fill_value=np.nan)
    )

    mask = heatmap_data.isna()

    num_queries, num_genomes = heatmap_data.shape

    if num_queries < 1 or num_genomes < 2:
        print("\033[91mInsufficient data to generate heatmap (need at least 1 query and 2 genomes)\033[0m")
        return

    # Dynamic sizing and font control
    if num_queries <= 10 and num_genomes <= 10:
        figsize = (10, 8)
        x_fontsize = 15
        y_fontsize = 15
        legend_fontsize = 15
        legend_title_size = 20
    elif num_queries <= 10 and 11 <= num_genomes <= 50:
        figsize = (10, 6)
        x_fontsize = 12
        y_fontsize = 15
        legend_fontsize = 15
        legend_title_size = 20
    elif 11 <= num_queries <= 100 and num_genomes <= 10:
        figsize = (10, 16)
        x_fontsize = 15
        y_fontsize = 5
        legend_fontsize = 15
        legend_title_size = 20
    elif 101 <= num_queries <= 150 and 11 <= num_genomes <= 40:
        figsize = (20, 24)
        x_fontsize = 5
        y_fontsize = 12
        legend_fontsize = 15
        legend_title_size = 20
    elif 101 <= num_queries <= 150 and 41 <= num_genomes <= 100:
        figsize = (30, 30)
        x_fontsize = 5
        y_fontsize = 5
        legend_fontsize = 20
        legend_title_size = 25
    else:
        figsize = (40, 40)
        x_fontsize = 5
        y_fontsize = 5
        legend_fontsize = 25
        legend_title_size = 30

    # Plotting
    plt.figure(figsize=figsize)
    cmap = sns.color_palette("plasma_r", as_cmap=True)
    cmap.set_bad(color='white')

    ax = sns.heatmap(
        heatmap_data,
        cmap=cmap,
        annot=False,
        fmt=".1f",
        linewidths=0.01,
        cbar_kws={'label': 'Average % Identity'},
        mask=mask
    )
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=legend_fontsize)
    cbar.set_label("Average % Identity", fontsize=legend_title_size)

    plt.xlabel("Genome", fontsize=x_fontsize)
    plt.ylabel("Query", fontsize=y_fontsize)
    plt.xticks(rotation=45, ha='right', fontsize=x_fontsize)
    plt.yticks(fontsize=y_fontsize)
    plt.title("BLAST Query Hits (Avg % Identity)")
    plt.tight_layout()

    # Output format
    output_dir = args.output_dir if hasattr(args, 'output_dir') else '.'
    out_ext = "pdf" if getattr(args, "pdf", False) else "png"
    out_file = f"query_identity_heatmap.{out_ext}"
    out_path = os.path.join(output_dir, out_file)

    plt.savefig(out_path, dpi=300 if out_ext == "png" else None)
    print(f"\n\033[92mIdentity heatmap saved to {out_path}\033[0m")


def run_sequence_logo(motif_df, args):
    import pandas as pd
    import logomaker
    import matplotlib.pyplot as plt
    from collections import Counter

    output_format = "pdf" if getattr(args, "pdf", False) else "png"
    output_path = f"motif_logo.{output_format}"
    dpi = 300

    df = motif_df.copy()

    if 'motif_match' not in df.columns:
        print("\033[91mError: 'motif_match' column not found in file.\033[0m")
        return

    motifs = df['motif_match'].dropna().astype(str).tolist()

    if len(motifs) == 0:
        print("\033[91mNo motif matches to visualize.\033[0m")
        return

    valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
    motifs = [m for m in motifs if all(aa in valid_aas for aa in m)]

    if len(motifs) == 0:
        print("\033[91mNo valid motifs remain after filtering.\033[0m")
        return

    motif_len = len(motifs[0])
    if any(len(m) != motif_len for m in motifs):
        print("\033[91mError: Not all motifs are the same length\033[0m")
        return

    position_counts = [Counter() for _ in range(motif_len)]
    for motif in motifs:
        for i, aa in enumerate(motif):
            position_counts[i][aa] += 1

    df_logo = pd.DataFrame.from_dict(
        {i: dict(counter) for i, counter in enumerate(position_counts)},
        orient='columns'
    ).fillna(0)

    # Normalize to frequencies
    df_logo = df_logo.div(df_logo.sum(axis=0), axis=1).fillna(0)

    # Transpose for logomaker (rows = positions, columns = amino acids)
    df_logo = df_logo.T
    df_logo.index.name = "amino_acid"
    df_logo.columns.name = "position"

    # Color mapping by amino acid category
    color_scheme = {
        'A': '#f222ff', 'C': '#ff2975', 'D': 'red', 'E': 'red',
        'F': 'blue', 'G': '#8c1eff', 'H': 'purple', 'I': '#bb496c',
        'K': 'orange', 'L': '#920075', 'M': 'green', 'N': 'pink',
        'P': '#ff901f', 'Q': '#540d6e', 'R': 'orange', 'S': 'cyan',
        'T': 'cyan', 'V': 'green', 'W': '#00e78a', 'Y': 'blue'
    }

    plt.figure(figsize=(motif_len, 3))
    logo = logomaker.Logo(df_logo, font_name='DejaVu Sans Mono', color_scheme=color_scheme)
    plt.title("Motif Match Sequence Logo")
    plt.xlabel("Motif Position")
    plt.ylabel("Amino Acid Frequency")
    plt.tight_layout()

    if output_format == "png":
      plt.savefig(output_path, dpi=dpi)
    else:
      plt.savefig(output_path)
  
    print(f"\n\033[92mSequence logo saved to {output_path}\033[0m\n")
