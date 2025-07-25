# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Elijah Bring Horvath

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
    
    def get_dynamic_plot_settings(num_queries, num_genomes):
        base_size = 5.0
        min_col_width = 1.5
        min_row_height = 0.4

        width = max(10, base_size + min_col_width * num_genomes)
        height = max(15, base_size + min_row_height * num_queries)
        figsize = (width, height)

        x_fontsize = max(5, 38 - 0.15 * num_genomes)
        x_title = max(20, 38 - 0.02 * num_queries)
        y_fontsize = max(5, 38 - 0.02 * num_queries)
        y_title = max(20, 38 - 0.02 * num_queries)

        legend_fontsize = max(15, min(30, 0.8 * (x_fontsize + y_fontsize)))
        legend_title_size = legend_fontsize + 8

        return figsize, x_fontsize, y_fontsize, legend_fontsize, legend_title_size, y_title, x_title

    strongest_hits = (
        filtered_df.sort_values(['query_file_name', 'database', 'pident'], ascending=[True, True, False])
        .drop_duplicates(subset=['query_file_name', 'database'], keep='first')
    )

    heatmap_data = (
        strongest_hits
        .pivot(index='query_file_name', columns='database', values='pident')
    )

    mask = heatmap_data.isna()

    num_queries, num_genomes = heatmap_data.shape

    if num_queries < 1 or num_genomes < 2:
        print("\033[91mInsufficient data to generate heatmap (need at least 1 query and 2 genomes)\033[0m")
        return

    figsize, x_fontsize, y_fontsize, legend_fontsize, legend_title_size, y_title, x_title = get_dynamic_plot_settings(num_queries, num_genomes)

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
    ax.set_aspect('auto')

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=legend_fontsize)
    cbar.set_label("Average % Identity", fontsize=legend_title_size)

    # plt.figtext(0.01, 0.01, "Strongest hit per query per genome",
    #             horizontalalignment='left', fontsize=y_fontsize, style='italic')
    plt.xlabel("Genome", fontsize=x_title)
    plt.ylabel("Query", fontsize=y_title)
    plt.xticks(rotation=45, ha='right', fontsize=x_fontsize)
    plt.yticks(fontsize=y_fontsize, rotation=0)
    plt.title("BLAST Query Hits (Avg % Identity)\nStrongest hit per query per genome", fontsize=x_fontsize)
    plt.tight_layout()

    # Output format
    output_dir = args.output if hasattr(args, 'output') else '.'
    out_ext = "pdf" if getattr(args, "pdf", False) else "png"
    out_file = f"query_identity_heatmap.{out_ext}"
    out_path = os.path.join(output_dir, out_file)

    plt.savefig(out_path, dpi=300 if out_ext == "png" else None, bbox_inches='tight', pad_inches=0.1)
    heatmap_data.to_csv(os.path.join(output_dir, "heatmap_matrix.csv"))
    print("\033[95mVisualization Files:\n\033[0m"
          f"\033[92mIdentity heatmap saved to {out_path}\033[0m\n"
          f"\033[92mHeatmap data matrix saved to {output_dir}heatmap_matrix.csv\033[0m")


def run_sequence_logo(motif_df, args):
    import os
    import pandas as pd
    import logomaker
    import matplotlib.pyplot as plt
    from collections import Counter

    output_dir = args.output if hasattr(args, 'output') else '.'
    out_ext = "pdf" if getattr(args, "pdf", False) else "png"

    # Detect motif columns
    motif_cols = [col for col in motif_df.columns if col.startswith("motif_") and not col.endswith("_pattern")]
    if not motif_cols:
        print("\033[91mNo motif columns found for visualization.\033[0m")
        return

    for motif_col in motif_cols:
        pattern_col = f"{motif_col}_pattern"
        if pattern_col not in motif_df.columns:
            print(f"\033[93mSkipping {motif_col} (no pattern info found)\033[0m")
            continue

        motifs = (
            motif_df[motif_col]
            .dropna()
            .astype(str)
            .apply(str.strip)
            .loc[lambda s: s != '']
            .tolist()
        )

        motifs = [m for m in motifs if set(m).issubset(set("ACDEFGHIKLMNPQRSTVWY"))]

        if not motifs:
            print(f"\033[93mNo valid motifs found in {motif_col}. Skipping.\033[0m")
            continue

        motif_len = len(motifs[0])
        if any(len(m) != motif_len for m in motifs):
            print(f"\033[91mSkipping {motif_col}: Not all motifs are the same length\033[0m")
            continue

        # Count amino acids at each position
        position_counts = [Counter() for _ in range(motif_len)]
        for motif in motifs:
            for i, aa in enumerate(motif):
                position_counts[i][aa] += 1

        df_logo = pd.DataFrame.from_dict(
            {i: dict(counter) for i, counter in enumerate(position_counts)},
            orient='columns'
        ).fillna(0)

        # Normalize and format
        df_logo = df_logo.div(df_logo.sum(axis=0), axis=1).fillna(0).T
        df_logo.index.name = "amino_acid"
        df_logo.columns.name = "position"

        # Color scheme
        color_scheme = {
            'A': '#f222ff', 'C': '#ff2975', 'D': 'red', 'E': 'red',
            'F': 'blue', 'G': '#8c1eff', 'H': 'purple', 'I': '#bb496c',
            'K': 'orange', 'L': '#920075', 'M': 'green', 'N': 'pink',
            'P': '#ff901f', 'Q': '#540d6e', 'R': 'orange', 'S': 'cyan',
            'T': 'cyan', 'V': 'green', 'W': '#00e78a', 'Y': 'blue'
        }

        plt.figure(figsize=(motif_len, 3))
        logomaker.Logo(df_logo, font_name='DejaVu Sans Mono', color_scheme=color_scheme)

        motif_pattern = motif_df[pattern_col].dropna().iloc[0]
        plt.title(f"Sequence Logo: {motif_col} ({motif_pattern})")
        plt.xlabel("Motif Position")
        plt.ylabel("Amino Acid Frequency")
        plt.tight_layout()

        out_file = f"{motif_col}_logo.{out_ext}"
        out_path = os.path.join(output_dir, out_file)
        plt.savefig(out_path, dpi=300 if out_ext == "png" else None, bbox_inches='tight', pad_inches=0.1)
        plt.close()

        print(f"\033[92m{motif_col} sequence logo saved to {out_path}\033[0m")
