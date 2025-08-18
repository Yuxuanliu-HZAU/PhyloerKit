import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import argparse
import pandas as pd
import chardet
from typing import Optional
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF

def svg_to_pdf(svg_path, pdf_path):
    """
    Convert SVG to PDF using svglib (suitable for simple SVG files)
    """
    drawing = svg2rlg(svg_path)
    renderPDF.drawToFile(drawing, pdf_path)


def generate_heatmap_svg(
    input_file: str,
    output_svg: str,
    reverse: bool = False,
    cell_size: int = 20,
    show_YN_label: bool = True,
    present_color: str = 'black',
    absent_color: str = 'red',
    show_grid: bool = True,
    grid_width: float = 0.6,
    grid_color: str = 'black',
    title: str = 'Occupancy Matrix',
    xlabel: str = 'X_label',
    ylabel: str = 'Y_label'
) -> None:
    """Generate a heatmap of occupancy matrix (Y: present, N: absent)"""
    try:
        # Detect file encoding and read data
        with open(input_file, 'rb') as f:
            encoding = chardet.detect(f.read(10000))['encoding']
        
        csv_data = pd.read_csv(input_file, encoding=encoding)
        
        # Row and column reversal processing
        if reverse:
            csv_data = csv_data.set_index(csv_data.columns[0]).T.reset_index()
            csv_data.columns.name = None

        # Calculate gene coverage
        y_coverage = {}
        for gene in csv_data.columns[1:]:  # Skip the first column (species name)
            present_count = (csv_data[gene] == 'Y').sum()
            coverage_rate = present_count / len(csv_data) * 100
            y_coverage[gene] = coverage_rate

        # Sort genes by coverage rate
        sorted_genes = sorted(y_coverage.keys(), key=lambda x: (-y_coverage[x], x))

        # Prepare plotting data
        n_species = len(csv_data)
        n_genes = len(sorted_genes)
        matrix = np.zeros((n_species, n_genes))

        for i in range(n_species):
            for j, gene in enumerate(sorted_genes):
                val = csv_data.iloc[i, csv_data.columns.get_loc(gene)]
                matrix[i, j] = 1 if str(val).upper() == 'Y' else 0

        # Calculate canvas size
        margin = 100
        label_space = 150
        canvas_width = margin + label_space + n_genes * cell_size + margin
        canvas_height = margin + label_space + n_species * cell_size + margin

        # Create figure
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['axes.unicode_minus'] = False
        fig, ax = plt.subplots(figsize=(canvas_width/100, canvas_height/100), dpi=100)

        # Draw heatmap
        cmap = plt.matplotlib.colors.ListedColormap([absent_color, present_color])
        im = ax.imshow(matrix, cmap=cmap, aspect='equal')

        # Set axes
        ax.set_xticks(range(n_genes))
        ax.set_yticks(range(n_species))
        ax.set_xticklabels(sorted_genes, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(csv_data.iloc[:, 0], fontsize=8)  # The first column is the species name

        # Add title and labels
        ax.set_title(title, fontsize=12, pad=20, fontweight='bold')
        ax.set_xlabel(xlabel, fontsize=10, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=10, fontweight='bold')

        # Add grid lines
        if show_grid:
            ax.set_xticks(np.arange(-0.5, n_genes, 1), minor=True)
            ax.set_yticks(np.arange(-0.5, n_species, 1), minor=True)
            ax.grid(which='minor', color=grid_color, linestyle='-', linewidth=grid_width)

        # Add Y/N labels
        if show_YN_label:
            for i in range(n_species):
                for j in range(n_genes):
                    color = 'white' if matrix[i, j] == 1 else 'black'
                    ax.text(j, i, 'Y' if matrix[i, j] == 1 else 'N',
                           ha='center', va='center', color=color,
                           fontsize=max(6, cell_size//4), fontweight='bold')

        # Add legend
        legend_elements = [
            plt.Rectangle((0,0), 1, 1, facecolor=present_color, label='Present (Y)'),
            plt.Rectangle((0,0), 1, 1, facecolor=absent_color, label='Absent (N)')
        ]
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))

        # Save result
        plt.tight_layout()
        plt.savefig(output_svg, format='svg', dpi=300, bbox_inches='tight')
        plt.close()

    except Exception as e:
        print(f"Error: {str(e)}")
        if "UnicodeDecodeError" in str(e):
            print("Tips: Try specifying the file encoding like --encoding gbk")

def main():
    parser = argparse.ArgumentParser(description="Generate a heatmap of occupancy matrix (Y: present, N: absent)")
    parser.add_argument('-i', '--input', required=True, help='The path of the input csv file')
    parser.add_argument('-o', '--output', required=True, help='The path of the output svgfile')
    parser.add_argument('--reverse', action='store_true', help='Reverse the rows and columns')
    parser.add_argument('--cell_size', type=int, default=20, help='The size of the cells')
    parser.add_argument('--show_YN_label', action='store_true', help='Show the Y/N labels')
    parser.add_argument('--present_color', default='black', help='The color of the present (Y)')
    parser.add_argument('--absent_color', default='red', help='The color of the absent (N)')
    parser.add_argument('--show_grid', action='store_true', help='Show the grid lines')
    parser.add_argument('--grid_width', type=float, default=0.6, help='The width of the grid lines')
    parser.add_argument('--grid_color', default='black', help='The color of the grid lines')
    parser.add_argument('--title', default='Occupancy Matrix', help='The title of the plot')
    parser.add_argument('--xlabel', default='X_label', help='The label of the x-axis')
    parser.add_argument('--ylabel', default='Y_label', help='The label of the y-axis')
    parser.add_argument('--pdf', action='store_true', help='Whether to convert the SVG to PDF')

    args = parser.parse_args()
    generate_heatmap_svg(args.input, args.output, args.reverse, args.cell_size, args.show_YN_label, args.present_color, args.absent_color, args.show_grid, args.grid_width, args.grid_color, args.title, args.xlabel, args.ylabel)
    if os.path.exists(args.output):
        print(f"[phyloerkit ✓]: The heatmap (svg format) has been saved to: {args.output}")
    else:
        print(f"[phyloerkit Error]: The heatmap (svg format) has not been saved to: {args.output}")
    
    if args.pdf:
        output_pdf = args.output.replace('.svg', '.pdf')
        svg_to_pdf(args.output, output_pdf)
        if os.path.exists(output_pdf):
            print(f"[phyloerkit ✓]: The heatmap (pdf format) has been saved to: {output_pdf}")
        else:
            print(f"[phyloerkit Error]: The heatmap (pdf format) has not been saved to: {output_pdf}")

if __name__ == "__main__":
    main()
