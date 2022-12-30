import pandas as pd
import numpy as np
import panel as pn
import altair as alt
from pathlib import Path

# Import plotting functions from plotting
from plotting import (
    bracken_raw,
    contig_quality,
    kaiju_raw,
    kaiju_megahit,
    cat_megahit,
    bowtie2_alignment_plot,
)
from utils import parse_bowtielog, parse_fastp_report

pn.extension("tabulator")
pn.extension("vega", sizing_mode="stretch_width", template="fast")
pn.widgets.Tabulator.theme = "modern"


# Header
def header(
    text: str,
    bg_color: str = "#04c273",
    height: int = 150,
    fontsize: str = "px20",
    textalign: str = "center",
):
    """
    Template for markdown header like block
    """
    return pn.pane.Markdown(
        f"""
        {text}
        """,
        background=bg_color,
        height=height,
        margin=10,
        style={
            "color": "white",
            "padding": "10px",
            "text-align": f"{textalign}",
            "font-size": f"{fontsize}",
        },
    )


# Generate the report
def panel_report(
    sample: str,
    coverage_plot_path: str,
    outfolder: str,
) -> None:
    """
    Generates Panel report
    """
    # --- IO --- #
    sample = Path(sample)
    outfolder = Path(outfolder)
    sample_name = sample.parts[-1]

    # --- Alignment and Read Statistics --- #

    # Number of reads and number of reads aligned to reference genome (from bowtie2logfile)
    bowtie2log = list(sample.rglob("*bowtie_raw.log"))[0]

    total_reads, percent_aligned = parse_bowtielog.parse_alignments(bowtie2log)
    number_aligned = int(total_reads * percent_aligned / 100)
    number_unaligned = total_reads - number_aligned

    # Markdown with above text
    alignment_stats = pn.pane.Markdown(
        f"""
        ### Total Number of Reads: 
        {total_reads}
        ### Reads aligned to Human Genome: 
        {number_aligned} ({percent_aligned}%)
        ### Reads NOT aligned to Human Genome:
        {number_unaligned} ({100 - percent_aligned}%)
        """,
        name="Alignment Stats",
    )

    # Bowtie2 alignment plot:
    bowtie_plot = bowtie2_alignment_plot.plot_alignment(bowtie2log).interactive()
    bowtie_plot_pane = pn.pane.Vega(
        bowtie_plot, sizing_mode="stretch_both", name="Alignment Plot"
    )

    # fastp report
    fastp_report = list(sample.rglob("*fastp/*.html"))[0]

    fastp_df = parse_fastp_report.parse_fastp(fastp_report)

    fastp_table = pn.widgets.Tabulator(
        fastp_df, layout="fit_columns", show_index=False, name="Read Summary from FASTP"
    )

    # Header for this section
    alignment_subheader = header(
        text=f"## Alignment and Read statistics",
        bg_color="#04c273",
        height=80,
        textalign="left",
    )

    bowtie_and_stats = pn.Column(
        alignment_stats, pn.layout.Divider(), bowtie_plot_pane, name="Alignment"
    )

    # Section
    alignment_tab = pn.Tabs(bowtie_and_stats, fastp_table)
    alignment_section = pn.Column(alignment_subheader, alignment_tab)

    # --- Raw Classification --- #

    number = 10
    # Raw bracken and kaiju report
    cleaned_bracken_report = list(sample.rglob("*bracken_raw.csv"))[0]
    cleaned_kaiju_report = list(sample.rglob("*kaiju_raw.csv"))[0]

    # Raw bracken and kaiju plots
    bracken_bar_plot = bracken_raw.bar_chart_bracken_raw(
        cleaned_bracken_report, number=number, virus_only=True
    ).interactive()

    bracken_domain_bar_plot = bracken_raw.bar_chart_bracken_raw(
        cleaned_bracken_report, level="domain", virus_only=False
    ).interactive()

    kaiju_raw_plot = kaiju_raw.bar_chart_kaiju_raw(
        file=cleaned_kaiju_report
    ).interactive()

    # Vega panes
    bracken_bar_plot_pane = pn.pane.Vega(
        bracken_bar_plot, sizing_mode="stretch_both", name="Kraken Virus Only"
    )
    bracken_domain_bar_plot_pane = pn.pane.Vega(
        bracken_domain_bar_plot, sizing_mode="stretch_both", name="Kraken All Domains"
    )
    kaiju_raw_plot_pane = pn.pane.Vega(
        kaiju_raw_plot, sizing_mode="stretch_both", name="Kaiju"
    )

    # Header for this section
    raw_header = header(
        text=f"## Classification of Raw Reads",
        bg_color="#04c273",
        height=80,
        textalign="left",
    )

    # Section
    raw_tab = pn.Tabs(
        bracken_bar_plot_pane, bracken_domain_bar_plot_pane, kaiju_raw_plot_pane
    )
    raw_section = pn.Column(raw_header, raw_tab)

    # --- Contig Classification --- #

    # Contigs (Megahit)
    megahit_csv = list(sample.rglob("megahit/*.csv"))[0]
    megahit_histogram = contig_quality.megahit_contig_histogram(
        file=megahit_csv
    ).interactive()

    # Contigs (CAT and Kaiju)
    kaiju_megahit_report = list(sample.rglob("*megahit.out"))[0]
    cat_megahit_out = list(sample.rglob("*contigs_names.txt"))[0]
    cat_kaiju_csv = list(sample.rglob("*cat_kaiju_merged.csv"))[0]

    # plots
    kaiju_bar_plot = kaiju_megahit.bar_chart_kaiju_megahit(
        file=kaiju_megahit_report
    ).interactive()
    cat_bar_plot = cat_megahit.bar_chart_cat_megahit(file=cat_megahit_out).interactive()

    # Vega panes
    megahit_histogram_pane = pn.pane.Vega(
        megahit_histogram, sizing_mode="stretch_both", name="Contig Histogram"
    )
    kaiju_bar_plot_pane = pn.pane.Vega(
        kaiju_bar_plot, sizing_mode="stretch_both", name="Kaiju"
    )
    cat_bar_plot_pane = pn.pane.Vega(
        cat_bar_plot, sizing_mode="stretch_both", name="CAT"
    )

    # cat and kaiju dataframe
    cat_kaiju_csv = list(sample.rglob("*cat_kaiju_merged.csv"))[0]
    cat_kaiju_df = pd.read_csv(cat_kaiju_csv)[
        ["name", "taxon_id", "length", "last_level_kaiju", "last_level_cat", "sequence"]
    ]
    cat_kaiju_table = pn.widgets.Tabulator(
        cat_kaiju_df,
        editors={"sequence": {"type": "editable", "value": False}},
        layout="fit_columns",
        pagination="local",
        page_size=15,
        show_index=False,
        name="Contig Table",
    )

    # Header for this section
    contig_header = header(
        text=f"## Classification of Contigs",
        bg_color="#04c273",
        height=80,
        textalign="left",
    )

    # Section
    contig_tab = pn.Tabs(
        megahit_histogram_pane, kaiju_bar_plot_pane, cat_bar_plot_pane, cat_kaiju_table
    )
    contig_section = pn.Column(contig_header, contig_tab)

    # --- Coverage plots --- #

    # IO
    coverage_plot_path = Path(coverage_plot_path)
    coverage_plots = [
        x
        for x in coverage_plot_path.rglob(f"{sample_name}/*.svg")
        if not "ipynb" in str(x)
    ]

    coverage_tab = pn.Tabs()
    if coverage_plots:
        for plot in coverage_plots:
            name = plot.stem.split("_")[-1].replace(".", " ")
            SVG_pane = pn.pane.SVG(plot, name=name)
            coverage_tab.append(SVG_pane)
    else:
        no_plots = pn.pane.Markdown(
            "## No Coverage plots Available", name="No Coverage Plots"
        )
        coverage_tab.append(no_plots)

    # Header for this section
    coverage_header = header(
        text=f"## Alignment Coverage", bg_color="#04c273", height=80, textalign="left"
    )

    # Section
    coverage_section = pn.Column(coverage_header, coverage_tab)

    # --- Information about programs used --- #

    kaiju_and_kraken_info = pn.pane.Markdown(
        f"""
        ## Kaiju and Kraken2
        Kaiju and Kraken2 are both programs for taxonomic classification of DNA sequences.
        They both use a database of known genetic markers and a sequence alignment algorithm to compare a query sequence to the reference 
        index and score the alignments between the two
        
        * Algorithms: Kaiju uses the MEGAN algorithm for taxonomic classification, while Kraken2 uses a custom algorithm called "k-mer counting." 
        The specific algorithm used can affect the speed and accuracy of the classification.
        * Reference database: Kaiju and Kraken2 use different reference databases for taxonomic classification. 
        Kaiju uses a database called "NINJA" (Non-redundant Improved and Normalized Just-in-time Annotations), 
        which includes annotations for a wide variety of organisms.
        Kraken2 uses a database called "Minikraken," 
        which includes annotations for a more limited set of organisms but is more comprehensive for those organisms.
        """,
        name="Kaiju and Kraken",
    )
    megahit_info = pn.pane.Markdown(
        f"""
        ## MEGAHIT
        MegaHit is a program for assembling DNA sequences, also known as "contigs," from raw sequencing data. 
        It is designed to take a large number of short DNA sequences, called reads, 
        and use them to reconstruct the full-length sequences from which they were derived.
        
        * MEGAHIT processes the raw sequencing data to filter out low-quality reads and remove contaminants.
        * It aligns the reads to a set of known sequences called "k-mers," which are short, fixed-length substrings of the genome.
        * Based on the alignments, MegaHit identifies overlaps between the reads and assembles them into longer contigs using a de Bruijn graph algorithm.
        * The final output of MegaHit is a set of assembled contigs, which can be used for further analysis such as gene prediction or phylogenetic analysis.
        """,
        name="MEGAHIT",
    )

    cat_info = pn.pane.Markdown(
        f"""
        ## CAT
        CAT is designed to identify the species or other taxonomic group of origin for a given DNA sequence, 
        based on the presence or absence of specific genetic markers.
        
        * CAT uses a database of known genetic markers for different taxonomic groups to create a reference index.
        * The DNA sequence to be classified is compared to the reference index using a sequence alignment algorithm, such as BLAST.
        * CAT scores each alignment between the amplicon and the reference markers based on how well they match.
        * The taxonomic group with the highest score is considered the most likely origin of the sequence.
        """,
        name="CAT",
    )

    # Header for this section
    information_header = header(
        text=f"## Information About Programs Used",
        bg_color="#04c273",
        height=80,
        textalign="left",
    )

    # Section
    information_tab = pn.Tabs(kaiju_and_kraken_info, megahit_info, cat_info)
    information_section = pn.Column(information_header, information_tab)

    # --- Create the report --- #

    # header
    head = header(
        text=f"""
        # Pandemic Preparedness Report
        ## Report of Sample {sample.parts[-1]}
        """,
        fontsize="20px",
        bg_color="#011a01",
        height=185,
    )

    all_tabs = pn.Tabs(
        ("Alignment Stats", alignment_section),
        ("Classification of Raw Reads", raw_section),
        ("Classification of Contigs", contig_section),
        ("Alignment Coverage", coverage_section),
        ("Information About Programs", information_section),
        tabs_location="left",
    )

    outfile = outfolder / f"{sample_name}_report.html"
    report = pn.Column(
        head,
        pn.layout.Divider(),
        all_tabs,
    ).save(outfile, title=f"Report {sample_name}")
