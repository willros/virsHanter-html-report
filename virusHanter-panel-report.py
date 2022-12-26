import pandas as pd
import numpy as np
import panel as pn
import altair as alt
from pathlib import Path
# Import plotting functions from plotting
from plotting import bracken_raw, contig_quality, kaiju_raw, kaiju_megahit, cat_megahit, bowtie2_alignment_plot
from utils import parse_bowtielog, parse_fastp_report

pn.extension("tabulator")
pn.extension("vega", sizing_mode="stretch_width", template="fast")
pn.widgets.Tabulator.theme = 'modern'


# function to create the Panel report
def panel_template_report(
    sample_name: str,
    out_path: str,
    total_reads: int,
    number_aligned: int,
    number_unaligned: int,
    bowtie_plot: str,
    fastp_df: str,
    megahit_histogram: str,
    kaiju_raw: str,
    kraken_raw: str,
    kaiju_and_cat: str,
    cat_kaiju_df: str,
    svg: str,
) -> None:
    """
    Creates Panel report
    """
    pass




# function that writes the report using the right samples
def create_panel_report(
    sample: str,
    out_path: str,
) -> None:
    """
    Generate the report
    """
    # Sample and sample name
    sample = Path(sample)
    sample_name = sample.parts[-1]
    
    # Number of bars to include in the figures:
    number = 10
    
    # Number of reads and number of reads aligned to reference genome (from bowtie2logfile)
    bowtie2log = list(sample.rglob("*bowtie_raw.log"))[0]
    
    total_reads, percent_aligned = parse_bowtielog.parse_alignments(bowtie2log)
    number_aligned = int(total_reads * percent_aligned / 100)
    number_unaligned = total_reads - number_aligned
    
    # Bowtie2 alignment plot:
    bowtie_plot = bowtie2_alignment_plot.plot_alignment(bowtie2log).interactive()
    
    # fastp dataframe
    fastp_report = list(sample.rglob("*fastp/*.html"))[0]
    fastp_df = parse_fastp_report.parse_fastp(fastp_report)
    fastp_table = pn.widgets.Tabulator(
        fastp_df, 
        layout='fit_columns',
        show_index=False,
        name="FastP Summary"
    )
    
    # Raw bracken and kaiju report
    cleaned_bracken_report = list(sample.rglob("*bracken_raw.csv"))[0]
    cleaned_kaiju_report = list(sample.rglob("*kaiju_raw.csv"))[0]
    
    # Raw bracken and kaiju plots
    bracken_bar_plot = bracken_raw.bar_chart_bracken_raw(
        cleaned_bracken_report, number=number,virus_only=True
    )

    bracken_domain_bar_plot = bracken_raw.bar_chart_bracken_raw(
        cleaned_bracken_report, level="domain", virus_only=False
    )

    kaiju_raw_plot = kaiju_raw.bar_chart_kaiju_raw(file=cleaned_kaiju_report).interactive()

    species_and_domain_bracken = (
        alt.hconcat(bracken_bar_plot, bracken_domain_bar_plot)
        .resolve_scale(color="independent")
    ).interactive()
    
    
    # Contigs (Megahit)
    megahit_csv = list(sample.rglob("megahit/*.csv"))[0]
    megahit_histogram = contig_quality.megahit_contig_histogram(file=megahit_csv).interactive()
    
    # Contigs (CAT and Kaiju)
    # files
    kaiju_megahit_report = list(sample.rglob("*megahit.out"))[0]
    cat_megahit_out = list(sample.rglob("*contigs_names.txt"))[0]
    cat_kaiju_csv = list(sample.rglob("*cat_kaiju_merged.csv"))[0]

    # plots
    kaiju_bar_plot = kaiju_megahit.bar_chart_kaiju_megahit(file=kaiju_megahit_report)
    cat_bar_plot = cat_megahit.bar_chart_cat_megahit(file=cat_megahit_out)
    kaiju_and_cat = (
        alt.hconcat(kaiju_bar_plot, cat_bar_plot)
        .resolve_scale(color="independent")
        .interactive()
    )

    # cat and kaiju dataframe
    cat_kaiju_df = pd.read_csv(cat_kaiju_csv)[["name", "taxon_id", "length", "last_level_kaiju", "last_level_cat", "sequence"]]
    cat_kaiju_table = pn.widgets.Tabulator(
        cat_kaiju_df, 
        editors={
            'sequence': {'type': 'editable', 'value': False}
        },
        layout='fit_columns',
        show_index=False,
        name="Contig Table"
    )


    # test svg
    svg = return_svg("visualization (1).svg")
    
    # generate the html report
    panel_template_report(
        sample_name=sample_name, 
        out_path=out_path, 
        total_reads=total_reads,
        number_aligned=number_aligned,
        number_unaligned=number_unaligned,
        kraken_raw=species_and_domain_bracken, 
        kaiju_raw=kaiju_raw_plot, 
        svg=svg,
        bowtie_plot=bowtie_plot,
        fastp_df=fastp_df,
        megahit_histogram=megahit_histogram,
        kaiju_and_cat=kaiju_and_cat,
        cat_kaiju_df=cat_kaiju_df,
    )