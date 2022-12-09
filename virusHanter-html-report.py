from pathlib import Path
import pandas as pd
import json
import altair as alt
# Import plotting functions from plotting
from plotting import bracken_raw, contig_quality, kaiju_raw, kaiju_megahit, cat_megahit
from utils import parse_bowtielog

# function to read in svg code
def return_svg(svg: str):
    with open(svg, "r") as f:
        return f.read()
    
# function to create the report
def html_template_report(
    sample_name: str,
    out_path: str,
    number_of_reads: int,
    percent_aligned: float,
    plot1: str,
    plot2: str,
    svg: str,
) -> None:
    """
    Creates html report
    """
  
    html = f"""
    <!DOCTYPE html>
    <html>
        <head>
            <title>Report of {sample_name} </title>
            <link href="https://cdn.jsdelivr.net/npm/@mdi/font@latest/css/materialdesignicons.min.css" rel="stylesheet">
            <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
            <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
            <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
            <style>
            header {{
            background-color: #333;
            color: white;
            padding: 20px;
            text-align: center;
            }}
            </style>
        </head>
        <body>
            <header>
                <h1>
                Report of {sample_name}
                Number of reads: {number_of_reads}
                Percent aligned to the human genome: {percent_aligned}%
                </h1>
            </header>
            <h1 style="text-align:center;">
                Welcome to The report of the sample! 
            </h1>
            
            <!-- PLOT1 -->
            <h3 style="font-style:italic; text-align:center;"> 
                This is a plot of a beautiful graph I just made! Enjoy!
            </h3>
            
            <div id="vis" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;"></div>
            <script type="text/javascript">
                vegaEmbed("#vis", {plot1});
            </script>
            
            <!-- PLOT2 -->
            <h3 style="font-style:italic; text-align:center;"> 
                This is a plot of a beautiful graph I just made! Enjoy!
            </h3>
            
            <div id="vis2" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;"></div>
            <script type="text/javascript">
                vegaEmbed("#vis2", {plot2});
            </script>
            
            <h1> 
                This is a PDF 
            </h1>
            <!-- <iframe src="{{pdf}}#toolbar=0" width="100%" height="900px" border="0"> -->
            </iframe>
            
            <h1 style="color:orange; text-align:center; font-style:italic;">
                This is an embedded SVG
            </h1>
            <div style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;">
                {svg}
            </div>
            
        </body>
    </html>
    """
    
    output = Path(out_path) / f"{sample_name}-report.html"
    with open(output, "w") as f:
        print(html, file=f)
        
        
# function that writes the report using the right samples
def create_report(
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
    number_of_reads, percent_aligned = parse_bowtielog.parse_alignments(bowtie2log)
    
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

    kaiju_bar_plot = kaiju_raw.bar_chart_kaiju_raw(file=cleaned_kaiju_report).to_json()

    species_and_domain_bracken = (
        alt.hconcat(bracken_bar_plot, bracken_domain_bar_plot)
        .resolve_scale(color="independent")
    ).to_json()


    # test svg
    svg = return_svg("visualization (1).svg")
    
    # generate the html report
    html_template_report(
        sample_name=sample_name, 
        out_path=out_path, 
        number_of_reads=number_of_reads,
        percent_aligned=percent_aligned,
        plot1=species_and_domain_bracken, 
        plot2=kaiju_bar_plot, 
        svg=svg,
    )

# read in the data (testing sample11 for now)
sample_folder = Path("../virusclassification_nextflow/results/sample11_S6/")


# create report
create_report(sample=sample_folder, out_path=".")