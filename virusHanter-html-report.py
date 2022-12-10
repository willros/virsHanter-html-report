from pathlib import Path
import pandas as pd
import json
import altair as alt
# Import plotting functions from plotting
from plotting import bracken_raw, contig_quality, kaiju_raw, kaiju_megahit, cat_megahit, bowtie2_alignment_plot
from utils import parse_bowtielog

# function to read in svg code
def return_svg(svg: str):
    with open(svg, "r") as f:
        return f.read()
    
# function to create the report
def html_template_report(
    sample_name: str,
    out_path: str,
    total_reads: int,
    number_aligned: int,
    number_unaligned: int,
    bowtie_plot: str,
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

                h1 {{
                color: #2ecc71;
                }}

                h2 {{
                color: #3498db;
                }}

                h3 {{
                color: #9b59b6;
                }}

                h4 {{
                color: #f1c40f;
                }}

                div {{
                color: #e67e22;
                }}

            </style>

        </head>
        <body>
            <header>
                <h1>
                Report of {sample_name}
                </h1>
            </header>
            <h1 style="text-align:center;">
                Welcome to The report of the sample! 
            </h1>
            
            <!-- PLOT ALIGNED -->
            <h3 style="font-style:italic; text-align:center;"> 
                This is the bowtie plot!
                Number of reads: {total_reads}
                Number aligned to the human genome: {number_aligned}
                Number unaligned: {number_unaligned}
            </h3>
            
            
            <div id="aligned" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;"></div>
            <script type="text/javascript">
                vegaEmbed("#aligned", {bowtie_plot});
            </script>
            
            <hr style="width: 100%; color: #FF0000">
            
            <!-- PLOT1 -->
            <h1> This is another plot </h1>
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
    
    total_reads, percent_aligned = parse_bowtielog.parse_alignments(bowtie2log)
    number_aligned = int(total_reads * percent_aligned / 100)
    number_unaligned = total_reads - number_aligned
    
    # Bowtie2 alignment plot:
    bowtie_plot = bowtie2_alignment_plot.plot_alignment(bowtie2log).to_json()
    
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
        total_reads=total_reads,
        number_aligned=number_aligned,
        number_unaligned=number_unaligned,
        plot1=species_and_domain_bracken, 
        plot2=kaiju_bar_plot, 
        svg=svg,
        bowtie_plot=bowtie_plot,
    )

# read in the data (testing sample11 for now)
sample_folder = Path("../virusclassification_nextflow/results/sample11_S6/")


# create report
create_report(sample=sample_folder, out_path=".")