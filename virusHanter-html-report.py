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
                body {{
                font-family: Arial, sans-serif;
                font-size: 16px;
                line-height: 1.5;
                }}
                header {{
                background: linear-gradient(to right, #000000, #333333);
                color: white;
                padding: 40px;
                text-align: center;
                box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0 6px 20px 0 rgba(0, 0, 0, 0.19);
                }}
                h1 {{
                font-size: 3em;
                margin: 0;
                text-align: center;
                }}
                h2 {{
                font-size: 2em;
                margin: 0;
                }}
                section {{
                margin: 80px 0;
                }}
                hr {{
                border: 0;
                height: 1px;
                background: #333;
                background-image: linear-gradient(to right, #ccc, #333, #ccc);
                }}
            </style>

        </head>
        <body>
            <header>
                <h1>
                Pandemic Prepardeness Report 
                </h1>
                <h2>
                Report of {sample_name}
                </h2>
            </header>
            <h1 style="text-align:center;">
                Welcome to The report of the sample! 
            </h1>
            
            

        <div style="border: 1px solid black; padding: 10px; margin: 10px; width: 600px; border-radius: 10px; box-shadow: 5px 5px 5px grey; background-color: #cccccc;">
            <h1>Kraken2 and Kaiju</h1>
            <ul>
            <li>Kraken2 and Kaiju are tools for identifying the taxonomic origin of <strong>DNA sequences</strong></li>
            <li>They use <strong>hierarchical classification systems</strong> and <strong>databases of pre-computed k-mer hashes</strong> to classify query sequences</li>
            <li>Both programs may be used to <strong>rapidly and accurately classify large numbers of sequences</strong>, even in the absence of significant sequence similarity to known reference sequences</li>
            </ul>
        </div>
            
            <!-- PLOT ALIGNED -->
            <h3 style="font-style:italic; text-align:center;"> 
                This is the bowtie plot!
                Number of reads: {total_reads}
                Number aligned to the human genome: {number_aligned}
                Number unaligned: {number_unaligned}
            </h3>
            
            
            <div id="aligned" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;margin:40px;"></div>
            <script type="text/javascript">
                vegaEmbed("#aligned", {bowtie_plot});
            </script>
            
            <hr />
            
            <!-- PLOT1 -->
            <h1> This is another plot </h1>
            <h3 style="font-style:italic; text-align:center;"> 
                This is a plot of a beautiful graph I just made! Enjoy!
            </h3>
            
            <div id="vis" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;margin:40px;"></div>
            <script type="text/javascript">
                vegaEmbed("#vis", {plot1});
            </script>
            
            <!-- PLOT2 -->
            <h3 style="font-style:italic; text-align:center;"> 
                This is a plot of a beautiful graph I just made! Enjoy!
            </h3>
            
            <div id="vis2" style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;margin:40px;"></div>
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
            <div style="display:flex;justify-content:center;align-items:center;width:100%;height:100%;margin:40px;">
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
sample_folder = Path("testdata/sample11_S6")


# create report
create_report(sample=sample_folder, out_path=".")